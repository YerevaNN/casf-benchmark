#!/usr/bin/env python3
"""Freeze and extract a generator-blind PLINDER ATC48 pharmacophore cohort.

The cohort is the ATC/DrugBank-annotated PLINDER 2024-06/v2 table after
removing every identity in the original PLINDER-23 shortlist and molecules
containing elements outside the conformer-method panel.  For each retained
identity, up to eight distinct PDB entries are selected by a seeded SHA-256
ordering.  Reference validity is assessed without loading generator output.

Ligand instance coordinates come from the RCSB ModelServer and source mmCIFs
from RCSB.  These are the deposited coordinates indexed by PLINDER; this route
avoids the retired anonymous HTTP endpoint for the historical PLINDER archives.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import pickle
import time
from pathlib import Path
from typing import Any

import gemmi
import pandas as pd
import requests
from rdkit import Chem

from casf_benchmark.analysis.pharmacophore import (
    canonical_isomeric_smiles,
    molecule_coordinate_hash,
    sha256,
)


DEFAULT_SOURCE = Path("/mnt/weka/hrant/plinder_frequency/2024-06-v2")
DEFAULT_OLD_AUDIT = Path("docs/plinder23_audit/reference_audit.csv")
DEFAULT_OUTPUT = Path("data/pharmacophore_atc48")
ALLOWED_ATOMIC_NUMBERS = {6, 7, 8, 9, 15, 16, 17, 35, 53}


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def stable_rank(seed: str, row: pd.Series) -> str:
    value = "|".join(
        [
            seed,
            str(row.get("inchikey", row.get("inchikeys"))),
            str(row["pdb_id"]),
            str(row["system_id"]),
            str(row["ligand_id"]),
        ]
    )
    return hashlib.sha256(value.encode()).hexdigest()


def download(url: str, destination: Path, session: requests.Session) -> None:
    if destination.is_file() and destination.stat().st_size:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    for attempt in range(5):
        try:
            with session.get(url, timeout=(30, 180), stream=True) as response:
                response.raise_for_status()
                with temporary.open("wb") as handle:
                    for chunk in response.iter_content(1024 * 1024):
                        if chunk:
                            handle.write(chunk)
            if temporary.stat().st_size == 0:
                raise RuntimeError("empty response")
            os.replace(temporary, destination)
            return
        except Exception:
            temporary.unlink(missing_ok=True)
            if attempt == 4:
                raise
            time.sleep(2**attempt)


def covalent_to_ligand(mmcif_path: Path, ligand_asym: str) -> bool:
    block = gemmi.cif.read(str(mmcif_path)).sole_block()
    category = block.get_mmcif_category("_struct_conn.")
    if not category:
        return False
    for values in zip(*category.values()):
        row = dict(zip(category, values, strict=True))
        connection_type = row.get("conn_type_id", "").lower()
        partners = {
            row.get("ptnr1_label_asym_id", ""),
            row.get("ptnr2_label_asym_id", ""),
        }
        if connection_type.startswith("covale") and ligand_asym in partners:
            return True
    return False


def load_ligand(
    path: Path,
    expected_smiles: str,
    *,
    allow_stereo_specialization: bool,
) -> tuple[Chem.Mol, str, str]:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
    molecules = [mol for mol in supplier if mol is not None]
    if len(molecules) != 1:
        raise ValueError(f"expected one readable ligand, found {len(molecules)}")
    mol = molecules[0]
    if mol.GetNumConformers() != 1:
        raise ValueError("ligand instance has no unique conformer")
    actual = canonical_isomeric_smiles(mol)
    if actual == expected_smiles:
        return mol, actual, "exact"
    expected = Chem.MolFromSmiles(expected_smiles)
    actual_nonstereo = Chem.Mol(mol)
    expected_nonstereo = Chem.Mol(expected)
    Chem.RemoveStereochemistry(actual_nonstereo)
    Chem.RemoveStereochemistry(expected_nonstereo)
    same_nonstereo_graph = canonical_isomeric_smiles(
        actual_nonstereo
    ) == canonical_isomeric_smiles(expected_nonstereo)
    expected_stereo_is_preserved = Chem.RemoveHs(mol).HasSubstructMatch(
        expected,
        useChirality=True,
    )
    if allow_stereo_specialization and same_nonstereo_graph and expected_stereo_is_preserved:
        return mol, actual, "coordinate_resolved_unspecified_stereo"
    raise ValueError(f"identity_or_stereo_mismatch:{actual}")


def choose_cohort(source: Path, old_audit: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    atc = pd.read_csv(source / "druglike_atc_candidates.tsv", sep="\t")
    old = pd.read_csv(old_audit, keep_default_na=False)
    old_inchikeys = set(old["inchikey"].astype(str))
    selection_rows = []
    retained = []
    for row in atc.itertuples(index=False):
        mol = Chem.MolFromSmiles(row.smiles)
        if mol is None:
            raise ValueError(f"ATC source contains unreadable SMILES: {row.ccd_resolved}")
        atomic_numbers = sorted({atom.GetAtomicNum() for atom in mol.GetAtoms()})
        reason = ""
        if row.inchikey in old_inchikeys:
            reason = "original_plinder23_identity"
        elif not set(atomic_numbers).issubset(ALLOWED_ATOMIC_NUMBERS):
            reason = "unsupported_element_for_method_panel"
        record = {
            "ccd": row.ccd_resolved,
            "name": row.preferred_name,
            "inchikey": row.inchikey,
            "source_smiles": row.smiles,
            "canonical_smiles": canonical_isomeric_smiles(mol),
            "n_pdbs": int(row.n_pdbs),
            "atomic_numbers": ";".join(map(str, atomic_numbers)),
            "retained": not reason,
            "exclusion_reason": reason,
        }
        selection_rows.append(record)
        if not reason:
            retained.append(record)
    cohort = pd.DataFrame(retained)
    if len(cohort) != 48:
        raise RuntimeError(f"Frozen selection rule produced {len(cohort)} rather than 48")
    if cohort["inchikey"].duplicated().any():
        raise RuntimeError("Frozen cohort contains duplicate identities")
    return cohort, pd.DataFrame(selection_rows)


def build(args: argparse.Namespace) -> None:
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    cohort, selection_audit = choose_cohort(args.source, args.old_audit)
    selection_audit.to_csv(output / "cohort_selection_audit.csv", index=False)

    index = pd.read_parquet(
        args.source / "ligands_per_system.parquet",
        columns=[
            "pdb_id",
            "system_id",
            "ligand_rdkit_canonical_smiles",
            "ligand_ccd_code",
            "ligand_id",
            "inchikeys",
        ],
    )
    index = index[index["inchikeys"].isin(set(cohort["inchikey"]))].copy()
    metadata = cohort.set_index("inchikey").to_dict(orient="index")
    index["selection_rank"] = index.apply(
        lambda row: stable_rank(args.seed, row), axis=1
    )
    index = index.sort_values("selection_rank").drop_duplicates(
        ["inchikeys", "pdb_id"], keep="first"
    )

    raw_dir = output / "raw"
    combined_dir = output / "references"
    pdb_dir = output / "pdb"
    raw_dir.mkdir(exist_ok=True)
    combined_dir.mkdir(exist_ok=True)
    pdb_dir.mkdir(exist_ok=True)
    session = requests.Session()
    session.headers["User-Agent"] = "casf-benchmark-pharmacophore-atc48/1.0"
    audit_rows: list[dict[str, Any]] = []
    attempt_rows: list[dict[str, Any]] = []
    input_mapping: dict[str, dict[str, str]] = {}

    for molecule_number, (inchikey, group) in enumerate(
        index.groupby("inchikeys", sort=True), start=1
    ):
        info = metadata[inchikey]
        source_smiles = str(info["canonical_smiles"])
        expected_smiles: str | None = None
        stereo_resolution = ""
        accepted: list[tuple[pd.Series, Chem.Mol, Path, Path, bool]] = []
        for _, row in group.sort_values("selection_rank").iterrows():
            if len(accepted) >= args.references_per_molecule:
                break
            pdb_id = str(row["pdb_id"]).lower()
            system_id = str(row["system_id"])
            ligand_id = str(row["ligand_id"])
            ligand_instance = ligand_id.rsplit("__", 1)[-1]
            ligand_asym = ligand_instance.split(".", 1)[-1]
            raw_path = raw_dir / str(info["ccd"]) / f"{pdb_id}_{ligand_instance}.sdf"
            mmcif_path = pdb_dir / f"{pdb_id}.cif.gz"
            attempt = {
                "inchikey": inchikey,
                "ccd": info["ccd"],
                "pdb_id": pdb_id,
                "system_id": system_id,
                "ligand_id": ligand_id,
                "ligand_asym": ligand_asym,
                "selection_rank": row["selection_rank"],
                "accepted": False,
                "rejection_reason": "",
            }
            try:
                ligand_url = (
                    f"https://models.rcsb.org/v1/{pdb_id}/ligand?"
                    f"label_asym_id={ligand_asym}&encoding=sdf"
                )
                download(ligand_url, raw_path, session)
                mol, actual_smiles, identity_status = load_ligand(
                    raw_path,
                    expected_smiles or source_smiles,
                    allow_stereo_specialization=expected_smiles is None,
                )
                download(
                    f"https://files.rcsb.org/download/{pdb_id.upper()}.cif.gz",
                    mmcif_path,
                    session,
                )
                covalent = covalent_to_ligand(mmcif_path, ligand_asym)
                if args.exclude_covalent and covalent:
                    raise ValueError("covalent_complex_excluded")
                if expected_smiles is None:
                    expected_smiles = actual_smiles
                    stereo_resolution = identity_status
                accepted.append((row, mol, raw_path, mmcif_path, covalent))
                attempt["accepted"] = True
            except Exception as exc:
                attempt["rejection_reason"] = f"{type(exc).__name__}:{exc}"
            attempt_rows.append(attempt)
        if not accepted:
            raise RuntimeError(f"No valid deposited reference for {info['ccd']} {inchikey}")
        assert expected_smiles is not None

        combined_path = combined_dir / f"{info['ccd']}.sdf"
        with Chem.SDWriter(str(combined_path)) as writer:
            for reference_index, (row, mol, raw_path, mmcif_path, covalent) in enumerate(
                accepted
            ):
                record = Chem.Mol(mol)
                record.SetProp("_Name", f"{info['ccd']}__{row['pdb_id']}")
                record.SetProp("molecule_name", str(info["name"]))
                record.SetProp("ccd", str(info["ccd"]))
                record.SetProp("pdb_id", str(row["pdb_id"]))
                record.SetProp("system_id", str(row["system_id"]))
                record.SetProp("ligand_id", str(row["ligand_id"]))
                writer.write(record)
                ligand_instance = str(row["ligand_id"]).rsplit("__", 1)[-1]
                ligand_asym = ligand_instance.split(".", 1)[-1]
                audit_rows.append(
                    {
                        "ccd": info["ccd"],
                        "pdb_id": str(row["pdb_id"]),
                        "system_id": str(row["system_id"]),
                        "ligand_id": str(row["ligand_id"]),
                        "expected_smiles": expected_smiles,
                        "actual_smiles": canonical_isomeric_smiles(mol),
                        "source_atc_smiles": source_smiles,
                        "stereo_resolution": stereo_resolution,
                        "coordinate_stereo_matches": True,
                        "identity_matches": True,
                        "inchikey": inchikey,
                        "assembly": str(row["system_id"]).split("__")[1],
                        "ligand_asym": ligand_asym,
                        "covalent_flag": covalent,
                        "geometry_pass": False,
                        "geometry_failure_reasons": "not_reference_posebusters_screened",
                        "name": info["name"],
                        "conf_index": reference_index,
                        "reference_index": reference_index,
                        "original_reference_index": reference_index,
                        "raw_sdf_paths": str(raw_path.resolve()),
                        "source_sdf": str(combined_path.resolve()),
                        "coordinate_sha256": molecule_coordinate_hash(mol),
                        "raw_sdf_sha256": sha256(raw_path),
                        "mmcif_path": str(mmcif_path.resolve()),
                        "mmcif_sha256": sha256(mmcif_path),
                        "selection_rank": row["selection_rank"],
                    }
                )
        input_mapping[expected_smiles] = {
            "name": str(info["name"]),
            "ccd": str(info["ccd"]),
        }
        print(
            f"references {molecule_number}/{len(cohort)}: {info['ccd']} "
            f"{len(accepted)}/{args.references_per_molecule}",
            flush=True,
        )

    audit = pd.DataFrame(audit_rows)
    attempts = pd.DataFrame(attempt_rows)
    missing = set(cohort["inchikey"]) - set(audit["inchikey"])
    if missing:
        raise RuntimeError(f"Missing selected identities: {sorted(missing)}")
    audit.to_csv(output / "reference_audit.csv", index=False)
    attempts.to_csv(output / "reference_extraction_attempts.csv", index=False)
    cohort.to_csv(output / "cohort.csv", index=False)
    with (output / "generation_input.pickle").open("wb") as handle:
        pickle.dump(input_mapping, handle, protocol=pickle.HIGHEST_PROTOCOL)
    manifest = {
        "status": "frozen",
        "selection_rule": (
            "all ATC/DrugBank-annotated PLINDER 2024-06/v2 identities excluding "
            "the original PLINDER-23 identities and molecules outside the declared "
            "method-panel element set"
        ),
        "selection_seed": args.seed,
        "pdb_selection": (
            "one indexed ligand instance per distinct PDB, ordered by seeded SHA-256; "
            "first reference-valid noncovalent entries retained"
        ),
        "generator_results_examined": False,
        "molecule_count": int(audit["inchikey"].nunique()),
        "reference_count": len(audit),
        "references_per_molecule_maximum": args.references_per_molecule,
        "molecules_below_reference_cap": int(
            (audit.groupby("inchikey").size() < args.references_per_molecule).sum()
        ),
        "source_files": {
            str(args.source / "druglike_atc_candidates.tsv"): sha256(
                args.source / "druglike_atc_candidates.tsv"
            ),
            str(args.source / "ligands_per_system.parquet"): sha256(
                args.source / "ligands_per_system.parquet"
            ),
            str(args.old_audit): sha256(args.old_audit),
        },
        "outputs": {
            name: sha256(output / name)
            for name in [
                "cohort.csv",
                "cohort_selection_audit.csv",
                "reference_audit.csv",
                "reference_extraction_attempts.csv",
                "generation_input.pickle",
            ]
        },
        "coordinate_source": "RCSB ModelServer ligand instance SDF and RCSB mmCIF",
    }
    atomic_json(output / "freeze_manifest.json", manifest)
    print(json.dumps(manifest, indent=2), flush=True)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    result.add_argument("--old-audit", type=Path, default=DEFAULT_OLD_AUDIT)
    result.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    result.add_argument("--references-per-molecule", type=int, default=8)
    result.add_argument("--seed", default="plinder-atc48-v1")
    result.add_argument("--include-covalent", action="store_false", dest="exclude_covalent")
    result.set_defaults(exclude_covalent=True)
    return result


if __name__ == "__main__":
    build(parser().parse_args())
