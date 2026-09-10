#!/usr/bin/env python3
"""Convert a molgen3D `generation_results.pickle` into the raw external-generation
input contract expected by `casf-materialize-generation-sets`
(see docs/generation_methods.md#external-learned-models):

    {run_root}/generation/manifest.tsv
    {run_root}/generation/{method}/{mol_id}.sdf

`generation_results.pickle` is a dict: SMILES -> list[rdkit.Chem.Mol] (one
conformer per Mol, heavy atoms only). SMILES keys are matched against the
CASF/ChEMBL3D intersection mapping CSV on `chembl3d_isomeric_smiles` (falling
back to `casf_heavy_isomeric_smiles`) to recover the CASF `ligand_id`(s) used as
`mol_id` everywhere downstream. One SMILES/topology can map to several CASF
`ligand_id`s (e.g. the same ChEMBL3D compound co-crystallized in multiple PDB
entries, common in the ref set) -- the same generated conformers are written
under every matching ligand_id.
"""

from __future__ import annotations

import argparse
import csv
import pickle
from collections import defaultdict
from pathlib import Path

from rdkit.Chem import SDWriter


def load_mapping(chembl_map_csv: Path) -> dict[str, list[str]]:
    smi_to_ligand_ids: dict[str, list[str]] = defaultdict(list)
    with open(chembl_map_csv, newline="") as f:
        for row in csv.DictReader(f):
            smi = row.get("chembl3d_isomeric_smiles") or row.get("casf_heavy_isomeric_smiles")
            if smi:
                smi_to_ligand_ids[smi].append(row["ligand_id"])
    return dict(smi_to_ligand_ids)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--pickle", type=Path, required=True)
    parser.add_argument("--chembl-map-csv", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--method", required=True, help="Untiered source method name, e.g. qwen_step29600_core")
    args = parser.parse_args()

    smi_to_ligand_ids = load_mapping(args.chembl_map_csv)
    with open(args.pickle, "rb") as f:
        data = pickle.load(f)

    generation_dir = args.run_root / "generation"
    method_dir = generation_dir / args.method
    method_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    unmatched = []
    for smi, mols in data.items():
        ligand_ids = smi_to_ligand_ids.get(smi)
        if not ligand_ids:
            unmatched.append(smi)
            continue
        for ligand_id in ligand_ids:
            sdf_path = method_dir / f"{ligand_id}.sdf"
            writer = SDWriter(str(sdf_path))
            for conf_id, mol in enumerate(mols):
                mol.SetProp("_Name", ligand_id)
                mol.SetProp("mol_id", ligand_id)
                mol.SetProp("input_smiles", smi)
                mol.SetProp("generation_method", args.method)
                mol.SetProp("conf_id", str(conf_id))
                writer.write(mol)
            writer.close()
            manifest_rows.append(
                {
                    "mol_id": ligand_id,
                    "generation_method": args.method,
                    "input_smiles": smi,
                    "kept_confs": len(mols),
                }
            )

    if unmatched:
        raise SystemExit(
            f"{len(unmatched)} pickle SMILES had no mapping-CSV match, e.g. {unmatched[:5]}"
        )

    manifest_path = generation_dir / "manifest.tsv"
    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["mol_id", "generation_method", "input_smiles", "kept_confs"], delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"Wrote {len(manifest_rows)} ligand SDF(s) under {method_dir}")
    print(f"Wrote manifest: {manifest_path}")


if __name__ == "__main__":
    main()
