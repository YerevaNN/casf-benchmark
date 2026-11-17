"""Run the staged PLINDER-22 pharmacophore-query matching pilot."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import pickle
import shutil
import subprocess
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
import yaml

from casf_benchmark.analysis.pharmacophore import (
    PharmitBackend,
    atomic_json_dump,
    backend_smoke_test,
    canonical_isomeric_smiles,
    distinct_centers,
    exact_feature_correspondence_transform,
    explicit_atom_mapping,
    exact_feature_membership,
    export_candidate_library,
    export_normalized_feature_library,
    hit_probability,
    load_pickle,
    maximum_triangle_area,
    mixture_hit_probability,
    molecule_coordinate_hash,
    paired_bootstrap_difference,
    read_match_titles,
    reconstruct_native_receptor,
    sha256,
    stable_query_id,
    stable_smiles_id,
    transform_query,
)


DEFAULT_CONFIG = Path(__file__).resolve().parents[1] / "config" / "pharmacophore_plinder22.yaml"


def load_config(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    config["_path"] = str(path.resolve())
    config["_sha256"] = sha256(path)
    return config


def output_dir(config: dict[str, Any]) -> Path:
    path = Path(config["output_dir"])
    path.mkdir(parents=True, exist_ok=True)
    return path


def backend(config: dict[str, Any]) -> PharmitBackend:
    settings = config["backend"]
    executable = Path(settings["executable"])
    if not executable.is_file():
        raise FileNotFoundError(f"Pharmit executable not found: {executable}")
    return PharmitBackend(
        executable=executable,
        timeout_seconds=int(settings["timeout_seconds"]),
        skip_triplet_fingerprints=bool(settings["skip_triplet_fingerprints"]),
        skip_geometric_pruning=bool(settings["skip_geometric_pruning"]),
    )


def bool_value(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def radius_tag(radius: float) -> str:
    return f"r{radius:g}".replace(".", "p")


def write_frame_atomic(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    frame.to_csv(temporary, index=False)
    os.replace(temporary, path)


def conda_package_inventory(environment: Path) -> list[dict[str, Any]]:
    """Read package provenance without requiring conda or mamba on the worker."""
    metadata_dir = environment / "conda-meta"
    metadata_paths = sorted(metadata_dir.glob("*.json"))
    if not metadata_paths:
        raise FileNotFoundError(f"No conda package metadata found in {metadata_dir}")
    fields = ("name", "version", "build", "build_number", "channel", "subdir", "url")
    packages = []
    for path in metadata_paths:
        record = json.loads(path.read_text(encoding="utf-8"))
        if not record.get("name") or not record.get("version"):
            raise ValueError(f"Incomplete conda package metadata: {path}")
        packages.append({field: record.get(field) for field in fields})
    return sorted(
        packages,
        key=lambda record: (record["name"], record["version"], record["build"] or ""),
    )


def stage_inventory(config: dict[str, Any]) -> None:
    out = output_dir(config)
    audit = pd.read_csv(config["reference_audit"], keep_default_na=False)
    expected_reference_count = config.get("expected_reference_count")
    if expected_reference_count is not None and len(audit) != int(expected_reference_count):
        raise ValueError(
            f"Expected {int(expected_reference_count):,} references, found {len(audit):,}"
        )
    records = []
    expected_smiles = set(audit["expected_smiles"])
    for method, settings in config["methods"].items():
        path = Path(settings["path"])
        pool = load_pickle(path)
        missing = expected_smiles - set(pool)
        if missing:
            raise ValueError(f"{method} is missing {len(missing)} PLINDER molecules")
        extra = set(pool) - expected_smiles
        for smiles in sorted(expected_smiles):
            records.append(
                {
                    "method": method,
                    "label": settings["label"],
                    "path": str(path),
                    "path_sha256": sha256(path),
                    "smiles": smiles,
                    "smiles_id": stable_smiles_id(smiles),
                    "stored_entries": len(pool[smiles]),
                    "extra_pool_identities": len(extra),
                }
            )
    inventory = pd.DataFrame(records)
    primary_counts = inventory.groupby("method")["stored_entries"].agg(["min", "max", "sum"])
    if (primary_counts["min"] < 1).any():
        raise ValueError(f"At least one supplied pool is empty:\n{primary_counts}")
    for method, settings in config["methods"].items():
        counts = inventory.loc[inventory["method"] == method, "stored_entries"]
        minimum = int(settings.get("minimum_pool_size", 1))
        maximum = settings.get("maximum_pool_size")
        if int(counts.min()) < minimum:
            raise ValueError(
                f"{method} has a pool smaller than its declared minimum {minimum}: "
                f"{int(counts.min())}"
            )
        if maximum is not None and int(counts.max()) > int(maximum):
            raise ValueError(
                f"{method} has a pool larger than its declared maximum {int(maximum)}: "
                f"{int(counts.max())}"
            )
    write_frame_atomic(inventory, out / "pool_inventory_completed.csv")
    original_hashes: dict[str, str] | None = None
    original_prefix_verified = False
    if "original_qwen" in config:
        original = config["original_qwen"]
        original_hashes = {
            "cache": sha256(Path(original["path"])),
            "metadata": sha256(Path(original["metadata_path"])),
        }
        if original_hashes != {
            "cache": original["sha256"],
            "metadata": original["metadata_sha256"],
        }:
            raise RuntimeError("Original Qwen artifacts no longer match their immutable hashes")
        original_pool = load_pickle(Path(original["path"]))
        qwen_method = config["analysis"]["mixture_sources"][0]
        completed_pool = load_pickle(Path(config["methods"][qwen_method]["path"]))
        for smiles, original_molecules in original_pool.items():
            before = [molecule_coordinate_hash(mol) for mol in original_molecules]
            after = [
                molecule_coordinate_hash(mol)
                for mol in completed_pool[smiles][: len(original_molecules)]
            ]
            if before != after:
                raise RuntimeError(f"Completed Qwen cache changed an original record: {smiles}")
        original_prefix_verified = True
    source_root = Path("/mnt/weka/mbedrosian/tool-sources/casf-pharmit-pilot")
    source_directories = {
        "pharmit": "pharmit",
        "openbabel": "openbabel",
        "smina": "smina-sourceforge",
        "fastcgi": "fcgi2",
        "bitmagic": "BitMagic",
    }
    source_commits = {}
    for name, directory in source_directories.items():
        source_commits[name] = subprocess.run(
            ["git", "-C", str(source_root / directory), "rev-parse", "HEAD"],
            check=True,
            text=True,
            capture_output=True,
        ).stdout.strip()
    packages = conda_package_inventory(Path(config["backend"]["environment"]))
    linked_libraries = subprocess.run(
        ["ldd", str(config["backend"]["executable"])],
        check=True,
        text=True,
        capture_output=True,
    ).stdout.splitlines()
    build_manifest = {
        "reproduction_command": "bash scripts/build_pharmit_pilot.sh",
        "isolated_environment": config["backend"]["environment"],
        "source_root": str(source_root),
        "source_commits": source_commits,
        "executable": config["backend"]["executable"],
        "executable_sha256": sha256(Path(config["backend"]["executable"])),
        "packages": packages,
        "linked_libraries": linked_libraries,
        "compatibility_patches": [
            "Open Babel obutil.h: include <ctime>",
            "Pharmit SpinLock.h: include <chrono>",
            "Pharmit CMake: conda jsoncpp_lib target",
            "LEMON umbrella header: expose every enabled solver API",
            "Pharmit BitMagic source replaced by pinned upstream commit",
            "Pharmit exact-search option: disable non-conservative higher-order triplet pruning",
        ],
        "feature_definitions_changed_by_patches": False,
        "exact_search_semantics_patch": True,
    }
    atomic_json_dump(out / "backend_build_manifest.json", build_manifest)
    manifest = {
        "status": "passed",
        "config_path": config["_path"],
        "config_sha256": config["_sha256"],
        "reference_count": len(audit),
        "molecule_count": len(expected_smiles),
        "method_count": len(config["methods"]),
        "candidate_count": int(inventory["stored_entries"].sum()),
        "pool_size_summary": primary_counts.reset_index().to_dict(orient="records"),
        "budget_semantics": (
            "accepted cached conformers, not attempts or compute; requested draws are "
            "capped at the supplied method-molecule pool size"
        ),
        "original_qwen_hashes": original_hashes,
        "original_qwen_prefix_topology_and_double_coordinates_verified": (
            original_prefix_verified
        ),
        "backend": {
            **config["backend"],
            "executable_sha256": sha256(Path(config["backend"]["executable"])),
        },
        "backend_build_manifest_sha256": sha256(out / "backend_build_manifest.json"),
    }
    atomic_json_dump(out / "inventory_manifest.json", manifest)
    print(json.dumps(manifest, indent=2), flush=True)


def stage_smoke(config: dict[str, Any]) -> None:
    report = backend_smoke_test(backend(config), output_dir(config) / "backend_smoke")
    print(json.dumps(report, indent=2), flush=True)


def _load_sdf_records(path: str, cache: dict[str, list[Any]]) -> list[Any]:
    from rdkit import Chem

    if path not in cache:
        # Keep physical SDF record positions intact: ``reference_index`` is a
        # row index, so filtering an earlier unreadable record would silently
        # attach a query to the wrong conformer.
        cache[path] = list(Chem.SDMolSupplier(path, removeHs=False, sanitize=True))
    return cache[path]


def _load_single_sdf(path: str) -> Any:
    from rdkit import Chem

    records = [
        mol
        for mol in Chem.SDMolSupplier(path, removeHs=False, sanitize=True)
        if mol is not None
    ]
    if len(records) != 1:
        raise ValueError(f"Expected one molecule in {path}, found {len(records)}")
    return records[0]


def _query_feature_rows(query_id: str, query: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for index, point in enumerate(query.get("points", [])):
        rows.append(
            {
                "query_id": query_id,
                "feature_index": index,
                "feature_type": point["name"],
                "x": point["x"],
                "y": point["y"],
                "z": point["z"],
                "radius": point["radius"],
                "vector_on": point["vector_on"],
                "enabled": point["enabled"],
                "size": point.get("size", ""),
            }
        )
    return rows


def _plot_overlay(
    *,
    ligand: Any,
    receptor_path: Path,
    query: dict[str, Any],
    title: str,
    output_path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ligand_xyz = np.asarray(ligand.GetConformer().GetPositions(), dtype=float)
    receptor_xyz = []
    for line in receptor_path.read_text().splitlines():
        if line.startswith("ATOM"):
            receptor_xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    receptor_xyz = np.asarray(receptor_xyz, dtype=float)
    points = query["points"]
    colors = {
        "Aromatic": "#8e44ad",
        "HydrogenAcceptor": "#e74c3c",
        "HydrogenDonor": "#3498db",
        "Hydrophobic": "#f39c12",
        "NegativeIon": "#c0392b",
        "PositiveIon": "#2980b9",
    }
    figure = plt.figure(figsize=(7, 6))
    axis = figure.add_subplot(111, projection="3d")
    if len(receptor_xyz):
        distance = np.linalg.norm(
            receptor_xyz[:, None, :] - ligand_xyz[None, :, :], axis=2
        ).min(axis=1)
        close = receptor_xyz[distance <= 8.0]
        axis.scatter(*close.T, s=5, c="#bdc3c7", alpha=0.18, label="protein ≤8 Å")
    axis.scatter(*ligand_xyz.T, s=24, c="#2c3e50", alpha=0.8, label="native ligand")
    seen = set()
    for point in points:
        label = point["name"] if point["name"] not in seen else None
        seen.add(point["name"])
        axis.scatter(
            point["x"], point["y"], point["z"],
            s=100, c=colors.get(point["name"], "#16a085"), marker="o", label=label,
            edgecolors="black", linewidths=0.4,
        )
    axis.set_title(title)
    axis.set_xlabel("x (Å)")
    axis.set_ylabel("y (Å)")
    axis.set_zlabel("z (Å)")
    axis.legend(fontsize=7, loc="best")
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=160)
    plt.close(figure)


def _run_native_self_match(
    *,
    config: dict[str, Any],
    manifest: pd.DataFrame,
    native_molecules: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    """Group native references by molecule and require self/reoriented recovery."""
    from rdkit import Chem

    out = output_dir(config)
    engine = backend(config)
    primary_radius = float(config["query"]["primary_radius"])
    tag = radius_tag(primary_radius)
    qa_root = out / "reference_qa"
    results: dict[str, dict[str, Any]] = {}
    angle = 0.613
    rotation = np.asarray(
        [
            [math.cos(angle), 0.0, math.sin(angle)],
            [0.0, 1.0, 0.0],
            [-math.sin(angle), 0.0, math.cos(angle)],
        ]
    )
    translation = np.asarray([11.0, -4.5, 2.25])
    candidates = manifest[manifest["reference_eligible_pre_self_match"]].copy()
    for smiles, group in candidates.groupby("expected_smiles", sort=True):
        sid = stable_smiles_id(smiles)
        group_dir = qa_root / sid
        group_dir.mkdir(parents=True, exist_ok=True)
        library = group_dir / "native_library.sdf"
        with Chem.SDWriter(str(library)) as writer:
            for query_id in group["query_id"]:
                mol = Chem.Mol(native_molecules[query_id])
                mol.SetProp("_Name", query_id)
                writer.write(mol)
        database = group_dir / "database"
        engine.create_database(library, database)
        query_paths = [out / "queries" / tag / f"{qid}.json" for qid in group["query_id"]]
        transformed_paths = []
        for query_id, query_path in zip(group["query_id"], query_paths, strict=True):
            transformed = transform_query(json.loads(query_path.read_text()), rotation, translation)
            transformed_path = group_dir / f"{query_id}.transformed.json"
            atomic_json_dump(transformed_path, transformed)
            transformed_paths.append(transformed_path)
        original_outputs = [group_dir / f"{qid}.self.sdf" for qid in group["query_id"]]
        transformed_outputs = [
            group_dir / f"{qid}.transformed.sdf" for qid in group["query_id"]
        ]
        batch_size = int(config["backend"]["query_batch_size"])
        for start in range(0, len(query_paths), batch_size):
            stop = start + batch_size
            engine.search_many(
                database_path=database,
                query_paths=query_paths[start:stop],
                output_paths=original_outputs[start:stop],
                max_hits=len(group) + 1,
                workers=int(config["backend"]["workers"]),
            )
            engine.search_many(
                database_path=database,
                query_paths=transformed_paths[start:stop],
                output_paths=transformed_outputs[start:stop],
                max_hits=len(group) + 1,
                workers=int(config["backend"]["workers"]),
            )
        for query_id, original_path, transformed_path in zip(
            group["query_id"], original_outputs, transformed_outputs, strict=True
        ):
            self_titles = read_match_titles(original_path)
            transformed_titles = read_match_titles(transformed_path)
            results[query_id] = {
                "self_match_pass": query_id in self_titles,
                "rigid_invariance_pass": query_id in transformed_titles,
                "self_match_count": len(self_titles),
                "transformed_match_count": len(transformed_titles),
            }
            original_path.unlink()
            transformed_path.unlink()
    return results


def stage_queries(config: dict[str, Any], *, limit: int | None = None) -> None:
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    out = output_dir(config)
    query_settings = config["query"]
    audit = pd.read_csv(config["reference_audit"], keep_default_na=False)
    if limit is not None:
        audit = audit.iloc[:limit].copy()
    engine = backend(config)
    radii = [
        float(query_settings["primary_radius"]),
        *map(float, query_settings["sensitivity_radii"]),
    ]
    cleaned_cache: dict[str, list[Any]] = {}
    native_molecules: dict[str, Any] = {}
    manifest_rows = []
    feature_rows = []
    maps_dir = out / "reference_maps"
    receptors_dir = out / "native_receptors"
    errors_dir = out / "query_errors"
    errors_dir.mkdir(parents=True, exist_ok=True)
    for stale_error in errors_dir.glob("*.txt"):
        stale_error.unlink()
    preview_names = {"Actinonin", "Imatinib", "Amprenavir", "Nicotine", "Fosmidomycin"}
    previewed: set[str] = set()

    for position, row in enumerate(audit.itertuples(index=False), start=1):
        query_id = stable_query_id(row)
        record = {
            "query_id": query_id,
            "name": row.name,
            "ccd": row.ccd,
            "conf_index": int(row.conf_index),
            "reference_index": int(row.reference_index),
            "original_reference_index": int(row.original_reference_index),
            "expected_smiles": row.expected_smiles,
            "smiles_id": stable_smiles_id(row.expected_smiles),
            "pdb_id": row.pdb_id,
            "system_id": row.system_id,
            "ligand_id": row.ligand_id,
            "ligand_asym": row.ligand_asym,
            "assembly": str(row.assembly),
            "raw_sdf_path": row.raw_sdf_paths,
            "cleaned_sdf_path": row.source_sdf,
            "covalent_flag": bool_value(row.covalent_flag),
            "reference_geometry_pass": bool_value(row.geometry_pass),
            "reference_geometry_failure_reasons": row.geometry_failure_reasons,
            "identity_matches_audit": bool_value(row.identity_matches),
            "stereo_matches_audit": bool_value(row.coordinate_stereo_matches),
            "reference_eligible_pre_self_match": False,
            "eligible": False,
            "exclusion_reason": "",
        }
        try:
            if record["covalent_flag"] and bool(query_settings["exclude_covalent"]):
                raise ValueError("covalent_complex_excluded")
            if not record["identity_matches_audit"]:
                raise ValueError("audit_identity_mismatch")
            if not record["stereo_matches_audit"]:
                raise ValueError("audit_stereo_mismatch")
            native = _load_single_sdf(row.raw_sdf_paths)
            if canonical_isomeric_smiles(native) != row.expected_smiles:
                raise ValueError("reverified_native_identity_mismatch")
            cleaned_records = _load_sdf_records(row.source_sdf, cleaned_cache)
            if int(row.reference_index) >= len(cleaned_records):
                raise ValueError("cleaned_conformer_index_out_of_range")
            cleaned = cleaned_records[int(row.reference_index)]
            if cleaned is None:
                raise ValueError("cleaned_conformer_unreadable")
            if canonical_isomeric_smiles(cleaned) != row.expected_smiles:
                raise ValueError("reverified_cleaned_identity_mismatch")
            atom_map = explicit_atom_mapping(cleaned, native)
            if atom_map["rmsd"] > float(query_settings["cleaned_native_max_rmsd"]):
                raise ValueError("cleaned_native_mapping_rmsd_exceeded")
            receptor_path = receptors_dir / f"{query_id}.pdb"
            native_map = reconstruct_native_receptor(
                mmcif_path=Path(config["pdb_cache"]) / f"{row.pdb_id}.cif.gz",
                assembly_name=str(row.assembly),
                ligand_asym=row.ligand_asym,
                native_mol=native,
                receptor_path=receptor_path,
                contact_cutoff=float(query_settings["contact_cutoff"]),
            )
            if native_map["rmsd_without_superposition"] > float(
                query_settings["native_copy_max_rmsd"]
            ):
                raise ValueError("native_assembly_copy_rmsd_exceeded")
            provenance = {"cleaned_native": atom_map, "native_assembly": native_map}
            atomic_json_dump(maps_dir / f"{query_id}.json", provenance)
            primary_path = out / "queries" / radius_tag(radii[0]) / f"{query_id}.json"
            query = engine.extract_interaction_query(
                Path(row.raw_sdf_paths),
                receptor_path,
                primary_path,
                radii[0],
            )
            # Pharmit serializes an empty interaction set as JSON null rather than
            # an empty array.  Normalize it before eligibility calculations and
            # before writing the radius-sensitivity variants.
            points = query.get("points") or []
            query["points"] = points
            centers = distinct_centers(
                points, float(query_settings["minimum_center_separation"])
            )
            triangle_area = maximum_triangle_area(centers)
            record.update(
                {
                    "native_copy_rmsd": native_map["rmsd_without_superposition"],
                    "cleaned_native_mapping_rmsd": atom_map["rmsd"],
                    "n_contacting_protein_chains": len(native_map["protein_chains"]),
                    "nearby_metals": ";".join(native_map["context"]["metals"]),
                    "nearby_waters": len(native_map["context"]["waters"]),
                    "nearby_cofactors": ";".join(native_map["context"]["cofactors"]),
                    "nearby_nucleic_acids": ";".join(
                        native_map["context"]["nucleic_acids"]
                    ),
                    "feature_count": len(points),
                    "distinct_center_count": len(centers),
                    "maximum_triangle_area": triangle_area,
                    "feature_types": ";".join(sorted({point["name"] for point in points})),
                }
            )
            query["pilot_provenance"] = {
                "query_id": query_id,
                "source_ligand_sha256": sha256(Path(row.raw_sdf_paths)),
                "receptor_sha256": sha256(receptor_path),
                "backend_commit": config["backend"]["commit"],
                "positional_only": True,
            }
            atomic_json_dump(primary_path, query)
            for radius in radii[1:]:
                sensitivity = json.loads(json.dumps(query))
                for point in sensitivity["points"]:
                    point["radius"] = radius
                atomic_json_dump(
                    out / "queries" / radius_tag(radius) / f"{query_id}.json",
                    sensitivity,
                )
            feature_rows.extend(_query_feature_rows(query_id, query))
            if len(centers) < int(query_settings["minimum_distinct_centers"]):
                raise ValueError("fewer_than_four_distinct_contact_centers")
            if triangle_area < float(query_settings["minimum_triangle_area"]):
                raise ValueError("contact_centers_collinear")
            record["reference_eligible_pre_self_match"] = True
            native_molecules[query_id] = native
            if row.name in preview_names and row.name not in previewed:
                _plot_overlay(
                    ligand=native,
                    receptor_path=receptor_path,
                    query=query,
                    title=f"{row.name}: {row.pdb_id} / {query_id}",
                    output_path=out / "overlays" / f"{row.name.replace(' ', '_')}_{query_id}.png",
                )
                previewed.add(row.name)
        except Exception as exc:
            record["exclusion_reason"] = str(exc)
            (errors_dir / f"{query_id}.txt").write_text(
                f"{type(exc).__name__}: {exc}\n", encoding="utf-8"
            )
        manifest_rows.append(record)
        if position % 50 == 0:
            print(f"queries: {position}/{len(audit)}", flush=True)
            write_frame_atomic(pd.DataFrame(manifest_rows), out / "query_manifest.partial.csv")

    manifest = pd.DataFrame(manifest_rows)
    if limit is None:
        qa = _run_native_self_match(
            config=config,
            manifest=manifest,
            native_molecules=native_molecules,
        )
        manifest["self_match_pass"] = manifest["query_id"].map(
            lambda query_id: qa.get(query_id, {}).get("self_match_pass", False)
        )
        manifest["rigid_invariance_pass"] = manifest["query_id"].map(
            lambda query_id: qa.get(query_id, {}).get("rigid_invariance_pass", False)
        )
        manifest["self_match_count"] = manifest["query_id"].map(
            lambda query_id: qa.get(query_id, {}).get("self_match_count", 0)
        )
        manifest["transformed_match_count"] = manifest["query_id"].map(
            lambda query_id: qa.get(query_id, {}).get("transformed_match_count", 0)
        )
        self_fail = manifest["reference_eligible_pre_self_match"] & ~(
            manifest["self_match_pass"] & manifest["rigid_invariance_pass"]
        )
        manifest.loc[self_fail, "exclusion_reason"] = "reference_self_match_or_invariance_failure"
        manifest["eligible"] = manifest["reference_eligible_pre_self_match"] & ~self_fail
    else:
        manifest["eligible"] = manifest["reference_eligible_pre_self_match"]
    write_frame_atomic(manifest, out / "query_manifest.csv")
    write_frame_atomic(pd.DataFrame(feature_rows), out / "query_features.csv")
    exclusions = (
        manifest.loc[~manifest["eligible"], "exclusion_reason"]
        .replace("", "unclassified")
        .value_counts()
        .rename_axis("exclusion_reason")
        .reset_index(name="count")
    )
    write_frame_atomic(exclusions, out / "query_exclusions.csv")
    freeze = {
        "status": "frozen" if limit is None else "limited_debug_run",
        "config_sha256": config["_sha256"],
        "total_references": len(manifest),
        "eligible_queries": int(manifest["eligible"].sum()),
        "eligible_molecules": int(manifest.loc[manifest["eligible"], "expected_smiles"].nunique()),
        "exclusions": dict(zip(exclusions["exclusion_reason"], exclusions["count"], strict=True)),
        "radii": radii,
        "positional_only": True,
        "generator_results_examined_for_eligibility": False,
        "query_manifest_sha256": sha256(out / "query_manifest.csv"),
    }
    atomic_json_dump(out / "query_freeze_manifest.json", freeze)
    (out / "query_manifest.partial.csv").unlink(missing_ok=True)
    print(json.dumps(freeze, indent=2), flush=True)


def _write_single_witness(source: Path, destination: Path) -> None:
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(source), removeHs=False)
    first = next((mol for mol in supplier if mol is not None), None)
    if first is None:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(destination)) as writer:
        writer.write(first)


def _audit_candidate_export(
    source_molecules: Sequence[Any],
    records: Sequence[dict[str, Any]],
    sdf_path: Path,
) -> dict[str, Any]:
    """Verify exported record identity, order, and SDF coordinate precision."""
    from rdkit import Chem

    exported_records = [record for record in records if record["exported"]]
    exported_molecules = list(
        Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
    )
    if len(exported_molecules) != len(exported_records):
        raise ValueError(f"Exported SDF record count changed: {sdf_path}")
    maximum_delta = 0.0
    for exported, record in zip(exported_molecules, exported_records, strict=True):
        if exported is None:
            raise ValueError(f"Unreadable exported candidate: {record['candidate_id']}")
        if exported.GetProp("_Name") != record["candidate_id"]:
            raise ValueError("Candidate title/order changed during SDF export")
        source = source_molecules[int(record["candidate_index"])]
        if source.GetNumAtoms() != exported.GetNumAtoms():
            raise ValueError("Candidate atom count changed during SDF export")
        source_atomic_numbers = [atom.GetAtomicNum() for atom in source.GetAtoms()]
        exported_atomic_numbers = [atom.GetAtomicNum() for atom in exported.GetAtoms()]
        if source_atomic_numbers != exported_atomic_numbers:
            raise ValueError("Candidate atom order changed during SDF export")
        delta = np.abs(
            np.asarray(source.GetConformer().GetPositions(), dtype=float)
            - np.asarray(exported.GetConformer().GetPositions(), dtype=float)
        )
        maximum_delta = max(maximum_delta, float(delta.max(initial=0.0)))
    if maximum_delta > 0.001:
        raise ValueError(
            f"SDF serialization changed a coordinate by {maximum_delta:.6g} Å"
        )
    return {
        "record_count": len(exported_records),
        "titles_and_atom_order_verified": True,
        "maximum_absolute_coordinate_serialization_delta_angstrom": maximum_delta,
        "coordinate_serialization_tolerance_angstrom": 0.001,
    }


_EXACT_FEATURE_COORDINATES: Sequence[np.ndarray] | np.ndarray | None = None
_EXACT_FEATURE_TYPES: Sequence[str] | Sequence[Sequence[str]] = ()


def _initialize_exact_feature_worker(
    coordinates: Sequence[np.ndarray] | np.ndarray,
    feature_types: Sequence[str] | Sequence[Sequence[str]],
) -> None:
    global _EXACT_FEATURE_COORDINATES, _EXACT_FEATURE_TYPES
    _EXACT_FEATURE_COORDINATES = coordinates
    _EXACT_FEATURE_TYPES = feature_types


def _exact_feature_worker(
    task: tuple[str, float, list[dict[str, Any]]]
) -> tuple[str, float, np.ndarray]:
    query_id, radius, points = task
    if _EXACT_FEATURE_COORDINATES is None:
        raise RuntimeError("Exact feature worker was not initialized")
    matches = exact_feature_membership(
        points,
        _EXACT_FEATURE_COORDINATES,
        _EXACT_FEATURE_TYPES,
        radius,
    )
    return query_id, radius, matches


def _prepare_candidate_feature_cache(
    *,
    config: dict[str, Any],
    engine: PharmitBackend,
    method: str,
    smiles: str,
    mols: Sequence[Any],
    candidate_records: Sequence[dict[str, Any]],
) -> tuple[
    Sequence[np.ndarray] | np.ndarray,
    Sequence[str] | Sequence[Sequence[str]],
    list[str],
    dict[str, Any],
]:
    """Recognize Pharmit features once and freeze them in a compact cache."""
    out = output_dir(config)
    sid = stable_smiles_id(smiles)
    exported = [record for record in candidate_records if record["exported"]]
    exported_ids = [str(record["candidate_id"]) for record in exported]
    cache_path = out / "candidate_features" / method / f"{sid}.npz"
    metadata_path = cache_path.with_suffix(".metadata.json")
    backend_hash = sha256(engine.executable)
    if cache_path.is_file() and metadata_path.is_file():
        metadata = json.loads(metadata_path.read_text())
        with np.load(cache_path, allow_pickle=True) as archive:
            coordinates = archive["coordinates"]
            stored_types = archive["feature_types"]
            if stored_types.ndim == 1 and all(
                isinstance(value, str) for value in stored_types.tolist()
            ):
                feature_types: Sequence[str] | Sequence[Sequence[str]] = list(
                    map(str, stored_types)
                )
            else:
                feature_types = [tuple(map(str, values)) for values in stored_types]
            stored_ids = list(map(str, archive["candidate_ids"]))
        if (
            metadata.get("backend_sha256") == backend_hash
            and stored_ids == exported_ids
            and len(coordinates) == len(exported_ids)
        ):
            return coordinates, feature_types, stored_ids, metadata

    feature_sdf = out / "feature_libraries" / method / f"{sid}.sdf"
    feature_molecules = [mols[int(record["candidate_index"])] for record in exported]
    export_audit = export_normalized_feature_library(
        smiles=smiles,
        mols=feature_molecules,
        candidate_ids=exported_ids,
        sdf_path=feature_sdf,
    )
    raw_path = out / "candidate_features" / method / f"{sid}.recognition.json"
    feature_sets = engine.recognize_features(feature_sdf, raw_path)
    if len(feature_sets) != len(exported_ids):
        raise RuntimeError("Pharmit feature recognition changed the candidate count")
    # Pharmit may merge nearby hydrophobic groups, so feature counts can vary
    # with conformation. Preserve every per-conformer signature and coordinate
    # set; the exact matcher already supports ragged typed feature collections.
    signatures = [tuple(str(point["name"]) for point in points) for points in feature_sets]
    if not signatures:
        raise RuntimeError("Candidate Pharmit feature recognition returned no records")
    coordinate_sets = [
        np.asarray([[point[axis] for axis in "xyz"] for point in points], dtype=float)
        for points in feature_sets
    ]
    stable_signature = all(signature == signatures[0] for signature in signatures[1:])
    if stable_signature:
        feature_types = list(signatures[0])
        coordinates: Sequence[np.ndarray] | np.ndarray = np.asarray(
            coordinate_sets, dtype=float
        )
        stored_types = np.asarray(feature_types, dtype=object)
    else:
        coordinates_array = np.empty(len(coordinate_sets), dtype=object)
        types_array = np.empty(len(signatures), dtype=object)
        for index, (coordinate_set, signature) in enumerate(
            zip(coordinate_sets, signatures, strict=True)
        ):
            coordinates_array[index] = coordinate_set
            types_array[index] = signature
        coordinates = coordinates_array
        feature_types = list(signatures)
        stored_types = types_array
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache_path,
        coordinates=coordinates,
        feature_types=stored_types,
        candidate_ids=np.asarray(exported_ids, dtype=object),
    )
    raw_path.unlink(missing_ok=True)
    metadata = {
        "method": method,
        "smiles": smiles,
        "smiles_id": sid,
        "candidate_count": len(exported_ids),
        "minimum_feature_count": min(map(len, signatures)),
        "maximum_feature_count": max(map(len, signatures)),
        "unique_feature_type_signatures": len(set(signatures)),
        "signature_stable_across_conformers": stable_signature,
        "recognition_backend": "Pharmit pharma command",
        "backend_commit": config["backend"]["commit"],
        "backend_sha256": backend_hash,
        "feature_definitions_changed": False,
        "normalized_export": export_audit,
        "cache_path": str(cache_path),
        "cache_sha256": sha256(cache_path),
    }
    atomic_json_dump(metadata_path, metadata)
    return coordinates, feature_types, exported_ids, metadata


def _match_group_exact_features(
    *,
    config: dict[str, Any],
    engine: PharmitBackend,
    method: str,
    smiles: str,
    mols: Sequence[Any],
    candidate_records: list[dict[str, Any]],
    query_group: pd.DataFrame,
    radii: Sequence[float],
) -> list[dict[str, Any]]:
    """Compute complete memberships without Pharmit's lossy triplet index."""
    out = output_dir(config)
    sid = stable_smiles_id(smiles)
    coordinates, feature_types, exported_ids, feature_metadata = (
        _prepare_candidate_feature_cache(
            config=config,
            engine=engine,
            method=method,
            smiles=smiles,
            mols=mols,
            candidate_records=candidate_records,
        )
    )
    all_ids = [str(record["candidate_id"]) for record in candidate_records]
    index_by_id = {candidate_id: index for index, candidate_id in enumerate(all_ids)}
    query_ids = query_group["query_id"].astype(str).tolist()
    tasks: list[tuple[str, float, list[dict[str, Any]]]] = []
    for query_id in query_ids:
        primary_query = json.loads(
            (out / "queries" / radius_tag(float(radii[0])) / f"{query_id}.json").read_text()
        )
        for radius in radii:
            tasks.append((query_id, float(radius), primary_query["points"]))
    started = time.perf_counter()
    workers = min(int(config["backend"]["workers"]), max(1, len(tasks)))
    if workers == 1:
        _initialize_exact_feature_worker(coordinates, feature_types)
        completed = list(map(_exact_feature_worker, tasks))
    else:
        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=_initialize_exact_feature_worker,
            initargs=(coordinates, feature_types),
        ) as executor:
            completed = list(executor.map(_exact_feature_worker, tasks, chunksize=1))
    elapsed = time.perf_counter() - started
    membership = {(query_id, radius): values for query_id, radius, values in completed}
    results = []
    for radius in radii:
        matches = np.zeros((len(query_ids), len(all_ids)), dtype=bool)
        for row_index, query_id in enumerate(query_ids):
            exported_matches = membership[(query_id, float(radius))]
            for candidate_id, is_match in zip(exported_ids, exported_matches, strict=True):
                if is_match:
                    matches[row_index, index_by_id[candidate_id]] = True
        matrix_path = out / "matches" / method / sid / f"{radius_tag(float(radius))}.npz"
        matrix_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            matrix_path,
            matches=matches,
            query_ids=np.asarray(query_ids, dtype=object),
            candidate_ids=np.asarray(all_ids, dtype=object),
        )
        result = {
            "method": method,
            "smiles": smiles,
            "smiles_id": sid,
            "radius": float(radius),
            "query_count": len(query_ids),
            "candidate_count": len(all_ids),
            "exported_candidate_count": len(exported_ids),
            "unexported_candidate_count": len(all_ids) - len(exported_ids),
            "minimum_matches": int(matches.sum(axis=1).min()) if len(query_ids) else 0,
            "maximum_matches": int(matches.sum(axis=1).max()) if len(query_ids) else 0,
            "total_query_candidate_hits": int(matches.sum()),
            "group_runtime_seconds_all_radii": elapsed,
            "search_complete": True,
            "matching_algorithm": "exhaustive_typed_feature_correspondence_then_proper_kabsch",
            "pair_distance_pruning_is_necessary_only": True,
            "per_feature_radius_check": True,
            "feature_recognition_backend": "Pharmit",
            "feature_cache_sha256": feature_metadata["cache_sha256"],
            "numerical_tolerance_angstrom": 0.0002,
            "matrix_path": str(matrix_path),
            "matrix_sha256": sha256(matrix_path),
        }
        atomic_json_dump(matrix_path.with_suffix(".metadata.json"), result)
        results.append(result)
    return results


def _match_group(
    *,
    config: dict[str, Any],
    engine: PharmitBackend,
    method: str,
    smiles: str,
    candidate_records: list[dict[str, Any]],
    query_group: pd.DataFrame,
    database: Path,
    radius: float,
) -> dict[str, Any]:
    out = output_dir(config)
    sid = stable_smiles_id(smiles)
    tag = radius_tag(radius)
    query_ids = query_group["query_id"].astype(str).tolist()
    query_paths = [out / "queries" / tag / f"{query_id}.json" for query_id in query_ids]
    search_dir = out / "search_work" / method / sid / tag
    search_dir.mkdir(parents=True, exist_ok=True)
    output_paths = [search_dir / f"{query_id}.sdf" for query_id in query_ids]
    for output_path in output_paths:
        output_path.unlink(missing_ok=True)
    exported_records = [record for record in candidate_records if record["exported"]]
    exported_ids = [record["candidate_id"] for record in exported_records]
    all_ids = [record["candidate_id"] for record in candidate_records]
    index_by_id = {candidate_id: index for index, candidate_id in enumerate(all_ids)}
    matches = np.zeros((len(query_ids), len(all_ids)), dtype=bool)
    batch_size = int(config["backend"]["query_batch_size"])
    logs = []
    started = time.perf_counter()
    for start in range(0, len(query_ids), batch_size):
        stop = min(len(query_ids), start + batch_size)
        logs.append(
            engine.search_many(
                database_path=database,
                query_paths=query_paths[start:stop],
                output_paths=output_paths[start:stop],
                max_hits=len(exported_ids) + 1,
                workers=int(config["backend"]["workers"]),
            )
        )
    match_counts = []
    witness_written = False
    for row_index, (query_id, result_path) in enumerate(
        zip(query_ids, output_paths, strict=True)
    ):
        titles = read_match_titles(result_path)
        unknown = set(titles) - set(exported_ids)
        if unknown:
            raise ValueError(f"Pharmit returned unknown candidate IDs: {sorted(unknown)[:3]}")
        if len(titles) > len(exported_ids):
            raise RuntimeError("Search returned more candidates than were indexed")
        for title in titles:
            matches[row_index, index_by_id[title]] = True
        match_counts.append(len(titles))
        if titles and not witness_written:
            _write_single_witness(
                result_path,
                out / "witnesses" / method / sid / f"{tag}_{query_id}.sdf",
            )
            witness_written = True
        result_path.unlink(missing_ok=True)
    elapsed = time.perf_counter() - started
    matrix_path = out / "matches" / method / sid / f"{tag}.npz"
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        matrix_path,
        matches=matches,
        query_ids=np.asarray(query_ids, dtype=object),
        candidate_ids=np.asarray(all_ids, dtype=object),
    )
    log_path = matrix_path.with_suffix(".log")
    log_path.write_text("\n".join(logs), encoding="utf-8")
    result = {
        "method": method,
        "smiles": smiles,
        "smiles_id": sid,
        "radius": radius,
        "query_count": len(query_ids),
        "candidate_count": len(all_ids),
        "exported_candidate_count": len(exported_ids),
        "unexported_candidate_count": len(all_ids) - len(exported_ids),
        "minimum_matches": min(match_counts, default=0),
        "maximum_matches": max(match_counts, default=0),
        "total_query_candidate_hits": int(matches.sum()),
        "runtime_seconds": elapsed,
        "search_complete": True,
        "max_hits": len(exported_ids) + 1,
        "max_hits_exceeds_database_size": True,
        "max_orientations_per_conformer": 1,
        "triplet_fingerprint_pruning_disabled": engine.skip_triplet_fingerprints,
        "higher_order_triplet_pruning_disabled": engine.skip_geometric_pruning,
        "matrix_path": str(matrix_path),
        "matrix_sha256": sha256(matrix_path),
    }
    atomic_json_dump(matrix_path.with_suffix(".metadata.json"), result)
    return result


def stage_match(
    config: dict[str, Any],
    *,
    methods: Sequence[str] | None = None,
    molecule_limit: int | None = None,
) -> None:
    out = output_dir(config)
    freeze = json.loads((out / "query_freeze_manifest.json").read_text())
    if freeze["status"] != "frozen":
        raise RuntimeError("Full matching requires a frozen, non-debug query manifest")
    if freeze["config_sha256"] != config["_sha256"]:
        raise RuntimeError("Config changed after reference-only query freeze")
    query_manifest = pd.read_csv(out / "query_manifest.csv", keep_default_na=False)
    eligible = query_manifest[query_manifest["eligible"].map(bool_value)].copy()
    if eligible.empty:
        raise RuntimeError("No eligible queries")
    method_names = list(methods or config["methods"])
    unknown = set(method_names) - set(config["methods"])
    if unknown:
        raise ValueError(f"Unknown methods: {sorted(unknown)}")
    smiles_order = sorted(eligible["expected_smiles"].unique())
    if molecule_limit is not None:
        smiles_order = smiles_order[:molecule_limit]
    radii = [
        float(config["query"]["primary_radius"]),
        *map(float, config["query"]["sensitivity_radii"]),
    ]
    engine = backend(config)
    candidate_rows = []
    completeness_rows = []
    for method in method_names:
        pool = load_pickle(Path(config["methods"][method]["path"]))
        for molecule_number, smiles in enumerate(smiles_order, start=1):
            sid = stable_smiles_id(smiles)
            print(
                f"match export/db: {method} molecule {molecule_number}/{len(smiles_order)}",
                flush=True,
            )
            sdf_path = out / "candidate_libraries" / method / f"{sid}.sdf"
            records_path = out / "candidate_libraries" / method / f"{sid}.records.json"
            database = out / "pharmit_databases" / method / sid
            marker = database / ".pilot_complete.json"
            if records_path.is_file() and marker.is_file():
                records = json.loads(records_path.read_text())
                marker_payload = json.loads(marker.read_text())
            else:
                records = export_candidate_library(
                    method=method,
                    smiles=smiles,
                    mols=pool[smiles],
                    sdf_path=sdf_path,
                )
                atomic_json_dump(records_path, records)
                create_log = engine.create_database(sdf_path, database)
                marker_payload = {
                    "method": method,
                    "smiles": smiles,
                    "candidate_count": len(records),
                    "exported_candidate_count": sum(r["exported"] for r in records),
                    "library_sha256": sha256(sdf_path),
                    "backend_sha256": sha256(engine.executable),
                    "create_log": create_log,
                }
                atomic_json_dump(marker, marker_payload)
            records_refreshed = False
            for record in records:
                exact_hash = molecule_coordinate_hash(
                    pool[smiles][int(record["candidate_index"])]
                )
                if record.get("coordinate_sha256") != exact_hash:
                    record["coordinate_sha256"] = exact_hash
                    records_refreshed = True
                semantics = "RDKit topology plus double-precision coordinates"
                if record.get("coordinate_hash_semantics") != semantics:
                    record["coordinate_hash_semantics"] = semantics
                    records_refreshed = True
            if records_refreshed:
                atomic_json_dump(records_path, records)
            if "export_audit" not in marker_payload:
                marker_payload["export_audit"] = _audit_candidate_export(
                    pool[smiles], records, sdf_path
                )
                atomic_json_dump(marker, marker_payload)
            candidate_rows.extend(records)
            group = eligible[eligible["expected_smiles"] == smiles].sort_values(
                "reference_index"
            )
            reusable_results = []
            for radius in radii:
                matrix_path = out / "matches" / method / sid / f"{radius_tag(radius)}.npz"
                metadata_path = matrix_path.with_suffix(".metadata.json")
                if matrix_path.is_file() and metadata_path.is_file():
                    with np.load(matrix_path, allow_pickle=True) as archive:
                        matrix = archive["matches"]
                        stored_queries = list(map(str, archive["query_ids"]))
                        stored_candidates = list(map(str, archive["candidate_ids"]))
                        expected_queries = group["query_id"].astype(str).tolist()
                        expected_candidates = [record["candidate_id"] for record in records]
                        metadata = json.loads(metadata_path.read_text())
                        if (
                            matrix.shape == (len(group), len(records))
                            and stored_queries == expected_queries
                            and stored_candidates == expected_candidates
                            and metadata.get("matching_algorithm")
                            == "exhaustive_typed_feature_correspondence_then_proper_kabsch"
                        ):
                            metadata["reused"] = True
                            reusable_results.append(metadata)
            if len(reusable_results) == len(radii):
                completeness_rows.extend(reusable_results)
            else:
                exact_results = _match_group_exact_features(
                    config=config,
                    engine=engine,
                    method=method,
                    smiles=smiles,
                    mols=pool[smiles],
                    candidate_records=records,
                    query_group=group,
                    radii=radii,
                )
                for result in exact_results:
                    result["reused"] = False
                    completeness_rows.append(result)
            write_frame_atomic(pd.DataFrame(candidate_rows), out / "candidate_manifest.partial.csv")
            write_frame_atomic(
                pd.DataFrame(completeness_rows), out / "search_completeness.partial.csv"
            )
    candidates = pd.DataFrame(candidate_rows).drop_duplicates(
        ["method", "smiles", "candidate_index"]
    )
    completeness = pd.DataFrame(completeness_rows).drop_duplicates(
        ["method", "smiles", "radius"], keep="last"
    )
    write_frame_atomic(candidates, out / "candidate_manifest.csv")
    write_frame_atomic(completeness, out / "search_completeness.csv")
    (out / "candidate_manifest.partial.csv").unlink(missing_ok=True)
    (out / "search_completeness.partial.csv").unlink(missing_ok=True)
    radius_membership_violations = 0
    for method in method_names:
        for smiles in smiles_order:
            loaded = [
                _load_match_matrix(out, method, smiles, radius) for radius in radii
            ]
            query_ids = [item[0] for item in loaded]
            candidate_ids = [item[1] for item in loaded]
            if any(ids != query_ids[0] for ids in query_ids[1:]) or any(
                ids != candidate_ids[0] for ids in candidate_ids[1:]
            ):
                raise ValueError("Radius matrices do not have identical indexed membership")
            for lower, upper in zip(loaded, loaded[1:]):
                radius_membership_violations += int(np.count_nonzero(lower[2] & ~upper[2]))
    if radius_membership_violations:
        raise RuntimeError(
            f"Larger query radii lost {radius_membership_violations} candidate memberships"
        )
    manifest = {
        "status": "complete",
        "config_sha256": config["_sha256"],
        "methods": method_names,
        "molecules": len(smiles_order),
        "eligible_queries": len(eligible),
        "radii": radii,
        "candidate_records": len(candidates),
        "all_searches_complete": bool(completeness["search_complete"].all()),
        "all_candidates_exported": bool(candidates["exported"].map(bool_value).all()),
        "radius_membership_monotonicity_violations": radius_membership_violations,
        "matching_algorithm": "exhaustive_typed_feature_correspondence_then_proper_kabsch",
        "pharmit_indexed_search_used_for_final_membership": False,
        "query_freeze_sha256": sha256(out / "query_freeze_manifest.json"),
        "candidate_manifest_sha256": sha256(out / "candidate_manifest.csv"),
        "search_completeness_sha256": sha256(out / "search_completeness.csv"),
    }
    atomic_json_dump(out / "match_manifest.json", manifest)
    print(json.dumps(manifest, indent=2), flush=True)


def stage_posebusters(
    config: dict[str, Any],
    *,
    methods: Sequence[str] | None = None,
    molecule_indices: Sequence[int] | None = None,
    finalize: bool = True,
) -> None:
    """Replay the common candidate-level PoseBusters gate without replenishment."""
    from rdkit import Chem, RDLogger

    from casf_benchmark.analysis.metrics import mol_no_h
    from casf_benchmark.generation.conformer_sets import posebusters_geometry_passes

    RDLogger.DisableLog("rdApp.*")
    out = output_dir(config)
    pb_dir = out / "posebusters"
    pb_dir.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    failure_rows = []
    settings = config["posebusters"]
    qwen_method = config["analysis"]["mixture_sources"][0]
    original_qwen = (
        load_pickle(Path(config["original_qwen"]["path"]))
        if "original_qwen" in config
        else None
    )
    method_names = list(methods or config["methods"])
    unknown = set(method_names) - set(config["methods"])
    if unknown:
        raise ValueError(f"Unknown methods: {sorted(unknown)}")
    all_expected_smiles = sorted(
        pd.read_csv(config["reference_audit"], keep_default_na=False)["expected_smiles"].unique()
    )
    if molecule_indices is None:
        expected_smiles = all_expected_smiles
    else:
        invalid = sorted(
            index
            for index in set(map(int, molecule_indices))
            if index < 0 or index >= len(all_expected_smiles)
        )
        if invalid:
            raise ValueError(
                f"PoseBusters molecule indices outside 0..{len(all_expected_smiles) - 1}: "
                f"{invalid}"
            )
        expected_smiles = [all_expected_smiles[index] for index in molecule_indices]
    for method in method_names:
        method_settings = config["methods"][method]
        pool = load_pickle(Path(method_settings["path"]))
        for number, smiles in enumerate(expected_smiles, start=1):
            sid = stable_smiles_id(smiles)
            mask_path = pb_dir / method / f"{sid}.npz"
            checks_path = pb_dir / method / f"{sid}.checks.jsonl"
            mask_path.parent.mkdir(parents=True, exist_ok=True)
            if mask_path.is_file() and checks_path.is_file():
                mask = np.load(mask_path)["passes"].astype(bool)
                check_rows = [json.loads(line) for line in checks_path.read_text().splitlines()]
            else:
                mols = [mol_no_h(mol) for mol in pool[smiles]]
                if any(mol.GetNumConformers() == 0 for mol in mols):
                    raise ValueError(f"Conformer missing during PB replay: {method} {smiles}")
                reference = Chem.MolFromSmiles(smiles)
                renumbered = []
                for candidate_index, mol in enumerate(mols):
                    match = (
                        mol.GetSubstructMatch(reference)
                        if mol.GetNumAtoms() == reference.GetNumAtoms()
                        else ()
                    )
                    if len(match) != reference.GetNumAtoms():
                        raise ValueError(
                            "Cannot renumber candidate into reference atom order: "
                            f"{method} {smiles} candidate {candidate_index}"
                        )
                    renumbered.append(Chem.RenumberAtoms(mol, list(match)))
                result = posebusters_geometry_passes(
                    renumbered,
                    reference,
                    max_workers=int(settings["workers"]),
                    energy_num_threads=int(settings["energy_num_threads"]),
                )
                mask = np.asarray(result.passes, dtype=bool)
                check_rows = result.check_rows
                np.savez_compressed(mask_path, passes=mask)
                with checks_path.open("w", encoding="utf-8") as handle:
                    for candidate_index, row in enumerate(check_rows):
                        handle.write(
                            json.dumps(
                                {
                                    "candidate_index": candidate_index,
                                    "candidate_id": f"{method}__{sid}__{candidate_index:04d}",
                                    **row,
                                },
                                sort_keys=True,
                            )
                            + "\n"
                        )
            if len(mask) != len(pool[smiles]):
                raise RuntimeError("PoseBusters mask length differs from supplied pool")
            failures = Counter()
            for passed, checks in zip(mask, check_rows, strict=True):
                if not passed:
                    failures.update(
                        key
                        for key, value in checks.items()
                        if key not in {"candidate_index", "candidate_id"} and not bool(value)
                    )
            original_count = (
                len(original_qwen[smiles])
                if original_qwen is not None and method == qwen_method
                else len(mask)
            )
            summary_rows.append(
                {
                    "method": method,
                    "smiles": smiles,
                    "smiles_id": sid,
                    "candidate_count": len(mask),
                    "pass_count": int(mask.sum()),
                    "fail_count": int((~mask).sum()),
                    "pass_rate": float(mask.mean()),
                    "archived_comparable_candidate_count": original_count,
                    "archived_comparable_pass_count": int(mask[:original_count].sum()),
                    "replacement_candidate_count": len(mask) - original_count,
                    "replacement_pass_count": int(mask[original_count:].sum()),
                    "mask_path": str(mask_path),
                    "mask_sha256": sha256(mask_path),
                }
            )
            failure_rows.extend(
                {
                    "method": method,
                    "smiles": smiles,
                    "smiles_id": sid,
                    "check": check,
                    "failure_count": count,
                }
                for check, count in failures.items()
            )
            print(
                f"PoseBusters {method} {number}/{len(expected_smiles)}: "
                f"{int(mask.sum())}/{len(mask)}",
                flush=True,
            )
            if finalize:
                write_frame_atomic(pd.DataFrame(summary_rows), pb_dir / "summary.partial.csv")
    if not finalize:
        print(
            json.dumps(
                {
                    "status": "shard_complete",
                    "methods": method_names,
                    "molecule_indices": (
                        list(map(int, molecule_indices))
                        if molecule_indices is not None
                        else list(range(len(all_expected_smiles)))
                    ),
                    "molecules": len(expected_smiles),
                },
                indent=2,
            ),
            flush=True,
        )
        return
    summary = pd.DataFrame(summary_rows)
    write_frame_atomic(summary, pb_dir / "summary.csv")
    write_frame_atomic(pd.DataFrame(failure_rows), pb_dir / "failure_counts.csv")
    (pb_dir / "summary.partial.csv").unlink(missing_ok=True)
    archived = {
        str(method): int(total)
        for method, total in settings.get("archived_pass_totals", {}).items()
    }
    replay = summary.groupby("method")["pass_count"].sum().astype(int).to_dict()
    comparable_replay = (
        summary.groupby("method")["archived_comparable_pass_count"]
        .sum()
        .astype(int)
        .to_dict()
    )
    replacement_replay = (
        summary.groupby("method")["replacement_pass_count"].sum().astype(int).to_dict()
    )
    comparison = [
        {
            "method": method,
            "archived_pass_count": archived[method],
            "replayed_comparable_pass_count": comparable_replay[method],
            "comparable_difference": comparable_replay[method] - archived[method],
            "completed_pool_pass_count": replay[method],
            "replacement_pass_count": replacement_replay[method],
        }
        for method in method_names
        if method in archived
    ]
    comparison_columns = [
        "method",
        "archived_pass_count",
        "replayed_comparable_pass_count",
        "comparable_difference",
        "completed_pool_pass_count",
        "replacement_pass_count",
    ]
    write_frame_atomic(
        pd.DataFrame(comparison, columns=comparison_columns),
        pb_dir / "archived_comparison.csv",
    )
    manifest = {
        "status": "complete",
        "config_sha256": config["_sha256"],
        "methods": method_names,
        "candidate_masks_saved": True,
        "failed_candidates_consume_budget": True,
        "post_gate_replenishment": False,
        "candidate_atom_order": (
            "heavy-atom candidates renumbered to the coordinate-free reference SMILES "
            "order, matching the archived evaluator"
        ),
        "settings": settings,
        "replayed_totals": replay,
        "archived_comparable_replayed_totals": comparable_replay,
        "replacement_pass_totals": replacement_replay,
        "archived_totals": archived,
        "archived_results_overwritten": False,
    }
    atomic_json_dump(pb_dir / "manifest.json", manifest)
    print(json.dumps(manifest, indent=2), flush=True)


def _load_match_matrix(
    out: Path,
    method: str,
    smiles: str,
    radius: float,
) -> tuple[list[str], list[str], np.ndarray]:
    path = out / "matches" / method / stable_smiles_id(smiles) / f"{radius_tag(radius)}.npz"
    with np.load(path, allow_pickle=True) as archive:
        query_ids = list(map(str, archive["query_ids"]))
        candidate_ids = list(map(str, archive["candidate_ids"]))
        matches = archive["matches"].astype(bool)
    if matches.shape != (len(query_ids), len(candidate_ids)):
        raise ValueError(f"Malformed match matrix: {path}")
    return query_ids, candidate_ids, matches


def _bootstrap_mean_interval(values: np.ndarray, replicates: int, seed: int) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(values), size=(replicates, len(values)))
    samples = values[indices].mean(axis=1)
    return float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))


def _build_probability_tables(config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = output_dir(config)
    manifest = pd.read_csv(out / "query_manifest.csv", keep_default_na=False)
    eligible = manifest[manifest["eligible"].map(bool_value)].copy()
    radii = [
        float(config["query"]["primary_radius"]),
        *map(float, config["query"]["sensitivity_radii"]),
    ]
    budgets = list(map(int, config["analysis"]["budgets"]))
    methods = list(config["methods"])
    qwen_method, loqi_method = config["analysis"]["mixture_sources"]
    original_qwen = (
        load_pickle(Path(config["original_qwen"]["path"]))
        if "original_qwen" in config
        else None
    )
    pb_available = (out / "posebusters" / "manifest.json").is_file()
    query_rows = []
    curve_rows = []
    for smiles, group in eligible.groupby("expected_smiles", sort=True):
        expected_query_ids = group.sort_values("reference_index")["query_id"].astype(str).tolist()
        matrices: dict[tuple[str, float], np.ndarray] = {}
        candidate_ids_by_method: dict[str, list[str]] = {}
        pb_masks: dict[str, np.ndarray] = {}
        for method in methods:
            for radius in radii:
                query_ids, candidate_ids, matches = _load_match_matrix(
                    out, method, smiles, radius
                )
                if query_ids != expected_query_ids:
                    raise ValueError(f"Query ID/order mismatch for {method} {smiles}")
                if not candidate_ids:
                    raise ValueError(f"Primary candidate pool is empty for {method} {smiles}")
                if method in candidate_ids_by_method and (
                    candidate_ids_by_method[method] != candidate_ids
                ):
                    raise ValueError(
                        f"Candidate ID/order differs across radii for {method} {smiles}"
                    )
                matrices[(method, radius)] = matches
                candidate_ids_by_method[method] = candidate_ids
            if pb_available:
                mask_path = (
                    out / "posebusters" / method / f"{stable_smiles_id(smiles)}.npz"
                )
                pb_masks[method] = np.load(mask_path)["passes"].astype(bool)
                if len(pb_masks[method]) != len(candidate_ids_by_method[method]):
                    raise ValueError("PB mask does not match candidate denominator")
        original_count = len(original_qwen[smiles]) if original_qwen is not None else 0
        for radius in radii:
            for query_index, query_id in enumerate(expected_query_ids):
                method_hits = {}
                method_pb_hits = {}
                for method in methods:
                    vector = matrices[(method, radius)][query_index]
                    population = len(candidate_ids_by_method[method])
                    hits = int(vector.sum())
                    method_hits[method] = hits
                    pb_hits = int((vector & pb_masks[method]).sum()) if pb_available else None
                    method_pb_hits[method] = pb_hits
                    query_rows.append(
                        {
                            "query_id": query_id,
                            "smiles": smiles,
                            "smiles_id": stable_smiles_id(smiles),
                            "method": method,
                            "radius": radius,
                            "candidate_count": population,
                            "hit_count": hits,
                            "hit_fraction": hits / population,
                            "pb_pass_hit_count": pb_hits,
                            "original_qwen_candidate_count": (
                                original_count if method == qwen_method else ""
                            ),
                            "original_qwen_hit_count": (
                                int(vector[:original_count].sum())
                                if method == qwen_method
                                else ""
                            ),
                        }
                    )
                    for budget in budgets:
                        realized_draws = min(budget, population)
                        curve_rows.append(
                            {
                                "query_id": query_id,
                                "smiles": smiles,
                                "smiles_id": stable_smiles_id(smiles),
                                "method": method,
                                "analysis": "primary",
                                "radius": radius,
                                "budget": budget,
                                "realized_draws": realized_draws,
                                "hit_probability": hit_probability(
                                    population, hits, realized_draws
                                ),
                            }
                        )
                        if pb_available:
                            curve_rows.append(
                                {
                                    "query_id": query_id,
                                    "smiles": smiles,
                                    "smiles_id": stable_smiles_id(smiles),
                                    "method": method,
                                    "analysis": "posebusters_gated_no_replenishment",
                                    "radius": radius,
                                    "budget": budget,
                                    "realized_draws": realized_draws,
                                    "hit_probability": hit_probability(
                                        population, int(pb_hits), realized_draws
                                    ),
                                }
                            )
                        if method == qwen_method and original_qwen is not None:
                            curve_rows.append(
                                {
                                    "query_id": query_id,
                                    "smiles": smiles,
                                    "smiles_id": stable_smiles_id(smiles),
                                    "method": "qwen_original_pool",
                                    "analysis": "original_qwen_sensitivity",
                                    "radius": radius,
                                    "budget": budget,
                                    "realized_draws": min(budget, original_count),
                                    "hit_probability": hit_probability(
                                        original_count,
                                        int(vector[:original_count].sum()),
                                        min(budget, original_count),
                                    ),
                                }
                            )
                for budget in budgets:
                    qwen_population = len(candidate_ids_by_method[qwen_method])
                    loqi_population = len(candidate_ids_by_method[loqi_method])
                    equal_draws = min(
                        budget // 2,
                        qwen_population,
                        loqi_population,
                    )
                    realized_mixture_draws = 2 * equal_draws
                    curve_rows.append(
                        {
                            "query_id": query_id,
                            "smiles": smiles,
                            "smiles_id": stable_smiles_id(smiles),
                            "method": "qwen_loqi_equal_mixture",
                            "analysis": "primary",
                            "radius": radius,
                            "budget": budget,
                            "realized_draws": realized_mixture_draws,
                            "hit_probability": mixture_hit_probability(
                                qwen_population,
                                method_hits[qwen_method],
                                loqi_population,
                                method_hits[loqi_method],
                                realized_mixture_draws,
                            ),
                        }
                    )
                    if pb_available:
                        curve_rows.append(
                            {
                                "query_id": query_id,
                                "smiles": smiles,
                                "smiles_id": stable_smiles_id(smiles),
                                "method": "qwen_loqi_equal_mixture",
                                "analysis": "posebusters_gated_no_replenishment",
                                "radius": radius,
                                "budget": budget,
                                "realized_draws": realized_mixture_draws,
                                "hit_probability": mixture_hit_probability(
                                    qwen_population,
                                    int(method_pb_hits[qwen_method]),
                                    loqi_population,
                                    int(method_pb_hits[loqi_method]),
                                    realized_mixture_draws,
                                ),
                            }
                        )
    return pd.DataFrame(query_rows), pd.DataFrame(curve_rows)


def _aggregate_curves(config: dict[str, Any], curves: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    molecule = (
        curves.groupby(
            ["analysis", "radius", "budget", "method", "smiles", "smiles_id"],
            as_index=False,
        )["hit_probability"]
        .mean()
        .rename(columns={"hit_probability": "molecule_coverage"})
    )
    replicates = int(config["analysis"]["bootstrap_replicates"])
    seed = int(config["analysis"]["bootstrap_seed"])
    rows = []
    for keys, group in molecule.groupby(
        ["analysis", "radius", "budget", "method"], sort=True
    ):
        analysis, radius, budget, method = keys
        values = group["molecule_coverage"].to_numpy(float)
        low, high = _bootstrap_mean_interval(
            values,
            replicates,
            seed + int(round(radius * 1000)) + int(budget),
        )
        rows.append(
            {
                "analysis": analysis,
                "radius": radius,
                "budget": budget,
                "method": method,
                "coverage": float(values.mean()),
                "ci_low": low,
                "ci_high": high,
                "molecule_count": len(values),
                "query_count": int(
                    curves[
                        (curves["analysis"] == analysis)
                        & (curves["radius"] == radius)
                        & (curves["budget"] == budget)
                        & (curves["method"] == method)
                    ]["query_id"].nunique()
                ),
            }
        )
    return molecule, pd.DataFrame(rows)


def _endpoint_comparisons(
    config: dict[str, Any],
    curves: pd.DataFrame,
    molecule: pd.DataFrame,
) -> pd.DataFrame:
    radius = float(config["query"]["primary_radius"])
    budget = max(map(int, config["analysis"]["budgets"]))
    subset = molecule[
        (molecule["analysis"] == "primary")
        & (molecule["radius"] == radius)
        & (molecule["budget"] == budget)
    ]
    pivot = subset.pivot(index="smiles", columns="method", values="molecule_coverage")
    qwen = config["analysis"]["mixture_sources"][0]
    baselines = list(config["analysis"].get("baseline_methods", []))
    choices = list(config["methods"]) + ["qwen_loqi_equal_mixture"]
    pairs = [("qwen_loqi_equal_mixture", qwen)]
    pairs.extend(
        (choice, baseline)
        for choice in choices
        for baseline in baselines
        if choice != baseline
    )
    rows = []
    replicates = int(config["analysis"]["bootstrap_replicates"])
    seed = int(config["analysis"]["bootstrap_seed"])
    query_subset = curves[
        (curves["analysis"] == "primary")
        & (curves["radius"] == radius)
        & (curves["budget"] == budget)
    ]
    query_pivot = query_subset.pivot(
        index="query_id", columns="method", values="hit_probability"
    )
    for pair_index, (left, right) in enumerate(dict.fromkeys(pairs)):
        if left not in pivot or right not in pivot:
            continue
        stats = paired_bootstrap_difference(
            pivot[left].to_numpy(float),
            pivot[right].to_numpy(float),
            replicates=replicates,
            seed=seed + pair_index,
        )
        query_difference = query_pivot[left] - query_pivot[right]
        rows.append(
            {
                "left_method": left,
                "right_method": right,
                "radius": radius,
                "budget": budget,
                **stats,
                "expected_queries_gained_net": float(query_difference.sum()),
                "expected_queries_gained_positive": float(
                    query_difference.clip(lower=0).sum()
                ),
                "expected_queries_lost": float(
                    (-query_difference.clip(upper=0)).sum()
                ),
            }
        )
    return pd.DataFrame(rows)


def _plot_curves(config: dict[str, Any], aggregate: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = output_dir(config)
    radius = float(config["query"]["primary_radius"])
    primary = aggregate[
        (aggregate["analysis"] == "primary") & (aggregate["radius"] == radius)
    ]
    labels = {method: settings["label"] for method, settings in config["methods"].items()}
    labels["qwen_loqi_equal_mixture"] = "Qwen/LoQI 50:50"
    colors = {
        "qwen_1p7b_v2_coordtoken_step41054": "#1f77b4",
        "loqi_druglike": "#ff7f0e",
        "rdkit_random_raw_druglike": "#2ca02c",
        "rdkit_random_minimized_22_druglike": "#d62728",
        "qwen_loqi_equal_mixture": "#9467bd",
    }
    figure, axis = plt.subplots(figsize=(8.2, 5.2))
    for method, group in primary.groupby("method"):
        group = group.sort_values("budget")
        axis.plot(
            group["budget"], group["coverage"], marker="o", linewidth=2,
            label=labels[method], color=colors.get(method),
        )
        axis.fill_between(
            group["budget"], group["ci_low"], group["ci_high"],
            alpha=0.12, color=colors.get(method),
        )
    axis.set_xscale("log")
    axis.set_xticks(config["analysis"]["budgets"])
    axis.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axis.set_ylim(0, 1.02)
    axis.set_xlabel("Accepted cached conformers sampled without replacement")
    axis.set_ylabel("Expected query coverage (equal molecule weight)")
    axis.set_title(f"PLINDER-22 positional pharmacophore coverage · radius {radius:g} Å")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8)
    figure.tight_layout()
    figures = out / "figures"
    figures.mkdir(exist_ok=True)
    figure.savefig(figures / "coverage_vs_budget_primary.png", dpi=180)
    figure.savefig(figures / "coverage_vs_budget_primary.pdf")
    plt.close(figure)

    gated = aggregate[
        (aggregate["analysis"] == "posebusters_gated_no_replenishment")
        & (aggregate["radius"] == radius)
    ]
    if not gated.empty:
        figure, axis = plt.subplots(figsize=(8.2, 5.2))
        for method, group in gated.groupby("method"):
            group = group.sort_values("budget")
            axis.plot(
                group["budget"], group["coverage"], marker="o", linewidth=2,
                label=labels[method], color=colors.get(method),
            )
        axis.set_xscale("log")
        axis.set_xticks(config["analysis"]["budgets"])
        axis.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        axis.set_ylim(0, 1.02)
        axis.set_xlabel("Supplied conformers sampled; failed candidates still consume budget")
        axis.set_ylabel("Expected PB-gated query coverage")
        axis.set_title(f"PoseBusters-gated sensitivity · radius {radius:g} Å")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
        figure.tight_layout()
        figure.savefig(figures / "coverage_vs_budget_posebusters.png", dpi=180)
        plt.close(figure)

    original_qwen = aggregate[
        (aggregate["analysis"] == "original_qwen_sensitivity")
        & (aggregate["radius"] == radius)
        & (aggregate["method"] == "qwen_original_pool")
    ]
    completed_qwen = primary[
        primary["method"] == config["analysis"]["mixture_sources"][0]
    ]
    if not original_qwen.empty:
        figure, axis = plt.subplots(figsize=(8.2, 5.2))
        for label, group, color in (
            ("Completed Qwen (1,000 each)", completed_qwen, "#1f77b4"),
            ("Original Qwen (958–1,000)", original_qwen, "#7f8c8d"),
        ):
            group = group.sort_values("budget")
            axis.plot(
                group["budget"], group["coverage"], marker="o", linewidth=2,
                label=label, color=color,
            )
        axis.set_xscale("log")
        axis.set_xticks(config["analysis"]["budgets"])
        axis.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        axis.set_ylim(0, 1.02)
        axis.set_xlabel("Cached conformers requested")
        axis.set_ylabel("Expected query coverage")
        axis.set_title("Qwen top-up sensitivity")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
        figure.tight_layout()
        figure.savefig(figures / "qwen_original_pool_sensitivity.png", dpi=180)
        plt.close(figure)

    endpoint = aggregate[
        (aggregate["analysis"] == "primary")
        & (aggregate["budget"] == max(config["analysis"]["budgets"]))
    ]
    figure, axis = plt.subplots(figsize=(8.2, 5.2))
    for method, group in endpoint.groupby("method"):
        group = group.sort_values("radius")
        axis.plot(
            group["radius"], group["coverage"], marker="o", linewidth=2,
            label=labels[method], color=colors.get(method),
        )
    axis.set_ylim(0, 1.02)
    axis.set_xlabel("Feature-sphere radius (Å)")
    axis.set_ylabel("Expected query coverage at budget 1,000")
    axis.set_title("Prespecified positional-tolerance sensitivity")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(figures / "tolerance_sensitivity.png", dpi=180)
    plt.close(figure)


def _publication_styles(config: dict[str, Any]) -> tuple[dict[str, str], dict[str, str]]:
    labels = {method: settings["label"] for method, settings in config["methods"].items()}
    labels["qwen_loqi_equal_mixture"] = "Qwen/LoQI 50:50"
    palette = [
        "#0072B2",
        "#E69F00",
        "#009E73",
        "#CC79A7",
        "#56B4E9",
        "#D55E00",
        "#F0E442",
        "#666666",
        "#8E44AD",
    ]
    colors = {
        method: palette[index % len(palette)]
        for index, method in enumerate(labels)
    }
    colors["qwen_loqi_equal_mixture"] = "#8E44AD"
    return labels, colors


def _plot_publication_coverage(
    config: dict[str, Any],
    aggregate: pd.DataFrame,
    breakdowns: pd.DataFrame,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = output_dir(config)
    figures = out / "figures"
    figures.mkdir(exist_ok=True)
    labels, colors = _publication_styles(config)
    radius = float(config["query"]["primary_radius"])
    maximum_budget = max(map(int, config["analysis"]["budgets"]))
    primary = aggregate[
        (aggregate["analysis"] == "primary") & (aggregate["radius"] == radius)
    ]
    order = list(config["methods"]) + ["qwen_loqi_equal_mixture"]

    figure, axes = plt.subplots(1, 2, figsize=(12.6, 5.1))
    for method in order:
        group = primary[primary["method"] == method].sort_values("budget")
        if group.empty:
            continue
        axes[0].plot(
            group["budget"],
            group["coverage"],
            marker="o",
            markersize=3.5,
            linewidth=1.7,
            label=labels[method],
            color=colors[method],
        )
    axes[0].set_xscale("log")
    axes[0].set_xticks(config["analysis"]["budgets"])
    axes[0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axes[0].set_ylim(0, 1.02)
    axes[0].set_xlabel("Requested conformer budget")
    axes[0].set_ylabel("Expected query coverage")
    axes[0].set_title("A  Coverage versus budget", loc="left", fontweight="bold")
    axes[0].grid(alpha=0.2)

    endpoint = primary[primary["budget"] == maximum_budget].set_index("method").loc[order]
    y = np.arange(len(order))
    values = endpoint["coverage"].to_numpy(float)
    axes[1].errorbar(
        values,
        y,
        xerr=np.vstack(
            [
                values - endpoint["ci_low"].to_numpy(float),
                endpoint["ci_high"].to_numpy(float) - values,
            ]
        ),
        fmt="none",
        ecolor="#7f8c8d",
        capsize=2,
        zorder=1,
    )
    for position, method in enumerate(order):
        axes[1].scatter(values[position], position, color=colors[method], s=35, zorder=2)
    axes[1].set_yticks(y, [labels[method] for method in order])
    axes[1].invert_yaxis()
    axes[1].set_xlim(0, 1.02)
    axes[1].set_xlabel("Coverage (95% molecule bootstrap CI)")
    axes[1].set_title(
        f"B  Endpoint at requested budget {maximum_budget:,}",
        loc="left",
        fontweight="bold",
    )
    axes[1].grid(axis="x", alpha=0.2)
    handles, legend_labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, legend_labels, loc="lower center", ncol=5, fontsize=8)
    figure.suptitle(
        f"{config.get('paper_label', config['pilot_id'])}: positional pharmacophore "
        f"coverage ({radius:g} Å)",
        fontsize=12,
    )
    figure.tight_layout(rect=(0, 0.12, 1, 0.95))
    figure.savefig(figures / "paper_pharmacophore_coverage.png", dpi=300)
    figure.savefig(figures / "paper_pharmacophore_coverage.pdf")
    plt.close(figure)

    focus = list(
        config.get("publication", {}).get(
            "difficulty_methods",
            [
                config["analysis"]["mixture_sources"][0],
                config["analysis"]["mixture_sources"][1],
                "qwen_loqi_equal_mixture",
            ],
        )
    )
    dimensions = [
        ("feature_count_group", ["4", "5", "6", "7", "8+"], "Contact features"),
        ("flexibility_group", ["0–2", "3–5", "6+"], "Rotatable bonds"),
    ]
    figure, axes = plt.subplots(1, 2, figsize=(11.4, 4.7), sharey=True)
    for axis, (dimension, groups, x_label) in zip(axes, dimensions, strict=True):
        subset = breakdowns[
            (breakdowns["dimension"] == dimension) & breakdowns["method"].isin(focus)
        ]
        for method in focus:
            method_rows = subset[subset["method"] == method].set_index("group").reindex(groups)
            x = np.arange(len(groups))
            values = method_rows["coverage"].to_numpy(float)
            axis.errorbar(
                x,
                values,
                yerr=np.vstack(
                    [
                        values - method_rows["ci_low"].to_numpy(float),
                        method_rows["ci_high"].to_numpy(float) - values,
                    ]
                ),
                marker="o",
                linewidth=1.8,
                capsize=2,
                label=labels[method],
                color=colors[method],
            )
        axis.set_xticks(np.arange(len(groups)), groups)
        axis.set_xlabel(x_label)
        axis.set_ylim(0, 1.02)
        axis.grid(alpha=0.2)
    axes[0].set_ylabel("Endpoint coverage (95% molecule bootstrap CI)")
    axes[0].set_title("A  Query complexity", loc="left", fontweight="bold")
    axes[1].set_title("B  Molecular flexibility", loc="left", fontweight="bold")
    axes[1].legend(fontsize=8, loc="lower left")
    figure.suptitle(
        f"{config.get('paper_label', config['pilot_id'])}: prespecified difficulty strata",
        fontsize=12,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    figure.savefig(figures / "paper_query_difficulty.png", dpi=300)
    figure.savefig(figures / "paper_query_difficulty.pdf")
    plt.close(figure)


def _plot_visual_examples(config: dict[str, Any], query_counts: pd.DataFrame) -> None:
    """Render deterministic Qwen-only, LoQI-only, and shared match overlays."""
    from rdkit import Chem

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = output_dir(config)
    radius = float(config["query"]["primary_radius"])
    qwen, loqi = config["analysis"]["mixture_sources"]
    subset = query_counts[query_counts["radius"] == radius]
    pivot = subset.pivot(index="query_id", columns="method", values="hit_count")
    if qwen not in pivot or loqi not in pivot:
        return
    categories = [
        ("Qwen only", qwen, pivot[(pivot[qwen] > 0) & (pivot[loqi] == 0)].index),
        ("LoQI only", loqi, pivot[(pivot[loqi] > 0) & (pivot[qwen] == 0)].index),
        ("Shared", qwen, pivot[(pivot[qwen] > 0) & (pivot[loqi] > 0)].index),
    ]
    manifest = pd.read_csv(out / "query_manifest.csv", keep_default_na=False).set_index(
        "query_id"
    )
    selected = []
    for category, method, query_ids in categories:
        available = manifest.loc[manifest.index.intersection(query_ids)].copy()
        if available.empty:
            continue
        available["feature_count_numeric"] = pd.to_numeric(
            available["feature_count"], errors="coerce"
        )
        available["distance_from_six"] = (
            available["feature_count_numeric"] - 6
        ).abs()
        row = available.sort_values(
            ["distance_from_six", "query_id"], kind="stable"
        ).iloc[0]
        selected.append((category, method, row.name, row))
    if not selected:
        return

    labels, colors = _publication_styles(config)
    feature_colors = {
        "Aromatic": "#8e44ad",
        "HydrogenAcceptor": "#e74c3c",
        "HydrogenDonor": "#3498db",
        "Hydrophobic": "#f39c12",
        "NegativeIon": "#c0392b",
        "PositiveIon": "#2980b9",
    }
    pools: dict[str, dict[str, list[Any]]] = {}
    figure = plt.figure(figsize=(5.2 * len(selected), 5.0))
    selection_rows = []

    def draw_bonds(axis: Any, molecule: Any, coordinates: np.ndarray, color: str, alpha: float) -> None:
        for bond in molecule.GetBonds():
            first = coordinates[bond.GetBeginAtomIdx()]
            second = coordinates[bond.GetEndAtomIdx()]
            axis.plot(
                [first[0], second[0]],
                [first[1], second[1]],
                [first[2], second[2]],
                color=color,
                alpha=alpha,
                linewidth=1.4,
            )

    for panel, (category, method, query_id, row) in enumerate(selected, start=1):
        smiles = str(row["expected_smiles"])
        query = json.loads(
            (out / "queries" / radius_tag(radius) / f"{query_id}.json").read_text()
        )
        query_coordinates = np.asarray(
            [[point[axis] for axis in "xyz"] for point in query["points"]], dtype=float
        )
        query_types = [str(point["name"]) for point in query["points"]]
        query_ids, candidate_ids, matches = _load_match_matrix(
            out, method, smiles, radius
        )
        match_row = query_ids.index(query_id)
        candidate_index = int(np.flatnonzero(matches[match_row])[0])
        sid = stable_smiles_id(smiles)
        with np.load(
            out / "candidate_features" / method / f"{sid}.npz",
            allow_pickle=True,
        ) as archive:
            feature_coordinates = np.asarray(
                archive["coordinates"][candidate_index], dtype=float
            )
            stored_types = archive["feature_types"]
            if stored_types.ndim == 1 and all(
                isinstance(value, str) for value in stored_types.tolist()
            ):
                candidate_types = list(map(str, stored_types))
            else:
                candidate_types = list(map(str, stored_types[candidate_index]))
        transform = exact_feature_correspondence_transform(
            query_coordinates,
            query_types,
            feature_coordinates,
            candidate_types,
            radius,
        )
        if transform is None:
            raise RuntimeError("Saved match lacks a reproducible correspondence transform")
        if method not in pools:
            pools[method] = load_pickle(Path(config["methods"][method]["path"]))
        candidate = Chem.RemoveHs(Chem.Mol(pools[method][smiles][candidate_index]))
        native = Chem.RemoveHs(_load_single_sdf(str(row["raw_sdf_path"])))
        candidate_xyz = np.asarray(
            candidate.GetConformer().GetPositions(), dtype=float
        )
        candidate_xyz = (
            candidate_xyz @ transform["rotation"] + transform["translation"]
        )
        native_xyz = np.asarray(native.GetConformer().GetPositions(), dtype=float)

        axis = figure.add_subplot(1, len(selected), panel, projection="3d")
        draw_bonds(axis, native, native_xyz, "#95a5a6", 0.65)
        draw_bonds(axis, candidate, candidate_xyz, colors[method], 0.9)
        axis.scatter(*native_xyz.T, s=10, c="#bdc3c7", alpha=0.6)
        axis.scatter(*candidate_xyz.T, s=12, c=colors[method], alpha=0.85)
        for point, coordinates in zip(query["points"], query_coordinates, strict=True):
            axis.scatter(
                *coordinates,
                s=85,
                c=feature_colors.get(point["name"], "#16a085"),
                edgecolors="black",
                linewidths=0.35,
            )
        combined = np.vstack([native_xyz, candidate_xyz, query_coordinates])
        center = combined.mean(axis=0)
        span = max(float(np.ptp(combined, axis=0).max()), 1.0) / 2
        axis.set_xlim(center[0] - span, center[0] + span)
        axis.set_ylim(center[1] - span, center[1] + span)
        axis.set_zlim(center[2] - span, center[2] + span)
        axis.set_title(
            f"{chr(64 + panel)}  {category}\n{row['name']} · {labels[method]}",
            loc="left",
            fontsize=10,
            fontweight="bold",
        )
        axis.set_xlabel("x (Å)")
        axis.set_ylabel("y (Å)")
        axis.set_zlabel("z (Å)")
        selection_rows.append(
            {
                "category": category,
                "method": method,
                "query_id": query_id,
                "name": row["name"],
                "pdb_id": row["pdb_id"],
                "smiles": smiles,
                "candidate_id": candidate_ids[candidate_index],
                "candidate_index": candidate_index,
                "feature_count": len(query_types),
                "maximum_feature_residual": float(np.max(transform["residuals"])),
            }
        )
    figure.suptitle(
        "Representative rigid pharmacophore matches: native (gray), generated "
        "conformer (colored), contact features (spheres)",
        fontsize=11,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    figures = out / "figures"
    figures.mkdir(exist_ok=True)
    figure.savefig(figures / "paper_match_examples.png", dpi=300)
    figure.savefig(figures / "paper_match_examples.pdf")
    plt.close(figure)
    write_frame_atomic(pd.DataFrame(selection_rows), out / "visual_example_selection.csv")


def _write_findings_report(
    config: dict[str, Any],
    aggregate: pd.DataFrame,
    comparisons: pd.DataFrame,
) -> None:
    out = output_dir(config)
    radius = float(config["query"]["primary_radius"])
    budget = max(config["analysis"]["budgets"])
    endpoint = aggregate[
        (aggregate["analysis"] == "primary")
        & (aggregate["radius"] == radius)
        & (aggregate["budget"] == budget)
    ].sort_values("coverage", ascending=False)
    manifest = pd.read_csv(out / "query_manifest.csv", keep_default_na=False)
    exclusions = manifest.loc[~manifest["eligible"].map(bool_value), "exclusion_reason"].value_counts()
    eligible = manifest[manifest["eligible"].map(bool_value)]
    labels = {method: settings["label"] for method, settings in config["methods"].items()}
    labels["qwen_loqi_equal_mixture"] = "Qwen/LoQI 50:50"
    labels["qwen_original_pool"] = "Original Qwen pool"
    topup_path = out / "qwen_topup" / "topup_manifest.json"
    topup = json.loads(topup_path.read_text()) if topup_path.is_file() else None
    lines = [
        f"# {config.get('paper_label', config['pilot_id'])} pharmacophore-query findings",
        "",
        "## Scope",
        "",
        "This pilot measures same-molecule structural coverage: whether cached conformer "
        "ensembles generated from molecular identity contain rigid arrangements of direct "
        "protein-contact pharmacophore features observed in experimental complexes. It is "
        "not a virtual-screening enrichment experiment.",
        "",
        "## Reference eligibility",
        "",
        f"- Retained references: {len(manifest):,}",
        f"- Eligible positional queries: {int(manifest['eligible'].map(bool_value).sum()):,}",
        f"- Eligible molecules: {eligible['expected_smiles'].nunique()}",
        "- Query eligibility was frozen from reference-only checks before generator matching.",
        "",
        "| Reference-only exclusion | Count |",
        "|---|---:|",
    ]
    lines.extend(f"| {reason} | {count:,} |" for reason, count in exclusions.items())
    lines.extend(
        [
            "",
            "The retained context includes nearby non-protein species for provenance, but the "
            "extractor received only directly contacting protein chains. Queries with nearby "
            f"metals: {(manifest.get('nearby_metals', '') != '').sum():,}; nearby nucleic-acid "
            f"chains: {(manifest.get('nearby_nucleic_acids', '') != '').sum():,}; nearby "
            f"cofactors: {(manifest.get('nearby_cofactors', '') != '').sum():,}.",
            "",
            "## Matching completeness",
            "",
            "Pharmit supplies the pinned feature definitions and recognizes every candidate "
            "feature set. Final membership does not use its indexed triplet search: the pilot "
            "exhaustively enumerates every injective type-compatible correspondence, applies "
            "only necessary pair-distance pruning, then requires a proper rigid Kabsch fit with "
            "every feature inside the stated radius. Candidate-level membership had zero "
            "set-inclusion violations as radii increased from 0.5 to 0.75 to 1.0 Å.",
            "",
        f"## Primary endpoint ({radius:g} Å, budget {budget:,})",
        "",
        "| Library | Coverage | 95% molecule-bootstrap interval |",
        "|---|---:|---:|",
        ]
    )
    if topup is not None:
        lines[lines.index("## Matching completeness"):lines.index("## Matching completeness")] = [
            "## Qwen completion",
            "",
            f"The immutable original pool contained {sum(topup['original_counts'].values()):,} "
            f"accepted conformers. The unbiased top-up made {topup['attempts']} attempts, "
            f"accepted {topup['accepted_additions']} additions, and rejected "
            f"{sum(topup['rejections_by_reason'].values())}; every molecule then had exactly "
            "1,000 accepted cached conformers. Original artifact hashes and original-prefix "
            "coordinate hashes were unchanged.",
            "",
        ]
    for row in endpoint.itertuples(index=False):
        lines.append(
            f"| {labels.get(row.method, row.method)} | {row.coverage:.4f} | "
            f"[{row.ci_low:.4f}, {row.ci_high:.4f}] |"
        )
    mix_qwen = comparisons[
        (comparisons["left_method"] == "qwen_loqi_equal_mixture")
        & (comparisons["right_method"] == config["analysis"]["mixture_sources"][0])
    ].iloc[0]
    lines.extend(
        [
            "",
            "## Principal paired comparison",
            "",
            f"The equal Qwen/LoQI mixture differs from Qwen alone by "
            f"{mix_qwen['difference']:+.4f} coverage points (95% interval "
            f"[{mix_qwen['ci_low']:+.4f}, {mix_qwen['ci_high']:+.4f}]).",
            "",
            "| Comparison | Difference | 95% paired-molecule interval | Expected queries gained | Expected queries lost |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in comparisons.itertuples(index=False):
        lines.append(
            f"| {labels.get(row.left_method, row.left_method)} − "
            f"{labels.get(row.right_method, row.right_method)} | {row.difference:+.4f} | "
            f"[{row.ci_low:+.4f}, {row.ci_high:+.4f}] | "
            f"{row.expected_queries_gained_positive:.2f} | {row.expected_queries_lost:.2f} |"
        )
    lines.extend(
        [
            "",
            "Each method uses its supplied pool without replacement; when a requested budget "
            "exceeds availability, all available conformers are used and `realized_draws` records "
            "the cap. The 1,000-conformer mixture allocates up to 500 Qwen and 500 LoQI draws. Expected "
            "gains/losses use the natural shared-Qwen coupling for mixture versus Qwen and "
            "deterministic full-pool endpoints for the single libraries.",
            "",
            "## Prespecified sensitivities",
            "",
        ]
    )
    tolerance = aggregate[
        (aggregate["analysis"] == "primary") & (aggregate["budget"] == budget)
    ].sort_values(["radius", "method"])
    lines.extend(["| Radius (Å) | Library | Coverage |", "|---:|---|---:|"])
    for row in tolerance.itertuples(index=False):
        lines.append(
            f"| {row.radius:g} | {labels.get(row.method, row.method)} | {row.coverage:.4f} |"
        )
    original_endpoint = aggregate[
        (aggregate["analysis"] == "original_qwen_sensitivity")
        & (aggregate["radius"] == radius)
        & (aggregate["budget"] == budget)
    ]
    if not original_endpoint.empty:
        completed = endpoint[endpoint["method"] == config["analysis"]["mixture_sources"][0]].iloc[0]
        original = original_endpoint.iloc[0]
        lines.extend(
            [
                "",
                f"Original Qwen-pool coverage was {original.coverage:.4f}; completed-Qwen "
                f"coverage was {completed.coverage:.4f} (marginal difference "
                f"{completed.coverage - original.coverage:+.4f}).",
            ]
        )
    gated = aggregate[
        (aggregate["analysis"] == "posebusters_gated_no_replenishment")
        & (aggregate["radius"] == radius)
        & (aggregate["budget"] == budget)
    ].sort_values("coverage", ascending=False)
    if not gated.empty:
        lines.extend(
            [
                "",
                "### PoseBusters-gated sensitivity",
                "",
                "Failed candidates remain in the sampling denominator and cannot match; no "
                "replacement or post-gate replenishment is performed.",
                "",
                "| Library | Gated coverage |",
                "|---|---:|",
            ]
        )
        for row in gated.itertuples(index=False):
            lines.append(f"| {labels.get(row.method, row.method)} | {row.coverage:.4f} |")
        comparison_path = out / "posebusters" / "archived_comparison.csv"
        if comparison_path.is_file():
            pb_comparison = pd.read_csv(comparison_path)
            if not pb_comparison.empty:
                lines.extend(
                    [
                        "",
                        "| Library | Archived passes | Comparable replay | Difference | Replacement passes |",
                        "|---|---:|---:|---:|---:|",
                    ]
                )
                for row in pb_comparison.itertuples(index=False):
                    lines.append(
                        f"| {labels.get(row.method, row.method)} | {row.archived_pass_count:,} | "
                        f"{row.replayed_comparable_pass_count:,} | {row.comparable_difference:+,} | "
                        f"{row.replacement_pass_count:,} |"
                    )

    rmsd_summary_path = out / "rmsd_agreement" / "summary.csv"
    if rmsd_summary_path.is_file():
        rmsd_summary = pd.read_csv(rmsd_summary_path)
        lines.extend(
            [
                "",
                "## Whole-molecule RMSD disagreement",
                "",
                f"Candidate correspondence was verified explicitly. Whole-molecule hits use "
                f"RMSD ≤ {config['analysis']['rmsd_hit_threshold']:g} Å; pharmacophore hits use "
                f"the {radius:g} Å positional query.",
                "",
                "| Library | Both | Pharmacophore only | RMSD only | Neither |",
                "|---|---:|---:|---:|---:|",
            ]
        )
        for row in rmsd_summary.itertuples(index=False):
            lines.append(
                f"| {labels.get(row.method, row.method)} | {row.both_hit_count:,} | "
                f"{row.pharmacophore_only_count:,} | {row.rmsd_only_count:,} | "
                f"{row.neither_count:,} |"
            )

    budget_note = (
        "All budgets count accepted cached conformers. Failed PoseBusters candidates, "
        "where that sensitivity is shown, still consume the draw and are not replenished."
    )
    if topup is not None:
        budget_note += " The extra Qwen generation effort is reported separately in `qwen_topup/`."
    else:
        budget_note += " Generation attempts and compute are not equated across methods."
    lines.extend(
        [
            "",
            "## Budget semantics and interpretation limits",
            "",
            budget_note,
            "",
            f"The {len(manifest):,} complexes are clustered within "
            f"{eligible['expected_smiles'].nunique()} selected molecules and are not "
            "independent samples. Uncertainty therefore resamples molecules. Feature definitions "
            "use the pinned local Pharmit backend and rigid matching uses the exhaustive positional "
            "correspondence check described above; no conformer was embedded, optimized, "
            "energy-ranked, or selected using the experimental references. Results "
            "measure experimental feature-arrangement coverage on this selected panel, not "
            "virtual-screening enrichment or independent model generalization.",
        ]
    )
    (out / "FINDINGS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _coverage_breakdowns(config: dict[str, Any], curves: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    from rdkit import Chem
    from rdkit.Chem import Lipinski

    out = output_dir(config)
    radius = float(config["query"]["primary_radius"])
    budget = max(map(int, config["analysis"]["budgets"]))
    endpoint = curves[
        (curves["analysis"] == "primary")
        & (curves["radius"] == radius)
        & (curves["budget"] == budget)
    ].copy()
    manifest = pd.read_csv(out / "query_manifest.csv", keep_default_na=False)
    eligible = manifest[manifest["eligible"].map(bool_value)].copy()
    descriptor_rows = []
    for smiles in sorted(eligible["expected_smiles"].unique()):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Cannot calculate flexibility descriptor: {smiles}")
        descriptor_rows.append(
            {"smiles": smiles, "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol))}
        )
    descriptors = pd.DataFrame(descriptor_rows)
    characteristics = eligible.merge(
        descriptors,
        left_on="expected_smiles",
        right_on="smiles",
        validate="many_to_one",
    )
    characteristics["feature_count"] = pd.to_numeric(
        characteristics["feature_count"], errors="raise"
    ).astype(int)
    characteristics["feature_count_group"] = pd.cut(
        characteristics["feature_count"],
        bins=[3, 4, 5, 6, 7, np.inf],
        labels=["4", "5", "6", "7", "8+"],
    ).astype(str)
    characteristics["flexibility_group"] = pd.cut(
        characteristics["rotatable_bonds"],
        bins=[-1, 2, 5, np.inf],
        labels=["0–2", "3–5", "6+"],
    ).astype(str)
    for column, label in (
        ("nearby_metals", "nearby metal"),
        ("nearby_nucleic_acids", "nearby nucleic acid"),
        ("nearby_cofactors", "nearby cofactor"),
    ):
        characteristics[label] = np.where(characteristics[column] != "", "present", "absent")
    characteristics["legacy_geometry_flag"] = np.where(
        characteristics["reference_geometry_pass"].map(bool_value), "pass", "flagged"
    )
    write_frame_atomic(characteristics, out / "eligible_query_characteristics.csv")

    joined = endpoint.merge(
        characteristics[
            [
                "query_id",
                "feature_count_group",
                "flexibility_group",
                "nearby metal",
                "nearby nucleic acid",
                "nearby cofactor",
                "legacy_geometry_flag",
            ]
        ],
        on="query_id",
        validate="many_to_one",
    )
    breakdown_rows = []
    replicates = int(config["analysis"]["bootstrap_replicates"])
    seed = int(config["analysis"]["bootstrap_seed"])
    dimensions = [
        "feature_count_group",
        "flexibility_group",
        "nearby metal",
        "nearby nucleic acid",
        "nearby cofactor",
        "legacy_geometry_flag",
    ]
    for dimension in dimensions:
        per_molecule = (
            joined.groupby(["method", dimension, "smiles"], as_index=False, observed=True)[
                "hit_probability"
            ]
            .mean()
            .rename(columns={"hit_probability": "molecule_coverage"})
        )
        query_counts = (
            joined.groupby(["method", dimension], observed=True)["query_id"]
            .nunique()
            .to_dict()
        )
        for (method, group_label), group in per_molecule.groupby(
            ["method", dimension], observed=True
        ):
            values = group["molecule_coverage"].to_numpy(float)
            group_seed = seed + sum(map(ord, f"{dimension}|{group_label}|{method}"))
            low, high = _bootstrap_mean_interval(values, replicates, group_seed)
            breakdown_rows.append(
                {
                    "dimension": dimension,
                    "group": group_label,
                    "method": method,
                    "radius": radius,
                    "budget": budget,
                    "coverage": float(values.mean()),
                    "ci_low": low,
                    "ci_high": high,
                    "molecule_count": int(group["smiles"].nunique()),
                    "query_count": int(query_counts[(method, group_label)]),
                    "small_group_fewer_than_five_molecules": bool(
                        group["smiles"].nunique() < 5
                    ),
                }
            )
    return characteristics, pd.DataFrame(breakdown_rows)


def stage_summarize(config: dict[str, Any]) -> None:
    out = output_dir(config)
    query_counts, curves = _build_probability_tables(config)
    if not curves["hit_probability"].between(0.0, 1.0).all():
        raise ValueError("A sampling probability is outside [0, 1]")
    budget_violations = 0
    for _, group in curves.groupby(
        ["query_id", "method", "analysis", "radius"], sort=False
    ):
        values = group.sort_values("budget")["hit_probability"].to_numpy(float)
        budget_violations += int(np.any(np.diff(values) < -1e-12))
    radius_violations = 0
    for _, group in curves.groupby(
        ["query_id", "method", "analysis", "budget"], sort=False
    ):
        values = group.sort_values("radius")["hit_probability"].to_numpy(float)
        radius_violations += int(np.any(np.diff(values) < -1e-12))
    if budget_violations or radius_violations:
        raise RuntimeError(
            "Expected coverage is not monotonic: "
            f"budget groups={budget_violations}, radius groups={radius_violations}"
        )
    completeness = pd.read_csv(out / "search_completeness.csv", keep_default_na=False)
    if not completeness["search_complete"].map(bool_value).all():
        raise RuntimeError("Summary cannot treat unresolved searches as negative matches")
    molecule, aggregate = _aggregate_curves(config, curves)
    comparisons = _endpoint_comparisons(config, curves, molecule)
    _, breakdowns = _coverage_breakdowns(config, curves)
    write_frame_atomic(query_counts, out / "per_query_match_counts.csv")
    write_frame_atomic(curves, out / "per_query_budget_curves.csv")
    write_frame_atomic(molecule, out / "per_molecule_budget_curves.csv")
    write_frame_atomic(aggregate, out / "aggregate_budget_curves.csv")
    write_frame_atomic(comparisons, out / "endpoint_comparisons.csv")
    write_frame_atomic(breakdowns, out / "coverage_breakdowns.csv")
    _plot_curves(config, aggregate)
    _plot_publication_coverage(config, aggregate, breakdowns)
    _plot_visual_examples(config, query_counts)
    _write_findings_report(config, aggregate, comparisons)
    manifest = {
        "status": "complete",
        "config_sha256": config["_sha256"],
        "aggregation": "mean queries within molecule, then equal mean across molecules",
        "sampling": "exact hypergeometric uniform subsets without replacement",
        "mixture": "independent equal draws from completed Qwen and LoQI pools",
        "bootstrap_unit": "molecule",
        "bootstrap_replicates": int(config["analysis"]["bootstrap_replicates"]),
        "bootstrap_seed": int(config["analysis"]["bootstrap_seed"]),
        "posebusters_sensitivity_included": bool(
            (out / "posebusters" / "manifest.json").is_file()
        ),
        "validation": {
            "probabilities_in_unit_interval": True,
            "budget_monotonicity_violations": budget_violations,
            "radius_monotonicity_violations": radius_violations,
            "all_searches_resolved": True,
            "availability_semantics": (
                "requested draws are capped at each supplied method-molecule pool; "
                "realized_draws is reported per query and budget"
            ),
        },
        "outputs": {
            name: sha256(out / name)
            for name in [
                "per_query_match_counts.csv",
                "per_query_budget_curves.csv",
                "per_molecule_budget_curves.csv",
                "aggregate_budget_curves.csv",
                "endpoint_comparisons.csv",
                "eligible_query_characteristics.csv",
                "coverage_breakdowns.csv",
                "FINDINGS.md",
            ]
        },
    }
    atomic_json_dump(out / "summary_manifest.json", manifest)
    print(json.dumps(manifest, indent=2), flush=True)


def _extend_qwen_rmsd_cache(
    *,
    original_matrices: dict[str, np.ndarray],
    original_pool: dict[str, list[Any]],
    completed_pool: dict[str, list[Any]],
    ground_truth: dict[str, dict[str, Any]],
    workers: int,
) -> dict[str, np.ndarray]:
    from concurrent.futures import ThreadPoolExecutor, as_completed

    from casf_benchmark.analysis.druglike_covmat import best_rmsd, remove_hydrogens

    extended = {}
    tasks = []
    for smiles, original in original_matrices.items():
        original_count = len(original_pool[smiles])
        if original.shape[1] != original_count:
            raise ValueError(f"Original Qwen RMSD columns do not match cache order: {smiles}")
        additions = remove_hydrogens(completed_pool[smiles][original_count:])
        references = remove_hydrogens(ground_truth[smiles]["confs"])
        appended = np.full((len(references), len(additions)), np.nan, dtype=np.float32)
        extended[smiles] = (original, appended)
        for reference_index, reference in enumerate(references):
            tasks.append((smiles, reference_index, reference, additions))

    def calculate(item: tuple[str, int, Any, list[Any]]) -> tuple[str, int, np.ndarray]:
        smiles, reference_index, reference, additions = item
        values = np.asarray(
            [best_rmsd(mol, reference) for mol in additions], dtype=np.float32
        )
        return smiles, reference_index, values

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(calculate, task) for task in tasks]
        for completed, future in enumerate(as_completed(futures), start=1):
            smiles, reference_index, values = future.result()
            extended[smiles][1][reference_index] = values
            if completed % 100 == 0 or completed == len(futures):
                print(f"replacement RMSD rows: {completed}/{len(futures)}", flush=True)
    return {
        smiles: np.concatenate([original, appended], axis=1)
        for smiles, (original, appended) in extended.items()
    }


def _molecule_weighted_agreement(
    agreement: pd.DataFrame,
    categories: Sequence[str],
    *,
    replicates: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate query-candidate agreement with equal weight per molecule."""
    fractions = agreement[["method", "smiles", "query_id"]].copy()
    denominator = agreement["candidate_count"].to_numpy(float)
    if np.any(denominator <= 0):
        raise ValueError("Agreement rows contain an empty candidate pool")
    for category in categories:
        fractions[category] = agreement[category].to_numpy(float) / denominator
    per_molecule = fractions.groupby(["method", "smiles"], as_index=False)[
        list(categories)
    ].mean()
    rows = []
    for method, group in per_molecule.groupby("method", sort=True):
        for category_index, category in enumerate(categories):
            values = group[category].to_numpy(float)
            low, high = _bootstrap_mean_interval(
                values,
                replicates,
                seed + category_index + sum(map(ord, str(method))),
            )
            rows.append(
                {
                    "method": method,
                    "category": category,
                    "fraction": float(values.mean()),
                    "ci_low": low,
                    "ci_high": high,
                    "molecule_count": len(values),
                }
            )
    return per_molecule, pd.DataFrame(rows)


def stage_rmsd(config: dict[str, Any], *, workers: int = 12) -> None:
    """Join pharmacophore membership to verified, candidate-aligned RMSD matrices."""
    from rdkit import Chem

    from casf_benchmark.analysis.druglike_covmat import (
        best_rmsd,
        compute_rmsd_matrices,
        load_rmsd_matrices,
        save_rmsd_matrices,
    )
    from casf_benchmark.analysis.plinder_panel import load_reference_sdf_dir

    out = output_dir(config)
    rmsd_dir = out / "rmsd_agreement"
    rmsd_dir.mkdir(parents=True, exist_ok=True)
    audit = pd.read_csv(config["reference_audit"], keep_default_na=False)
    cleaned_dirs = {str(Path(path).parent) for path in audit["source_sdf"].unique()}
    if len(cleaned_dirs) != 1:
        raise ValueError("Reference audit points to multiple cleaned SDF directories")
    ground_truth = load_reference_sdf_dir(Path(cleaned_dirs.pop()))
    methods = list(config["methods"])
    matrices: dict[str, dict[str, np.ndarray]] = {}
    cache_manifest = []
    correspondence_rows = []
    configured_caches = config["analysis"].get("rmsd_caches", {})
    qwen_method = config["analysis"]["mixture_sources"][0]

    def correspondence_checks(
        method: str,
        method_matrices: dict[str, np.ndarray],
        method_pool: dict[str, list[Any]],
    ) -> tuple[list[dict[str, Any]], list[str]]:
        rows: list[dict[str, Any]] = []
        failures: list[str] = []
        for smiles, reference_row in ground_truth.items():
            if smiles not in method_matrices:
                failures.append(f"missing matrix for {smiles}")
                continue
            matrix = method_matrices[smiles]
            expected_shape = (
                len(reference_row["confs"]),
                len(method_pool[smiles]),
            )
            if matrix.shape != expected_shape:
                failures.append(
                    f"shape mismatch for {smiles}: {matrix.shape} != {expected_shape}"
                )
                continue
            references = reference_row["confs"]
            candidates = method_pool[smiles]
            reference_anchors = sorted({0, len(references) // 2, len(references) - 1})
            candidate_anchors = sorted({0, len(candidates) // 2, len(candidates) - 1})
            for reference_index in reference_anchors:
                for candidate_index in candidate_anchors:
                    recalculated = best_rmsd(
                        Chem.RemoveHs(Chem.Mol(candidates[candidate_index])),
                        Chem.RemoveHs(Chem.Mol(references[reference_index])),
                    )
                    cached = float(matrix[reference_index, candidate_index])
                    delta = abs(recalculated - cached)
                    passed = (
                        (math.isnan(recalculated) and math.isnan(cached))
                        or (
                            math.isfinite(recalculated)
                            and math.isfinite(cached)
                            and delta <= 1e-3
                        )
                    )
                    rows.append(
                        {
                            "method": method,
                            "smiles": smiles,
                            "reference_index": reference_index,
                            "candidate_index": candidate_index,
                            "cached_rmsd": cached,
                            "recalculated_rmsd": recalculated,
                            "absolute_difference": delta,
                            "passed": passed,
                        }
                    )
                    if not passed:
                        failures.append(
                            f"{smiles} reference {reference_index} candidate "
                            f"{candidate_index}: cached={cached}, "
                            f"recalculated={recalculated}"
                        )
        return rows, failures

    for method in methods:
        method_pool = load_pickle(Path(config["methods"][method]["path"]))
        destination = rmsd_dir / f"{method}.npz"
        if destination.is_file():
            method_matrices = load_rmsd_matrices(destination)
            source_description = "pilot cache"
        elif (
            method == qwen_method
            and "original_qwen" in config
            and method in configured_caches
        ):
            original_cache = Path(configured_caches[method])
            original_matrices = load_rmsd_matrices(original_cache)
            original_pool = load_pickle(Path(config["original_qwen"]["path"]))
            method_matrices = _extend_qwen_rmsd_cache(
                original_matrices=original_matrices,
                original_pool=original_pool,
                completed_pool=method_pool,
                ground_truth=ground_truth,
                workers=workers,
            )
            save_rmsd_matrices(destination, method_matrices)
            source_description = f"{original_cache} plus replacement-only calculation"
        elif method in configured_caches and Path(configured_caches[method]).is_file():
            source = Path(configured_caches[method])
            method_matrices = load_rmsd_matrices(source)
            for smiles, matrix in method_matrices.items():
                if matrix.shape[1] != len(method_pool[smiles]):
                    raise ValueError(f"RMSD candidate count mismatch: {method} {smiles}")
            save_rmsd_matrices(destination, method_matrices)
            source_description = str(source)
        else:
            method_matrices, missing, all_nan = compute_rmsd_matrices(
                ground_truth, method_pool, workers
            )
            if missing or all_nan:
                raise RuntimeError(
                    f"RMSD calculation incomplete for {method}: missing={missing}, all_nan={all_nan}"
                )
            save_rmsd_matrices(destination, method_matrices)
            source_description = "calculated consistently for pilot"
        method_checks, failures = correspondence_checks(
            method, method_matrices, method_pool
        )
        rejected_source = ""
        if failures:
            if source_description == "calculated consistently for pilot":
                raise ValueError(
                    f"Fresh RMSD calculation failed correspondence checks for {method}: "
                    f"{failures[:3]}"
                )
            rejected_source = source_description
            print(
                f"RMSD cache rejected for {method}; recalculating after "
                f"{len(failures)} failed anchor/shape checks",
                flush=True,
            )
            method_matrices, missing, all_nan = compute_rmsd_matrices(
                ground_truth, method_pool, workers
            )
            if missing or all_nan:
                raise RuntimeError(
                    f"RMSD recalculation incomplete for {method}: "
                    f"missing={missing}, all_nan={all_nan}"
                )
            save_rmsd_matrices(destination, method_matrices)
            source_description = (
                "calculated consistently for pilot after rejecting cache: "
                f"{rejected_source}"
            )
            method_checks, failures = correspondence_checks(
                method, method_matrices, method_pool
            )
            if failures:
                raise ValueError(
                    f"Recalculated RMSD matrix failed correspondence checks for {method}: "
                    f"{failures[:3]}"
                )
        correspondence_rows.extend(method_checks)
        matrices[method] = method_matrices
        cache_manifest.append(
            {
                "method": method,
                "pilot_cache": str(destination),
                "pilot_cache_sha256": sha256(destination),
                "source": source_description,
                "candidate_pickle": config["methods"][method]["path"],
                "candidate_pickle_sha256": sha256(Path(config["methods"][method]["path"])),
                "candidate_order_verified_by_recalculated_anchors": True,
                "rejected_cache_source": rejected_source,
            }
        )

    write_frame_atomic(
        pd.DataFrame(correspondence_rows),
        rmsd_dir / "candidate_correspondence_checks.csv",
    )

    query_manifest = pd.read_csv(out / "query_manifest.csv", keep_default_na=False)
    eligible = query_manifest[query_manifest["eligible"].map(bool_value)].copy()
    radius = float(config["query"]["primary_radius"])
    threshold = float(config["analysis"]["rmsd_hit_threshold"])
    rows = []
    examples = []
    for smiles, group in eligible.groupby("expected_smiles", sort=True):
        ordered = group.sort_values("reference_index")
        for method in methods:
            query_ids, candidate_ids, pharma = _load_match_matrix(
                out, method, smiles, radius
            )
            if query_ids != ordered["query_id"].astype(str).tolist():
                raise ValueError("Pharmacophore/RMSD query correspondence mismatch")
            rmsd_matrix = matrices[method][smiles]
            for match_row, row in enumerate(ordered.itertuples(index=False)):
                rmsd = rmsd_matrix[int(row.reference_index)]
                rmsd_hit = np.isfinite(rmsd) & (rmsd <= threshold)
                pharma_hit = pharma[match_row]
                both = pharma_hit & rmsd_hit
                pharma_only = pharma_hit & ~rmsd_hit
                rmsd_only = ~pharma_hit & rmsd_hit
                neither = ~pharma_hit & ~rmsd_hit
                rows.append(
                    {
                        "query_id": row.query_id,
                        "smiles": smiles,
                        "smiles_id": stable_smiles_id(smiles),
                        "method": method,
                        "radius": radius,
                        "rmsd_threshold": threshold,
                        "candidate_count": len(candidate_ids),
                        "both_hit_count": int(both.sum()),
                        "pharmacophore_only_count": int(pharma_only.sum()),
                        "rmsd_only_count": int(rmsd_only.sum()),
                        "neither_count": int(neither.sum()),
                        "pharmacophore_hit_count": int(pharma_hit.sum()),
                        "rmsd_hit_count": int(rmsd_hit.sum()),
                    }
                )
                for category, mask in (
                    ("pharmacophore_without_rmsd", pharma_only),
                    ("rmsd_without_pharmacophore", rmsd_only),
                ):
                    for index in np.flatnonzero(mask)[:3]:
                        examples.append(
                            {
                                "query_id": row.query_id,
                                "method": method,
                                "candidate_id": candidate_ids[int(index)],
                                "candidate_index": int(index),
                                "category": category,
                                "pharmacophore_match": bool(pharma_hit[index]),
                                "rmsd": float(rmsd[index]),
                            }
                        )
    agreement = pd.DataFrame(rows)
    summary = (
        agreement.groupby("method", as_index=False)[
            [
                "both_hit_count",
                "pharmacophore_only_count",
                "rmsd_only_count",
                "neither_count",
            ]
        ]
        .sum()
    )
    categories = [
        "both_hit_count",
        "pharmacophore_only_count",
        "rmsd_only_count",
        "neither_count",
    ]
    molecule_agreement, molecule_summary = _molecule_weighted_agreement(
        agreement,
        categories,
        replicates=int(config["analysis"]["bootstrap_replicates"]),
        seed=int(config["analysis"]["bootstrap_seed"]),
    )
    write_frame_atomic(agreement, rmsd_dir / "per_query_agreement.csv")
    write_frame_atomic(summary, rmsd_dir / "summary.csv")
    write_frame_atomic(molecule_agreement, rmsd_dir / "per_molecule_agreement.csv")
    write_frame_atomic(molecule_summary, rmsd_dir / "molecule_weighted_summary.csv")
    write_frame_atomic(pd.DataFrame(examples), rmsd_dir / "disagreement_examples.csv")
    write_frame_atomic(pd.DataFrame(cache_manifest), rmsd_dir / "cache_manifest.csv")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    category_labels = ["both", "pharmacophore only", "RMSD only", "neither"]
    colors = ["#27ae60", "#8e44ad", "#e67e22", "#d5d8dc"]
    labels = {
        method: settings["label"] for method, settings in config["methods"].items()
    }
    totals = summary[categories].sum(axis=1).to_numpy(float)
    left = np.zeros(len(summary), dtype=float)
    figure, axis = plt.subplots(figsize=(8.4, 4.8))
    for category, label, color in zip(categories, category_labels, colors, strict=True):
        values = summary[category].to_numpy(float) / totals
        axis.barh(
            [labels[method] for method in summary["method"]],
            values,
            left=left,
            label=label,
            color=color,
        )
        left += values
    axis.set_xlim(0, 1)
    axis.set_xlabel("Fraction of eligible query–candidate pairs")
    axis.set_title(
        f"Pharmacophore ({radius:g} Å) vs whole-molecule RMSD (≤{threshold:g} Å)"
    )
    axis.legend(fontsize=8, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.28))
    figure.tight_layout()
    figure.savefig(rmsd_dir / "pharmacophore_rmsd_agreement.png", dpi=180)
    plt.close(figure)

    method_order = list(config["methods"])
    category_label_by_id = dict(zip(categories, category_labels, strict=True))
    figure, axes = plt.subplots(1, len(categories), figsize=(13.2, 4.8), sharey=True)
    y = np.arange(len(method_order))
    for axis, category in zip(axes, categories, strict=True):
        subset = molecule_summary[molecule_summary["category"] == category].set_index(
            "method"
        ).loc[method_order]
        values = subset["fraction"].to_numpy(float)
        axis.errorbar(
            values,
            y,
            xerr=np.vstack(
                [
                    values - subset["ci_low"].to_numpy(float),
                    subset["ci_high"].to_numpy(float) - values,
                ]
            ),
            fmt="o",
            color="#2c3e50",
            ecolor="#95a5a6",
            capsize=2,
        )
        axis.set_xlim(0, 1)
        axis.set_title(category_label_by_id[category], fontsize=9)
        axis.grid(axis="x", alpha=0.2)
    axes[0].set_yticks(y, [labels[method] for method in method_order])
    axes[0].set_ylabel("Method")
    figure.supxlabel(
        "Mean fraction within query, then molecule (95% molecule bootstrap CI)"
    )
    figure.suptitle(
        f"Molecule-weighted pharmacophore ({radius:g} Å) / RMSD "
        f"(≤{threshold:g} Å) agreement",
        fontsize=11,
    )
    figure.tight_layout()
    figure.savefig(
        rmsd_dir / "pharmacophore_rmsd_agreement_molecule_weighted.png",
        dpi=180,
    )
    figure.savefig(rmsd_dir / "pharmacophore_rmsd_agreement_molecule_weighted.pdf")
    plt.close(figure)
    manifest = {
        "status": "complete",
        "config_sha256": config["_sha256"],
        "workers": workers,
        "pharmacophore_radius": radius,
        "whole_molecule_rmsd_threshold": threshold,
        "candidate_correspondence": (
            "Matrix columns use indexed generation-pickle records and were checked by "
            "recalculating deterministic reference/candidate anchors; Qwen replacement "
            "columns were appended in combined-cache order."
        ),
        "query_correspondence": (
            "reference audit reference_index into the cleaned SDF and RMSD matrix row; "
            "conf_index is retained only as the pre-dedup source key"
        ),
        "candidate_anchor_checks": len(correspondence_rows),
        "candidate_anchor_checks_sha256": sha256(
            rmsd_dir / "candidate_correspondence_checks.csv"
        ),
        "agreement_aggregation": (
            "candidate category fractions averaged over queries within molecule, then "
            "equally over molecules; uncertainty bootstraps molecules"
        ),
        "molecule_weighted_summary_sha256": sha256(
            rmsd_dir / "molecule_weighted_summary.csv"
        ),
        "cache_manifest_sha256": sha256(rmsd_dir / "cache_manifest.csv"),
    }
    atomic_json_dump(rmsd_dir / "manifest.json", manifest)
    print(json.dumps(manifest, indent=2), flush=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    subparsers = parser.add_subparsers(dest="stage", required=True)
    subparsers.add_parser("inventory")
    subparsers.add_parser("smoke")
    queries = subparsers.add_parser("queries")
    queries.add_argument("--limit", type=int)
    match = subparsers.add_parser("match")
    match.add_argument("--method", action="append", dest="methods")
    match.add_argument("--molecule-limit", type=int)
    posebusters = subparsers.add_parser("posebusters")
    posebusters.add_argument("--method", action="append", dest="methods")
    posebusters.add_argument(
        "--molecule-index",
        action="append",
        type=int,
        dest="molecule_indices",
        help="Zero-based index in the sorted reference-audit molecule list; repeatable",
    )
    posebusters.add_argument(
        "--no-finalize",
        action="store_true",
        help="Write candidate masks/checks only; aggregate in a later unsharded replay",
    )
    subparsers.add_parser("summarize")
    rmsd = subparsers.add_parser("rmsd")
    rmsd.add_argument("--workers", type=int, default=12)
    all_stages = subparsers.add_parser("all")
    all_stages.add_argument("--skip-posebusters", action="store_true")
    all_stages.add_argument("--rmsd-workers", type=int, default=12)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    config = load_config(args.config)
    if args.stage == "inventory":
        stage_inventory(config)
    elif args.stage == "smoke":
        stage_smoke(config)
    elif args.stage == "queries":
        stage_queries(config, limit=args.limit)
    elif args.stage == "match":
        stage_match(config, methods=args.methods, molecule_limit=args.molecule_limit)
    elif args.stage == "posebusters":
        stage_posebusters(
            config,
            methods=args.methods,
            molecule_indices=args.molecule_indices,
            finalize=not args.no_finalize,
        )
    elif args.stage == "summarize":
        stage_summarize(config)
    elif args.stage == "rmsd":
        stage_rmsd(config, workers=args.workers)
    elif args.stage == "all":
        stage_inventory(config)
        stage_smoke(config)
        stage_queries(config)
        stage_match(config)
        if not args.skip_posebusters:
            stage_posebusters(config)
        stage_rmsd(config, workers=args.rmsd_workers)
        stage_summarize(config)
    else:  # pragma: no cover - argparse enforces the stage
        raise ValueError(args.stage)


if __name__ == "__main__":
    main()





