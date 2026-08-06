#!/usr/bin/env python3
"""Normalize external CASF generation pools into fixed/dynamic/chembl-count sets."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from casf_benchmark.generation.normalizer import (
    MANIFEST_PARTS_DIRNAME,
    load_method_list_from_csv,
    materialize_generation_root,
    merge_manifest_parts,
    qwen_wrong_artifact_paths,
    quarantine_paths,
    restore_manifest_backup,
    validate_sampled_generation_root,
    write_manifest,
    write_manifest_part,
)
from casf_benchmark.paths import (
    DEFAULT_CHEMBL_DATASET_ROOT,
    DEFAULT_CHEMBL_MAP_CSV,
    DEFAULT_CORE_LIGAND_DIR,
    DEFAULT_QWEN_GENERATION_ROOT,
)
from casf_benchmark.generation.conformer_sets import positive_int


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_QWEN_GENERATION_ROOT)
    parser.add_argument("--generation-dir", type=Path, default=None)
    parser.add_argument("--source-manifest", type=Path, default=None)
    parser.add_argument("--output-manifest", type=Path, default=None)
    parser.add_argument("--chembl-map-csv", type=Path, default=DEFAULT_CHEMBL_MAP_CSV)
    parser.add_argument("--ligand-dir", type=Path, default=DEFAULT_CORE_LIGAND_DIR)
    parser.add_argument("--chembl-dataset-root", type=Path, default=DEFAULT_CHEMBL_DATASET_ROOT)
    parser.add_argument("--method", action="append", dest="methods")
    parser.add_argument("--methods-csv", type=Path, default=None)
    parser.add_argument("--limit-molecules", type=int, default=None)
    parser.add_argument("--molecule-offset", type=int, default=0)
    parser.add_argument("--fixed-set-size", type=positive_int, default=1000)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--posebusters-workers", type=positive_int, default=1)
    parser.add_argument("--posebusters-energy-threads", type=positive_int, default=1)
    parser.add_argument("--manifest-parts-dir", type=Path, default=None)
    parser.add_argument("--write-manifest-part", action="store_true")
    parser.add_argument("--merge-manifest-parts", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--restore-manifest-backup", default="")
    parser.add_argument("--quarantine-wrong-qwen-artifacts", action="store_true")
    parser.add_argument("--quarantine-root", type=Path, default=None)
    return parser


def _methods_from_args(args: argparse.Namespace) -> list[str] | None:
    methods = list(args.methods or [])
    if args.methods_csv is not None:
        methods.extend(load_method_list_from_csv(args.methods_csv))
    return methods or None


def main() -> None:
    args = build_parser().parse_args()
    generation_dir = (args.generation_dir or args.root / "generation").resolve()
    source_manifest = (args.source_manifest or generation_dir / "manifest.tsv").resolve()
    output_manifest = (args.output_manifest or generation_dir / "manifest.tsv").resolve()
    parts_dir = (args.manifest_parts_dir or generation_dir / MANIFEST_PARTS_DIRNAME).resolve()

    if args.quarantine_wrong_qwen_artifacts:
        quarantine_root = args.quarantine_root
        if quarantine_root is None:
            quarantine_root = generation_dir / f"quarantine_wrong_qwen_{time.strftime('%Y%m%d_%H%M%S')}"
        moved = quarantine_paths(qwen_wrong_artifact_paths(generation_dir), quarantine_root.resolve())
        print(f"Quarantined {len(moved)} path(s) under {quarantine_root}", flush=True)

    if args.restore_manifest_backup:
        source_manifest = restore_manifest_backup(generation_dir, args.restore_manifest_backup).resolve()
        output_manifest = source_manifest
        print(f"Restored manifest from {args.restore_manifest_backup}: {source_manifest}", flush=True)

    if args.merge_manifest_parts:
        row_count = merge_manifest_parts(parts_dir, output_manifest)
        print(f"Merged {row_count} manifest row(s): {output_manifest}", flush=True)
        return

    if args.validate_only:
        issues = validate_sampled_generation_root(generation_dir, source_manifest)
        errors = [issue for issue in issues if issue.severity == "error"]
        for issue in issues:
            print(f"{issue.severity}: {issue.message}", flush=True)
        if errors:
            raise SystemExit(f"Validation failed with {len(errors)} error(s)")
        print(f"Validation passed: {source_manifest}", flush=True)
        return

    rows = materialize_generation_root(
        generation_dir=generation_dir,
        manifest_path=source_manifest,
        chembl_map_csv=args.chembl_map_csv.resolve(),
        ligand_dir=args.ligand_dir.resolve(),
        chembl_dataset_root=args.chembl_dataset_root.resolve(),
        methods=_methods_from_args(args),
        limit=args.limit_molecules,
        offset=args.molecule_offset,
        fixed_set_size=args.fixed_set_size,
        seed=args.seed,
        posebusters_workers=args.posebusters_workers,
        posebusters_energy_threads=args.posebusters_energy_threads,
    )
    if args.write_manifest_part:
        part_path = write_manifest_part(parts_dir, rows)
        print(f"Wrote {len(rows)} manifest row(s): {part_path}", flush=True)
    else:
        write_manifest(output_manifest, list(rows))
        print(f"Wrote {len(rows)} manifest row(s): {output_manifest}", flush=True)


if __name__ == "__main__":
    main()
