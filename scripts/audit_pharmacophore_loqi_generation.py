#!/usr/bin/env python3
"""Validate sharded LoQI generation and write a compact provenance manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


SHARD_RE = re.compile(r"(?P<molecule>\d{3})\.shard-(?P<shard>\d{3})-of-(?P<count>\d{3})\.pickle")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def coordinate_hash(mol: object) -> str:
    conformer = mol.GetConformer()
    coordinates = np.asarray(conformer.GetPositions(), dtype="<f8")
    atomic_numbers = np.asarray(
        [atom.GetAtomicNum() for atom in mol.GetAtoms()], dtype="<i2"
    )
    digest = hashlib.sha256()
    digest.update(atomic_numbers.tobytes())
    digest.update(coordinates.tobytes())
    return digest.hexdigest()


def load_smiles(path: Path) -> list[str]:
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    if not isinstance(payload, dict) or not payload:
        raise TypeError(f"Expected a non-empty SMILES-keyed mapping in {path}")
    return list(payload)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generation-input", type=Path, required=True)
    parser.add_argument("--generation-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--base-seed", type=int, required=True)
    parser.add_argument("--expected-conformers", type=int, default=1000)
    parser.add_argument("--slurm-jobs", default="")
    parser.add_argument("--recovery-note", action="append", default=[])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    smiles = load_smiles(args.generation_input)
    parts_dir = args.generation_dir / "parts"
    paths = sorted(parts_dir.glob("*.shard-*-of-*.pickle"))
    ignored_unsharded_parts = sorted(parts_dir.glob("[0-9][0-9][0-9].pickle"))
    parsed: dict[int, list[tuple[int, int, Path]]] = {}
    for path in paths:
        match = SHARD_RE.fullmatch(path.name)
        if match is None:
            raise ValueError(f"Malformed shard filename: {path}")
        molecule_index = int(match.group("molecule"))
        parsed.setdefault(molecule_index, []).append(
            (int(match.group("shard")), int(match.group("count")), path)
        )

    molecules = []
    all_coordinate_hashes: set[str] = set()
    coordinate_hashes_by_molecule: dict[int, list[str]] = {}
    layout: dict[int, int] = {}
    for molecule_index, expected_smiles in enumerate(smiles):
        entries = sorted(parsed.get(molecule_index, []))
        if not entries:
            raise FileNotFoundError(f"No shards found for molecule {molecule_index}")
        counts = {count for _, count, _ in entries}
        if len(counts) != 1:
            raise ValueError(f"Mixed shard counts for molecule {molecule_index}: {counts}")
        shard_count = counts.pop()
        shard_indices = [shard_index for shard_index, _, _ in entries]
        if shard_indices != list(range(shard_count)):
            raise FileNotFoundError(
                f"Incomplete shards for molecule {molecule_index}: {shard_indices}"
            )

        hashes = []
        conformer_count = 0
        for _, _, path in entries:
            with path.open("rb") as handle:
                payload = pickle.load(handle)
            if not isinstance(payload, dict) or list(payload) != [expected_smiles]:
                raise ValueError(f"Unexpected payload key(s) in {path}")
            conformers = payload[expected_smiles]
            if not isinstance(conformers, list):
                raise TypeError(f"Expected a conformer list in {path}")
            conformer_count += len(conformers)
            hashes.extend(coordinate_hash(mol) for mol in conformers)

        if conformer_count != args.expected_conformers:
            raise ValueError(
                f"Molecule {molecule_index} has {conformer_count} conformers; "
                f"expected {args.expected_conformers}"
            )
        unique_hashes = set(hashes)
        duplicate_count = len(hashes) - len(unique_hashes)
        if duplicate_count:
            raise ValueError(
                f"Molecule {molecule_index} has {duplicate_count} duplicate coordinate hashes"
            )
        overlap = all_coordinate_hashes.intersection(unique_hashes)
        if overlap:
            raise ValueError(
                f"Molecule {molecule_index} reuses {len(overlap)} coordinate hashes"
            )
        all_coordinate_hashes.update(unique_hashes)
        coordinate_hashes_by_molecule[molecule_index] = hashes
        layout[molecule_index] = shard_count
        molecules.append(
            {
                "molecule_index": molecule_index,
                "smiles": expected_smiles,
                "shard_count": shard_count,
                "conformer_count": conformer_count,
                "unique_coordinate_hashes": len(unique_hashes),
                "duplicate_coordinate_hashes": duplicate_count,
                "effective_seeds": [
                    args.base_seed + molecule_index * 10_000 + shard_index
                    for shard_index in range(shard_count)
                ],
            }
        )

    unexpected = sorted(set(parsed).difference(range(len(smiles))))
    if unexpected:
        raise ValueError(f"Unexpected molecule indices in shards: {unexpected}")
    merged_path = args.generation_dir / "generation_results.pickle"
    summary_path = args.generation_dir / "generation_summary.json"
    if not merged_path.is_file() or not summary_path.is_file():
        raise FileNotFoundError("Merged generation output and summary are required")
    with merged_path.open("rb") as handle:
        merged = pickle.load(handle)
    if not isinstance(merged, dict) or list(merged) != smiles:
        raise ValueError("Merged generation output does not preserve the frozen SMILES order")
    for molecule_index, expected_smiles in enumerate(smiles):
        merged_conformers = merged[expected_smiles]
        merged_hashes = [coordinate_hash(mol) for mol in merged_conformers]
        if merged_hashes != coordinate_hashes_by_molecule[molecule_index]:
            raise ValueError(
                "Merged coordinates or ordering differ from accepted shards for "
                f"molecule {molecule_index}"
            )

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("n_molecules") != len(smiles):
        raise ValueError("Generation summary has an unexpected molecule count")
    if summary.get("total_conformers") != len(all_coordinate_hashes):
        raise ValueError("Generation summary has an unexpected conformer count")
    summarized_layout = {
        int(row["molecule_index"]): int(row["part_count"])
        for row in summary.get("molecules", [])
    }
    if summarized_layout != layout:
        raise ValueError("Generation summary shard layout differs from accepted shards")

    layout_groups: dict[str, list[int]] = {}
    for molecule_index, shard_count in layout.items():
        layout_groups.setdefault(str(shard_count), []).append(molecule_index)
    job_ids = [value for value in args.slurm_jobs.split(",") if value]
    manifest = {
        "schema_version": 1,
        "status": "complete",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "method": "LoQI",
        "generation_input": str(args.generation_input.resolve()),
        "generation_input_sha256": sha256(args.generation_input),
        "generation_output": str(merged_path.resolve()),
        "generation_output_sha256": sha256(merged_path),
        "generation_summary_sha256": sha256(summary_path),
        "merged_shard_order_validated": True,
        "base_seed": args.base_seed,
        "effective_seed_formula": "base_seed + molecule_index * 10000 + shard_index",
        "n_molecules": len(smiles),
        "conformers_per_molecule": args.expected_conformers,
        "total_conformers": len(all_coordinate_hashes),
        "unique_coordinate_hashes": len(all_coordinate_hashes),
        "duplicate_coordinate_hashes": 0,
        "shard_file_count": len(paths),
        "shard_layout": layout_groups,
        "ignored_unsharded_part_count": len(ignored_unsharded_parts),
        "ignored_unsharded_parts": [
            {"path": str(path.resolve()), "sha256": sha256(path)}
            for path in ignored_unsharded_parts
        ],
        "slurm_job_ids": job_ids,
        "recovery_notes": args.recovery_note,
        "molecules": molecules,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(
        f"Validated {len(smiles)} molecules, {len(paths)} shards, "
        f"{len(all_coordinate_hashes)} unique conformers -> {args.output}"
    )


if __name__ == "__main__":
    main()
