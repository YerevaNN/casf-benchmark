#!/usr/bin/env python3
"""Audit exact molecular-identity overlap with the Qwen training CSV."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

import pandas as pd
from rdkit import Chem

from casf_benchmark.analysis.pharmacophore import sha256


DEFAULT_TRAINING = Path(
    "/mnt/weka/vtarasov/dedup_casf285_plinder23_20261007/"
    "protomer/final/merged_train_clean_shuffled.csv"
)


def run(args: argparse.Namespace) -> None:
    audit = pd.read_csv(args.reference_audit, keep_default_na=False)
    identities = audit[["ccd", "name", "expected_smiles"]].drop_duplicates()
    identities["inchikey"] = identities["expected_smiles"].map(
        lambda smiles: Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))
    )
    key_to_row = identities.set_index("inchikey").to_dict(orient="index")
    args.output.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8") as patterns:
        for key in sorted(key_to_row):
            patterns.write(key + "\n")
        patterns.flush()
        process = subprocess.Popen(
            ["rg", "-F", "-f", patterns.name, str(args.training_csv)],
            stdout=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        assert process.stdout is not None
        counts: Counter[str] = Counter()
        source_counts: Counter[tuple[str, str]] = Counter()
        examples: dict[str, dict[str, str]] = {}
        with args.training_csv.open() as training_handle:
            header = next(csv.reader([training_handle.readline()]))
        positions = {name: header.index(name) for name in ("inchikey", "dataset", "source_name")}
        for line in process.stdout:
            row = next(csv.reader([line]))
            key = row[positions["inchikey"]]
            if key not in key_to_row:
                continue
            counts[key] += 1
            source = row[positions["source_name"]]
            source_counts[(key, source)] += 1
            examples.setdefault(
                key,
                {
                    "dataset": row[positions["dataset"]],
                    "source_name": source,
                },
            )
        return_code = process.wait()
        if return_code not in {0, 1}:
            raise RuntimeError(f"rg failed with exit {return_code}")

    rows = []
    for key, identity in key_to_row.items():
        row = {
            **identity,
            "inchikey": key,
            "training_row_count": counts[key],
            "exact_identity_overlap": counts[key] > 0,
            "example_dataset": examples.get(key, {}).get("dataset", ""),
            "example_source_name": examples.get(key, {}).get("source_name", ""),
        }
        rows.append(row)
    pd.DataFrame(rows).sort_values(["exact_identity_overlap", "ccd"], ascending=[False, True]).to_csv(
        args.output / "exact_identity_overlap.csv", index=False
    )
    pd.DataFrame(
        [
            {"inchikey": key, "source_name": source, "training_row_count": count}
            for (key, source), count in source_counts.items()
        ]
    ).to_csv(args.output / "exact_identity_overlap_by_source.csv", index=False)
    manifest = {
        "reference_audit": str(args.reference_audit.resolve()),
        "reference_audit_sha256": sha256(args.reference_audit),
        "training_csv": str(args.training_csv.resolve()),
        "training_csv_size_bytes": args.training_csv.stat().st_size,
        "training_csv_mtime_ns": args.training_csv.stat().st_mtime_ns,
        "molecule_count": len(key_to_row),
        "molecules_with_exact_identity_overlap": sum(counts[key] > 0 for key in key_to_row),
        "total_matching_training_rows": sum(counts.values()),
        "scope": (
            "exact isomeric InChIKey only; this does not establish absence of near-neighbor "
            "or scaffold overlap"
        ),
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2), flush=True)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--reference-audit", type=Path, required=True)
    result.add_argument("--training-csv", type=Path, default=DEFAULT_TRAINING)
    result.add_argument("--output", type=Path, required=True)
    return result


if __name__ == "__main__":
    run(parser().parse_args())
