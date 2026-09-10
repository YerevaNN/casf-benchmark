#!/usr/bin/env python3
"""Combine per-checkpoint druglike-eval CSVs (from eval_druglike_conformers.py)
into a single SQLite database.

Tables:
  summary       -- one row per checkpoint (overall PB pass rate, fail-mode counts)
  per_molecule  -- one row per (checkpoint, molecule): PB pass rate, diversity, energy

Re-runnable: drops and recreates both tables from the current CSVs each time.
"""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--eval-dir",
        type=Path,
        required=True,
        help="Directory holding the *_summary.csv / *_per_molecule.csv written by "
        "eval_druglike_conformers.py",
    )
    parser.add_argument("--db-path", type=Path, default=None, help="Default: {eval-dir}/druglike_eval.sqlite")
    args = parser.parse_args()

    db_path = args.db_path or (args.eval_dir / "druglike_eval.sqlite")

    summary_frames = [
        pd.read_csv(f) for f in sorted(args.eval_dir.glob("*_summary.csv")) if not f.name.startswith("_")
    ]

    per_molecule_frames = []
    for f in sorted(args.eval_dir.glob("*_per_molecule.csv")):
        if f.name.startswith("_"):
            continue
        frame = pd.read_csv(f)
        if "label" not in frame.columns:
            # Older eval_druglike_conformers.py runs didn't stamp the checkpoint
            # label onto each per-molecule row; recover it from the filename.
            frame.insert(0, "label", f.name.removesuffix("_per_molecule.csv"))
        per_molecule_frames.append(frame)

    if not summary_frames:
        raise SystemExit(f"No *_summary.csv files found under {args.eval_dir}")

    summary_df = pd.concat(summary_frames, ignore_index=True).sort_values(
        "overall_pb_pass_rate", ascending=False
    )
    per_molecule_df = pd.concat(per_molecule_frames, ignore_index=True) if per_molecule_frames else pd.DataFrame()

    con = sqlite3.connect(str(db_path))
    summary_df.to_sql("summary", con, if_exists="replace", index=False)
    if not per_molecule_df.empty:
        per_molecule_df.to_sql("per_molecule", con, if_exists="replace", index=False)
    con.execute("CREATE INDEX IF NOT EXISTS idx_per_molecule_label ON per_molecule(label)")
    con.execute("CREATE INDEX IF NOT EXISTS idx_per_molecule_smiles ON per_molecule(smiles)")
    con.commit()
    con.close()

    print(f"Wrote {len(summary_df)} summary row(s), {len(per_molecule_df)} per-molecule row(s) -> {db_path}")


if __name__ == "__main__":
    main()
