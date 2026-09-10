#!/usr/bin/env python3
"""Merge COV/MAT eval outputs into the druglike-eval tables and publish them for
the dashboard.

Writes `extended_druglike_summary` / `extended_druglike_per_molecule` into the
extended-analysis sidecar DB, which surfaces them as tabs in the Streamlit
dashboard's "Extended Analysis" section (see streamlit_app.py::EXTENDED_TABLES).
This targets the sidecar rather than the main `casf_analysis_dashboard.sqlite` on
purpose: extended_* tables in the main DB would flip the dashboard's auto-detected
"Extended DB" default away from the sidecar (see default_extended_db_path()),
hiding the extended_* tables that already live there.

Inputs:
  - the druglike_eval.sqlite built by build_druglike_eval_db.py
  - {run_dir}/eval_*/covmat_results.txt and rmsd_matrix.csv, one run_dir per
    run, resolved from the druglike entries in generation_runs.yaml

Generator-agnostic: `label` is whatever the catalog says, so pushing RDKit, loqi or
any other method into the same two tables is a new YAML entry plus its eval outputs,
never a change here.

Outputs:
  - druglike_eval.sqlite summary/per_molecule enriched with COV-R/COV-P/MAT-R/MAT-P
    (and censored variants)
  - extended_druglike_summary / extended_druglike_per_molecule in the sidecar DB
"""

from __future__ import annotations

import argparse
import re
import sqlite3
from pathlib import Path

import pandas as pd

from casf_benchmark.catalog import (
    RUN_DESCRIPTOR_FIELDS as DESCRIPTOR_COLUMNS,
    describe_run,
    load_generation_run_entries,
)
from casf_benchmark.paths import DEFAULT_EXTENDED_DB

AGG_PATTERN = re.compile(
    r"\((?P<abbr>[A-Za-z-]+)\):\s*\n\s*Mean:\s*(?P<mean>[-\d.]+)\s*\n\s*Median:\s*(?P<median>[-\d.]+)"
)


def find_eval_dir(run_dir: Path) -> Path | None:
    candidates = sorted(run_dir.glob("eval_*"))
    return candidates[-1] if candidates else None


def parse_covmat_summary(text: str) -> dict:
    out: dict[str, float | int | None] = {}
    for m in AGG_PATTERN.finditer(text):
        abbr = m.group("abbr").lower().replace("-", "_")
        out[f"{abbr}_mean"] = float(m.group("mean"))
        out[f"{abbr}_median"] = float(m.group("median"))

    def grab_int(pattern: str) -> int | None:
        mm = re.search(pattern, text)
        return int(mm.group(1)) if mm else None

    def grab_float(pattern: str) -> float | None:
        mm = re.search(pattern, text)
        return float(mm.group(1)) if mm else None

    out["threshold"] = grab_float(r"Threshold:\s*([\d.]+)")
    out["molecule_success_rate"] = grab_float(r"Molecule success rate:\s*([\d.]+)")
    out["total_true_confs"] = grab_int(r"Total conformers in ground truth:\s*(\d+)")
    out["total_molecules_ground_truth"] = grab_int(r"Total molecules in ground truth:\s*(\d+)")
    out["n_missing_molecules"] = grab_int(r"Missing molecules \(no conformers\):\s*(\d+)")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--generation-results-root",
        type=Path,
        required=True,
        help="Root holding the per-run inference output directories named in "
        "generation_runs.yaml",
    )
    parser.add_argument(
        "--druglike-db",
        type=Path,
        required=True,
        help="druglike_eval.sqlite written by build_druglike_eval_db.py",
    )
    parser.add_argument(
        "--extended-db",
        type=Path,
        default=DEFAULT_EXTENDED_DB,
        help="Extended-analysis sidecar DB, alongside the other extended_* tables the "
        "dashboard reads (default: %(default)s)",
    )
    args = parser.parse_args()

    con = sqlite3.connect(str(args.druglike_db))
    summary = pd.read_sql("select * from summary", con)
    per_molecule = pd.read_sql("select * from per_molecule", con)
    con.close()

    # Idempotent re-runs: drop any columns this script previously added, so
    # merging never collides with a prior run's output already sitting in
    # druglike_eval.sqlite.
    parsed_cols = list(DESCRIPTOR_COLUMNS)
    covmat_summary_cols = [
        "threshold", "molecule_success_rate", "total_true_confs",
        "total_molecules_ground_truth", "n_missing_molecules",
        "cov_r_mean", "cov_r_median", "cov_p_mean", "cov_p_median",
        "mat_r_mean", "mat_r_median", "mat_p_mean", "mat_p_median",
        "cmat_r_mean", "cmat_r_median", "cmat_p_mean", "cmat_p_median",
    ]
    covmat_per_mol_cols = [
        "min_rmsd", "max_rmsd", "avg_rmsd", "cov_r_075", "cov_p_075",
        "mat_r", "mat_p", "num_true_confs", "num_gen_confs", "num_valid_rmsd_pairs",
    ]
    summary = summary.drop(columns=parsed_cols + covmat_summary_cols, errors="ignore")
    per_molecule = per_molecule.drop(columns=covmat_per_mol_cols, errors="ignore")

    summary_rows = []
    per_mol_rows = []
    descriptor_rows = []
    missing = []
    for entry in load_generation_run_entries():
        dirname = entry.cohorts.get("druglike")
        if not dirname:
            continue
        label = entry.label
        run_dir = args.generation_results_root / dirname
        eval_dir = find_eval_dir(run_dir)
        if eval_dir is None:
            missing.append(label)
            continue

        agg_text = (eval_dir / "covmat_results.txt").read_text()
        agg = parse_covmat_summary(agg_text)
        agg["label"] = label
        summary_rows.append(agg)
        descriptor_rows.append({"label": label, **describe_run(label, entry.descriptors)})

        rmsd = pd.read_csv(eval_dir / "rmsd_matrix.csv")
        rmsd = rmsd.rename(columns={"geom_smiles": "smiles"})
        rmsd.insert(0, "label", label)
        per_mol_rows.append(rmsd)

    if missing:
        print(f"WARNING: no eval_* dir found for: {missing} (skipped)")

    covmat_summary = pd.DataFrame(summary_rows)
    covmat_per_mol = pd.concat(per_mol_rows, ignore_index=True) if per_mol_rows else pd.DataFrame()

    summary_enriched = summary.merge(covmat_summary, on="label", how="left")
    descriptors = pd.DataFrame(descriptor_rows, columns=["label", *DESCRIPTOR_COLUMNS])
    summary_enriched = summary_enriched.merge(descriptors, on="label", how="left")
    # Descriptors read better as the leading columns of the table.
    summary_enriched = summary_enriched[
        [*DESCRIPTOR_COLUMNS, *[c for c in summary_enriched.columns if c not in DESCRIPTOR_COLUMNS]]
    ]
    per_molecule_enriched = per_molecule.merge(
        covmat_per_mol.drop(columns=["sub_smiles"], errors="ignore"),
        on=["label", "smiles"],
        how="left",
    )

    # Write back into the standalone druglike_eval.sqlite
    con = sqlite3.connect(str(args.druglike_db))
    summary_enriched.to_sql("summary", con, if_exists="replace", index=False)
    per_molecule_enriched.to_sql("per_molecule", con, if_exists="replace", index=False)
    con.execute("CREATE INDEX IF NOT EXISTS idx_per_molecule_label ON per_molecule(label)")
    con.execute("CREATE INDEX IF NOT EXISTS idx_per_molecule_smiles ON per_molecule(smiles)")
    con.commit()
    con.close()
    print(f"Enriched {args.druglike_db} (summary={len(summary_enriched)}, per_molecule={len(per_molecule_enriched)})")

    # Publish into the extended-analysis sidecar DB as extended_* tables
    con = sqlite3.connect(str(args.extended_db))
    summary_enriched.to_sql("extended_druglike_summary", con, if_exists="replace", index=False)
    per_molecule_enriched.to_sql("extended_druglike_per_molecule", con, if_exists="replace", index=False)
    con.commit()
    con.close()
    print(f"Published extended_druglike_summary / extended_druglike_per_molecule -> {args.extended_db}")


if __name__ == "__main__":
    main()
