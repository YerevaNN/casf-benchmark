#!/usr/bin/env python3
"""Replace the public PLINDER torsion rows with the ETKDG-independent pools.

Copies the current public dashboard databases and overwrites
``torsion_raw_druglike`` and ``torsion_minimized_22_druglike`` in the Plinder
tables. The row labels stay the same, so the existing dashboard still finds
them. CASF tables are copied unchanged. The reference-seeded values are removed
from these rows.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import sqlite3
from pathlib import Path

import pandas as pd

from casf_benchmark.paths import REPO_ROOT

BASELINE = REPO_ROOT / "data/results/staging/v2-1p7b-coordtoken-v1"
OUTPUT = REPO_ROOT / "data/results/staging/independent-torsion-v1"
PANEL = REPO_ROOT / "data/results/plinder_panel"
LABELS = {
    "torsion_independent_raw_druglike": "torsion_raw_druglike",
    "torsion_independent_minimized_druglike": "torsion_minimized_22_druglike",
}
DISPLAY = {
    "torsion_raw_druglike": "Torsion perturb (raw)",
    "torsion_minimized_22_druglike": "Torsion perturb (minimized, PLINDER-22)",
}
RECIPE = {
    "torsion_raw_druglike": (
        "five ETKDGv3 seeds from isomeric SMILES, 200 x 5 accepted conformers, "
        "+/-120 deg, DG clash 0.7; no experimental coordinates"
    ),
    "torsion_minimized_22_druglike": (
        "MMFF94s of each independent raw conformer, one-to-one, status 0 or 1; "
        "scored on PLINDER-22"
    ),
}
PLINDER_FILES = {
    "plinder_summary": "summary.csv",
    "plinder_per_molecule": "per_molecule.csv",
    "plinder_recovery_curve": "recovery.csv",
    "plinder_k_efficiency": "k_efficiency.csv",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_replacement(source_label: str, public_label: str, filename: str) -> pd.DataFrame:
    frame = pd.read_csv(PANEL / source_label / filename)
    frame["label"] = public_label
    frame["display_label"] = DISPLAY[public_label]
    frame["generator"] = "Torsion perturb"
    frame["recipe"] = RECIPE[public_label]
    return frame


def replace_rows(connection: sqlite3.Connection, table: str, frame: pd.DataFrame, label: str) -> None:
    dest_cols = [row[1] for row in connection.execute(f'PRAGMA table_info("{table}")')]
    extra = [column for column in frame.columns if column not in dest_cols]
    if extra:
        raise AssertionError(f"{table} has columns the published schema lacks: {extra}")
    before_other = connection.execute(
        f'SELECT COUNT(*) FROM "{table}" WHERE label != ?',
        (label,),
    ).fetchone()[0]
    connection.execute(f'DELETE FROM "{table}" WHERE label = ?', (label,))
    out = frame.copy()
    for column in dest_cols:
        if column not in out.columns:
            out[column] = pd.NA
    out[dest_cols].to_sql(table, connection, if_exists="append", index=False)
    after_other = connection.execute(
        f'SELECT COUNT(*) FROM "{table}" WHERE label != ?',
        (label,),
    ).fetchone()[0]
    if before_other != after_other:
        raise AssertionError(f"{table} changed rows for another method")
    got = connection.execute(
        f'SELECT COUNT(*) FROM "{table}" WHERE label = ?',
        (label,),
    ).fetchone()[0]
    if got != len(frame):
        raise AssertionError(f"{table} {label} has {got} rows, expected {len(frame)}")


def _median(series: pd.Series) -> float:
    return float(pd.to_numeric(series, errors="coerce").median())


def _mean(series: pd.Series) -> float:
    return float(pd.to_numeric(series, errors="coerce").mean())


def update_druglike_summary(connection: sqlite3.Connection, public_label: str, molecules: pd.DataFrame) -> None:
    total_confs = int(molecules["n_generated"].sum())
    total_pb = int(molecules["n_pb_pass"].sum())
    assignments = {
        "recipe": RECIPE[public_label],
        "display_label": DISPLAY[public_label],
        "n_molecules": int(len(molecules)),
        "total_confs": total_confs,
        "total_pb_pass": total_pb,
        "overall_pb_pass_rate": total_pb / total_confs,
        "mean_per_molecule_pb_pass_rate": _mean(molecules["pb_pass_rate"]),
        "pb_check_fail_counts_json": "{}",
        "cov_r_mean": _mean(molecules["cov_r_0p75"]),
        "cov_r_median": _median(molecules["cov_r_0p75"]),
        "cov_p_mean": _mean(molecules["cov_p_0p75"]),
        "cov_p_median": _median(molecules["cov_p_0p75"]),
        "mat_r_mean": _mean(molecules["mat_r"]),
        "mat_r_median": _median(molecules["mat_r"]),
        "mat_p_mean": _mean(molecules["mat_p"]),
        "mat_p_median": _median(molecules["mat_p"]),
        "cmat_r_mean": _mean(molecules["mat_r"]),
        "cmat_r_median": _median(molecules["mat_r"]),
        "cmat_p_mean": _mean(molecules["mat_p"]),
        "cmat_p_median": _median(molecules["mat_p"]),
        "total_true_confs": int(molecules["num_true_confs"].sum()),
        "n_missing_molecules": 0,
        "molecule_success_rate": 1.0,
    }
    sets = ", ".join(f'"{column}" = ?' for column in assignments)
    connection.execute(
        f'UPDATE extended_druglike_summary SET {sets} WHERE label = ?',
        (*assignments.values(), public_label),
    )


def update_druglike_molecules(connection: sqlite3.Connection, public_label: str, molecules: pd.DataFrame) -> None:
    current = pd.read_sql(
        "SELECT smiles FROM extended_druglike_per_molecule WHERE label = ?",
        connection,
        params=(public_label,),
    )
    if set(current["smiles"].astype(str)) != set(molecules["smiles"].astype(str)):
        raise AssertionError(f"{public_label} druglike molecules do not match the new pool")
    for record in molecules.to_dict(orient="records"):
        connection.execute(
            """
            UPDATE extended_druglike_per_molecule
            SET n_generated = ?, n_pb_pass = ?, pb_pass_rate = ?,
                cov_r_075 = ?, cov_p_075 = ?, mat_r = ?, mat_p = ?,
                num_true_confs = ?, num_gen_confs = ?, num_valid_rmsd_pairs = ?,
                pb_check_fail_counts_json = '{}'
            WHERE label = ? AND smiles = ?
            """,
            (
                int(record["n_generated"]),
                int(record["n_pb_pass"]),
                float(record["pb_pass_rate"]),
                float(record["cov_r_0p75"]),
                float(record["cov_p_0p75"]),
                float(record["mat_r"]),
                float(record["mat_p"]),
                int(record["num_true_confs"]),
                int(record["num_gen_confs"]),
                int(record["num_true_confs"]) * int(record["num_gen_confs"]),
                public_label,
                str(record["smiles"]),
            ),
        )


def update_druglike_k(connection: sqlite3.Connection, public_label: str, curves: pd.DataFrame) -> None:
    # The druglike K table has no threshold column. Its curves are the 0.75 Å draw.
    selected = curves[pd.to_numeric(curves["threshold"], errors="coerce").round(2).eq(0.75)].copy()
    if selected.empty:
        raise AssertionError(f"{public_label} has no 0.75 Å K curve")
    existing = connection.execute(
        "SELECT COUNT(*) FROM extended_druglike_k_efficiency WHERE label = ?",
        (public_label,),
    ).fetchone()[0]
    if existing == 0:
        return
    connection.execute("DELETE FROM extended_druglike_k_efficiency WHERE label = ?", (public_label,))
    rows = []
    for record in selected.to_dict(orient="records"):
        rows.append(
            {
                "label": public_label,
                "generator": "Torsion perturb",
                "display_label": DISPLAY[public_label],
                "model_size": None,
                "tokenizer": None,
                "recipe": RECIPE[public_label],
                "step": None,
                "k": int(record["k"]),
                "n_seeds": int(record["n_seeds"]),
                "n_molecules": int(record["n_molecules"]),
                "cov_r_mean": float(record["cov_r_mean"]),
                "cov_p_mean": float(record["cov_p_mean"]),
                "mat_r_mean": float(record["mat_r_mean"]),
                "mat_p_mean": float(record["mat_p_mean"]),
                "cov_r_median": float(record["cov_r_median"]),
                "cov_p_median": float(record["cov_p_median"]),
                "mat_r_median": float(record["mat_r_median"]),
                "mat_p_median": float(record["mat_p_median"]),
            }
        )
    frame = pd.DataFrame(rows)
    replace_rows(connection, "extended_druglike_k_efficiency", frame, public_label)


def publish(baseline: Path, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    casf_src = baseline / "casf_analysis_dashboard.sqlite"
    extended_src = baseline / "extended_casf_analysis.sqlite"
    casf_dest = output / casf_src.name
    extended_dest = output / extended_src.name
    shutil.copy2(casf_src, casf_dest)
    shutil.copy2(extended_src, extended_dest)
    if sha256_file(casf_dest) != sha256_file(casf_src):
        raise AssertionError("CASF dashboard database changed while copying")

    with sqlite3.connect(extended_dest) as connection:
        casf_families = connection.execute(
            "SELECT COUNT(*) FROM extended_frontier_summary WHERE family IN ('torsion_raw', 'torsion_minimized')"
        ).fetchone()[0]
        for source_label, public_label in LABELS.items():
            molecules = load_replacement(source_label, public_label, "per_molecule.csv")
            if len(molecules) != 22 or int(molecules["num_true_confs"].sum()) != 1067:
                raise AssertionError(f"{public_label} is not the 22-molecule, 1,067-reference panel")
            if int(molecules["n_generated"].min()) != 1000 or int(molecules["n_generated"].max()) != 1000:
                raise AssertionError(f"{public_label} is not 1,000 conformers per molecule")
            for table, filename in PLINDER_FILES.items():
                replace_rows(connection, table, load_replacement(source_label, public_label, filename), public_label)
            update_druglike_summary(connection, public_label, molecules)
            update_druglike_molecules(connection, public_label, molecules)
            update_druglike_k(connection, public_label, load_replacement(source_label, public_label, "k_efficiency.csv"))
        summary = pd.read_sql(
            "SELECT label, cov_r_0p75_mean, cov_p_0p75_mean, mat_r_mean, mat_p_mean, mean_pb_pass_rate, recipe "
            "FROM plinder_summary WHERE label LIKE 'torsion%' AND label NOT LIKE 'torsional%'",
            connection,
        )
        print(summary.to_string(index=False))
        after_families = connection.execute(
            "SELECT COUNT(*) FROM extended_frontier_summary WHERE family IN ('torsion_raw', 'torsion_minimized')"
        ).fetchone()[0]
        if after_families != casf_families:
            raise AssertionError("CASF torsion rows in the extended database changed")
    if sha256_file(casf_dest) != sha256_file(casf_src):
        raise AssertionError("CASF dashboard database was modified")
    print(f"CASF sqlite unchanged: {sha256_file(casf_dest)}")
    print(f"Wrote {extended_dest}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, default=BASELINE)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    args = parser.parse_args()
    publish(args.baseline_dir.resolve(), args.output_dir.resolve())


if __name__ == "__main__":
    main()
