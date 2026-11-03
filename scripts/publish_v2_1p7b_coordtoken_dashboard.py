#!/usr/bin/env python3
"""Append the Qwen 1.7B v2 CoordToken run to the public dashboard databases.

The core generation is the usual 94-ligand CASF core. The ref generation scored
the frozen 1,044-crystal selection, so those rows are stored as ``ref_dedup``
and are not added to the 1,236-entry ref cohort. Plinder panel tables are
appended. The strict druglike tables stay as they are: this campaign did not
score that protocol.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import sqlite3
from pathlib import Path

import pandas as pd

from casf_benchmark.analysis.dataset import build_comparison_master, build_global_per_ligand_long
from casf_benchmark.analysis.ref_dedup import DEFAULT_MANIFEST, load_selected_ids
from casf_benchmark.cli.build_dashboard_db import build_comparison_strata, validate_strata_ligands
from casf_benchmark.paths import REPO_ROOT

FAMILY = "qwen_1p7b_v2_coordtoken_step41054"
WORK_ROOT = Path("/mnt/weka/mbedrosian/pharma_generation_analysis/v2_1p7b_coordtoken_step41054")
RUNS = WORK_ROOT / "casf_benchmark_runs"
PANEL = WORK_ROOT / "plinder_panel" / FAMILY
BASELINE_TAG = "dashboard-data-v2-1p7b-binned-pm12-v1"
OUTPUT = REPO_ROOT / "data/results/staging/v2-1p7b-coordtoken-v1"
BASELINE_SHA256 = {
    "casf_analysis_dashboard.sqlite": "d823a25ae12b3ce4edc1d85846ad59f431b823d66147e5f304740df6d5a04f31",
    "extended_casf_analysis.sqlite": "f7cbb130296f8bc07691d9888a8857372a00c058860fe0bca82fc743b34896bf",
}
PLINDER_TABLES = {
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


def append_aligned(connection: sqlite3.Connection, table: str, frame: pd.DataFrame) -> int:
    if frame.empty:
        return 0
    dest_cols = [row[1] for row in connection.execute(f'PRAGMA table_info("{table}")')]
    extra = [column for column in frame.columns if column not in dest_cols]
    if extra:
        raise ValueError(f"{table} would add columns that the published schema does not have: {extra}")
    out = frame.copy()
    for column in dest_cols:
        if column not in out.columns:
            out[column] = pd.NA
    out[dest_cols].to_sql(table, connection, if_exists="append", index=False)
    return len(out)


def source_config(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "sources:",
                f"  - run_id: {FAMILY}_core",
                '    label: "v2 1.7B CoordToken (core)"',
                "    ligand_set: core",
                f"    root: {RUNS / f'{FAMILY}_core'}",
                f"  - run_id: {FAMILY}_ref",
                '    label: "v2 1.7B CoordToken (ref dedup)"',
                "    ligand_set: ref_dedup",
                f"    root: {RUNS / f'{FAMILY}_ref'}",
                "",
            ]
        ),
        encoding="utf-8",
    )


def analysis_frames(config: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    selected = set(load_selected_ids(DEFAULT_MANIFEST))
    if len(selected) != 1044:
        raise AssertionError(f"selection manifest has {len(selected)} ligands, expected 1044")
    long = build_global_per_ligand_long(config, output_csv=None)
    if set(long["ligand_set"].astype(str)) != {"core", "ref_dedup"}:
        raise AssertionError(f"unexpected ligand sets: {sorted(set(long['ligand_set'].astype(str)))}")
    if set(long["family"].astype(str)) != {FAMILY}:
        raise AssertionError("analysis tables contain another family")
    ref = long["ligand_set"].astype(str).eq("ref_dedup")
    ref_ids = set(long.loc[ref, "mol_id"].astype(str))
    if ref_ids != selected:
        raise AssertionError(
            f"ref generation is not the frozen selection: missing {len(selected - ref_ids)}, extra {len(ref_ids - selected)}"
        )
    long.loc[ref, "run_id"] = f"{FAMILY}_ref__ref_dedup"
    counts = long.groupby(["ligand_set", "tier"])["mol_id"].nunique().to_dict()
    expected = {
        ("core", "fixed"): 94,
        ("core", "dynamic"): 94,
        ("core", "chembl_count"): 94,
        ("ref_dedup", "fixed"): 1044,
        ("ref_dedup", "dynamic"): 1044,
        ("ref_dedup", "chembl_count"): 1044,
    }
    if counts != expected:
        raise AssertionError(f"ligand counts {counts} != {expected}")
    fixed = long[long["tier"].astype(str).eq("fixed")]
    if int((pd.to_numeric(fixed["conformer_count"], errors="coerce") <= 0).sum()):
        raise AssertionError("a fixed-tier ligand has no conformers")
    if int(pd.to_numeric(fixed["casf_best_rmsd"], errors="coerce").isna().sum()):
        raise AssertionError("a fixed-tier ligand has no crystal RMSD")
    comparison = build_comparison_master(long, missing_as_miss_ligand_sets=frozenset({"ref_dedup"}))
    strata = build_comparison_strata(long, missing_as_miss_ligand_sets=frozenset({"ref_dedup"}))
    validate_strata_ligands(comparison, strata)
    fixed_rows = comparison[comparison["tier"].astype(str).eq("fixed")].set_index("ligand_set")
    if int(fixed_rows.loc["ref_dedup", "defined_rmsd_ligands"]) != 1044:
        raise AssertionError("ref dedup does not have a defined RMSD for every selected ligand")
    if int(fixed_rows.loc["core", "ligands"]) != 94 or int(fixed_rows.loc["ref_dedup", "ligands"]) != 1044:
        raise AssertionError("fixed-tier comparison rows have the wrong ligand counts")
    if comparison["ligand_set"].astype(str).eq("ref").any():
        raise AssertionError("ref rows leaked into the new comparison table")
    return long, comparison, strata


def plinder_frames() -> dict[str, pd.DataFrame]:
    frames = {table: pd.read_csv(PANEL / filename) for table, filename in PLINDER_TABLES.items()}
    summary = frames["plinder_summary"]
    molecules = frames["plinder_per_molecule"]
    if list(summary["label"]) != [FAMILY] or int(summary["n_molecules"].iloc[0]) != 22:
        raise AssertionError("Plinder summary is not the 22-molecule 1.7B CoordToken run")
    if len(molecules) != 22 or int(molecules["num_true_confs"].sum()) != 1067:
        raise AssertionError("Plinder references are not the cleaned 1,067-conformer set")
    if molecules["name"].astype(str).str.contains("ortezomib", case=False).any():
        raise AssertionError("Bortezomib is still in the Plinder table")
    return frames


def assert_baseline(path: Path) -> None:
    digest = sha256_file(path)
    expected = BASELINE_SHA256[path.name]
    if digest != expected:
        raise AssertionError(f"{path.name} sha256 {digest} != {BASELINE_TAG} {expected}")


def publish(baseline_dir: Path, output_dir: Path) -> None:
    config = output_dir / "sources.yaml"
    output_dir.mkdir(parents=True, exist_ok=True)
    source_config(config)
    long, comparison, strata = analysis_frames(config)
    plinder = plinder_frames()
    for name in ("casf_analysis_dashboard.sqlite", "extended_casf_analysis.sqlite"):
        source = baseline_dir / name
        assert_baseline(source)
        destination = output_dir / name
        shutil.copy2(source, destination)
    with sqlite3.connect(output_dir / "casf_analysis_dashboard.sqlite") as connection:
        already = connection.execute(
            "SELECT COUNT(*) FROM comparison_rows WHERE family = ?",
            (FAMILY,),
        ).fetchone()[0]
        if already:
            raise AssertionError(f"baseline already has {already} comparison rows for {FAMILY}")
        before = {
            table: connection.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone()[0]
            for table in ("comparison_rows", "paper_sampling_curves", "paper_energy_summary")
        }
        sentinel = connection.execute(
            "SELECT ligands FROM comparison_rows WHERE family = ? AND ligand_set = 'ref' AND tier = 'fixed' AND row_type = 'generation'",
            ("qwen_0p6b_v2_binned_pm12_step66057",),
        ).fetchone()
        if sentinel is None or int(sentinel[0]) != 1236:
            raise AssertionError("the 0.6B v2 ref row is not the 1,236-entry cohort")
        print("comparison_rows", append_aligned(connection, "comparison_rows", comparison))
        print("comparison_strata", append_aligned(connection, "comparison_strata", strata))
        print("per_ligand_long", append_aligned(connection, "per_ligand_long", long))
        after_sentinel = connection.execute(
            "SELECT ligands FROM comparison_rows WHERE family = ? AND ligand_set = 'ref' AND tier = 'fixed' AND row_type = 'generation'",
            ("qwen_0p6b_v2_binned_pm12_step66057",),
        ).fetchone()
        if int(after_sentinel[0]) != 1236:
            raise AssertionError("publishing changed an existing ref row")
        binned = connection.execute(
            "SELECT ligands FROM comparison_rows WHERE family = ? AND ligand_set = 'ref_dedup' AND tier = 'fixed' AND row_type = 'generation'",
            ("qwen_1p7b_v2_binned_pm12_step66057",),
        ).fetchone()
        if binned is None or int(binned[0]) != 1044:
            raise AssertionError("publishing changed the 1.7B binned ref-dedup row")
        if connection.execute(f'SELECT COUNT(*) FROM paper_sampling_curves').fetchone()[0] != before["paper_sampling_curves"]:
            raise AssertionError("paper sampling curves changed")
        if connection.execute(f'SELECT COUNT(*) FROM paper_energy_summary').fetchone()[0] != before["paper_energy_summary"]:
            raise AssertionError("paper energy summary changed")
        added = connection.execute(
            "SELECT ligand_set, tier, ligands FROM comparison_rows WHERE family = ? AND row_type = 'generation' ORDER BY ligand_set, tier",
            (FAMILY,),
        ).fetchall()
        print("added", added)
    with sqlite3.connect(output_dir / "extended_casf_analysis.sqlite") as connection:
        druglike_before = connection.execute("SELECT COUNT(*) FROM extended_druglike_summary").fetchone()[0]
        for table, frame in plinder.items():
            present = connection.execute(f'SELECT COUNT(*) FROM "{table}" WHERE label = ?', (FAMILY,)).fetchone()[0]
            if present:
                raise AssertionError(f"{table} already contains {FAMILY}")
            print(table, append_aligned(connection, table, frame))
        druglike_after = connection.execute("SELECT COUNT(*) FROM extended_druglike_summary").fetchone()[0]
        if druglike_before != druglike_after:
            raise AssertionError("druglike summary changed")
        molecules = connection.execute(
            "SELECT COUNT(*), SUM(num_true_confs) FROM plinder_per_molecule WHERE label = ?",
            (FAMILY,),
        ).fetchone()
        if int(molecules[0]) != 22 or int(molecules[1]) != 1067:
            raise AssertionError(f"published Plinder molecules are {molecules}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    args = parser.parse_args()
    publish(args.baseline_dir.resolve(), args.output_dir.resolve())


if __name__ == "__main__":
    main()
