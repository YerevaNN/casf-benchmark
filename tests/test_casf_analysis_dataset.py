from __future__ import annotations

from pathlib import Path
import sqlite3

import pandas as pd

from casf_benchmark.analysis.dataset import (
    build_comparison_master,
    build_global_per_ligand_long,
    build_master_frame,
)
from casf_benchmark.cli.build_dashboard_db import build_dashboard_db
from casf_benchmark.stratum_bins import assign_stratum


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def write_config(path: Path, sources: list[tuple[str, str, str, Path]], refs: list[tuple[str, str, str, Path]] | None = None) -> None:
    lines = ["sources:"]
    for run_id, label, ligand_set, root in sources:
        lines.extend(
            [
                f"  - run_id: {run_id}",
                f"    label: {label}",
                f"    ligand_set: {ligand_set}",
                f"    root: {root}",
            ]
        )
    if refs:
        lines.append("reference_sources:")
        for run_id, label, ligand_set, root in refs:
            lines.extend(
                [
                    f"  - run_id: {run_id}",
                    f"    label: {label}",
                    f"    ligand_set: {ligand_set}",
                    f"    root: {root}",
                ]
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def gen_row(mol_id: str, method: str, ligand_set: str, rot: int, heavy: int, clusters: int = 2) -> dict[str, object]:
    family = method.removesuffix("_fixed").removesuffix("_dynamic").removesuffix("_chembl_count")
    tier = method.rsplit("_", 1)[-1]
    if method.endswith("_chembl_count"):
        tier = "chembl_count"
        family = method[: -len("_chembl_count")]
    return {
        "mol_id": mol_id,
        "ligand_set": ligand_set,
        "row_type": "generation",
        "family": family,
        "tier": tier,
        "variant": "4b_revisited" if family == "qwen_4b_revisited" else "",
        "method": method,
        "source": method,
        "display_label": family,
        "rotatable_bonds": rot,
        "heavy_atoms": heavy,
        "conformer_count": 10,
        "pb_input_confs": 10,
        "pb_pass_confs": 9,
        "pb_fail_confs": 1,
        "pb_pass_rate": 0.9,
        "generated_candidates": 12,
        "kept_confs": 9,
        "target_confs": 10,
        "greedy_clusters_1p0": clusters,
        "clusters_per_100_1p0": clusters * 10,
        "energy_median": 2.0,
        "casf_hit_0p5": 0.5,
        "casf_opt_hit_0p5": 0.25 if ligand_set == "core" else float("nan"),
    }


def ref_row(mol_id: str, source: str, ligand_set: str, rot: int, heavy: int) -> dict[str, object]:
    return {
        "mol_id": mol_id,
        "ligand_set": ligand_set,
        "row_type": "reference",
        "family": source,
        "tier": "reference",
        "variant": "",
        "method": source,
        "source": source,
        "display_label": source,
        "rotatable_bonds": rot,
        "heavy_atoms": heavy,
        "conformer_count": 1,
        "clash_input_confs": 1,
        "clash_fail_confs": 0,
        "pb_input_confs": 1,
        "pb_pass_confs": 1,
        "pb_fail_confs": 0,
        "pb_pass_rate": 1.0,
        "greedy_clusters_1p0": 1,
        "clusters_per_100_1p0": 100.0,
        "energy_median": 1.0,
        "casf_hit_0p5": 1.0,
        "casf_opt_hit_0p5": 1.0 if ligand_set == "core" else float("nan"),
    }


def test_build_master_from_per_ligand_long(tmp_path):
    root = tmp_path / "qwen_core"
    rows = [
        gen_row("lig_a", "qwen_4b_revisited_fixed", "core", 2, 18),
        gen_row("lig_b", "qwen_4b_revisited_fixed", "core", 7, 25, clusters=4),
    ]
    write_csv(root / "analysis" / "tables" / "geometric_per_ligand_long.csv", rows)
    config = tmp_path / "sources.yaml"
    write_config(config, [("qwen_core", "Qwen Core", "core", root)])

    global_long = build_global_per_ligand_long(config, tmp_path / "long.csv")
    master = build_comparison_master(global_long)

    assert len(global_long) == 2
    assert len(master) == 1
    row = master.iloc[0]
    assert row["family"] == "qwen_4b_revisited"
    assert row["variant"] == "4b_revisited"
    assert row["ligands"] == 2
    assert row["mean_clusters_1p0"] == 3.0
    assert row["selected_pool_total"] == 20.0


def test_build_master_frame_keeps_reference_rows_separate(tmp_path):
    gen_root = tmp_path / "loqi_core"
    ref_root = tmp_path / "reference_core"
    write_csv(
        gen_root / "analysis" / "tables" / "geometric_per_ligand_long.csv",
        [gen_row("lig_a", "loqi_raw_fixed", "core", 2, 18)],
    )
    write_csv(
        ref_root / "analysis" / "tables" / "geometric_per_ligand_long.csv",
        [ref_row("lig_a", source, "core", 2, 18) for source in ["casf_crystal", "chembl3d_gt"]],
    )
    config = tmp_path / "sources.yaml"
    write_config(
        config,
        [("loqi_core", "LOQI Core", "core", gen_root)],
        refs=[("reference_core", "Reference Core", "core", ref_root)],
    )

    _run_count, master = build_master_frame(config)

    assert set(master["row_type"]) == {"generation", "reference"}
    assert set(master["source"]) == {"loqi_raw_fixed", "casf_crystal", "chembl3d_gt"}
    assert len(master) == 3


def test_dashboard_db_writes_strata_and_omits_empty_bins(tmp_path):
    root = tmp_path / "loqi_core"
    ref_root = tmp_path / "reference_core"
    write_csv(
        root / "analysis" / "tables" / "geometric_per_ligand_long.csv",
        [
            gen_row("lig_a", "loqi_raw_fixed", "core", 2, 18),
            gen_row("lig_b", "loqi_raw_fixed", "core", 7, 25),
        ],
    )
    write_csv(
        ref_root / "analysis" / "tables" / "geometric_per_ligand_long.csv",
        [
            ref_row("lig_a", "chembl3d_gt", "core", 2, 18),
            ref_row("lig_b", "chembl3d_gt", "core", 7, 25),
        ],
    )
    config = tmp_path / "sources.yaml"
    write_config(
        config,
        [("loqi_core", "LOQI Core", "core", root)],
        refs=[("reference_core", "Reference Core", "core", ref_root)],
    )

    db_path = tmp_path / "dashboard.sqlite"
    build_dashboard_db(config, db_path)

    with sqlite3.connect(db_path) as connection:
        rows = pd.read_sql_query("SELECT * FROM comparison_rows", connection)
        strata = pd.read_sql_query("SELECT * FROM comparison_strata", connection)
        user_version = connection.execute("PRAGMA user_version").fetchone()[0]

    assert user_version == 2
    assert set(rows["row_type"]) == {"generation", "reference"}
    rb = strata[(strata["breakdown"] == "rotatable_bonds") & (strata["view_tier"] == "fixed")]
    assert "15+" not in set(rb["stratum"])
    loqi_total = rows[rows["method"] == "loqi_raw_fixed"].iloc[0]["ligands"]
    loqi_strata = rb[rb["method"] == "loqi_raw_fixed"]["ligands"].sum()
    assert loqi_strata == loqi_total


def test_assign_stratum_preserves_series_index():
    values = pd.Series([2, 7], index=[100, 200])

    strata = assign_stratum(values, "rotatable_bonds")

    assert list(strata.index) == [100, 200]
    assert list(strata.astype(str)) == ["0-4", "5-9"]
