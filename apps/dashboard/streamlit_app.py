#!/usr/bin/env python3
"""Streamlit dashboard for PB-once CASF analysis outputs."""

from __future__ import annotations

import os
import sqlite3
from pathlib import Path

import pandas as pd
import streamlit as st

from casf_benchmark.dashboard_table_help import render_table_help
from casf_benchmark.paths import DEFAULT_DASHBOARD_DB, DEFAULT_EXTENDED_DB as _DEFAULT_EXTENDED_DB

_DEFAULT_DB = DEFAULT_DASHBOARD_DB
DEFAULT_DB = Path(os.environ.get("CASF_DASHBOARD_DB", str(_DEFAULT_DB)))
DEFAULT_EXTENDED_DB = Path(os.environ.get("CASF_EXTENDED_DB", str(_DEFAULT_EXTENDED_DB)))
TIERS = ("fixed", "dynamic", "chembl_count")
BREAKDOWN_LABELS = {
    "Total": "total",
    "Rotatable bonds": "rotatable_bonds",
    "Heavy atoms": "heavy_atoms",
}
EXTENDED_TABLES = {
    "Paired vs ChEMBL3D-PB": {
        "table": "extended_paired_delta_vs_chembl3d_pb",
        "columns": [
            "ligand_set",
            "source",
            "family",
            "tier",
            "method",
            "stratum_type",
            "stratum",
            "n_common_ligands",
            "mean_delta_casf_best_rmsd",
            "median_delta_casf_best_rmsd",
            "win_rate_0p1A",
            "loss_rate_0p1A",
            "delta_hit_0p5",
            "delta_hit_0p75",
            "delta_hit_2p0",
            "catastrophic_loss_count_gt1A",
            "large_rescue_count_gt1A",
        ],
    },
    "Extremes vs ChEMBL3D-PB": {
        "table": "extended_paired_delta_extreme_ligands_vs_chembl3d_pb",
        "columns": [
            "ligand_set",
            "source",
            "family",
            "tier",
            "method",
            "baseline_source",
            "direction",
            "mol_id",
            "baseline_casf_best_rmsd",
            "method_casf_best_rmsd",
            "delta_casf_best_rmsd",
            "baseline_casf_hit_0p75",
            "method_casf_hit_0p75",
        ],
    },
    "Paired vs RDKit raw": {
        "table": "extended_paired_delta_vs_rdkit_raw",
        "columns": [
            "ligand_set",
            "source",
            "family",
            "tier",
            "method",
            "stratum_type",
            "stratum",
            "n_common_ligands",
            "mean_delta_casf_best_rmsd",
            "median_delta_casf_best_rmsd",
            "win_rate_0p1A",
            "loss_rate_0p1A",
            "delta_hit_0p5",
            "delta_hit_0p75",
            "delta_hit_2p0",
            "catastrophic_loss_count_gt1A",
            "large_rescue_count_gt1A",
        ],
    },
    "Extremes vs RDKit raw": {
        "table": "extended_paired_delta_extreme_ligands_vs_rdkit_raw",
        "columns": [
            "ligand_set",
            "source",
            "family",
            "tier",
            "method",
            "baseline_source",
            "direction",
            "mol_id",
            "baseline_casf_best_rmsd",
            "method_casf_best_rmsd",
            "delta_casf_best_rmsd",
            "baseline_casf_hit_0p75",
            "method_casf_hit_0p75",
        ],
    },
    "Frontier": {
        "table": "extended_frontier_summary",
        "columns": [
            "ligand_set",
            "source",
            "family",
            "tier",
            "method",
            "n_ligands",
            "mean_casf_hit_0p75",
            "mean_casf_best_rmsd",
            "mean_pb_pass_rate",
            "mean_conformer_count",
            "mean_greedy_clusters_1p0",
            "useful_clusters_per_100_confs",
            "hit_0p75_per_100_confs",
            "valid_hit_score",
            "frontier_rank",
        ],
    },
    "PB failures": {
        "table": "extended_pb_failure_mechanisms",
        "columns": [
            "ligand_set",
            "source",
            "family",
            "tier",
            "method",
            "n_ligands",
            "total_pb_fail_rate",
            "median_per_ligand_pb_pass_rate",
            "ligands_pb_pass_rate_lt_0p9",
            "ligands_pb_pass_rate_lt_0p5",
            "dominant_pb_failure_test",
            "second_dominant_pb_failure_test",
            "fraction_failures_tetrahedral_chirality",
            "fraction_failures_internal_steric_clash",
            "fraction_failures_energy_ratio",
            "fraction_failures_bond_geometry_tests",
        ],
    },
    "PB tail ligands": {
        "table": "extended_pb_failure_tail_ligands",
        "columns": [
            "ligand_set",
            "source",
            "family",
            "tier",
            "method",
            "mol_id",
            "rotatable_bonds",
            "heavy_atoms",
            "pb_pass_rate",
            "pb_input_confs",
            "top_three_failed_pb_tests",
            "casf_best_rmsd",
            "casf_hit_0p75",
        ],
    },
    "Difficulty": {
        "table": "extended_bound_pose_difficulty",
        "columns": [
            "ligand_set",
            "source",
            "family",
            "tier",
            "method",
            "baseline_source",
            "difficulty_bin",
            "n_common_ligands",
            "mean_delta_casf_best_rmsd",
            "hit_0p75",
            "rescue_rate",
            "catastrophic_loss_rate",
        ],
    },
    "Rescue cases": {
        "table": "extended_rescue_cases",
        "columns": [
            "ligand_set",
            "source",
            "family",
            "tier",
            "method",
            "mol_id",
            "rotatable_bonds",
            "heavy_atoms",
            "baseline_casf_best_rmsd",
            "method_casf_best_rmsd",
            "delta_casf_best_rmsd",
            "pb_pass_rate",
            "greedy_clusters_1p0",
        ],
    },
    "K efficiency": {
        "table": "extended_k_efficiency",
        "columns": [
            "ligand_set",
            "source",
            "family",
            "tier",
            "method",
            "sampling_mode",
            "summary_mode",
            "k",
            "status",
            "skip_reason",
            "n_ligands",
            "n_metric_ligands",
            "n_ligands_with_pb",
            "mean_n_used",
            "mean_best_rmsd_at_k",
            "hit_0p25_at_k",
            "hit_0p75_at_k",
            "hit_2p0_at_k",
            "ci_low",
            "ci_high",
            "mean_pb_pass_rate_at_k",
            "mean_conformers_available",
            "n_ligands_with_at_least_k",
        ],
    },
    "K per ligand": {
        "table": "extended_chembl_k_efficiency_per_ligand",
        "columns": [
            "ligand_set",
            "source",
            "family",
            "tier",
            "method",
            "sampling_mode",
            "seed",
            "k",
            "mol_id",
            "n_available_pb",
            "n_used",
            "casf_best_rmsd_at_k",
            "casf_hit_0p25_at_k",
            "casf_hit_0p5_at_k",
            "casf_hit_0p75_at_k",
            "casf_hit_2p0_at_k",
        ],
    },
    "K comparisons": {
        "table": "extended_k_efficiency_comparisons",
        "columns": [
            "ligand_set",
            "source",
            "family",
            "tier",
            "method",
            "comparison",
            "baseline_source",
            "baseline_sampling_mode",
            "k",
            "n_common_ligands",
            "mean_delta_best_rmsd_at_k",
            "median_delta_best_rmsd_at_k",
            "delta_hit_0p25_at_k",
            "delta_hit_0p5_at_k",
            "delta_hit_0p75_at_k",
            "delta_hit_2p0_at_k",
        ],
    },
    "Energy windows": {
        "table": "extended_energy_window_summary",
        "columns": [
            "ligand_set",
            "source",
            "family",
            "tier",
            "method",
            "energy_window",
            "n_ligands",
            "n_ligands_with_window_confs",
            "mean_fraction_window_confs",
            "mean_n_window_confs",
            "mean_best_rmsd_window",
            "hit_0p25_window",
            "hit_0p5_window",
            "hit_0p75_window",
            "hit_2p0_window",
            "mean_clusters_1p0_window",
            "mean_clusters_per_100_1p0_window",
            "mean_entropy_1p0_window",
            "mean_largest_cluster_fraction_1p0_window",
            "hit_retention_deltaE_10",
            "cluster_retention_deltaE_10",
            "mean_useful_low_energy_clusters",
        ],
    },
    "Energy per ligand": {
        "table": "extended_energy_window_per_ligand",
        "columns": [
            "ligand_set",
            "source",
            "family",
            "tier",
            "method",
            "mol_id",
            "energy_window",
            "n_pb_confs",
            "n_window_confs",
            "fraction_window_confs",
            "window_energy_min",
            "window_energy_max_delta",
            "casf_best_rmsd_window",
            "casf_hit_0p25_window",
            "casf_hit_0p5_window",
            "casf_hit_0p75_window",
            "casf_hit_2p0_window",
            "greedy_clusters_1p0_window",
            "clusters_per_100_1p0_window",
            "cluster_entropy_1p0_window",
            "largest_cluster_fraction_1p0_window",
            "useful_low_energy_clusters",
        ],
    },
    "Energy comparisons": {
        "table": "extended_energy_window_comparisons",
        "columns": [
            "ligand_set",
            "source",
            "family",
            "tier",
            "method",
            "comparison",
            "n_ligands",
            "delta_hit_0p75_window",
            "high_energy_hit_fraction",
        ],
    },
    "K RMSD manifest": {
        "table": "extended_per_conformer_casf_rmsd_manifest",
        "columns": [
            "ligand_set",
            "run_id",
            "run_label",
            "source",
            "family",
            "tier",
            "method",
            "mol_id",
            "conformers_available",
            "rotatable_bonds",
            "heavy_atoms",
        ],
    },
    "Sanity": {
        "table": "extended_sanity_checks",
        "columns": ["check", "severity", "status", "affected_rows", "details"],
    },
}

TABLE_HEIGHT_PX = 520


@st.cache_data(show_spinner=False)
def load_table(db_path: str, table: str, db_mtime_ns: int) -> pd.DataFrame:
    del db_mtime_ns
    with sqlite3.connect(db_path) as connection:
        return pd.read_sql_query(f'SELECT * FROM "{table}"', connection)


@st.cache_data(show_spinner=False)
def load_table_names(db_path: str, db_mtime_ns: int) -> set[str]:
    del db_mtime_ns
    with sqlite3.connect(db_path) as connection:
        rows = connection.execute(
            "SELECT name FROM sqlite_master WHERE type IN ('table', 'view')"
        ).fetchall()
    return {str(row[0]) for row in rows}


def select_view(frame: pd.DataFrame, ligand_set: str, tier: str, family: str) -> pd.DataFrame:
    out = frame[frame["ligand_set"].astype(str) == ligand_set].copy()
    if "view_tier" in out.columns:
        out = out[out["view_tier"].astype(str) == tier]
    else:
        out = out[
            ((out["row_type"].astype(str) == "generation") & (out["tier"].astype(str) == tier))
            | (out["row_type"].astype(str) == "reference")
        ]
    if family != "All":
        out = out[out["family"].astype(str) == family]
    return sort_rows(out)


def sort_rows(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame
    out = frame.copy()
    out["_row_type_order"] = out["row_type"].map({"generation": 0, "reference": 1}).fillna(9)
    out["_tier_order"] = out["tier"].map({"fixed": 0, "dynamic": 1, "chembl_count": 2, "reference": 3}).fillna(9)
    sort_cols = [col for col in ("stratum", "_row_type_order", "family", "_tier_order", "method") if col in out.columns]
    out = out.sort_values(sort_cols, kind="stable")
    return out.drop(columns=["_row_type_order", "_tier_order"]).reset_index(drop=True)


def _build_column_config(frame: pd.DataFrame, columns: list[str]) -> dict[str, st.column_config.Column]:
    """Column display config: pin method, format numbers to three decimals."""
    config: dict[str, st.column_config.Column] = {}
    for column in columns:
        if column not in frame.columns:
            continue
        series = frame[column]
        if column == "method":
            config[column] = st.column_config.TextColumn(column, pinned=True, width="large")
        elif pd.api.types.is_bool_dtype(series):
            config[column] = st.column_config.CheckboxColumn(column)
        elif pd.api.types.is_integer_dtype(series):
            config[column] = st.column_config.NumberColumn(column, format="%d")
        elif pd.api.types.is_float_dtype(series):
            config[column] = st.column_config.NumberColumn(column, format="%.3f")
    return config


def _numeric_columns(frame: pd.DataFrame, columns: list[str]) -> list[str]:
    return [
        column
        for column in columns
        if column in frame.columns
        and pd.api.types.is_numeric_dtype(frame[column])
        and not pd.api.types.is_bool_dtype(frame[column])
    ]


def _apply_color_scale(frame: pd.DataFrame, columns: list[str]):
    numeric = _numeric_columns(frame, columns)
    if not numeric:
        return frame
    return frame.style.background_gradient(cmap="viridis", subset=numeric, axis=0)


def _apply_table_view_controls(data: pd.DataFrame, cols: list[str], name: str) -> tuple[pd.DataFrame, list[str], bool]:
    """Let users hide columns and filter rows before rendering."""
    with st.expander("Columns & rows", expanded=False):
        visible = st.multiselect(
            "Visible columns",
            cols,
            default=cols,
            key=f"{name}_visible_cols",
        )
        visible = [column for column in cols if column in visible]
        if not visible:
            st.warning("Select at least one column.")
            visible = cols

        filtered = data
        if "method" in filtered.columns:
            methods = sorted(filtered["method"].dropna().astype(str).unique())
            if methods:
                selected_methods = st.multiselect(
                    "Methods (rows)",
                    methods,
                    default=methods,
                    key=f"{name}_methods",
                )
                filtered = filtered[filtered["method"].astype(str).isin(selected_methods)]
        elif "display_label" in filtered.columns:
            labels = sorted(filtered["display_label"].dropna().astype(str).unique())
            if labels:
                selected_labels = st.multiselect(
                    "Labels (rows)",
                    labels,
                    default=labels,
                    key=f"{name}_display_labels",
                )
                filtered = filtered[filtered["display_label"].astype(str).isin(selected_labels)]

        if "stratum" in filtered.columns and filtered["stratum"].nunique() > 1:
            strata = sorted(filtered["stratum"].dropna().astype(str).unique())
            selected_strata = st.multiselect(
                "Strata (rows)",
                strata,
                default=strata,
                key=f"{name}_strata",
            )
            filtered = filtered[filtered["stratum"].astype(str).isin(selected_strata)]

        if "mol_id" in filtered.columns:
            mol_filter = st.text_input(
                "Filter mol_id (substring)",
                key=f"{name}_mol_id_filter",
            )
            if mol_filter.strip():
                filtered = filtered[
                    filtered["mol_id"].astype(str).str.contains(mol_filter.strip(), case=False, na=False)
                ]

        color_by_value = st.checkbox(
            "Color numeric columns by value",
            value=False,
            key=f"{name}_color_scale",
            help="Heatmap per column (low → high). Sorting, pinning, and CSV export are unchanged.",
        )

    return filtered.reset_index(drop=True), visible, color_by_value


def display_table(frame: pd.DataFrame, columns: list[str], name: str) -> None:
    cols = [column for column in columns if column in frame.columns]
    color_by_value = False
    if not cols:
        data = frame
        visible_cols = list(data.columns)
    else:
        data, visible_cols, color_by_value = _apply_table_view_controls(frame[cols], cols, name)

    display_data = _apply_color_scale(data, visible_cols) if color_by_value else data
    st.dataframe(
        display_data,
        column_config=_build_column_config(data, visible_cols),
        column_order=visible_cols,
        use_container_width=True,
        height=TABLE_HEIGHT_PX,
        hide_index=True,
    )
    st.download_button(
        f"Download {name} CSV",
        data=data.to_csv(index=False).encode("utf-8"),
        file_name=f"{name}.csv",
        mime="text/csv",
    )
    render_table_help(name)


def filter_extended_table(
    frame: pd.DataFrame,
    ligand_set: str,
    tier: str,
    family: str,
    stratum_type: str,
) -> pd.DataFrame:
    out = frame.copy()
    if ligand_set != "All" and "ligand_set" in out.columns:
        out = out[out["ligand_set"].astype(str) == ligand_set]
    if tier != "All" and "tier" in out.columns:
        out = out[out["tier"].astype(str) == tier]
    if family != "All" and "family" in out.columns:
        out = out[out["family"].astype(str) == family]
    if stratum_type != "All" and "stratum_type" in out.columns:
        out = out[out["stratum_type"].astype(str) == stratum_type]
    return out.reset_index(drop=True)


def default_extended_db_path(db_path: Path, table_names: set[str]) -> Path:
    if any(name.startswith("extended_") for name in table_names):
        return db_path
    override = os.environ.get("CASF_EXTENDED_DB") or os.environ.get("CASF_EXTENDED_ANALYSIS_DB")
    return Path(override if override else str(DEFAULT_EXTENDED_DB))


def render_extended_analysis(db_path: Path, table_names: set[str]) -> None:
    st.divider()
    st.header("Extended Analysis")

    default_path = default_extended_db_path(db_path, table_names)
    with st.sidebar.expander("Extended analysis", expanded=False):
        extended_db_path = Path(
            st.text_input("Extended DB", str(default_path), key="extended_db_path")
        ).expanduser()

    if not extended_db_path.exists():
        st.info(f"Extended DB not found: {extended_db_path}")
        return

    extended_mtime_ns = extended_db_path.stat().st_mtime_ns
    extended_table_names = load_table_names(str(extended_db_path), extended_mtime_ns)
    available = [
        label
        for label, spec in EXTENDED_TABLES.items()
        if spec["table"] in extended_table_names
    ]
    if not available:
        st.info(f"No extended_* tables found in {extended_db_path}")
        return

    first_available = load_table(
        str(extended_db_path),
        EXTENDED_TABLES[available[0]]["table"],
        extended_mtime_ns,
    )
    ligand_sets = ["All"]
    if "ligand_set" in first_available.columns:
        ligand_sets.extend(sorted(first_available["ligand_set"].dropna().astype(str).unique()))

    filter_cols = st.columns(4)
    with filter_cols[0]:
        ligand_set = st.selectbox("Extended ligand set", ligand_sets, key="extended_ligand_set")
    with filter_cols[1]:
        tier = st.selectbox("Extended tier", ["All", *TIERS], key="extended_tier")
    with filter_cols[2]:
        families = ["All"]
        if "family" in first_available.columns:
            families.extend(sorted(first_available["family"].dropna().astype(str).unique()))
        family = st.selectbox("Extended family", families, key="extended_family")
    with filter_cols[3]:
        stratum_type = st.selectbox(
            "Extended stratum",
            ["All", "total", "rotatable_bond_bin", "heavy_atom_bin", "chembl_conformer_count_bin", "baseline_casf_difficulty_bin"],
            key="extended_stratum_type",
        )

    tabs = st.tabs(available)
    for tab, label in zip(tabs, available):
        spec = EXTENDED_TABLES[label]
        with tab:
            frame = load_table(str(extended_db_path), spec["table"], extended_mtime_ns)
            view = filter_extended_table(frame, ligand_set, tier, family, stratum_type)
            display_table(view, spec["columns"], spec["table"])


def main() -> None:
    st.set_page_config(page_title="CASF Analysis Dashboard", layout="wide")
    st.title("CASF Analysis Dashboard")

    db_default = Path(os.environ.get("CASF_DASHBOARD_DB", str(DEFAULT_DB)))
    db_path = Path(st.sidebar.text_input("Dashboard DB", str(db_default))).expanduser()
    if not db_path.exists():
        st.error(f"Dashboard DB not found: {db_path}")
        st.stop()

    db_mtime_ns = db_path.stat().st_mtime_ns
    table_names = load_table_names(str(db_path), db_mtime_ns)
    comparison_rows = load_table(str(db_path), "comparison_rows", db_mtime_ns)
    comparison_strata = load_table(str(db_path), "comparison_strata", db_mtime_ns)

    ligand_sets = sorted(comparison_rows["ligand_set"].dropna().astype(str).unique())
    ligand_set = st.sidebar.selectbox("Ligand set", ligand_sets, index=0 if "core" not in ligand_sets else ligand_sets.index("core"))
    tier = st.sidebar.selectbox("Tier", ["fixed", "dynamic", "chembl_count"], index=0)
    families = ["All", *sorted(comparison_rows["family"].dropna().astype(str).unique())]
    family = st.sidebar.selectbox("Family", families)

    with st.sidebar.expander("Break down aggregates", expanded=False):
        breakdown_label = st.radio("View", list(BREAKDOWN_LABELS), index=0)
    breakdown = BREAKDOWN_LABELS[breakdown_label]

    if breakdown == "total":
        view = select_view(comparison_rows, ligand_set, tier, family)
    else:
        strata = comparison_strata[comparison_strata["breakdown"].astype(str) == breakdown]
        view = select_view(strata, ligand_set, tier, family)

    tabs = st.tabs(["Overview", "Clustering", "Energy", "CASF hits", "CASF opt hits", "Funnel"])
    identity = [
        "stratum",
        "row_type",
        "display_label",
        "family",
        "tier",
        "method",
        "ligands",
        "ligands_scope",
    ]

    with tabs[0]:
        display_table(
            view,
            [
                *identity,
                "total_confs",
                "mean_confs_per_ligand",
                "pairwise_mean",
                "pairwise_p90",
                "mean_torsion_std_deg",
            ],
            "overview",
        )

    with tabs[1]:
        display_table(
            view,
            [
                *identity,
                "mean_clusters_0p5",
                "mean_clusters_1p0",
                "mean_clusters_2p0",
                "mean_clusters_3p0",
                "clusters_per_100_1p0",
                "cluster_entropy_1p0",
                "effective_clusters_1p0",
                "largest_cluster_fraction_1p0",
                "singleton_fraction_1p0",
            ],
            "clustering",
        )

    with tabs[2]:
        display_table(
            view,
            [*identity, "energy_min", "energy_max", "energy_median", "energy_std"],
            "energy",
        )

    with tabs[3]:
        display_table(
            view,
            [
                *identity,
                "casf_best_rmsd",
                "casf_median_rmsd",
                "casf_hit_0p25",
                "casf_hit_0p5",
                "casf_hit_0p75",
                "casf_hit_2p0",
            ],
            "casf_hits",
        )

    with tabs[4]:
        display_table(
            view,
            [
                *identity,
                "casf_opt_best_rmsd",
                "casf_opt_median_rmsd",
                "casf_opt_hit_0p25",
                "casf_opt_hit_0p5",
                "casf_opt_hit_0p75",
                "casf_opt_hit_2p0",
            ],
            "casf_opt_hits",
        )

    with tabs[5]:
        funnel = select_view(comparison_rows, ligand_set, tier, family)
        funnel = funnel[funnel["row_type"].astype(str) == "generation"]
        display_table(
            funnel,
            [
                "row_type",
                "display_label",
                "family",
                "tier",
                "method",
                "ligands",
                "generated_ligands",
                "selected_ligands",
                "kept_ligands",
                "target_confs_mean",
                "selected_pool_total",
                "pb_fail_total",
                "kept_confs_total",
                "pb_fail_rate_mean",
                "kept_vs_target_rate_mean",
            ],
            "funnel",
        )

    render_extended_analysis(db_path, table_names)


if __name__ == "__main__":
    main()
