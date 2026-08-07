from __future__ import annotations

import pandas as pd

from casf_benchmark.cli import extended_analysis as extended


def test_common_ligand_pairing_excludes_non_overlapping_ligands() -> None:
    baseline = pd.DataFrame(
        {
            "mol_id": ["a", "b", "baseline_only"],
            "casf_best_rmsd": [1.0, 2.0, 3.0],
            "casf_hit_0p5": [0, 0, 0],
            "casf_hit_0p75": [0, 0, 0],
            "casf_hit_2p0": [1, 1, 0],
            "conformer_count": [10, 20, 30],
            "rotatable_bonds": [1, 2, 3],
            "heavy_atoms": [10, 20, 30],
        }
    )
    method = pd.DataFrame(
        {
            "mol_id": ["a", "b", "method_only"],
            "casf_best_rmsd": [0.8, 2.2, 0.1],
            "casf_hit_0p5": [0, 0, 1],
            "casf_hit_0p75": [0, 0, 1],
            "casf_hit_2p0": [1, 0, 1],
            "pb_pass_rate": [1.0, 0.9, 1.0],
            "greedy_clusters_1p0": [2, 3, 4],
            "rotatable_bonds": [1, 2, 3],
            "heavy_atoms": [10, 20, 30],
        }
    )

    paired = extended.add_common_pair_columns(method, baseline)

    assert paired["mol_id"].tolist() == ["a", "b"]


def test_paired_delta_signs_and_win_loss_thresholds() -> None:
    paired = pd.DataFrame(
        {
            "baseline_casf_best_rmsd": [1.0, 1.0, 1.0],
            "method_casf_best_rmsd": [0.8, 1.3, 1.05],
            "baseline_casf_hit_0p5": [0, 0, 0],
            "method_casf_hit_0p5": [1, 0, 0],
            "baseline_casf_hit_0p75": [0, 1, 0],
            "method_casf_hit_0p75": [1, 0, 0],
            "baseline_casf_hit_2p0": [1, 1, 1],
            "method_casf_hit_2p0": [1, 1, 1],
        }
    )

    summary = extended.paired_metric_summary(paired, threshold=0.10)

    assert round(summary["mean_delta_casf_best_rmsd"], 6) == round((-.2 + .3 + .05) / 3, 6)
    assert summary["win_rate_0p1A"] == 1 / 3
    assert summary["loss_rate_0p1A"] == 1 / 3
    assert summary["delta_hit_0p75"] == 0


def test_rescue_and_catastrophic_loss_definitions() -> None:
    per_ligand = pd.DataFrame(
        [
            _row("chembl3d_gt_pb", "reference", "reference", "a", 1.0, 0, 1),
            _row("chembl3d_gt_pb", "reference", "reference", "b", 0.4, 1, 1),
            _row("rdkit_random_raw_fixed", "generation", "fixed", "a", 0.5, 1, 1),
            _row("rdkit_random_raw_fixed", "generation", "fixed", "b", 2.5, 0, 0),
        ]
    )

    difficulty, rescue_cases = extended.build_bound_pose_difficulty(per_ligand)

    assert difficulty["rescue_rate"].max() == 1.0
    assert difficulty["catastrophic_loss_rate"].max() == 1.0
    assert rescue_cases["mol_id"].tolist() == ["a"]


def test_non_dominated_frontier_ranking() -> None:
    frame = pd.DataFrame(
        {
            "hit": [0.9, 0.8, 0.7],
            "valid": [0.9, 0.8, 0.95],
            "count": [100, 200, 80],
        }
    )

    ranks = extended.non_dominated_frontier_ranks(frame, ["hit", "valid", "count"], [True, True, False])

    assert ranks.tolist() == [1, 2, 1]


def test_pb_failure_json_parsing_handles_empty_and_malformed_values() -> None:
    assert extended.parse_fail_counts_json(None) == {}
    assert extended.parse_fail_counts_json("") == {}
    assert extended.parse_fail_counts_json("{bad json") == {}
    assert extended.parse_fail_counts_json('{"energy_ratio": 2, "bad": "x"}') == {"energy_ratio": 2.0}


def test_k_efficiency_ligand_rows_caps_at_available_conformers() -> None:
    rows = extended.k_efficiency_ligand_rows([1.2, 0.4], [1, 5], seed=7, subsamples=4)

    by_k = {row["k"]: row for row in rows}
    assert by_k[5]["ligand_best_rmsd_at_k"] == 0.4
    assert by_k[5]["ligand_hit_0p5_at_k"] == 1.0
    assert by_k[5]["conformers_available"] == 2
    assert by_k[5]["ligand_has_at_least_k"] == 0
    assert by_k[1]["ligand_has_at_least_k"] == 1


def test_chembl_k_aggregate_keeps_capped_and_strict_denominators() -> None:
    rows = pd.DataFrame(
        [
            _chembl_k_row("a", 5, 5, 0.4),
            _chembl_k_row("b", 2, 2, 0.6),
            _chembl_k_row("c", 0, 0, float("nan")),
        ]
    )

    summary = extended.aggregate_chembl_k_rows(rows)
    capped = summary[summary["summary_mode"] == "capped_at_available"].iloc[0]
    strict = summary[summary["summary_mode"] == "strict_at_least_k"].iloc[0]

    assert capped["n_ligands"] == 3
    assert capped["n_metric_ligands"] == 3
    assert capped["n_ligands_with_pb"] == 2
    assert capped["hit_0p75_at_k"] == 2 / 3
    assert strict["n_ligands"] == 3
    assert strict["n_metric_ligands"] == 1
    assert strict["hit_0p75_at_k"] == 1.0


def test_energy_window_aggregate_counts_window_empty_ligands_as_misses() -> None:
    rows = pd.DataFrame(
        [
            _energy_window_row("a", "deltaE_10", 10, 2, 0.5, 1),
            _energy_window_row("b", "deltaE_10", 10, 0, float("nan"), 0),
        ]
    )

    summary = extended.aggregate_energy_window_rows(rows)
    row = summary.iloc[0]

    assert row["n_ligands"] == 2
    assert row["n_ligands_with_window_confs"] == 1
    assert row["hit_0p75_window"] == 0.5


def _chembl_k_row(mol_id: str, n_available: int, n_used: int, best: float) -> dict[str, object]:
    hit = float(best <= 0.75) if pd.notna(best) else 0.0
    return {
        "ligand_set": "core",
        "run_id": "reference_core",
        "run_label": "Core Reference",
        "source": "chembl3d_gt_pb",
        "family": "chembl3d_gt_pb",
        "tier": "reference",
        "method": "chembl3d_gt_pb",
        "sampling_mode": "first_k",
        "seed": float("nan"),
        "k": 5,
        "mol_id": mol_id,
        "n_available_pb": n_available,
        "n_used": n_used,
        "casf_best_rmsd_at_k": best,
        "casf_hit_0p25_at_k": 0.0,
        "casf_hit_0p5_at_k": hit,
        "casf_hit_0p75_at_k": hit,
        "casf_hit_2p0_at_k": hit,
    }


def _energy_window_row(
    mol_id: str,
    window: str,
    n_pb: int,
    n_window: int,
    best: float,
    hit_0p75: int,
) -> dict[str, object]:
    return {
        "ligand_set": "core",
        "run_id": "run",
        "run_label": "Run",
        "source": "method_fixed",
        "family": "method",
        "tier": "fixed",
        "method": "method_fixed",
        "mol_id": mol_id,
        "energy_window": window,
        "n_pb_confs": n_pb,
        "n_window_confs": n_window,
        "fraction_window_confs": n_window / n_pb,
        "casf_best_rmsd_window": best,
        "casf_hit_0p25_window": 0,
        "casf_hit_0p5_window": 0,
        "casf_hit_0p75_window": hit_0p75,
        "casf_hit_2p0_window": hit_0p75,
        "greedy_clusters_1p0_window": float(n_window),
        "clusters_per_100_1p0_window": float(n_window),
        "cluster_entropy_1p0_window": 0.0,
        "largest_cluster_fraction_1p0_window": 1.0 if n_window else float("nan"),
        "useful_low_energy_clusters": float("nan"),
    }


def _row(
    source: str,
    row_type: str,
    tier: str,
    mol_id: str,
    best_rmsd: float,
    hit_0p75: int,
    hit_2p0: int,
) -> dict[str, object]:
    family = source.removesuffix("_fixed")
    return {
        "ligand_set": "ref",
        "run_id": "run",
        "run_label": "Run",
        "source": source,
        "family": family,
        "tier": tier,
        "method": source,
        "row_type": row_type,
        "mol_id": mol_id,
        "casf_best_rmsd": best_rmsd,
        "casf_hit_0p5": int(best_rmsd <= 0.5),
        "casf_hit_0p75": hit_0p75,
        "casf_hit_2p0": hit_2p0,
        "conformer_count": 10,
        "pb_pass_rate": 0.9,
        "greedy_clusters_1p0": 2,
        "rotatable_bonds": 4,
        "heavy_atoms": 25,
    }
