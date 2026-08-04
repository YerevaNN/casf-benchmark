from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "extended_casf_analysis.py"
SPEC = importlib.util.spec_from_file_location("extended_casf_analysis", SCRIPT_PATH)
extended = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = extended
SPEC.loader.exec_module(extended)


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
