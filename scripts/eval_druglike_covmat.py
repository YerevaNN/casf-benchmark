#!/usr/bin/env python3
"""Compute GEOM-style COV/MAT metrics for a generated druglike conformer pool.

Both inputs use the established SMILES-keyed pickle contract. Ground-truth
values are metadata dictionaries containing ``confs``; generated values are
lists of RDKit molecules. Outputs match the files consumed by
``build_druglike_covmat.py``.
"""

from __future__ import annotations

import argparse
import pickle
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign


def best_rmsd(generated: Chem.Mol, reference: Chem.Mol) -> float:
    try:
        return float(rdMolAlign.GetBestRMS(generated, reference))
    except Exception:
        return float("nan")


def remove_hydrogens(mols: list[Chem.Mol]) -> list[Chem.Mol]:
    cleaned = []
    for mol in mols:
        if mol is None or mol.GetNumConformers() == 0:
            continue
        try:
            cleaned.append(Chem.RemoveHs(mol))
        except Exception:
            continue
    return cleaned


def compute_rmsd_matrices(
    ground_truth: dict, generated: dict, workers: int
) -> tuple[dict[str, np.ndarray], list[str], list[str]]:
    matrices: dict[str, np.ndarray] = {}
    missing = []
    work = []
    for smiles, row in ground_truth.items():
        generated_mols = remove_hydrogens(generated.get(smiles, []))
        if not generated_mols:
            missing.append(smiles)
            continue
        references = remove_hydrogens(row.get("confs", []))
        matrix = np.full((len(references), len(generated_mols)), np.nan, dtype=np.float32)
        matrices[smiles] = matrix
        for reference_index, reference in enumerate(references):
            work.append((smiles, reference_index, reference, generated_mols))

    def compute_row(item):
        smiles, reference_index, reference, generated_mols = item
        values = np.asarray(
            [best_rmsd(mol, reference) for mol in generated_mols], dtype=np.float32
        )
        return smiles, reference_index, values

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(compute_row, item) for item in work]
        completed = 0
        for future in as_completed(futures):
            smiles, reference_index, values = future.result()
            matrices[smiles][reference_index] = values
            completed += 1
            if completed % 100 == 0 or completed == len(futures):
                print(f"RMSD rows: {completed}/{len(futures)}", flush=True)

    all_nan = [smiles for smiles, matrix in matrices.items() if matrix.size == 0 or np.isnan(matrix).all()]
    return matrices, missing, all_nan


def finite_axis_min(matrix: np.ndarray, axis: int) -> np.ndarray:
    if matrix.size == 0:
        return np.asarray([], dtype=float)
    finite = np.isfinite(matrix)
    replaced = np.where(finite, matrix, np.inf)
    values = replaced.min(axis=axis)
    values[~np.isfinite(values)] = np.nan
    return values


def per_molecule_metrics(
    matrix: np.ndarray, threshold: float, dmax: float
) -> dict[str, float | int]:
    valid = matrix[np.isfinite(matrix)]
    min_true = finite_axis_min(matrix, axis=1)
    min_generated = finite_axis_min(matrix, axis=0)
    finite_true = min_true[np.isfinite(min_true)]
    finite_generated = min_generated[np.isfinite(min_generated)]

    def coverage(values: np.ndarray) -> float:
        return float(np.mean(values < threshold)) if len(values) else float("nan")

    def matching(values: np.ndarray) -> float:
        return float(np.mean(values)) if len(values) else float("nan")

    def censored_matching(values: np.ndarray, expected_size: int) -> float:
        censored = np.full(expected_size, dmax, dtype=float)
        if expected_size:
            finite_mask = np.isfinite(values)
            censored[finite_mask] = np.minimum(values[finite_mask], dmax)
        return float(np.mean(censored)) if expected_size else float("nan")

    return {
        "min_rmsd": float(valid.min()) if len(valid) else float("nan"),
        "max_rmsd": float(valid.max()) if len(valid) else float("nan"),
        "avg_rmsd": float(valid.mean()) if len(valid) else float("nan"),
        "cov_r_075": coverage(finite_true),
        "cov_p_075": coverage(finite_generated),
        "mat_r": matching(finite_true),
        "mat_p": matching(finite_generated),
        "cmat_r": censored_matching(min_true, matrix.shape[0]),
        "cmat_p": censored_matching(min_generated, matrix.shape[1]),
        "num_true_confs": int(matrix.shape[0]),
        "num_gen_confs": int(matrix.shape[1]),
        "num_valid_rmsd_pairs": int(len(valid)),
    }


def aggregate(rows: list[dict]) -> dict[str, float]:
    def stat(column: str, operation) -> float:
        values = pd.to_numeric(pd.Series([row[column] for row in rows]), errors="coerce")
        return float(operation(values))

    return {
        "cov_r_mean": stat("cov_r_075", np.nanmean),
        "cov_r_median": stat("cov_r_075", np.nanmedian),
        "cov_p_mean": stat("cov_p_075", np.nanmean),
        "cov_p_median": stat("cov_p_075", np.nanmedian),
        "mat_r_mean": stat("mat_r", np.nanmean),
        "mat_r_median": stat("mat_r", np.nanmedian),
        "mat_p_mean": stat("mat_p", np.nanmean),
        "mat_p_median": stat("mat_p", np.nanmedian),
        "cmat_r_mean": stat("cmat_r", np.nanmean),
        "cmat_r_median": stat("cmat_r", np.nanmedian),
        "cmat_p_mean": stat("cmat_p", np.nanmean),
        "cmat_p_median": stat("cmat_p", np.nanmedian),
    }


def write_report(
    path: Path,
    *,
    gen_pickle: Path,
    gt_pickle: Path,
    ground_truth: dict,
    generated: dict,
    missing: list[str],
    all_nan: list[str],
    metrics: dict[str, float],
    threshold: float,
    dmax: float,
    elapsed: float,
) -> None:
    successful = len(ground_truth) - len(set(missing) | set(all_nan))
    success_rate = successful / len(ground_truth) if ground_truth else float("nan")
    text = f"""\
================================================================================
COVMAT EVALUATION RESULTS
================================================================================

EVALUATION SUMMARY
----------------------------------------
Processed file: {gen_pickle}
Ground truth file: {gt_pickle}
Total molecules generated: {len(generated)}
Total conformers generated: {sum(len(value) for value in generated.values())}
Total molecules in ground truth: {len(ground_truth)}
Total conformers in ground truth: {sum(len(value.get("confs", [])) for value in ground_truth.values())}
Missing molecules (no conformers): {len(missing)}
All-NaN RMSD keys: {len(all_nan)}

EXECUTION RUNTIME
----------------------------------------
Total processing: {elapsed / 60:.4f} min
CovMat processing: {elapsed / 60:.4f} min
PoseBusters processing: 0 min

COVERAGE AND RECALL METRICS
----------------------------------------
Threshold: {threshold:g}
Coverage-Recall (COV-R):
  Mean:   {metrics["cov_r_mean"]:.4f}
  Median: {metrics["cov_r_median"]:.4f}
Coverage-Precision (COV-P):
  Mean:   {metrics["cov_p_mean"]:.4f}
  Median: {metrics["cov_p_median"]:.4f}
Matching-Recall (MAT-R):
  Mean:   {metrics["mat_r_mean"]:.4f}
  Median: {metrics["mat_r_median"]:.4f}
Matching-Precision (MAT-P):
  Mean:   {metrics["mat_p_mean"]:.4f}
  Median: {metrics["mat_p_median"]:.4f}

Censored matching (restricted mean, dmax={dmax:g} A)
Molecule success rate: {success_rate:.4f}
Censored Matching-Recall (cMAT-R):
  Mean:   {metrics["cmat_r_mean"]:.4f}
  Median: {metrics["cmat_r_median"]:.4f}
Censored Matching-Precision (cMAT-P):
  Mean:   {metrics["cmat_p_mean"]:.4f}
  Median: {metrics["cmat_p_median"]:.4f}
"""
    path.write_text(text, encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ground-truth-pickle", type=Path, required=True)
    parser.add_argument("--gen-pickle", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--num-workers", type=int, default=80)
    parser.add_argument("--threshold", type=float, default=0.75)
    parser.add_argument("--dmax", type=float, default=3.0)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    RDLogger.DisableLog("rdApp.*")
    with args.ground_truth_pickle.open("rb") as handle:
        ground_truth = pickle.load(handle)
    with args.gen_pickle.open("rb") as handle:
        generated = pickle.load(handle)

    started = time.perf_counter()
    matrices, missing, all_nan = compute_rmsd_matrices(
        ground_truth, generated, workers=args.num_workers
    )
    rows = []
    for smiles, row in ground_truth.items():
        if smiles not in matrices:
            continue
        metrics = per_molecule_metrics(matrices[smiles], args.threshold, args.dmax)
        rows.append({"geom_smiles": smiles, **metrics, "sub_smiles": ""})
    summary = aggregate(rows)
    elapsed = time.perf_counter() - started

    args.out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out_dir / "rmsd_matrix.csv", index=False)
    with (args.out_dir / "rmsd_matrix.pickle").open("wb") as handle:
        pickle.dump(
            {"threshold": args.threshold, "dmax": args.dmax, "summary": summary},
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
    write_report(
        args.out_dir / "covmat_results.txt",
        gen_pickle=args.gen_pickle,
        gt_pickle=args.ground_truth_pickle,
        ground_truth=ground_truth,
        generated=generated,
        missing=missing,
        all_nan=all_nan,
        metrics=summary,
        threshold=args.threshold,
        dmax=args.dmax,
        elapsed=elapsed,
    )
    print(f"Wrote COV/MAT outputs to {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
