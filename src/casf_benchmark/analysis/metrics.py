"""Geometric conformer evaluation metrics for CASF benchmark analysis."""

from __future__ import annotations

import itertools
import math
import random
import statistics
from typing import Sequence

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolTransforms

PAIRWISE_SAMPLE_CAP = 120
RANDOM_SEED = 42

CASF_GEOMETRIC_CLUSTER_THRESHOLDS = (0.5, 1.0, 2.0, 3.0)


def threshold_tag(threshold: float) -> str:
    return str(threshold).replace(".", "p")


def mol_no_h(mol: Chem.Mol) -> Chem.Mol:
    try:
        return Chem.RemoveHs(mol)
    except Exception:
        try:
            return Chem.RemoveHs(mol, sanitize=False)
        except Exception:
            return mol


def best_aligned_rmsd(gen: Chem.Mol, ref: Chem.Mol) -> float:
    """Heavy-atom symmetry-aware aligned RMSD between two single-conformer mols."""
    a = mol_no_h(gen)
    b = mol_no_h(ref)
    if a.GetNumAtoms() != b.GetNumAtoms():
        return math.nan
    try:
        return float(rdMolAlign.GetBestRMS(a, b, prbId=0, refId=0))
    except Exception:
        try:
            atom_map = [(idx, idx) for idx in range(a.GetNumAtoms())]
            return float(rdMolAlign.AlignMol(a, b, prbCid=0, refCid=0, atomMap=atom_map))
        except Exception:
            return math.nan


def sample_mols(mols: Sequence[Chem.Mol], cap: int = PAIRWISE_SAMPLE_CAP) -> list[Chem.Mol]:
    if len(mols) <= cap:
        return list(mols)
    rng = random.Random(RANDOM_SEED)
    indices = sorted(rng.sample(range(len(mols)), cap))
    return [mols[i] for i in indices]


def heavy_rmsd_same_topology(mol_a: Chem.Mol, mol_b: Chem.Mol) -> float:
    return best_aligned_rmsd(mol_a, mol_b)


def greedy_cluster_assignments(mols: Sequence[Chem.Mol], threshold: float) -> list[int]:
    assignments: list[int] = []
    centers: list[Chem.Mol] = []
    for mol in mols:
        if not centers:
            centers.append(mol)
            assignments.append(0)
            continue
        distances = [heavy_rmsd_same_topology(mol, center) for center in centers]
        finite = [(idx, value) for idx, value in enumerate(distances) if math.isfinite(value)]
        if not finite:
            centers.append(mol)
            assignments.append(len(centers) - 1)
            continue
        nearest_idx, nearest_dist = min(finite, key=lambda item: item[1])
        if nearest_dist < threshold:
            assignments.append(nearest_idx)
        else:
            centers.append(mol)
            assignments.append(len(centers) - 1)
    return assignments


def cluster_distribution_metrics(assignments: Sequence[int]) -> dict[str, float]:
    if len(assignments) == 0:
        return {
            "cluster_count": 0.0,
            "largest_cluster_fraction": math.nan,
            "cluster_entropy": math.nan,
            "effective_cluster_count": math.nan,
            "simpson_concentration": math.nan,
            "singleton_fraction": math.nan,
        }
    counts = np.bincount(np.asarray(assignments, dtype=int))
    counts = counts[counts > 0]
    total = float(counts.sum())
    probs = counts / total
    entropy = -float(np.sum(probs * np.log(probs))) if len(probs) > 1 else 0.0
    normalized_entropy = entropy / math.log(len(probs)) if len(probs) > 1 else 0.0
    simpson = float(np.sum(probs * probs))
    return {
        "cluster_count": float(len(counts)),
        "largest_cluster_fraction": float(np.max(probs)),
        "cluster_entropy": normalized_entropy,
        "effective_cluster_count": float(math.exp(entropy)),
        "simpson_concentration": simpson,
        "singleton_fraction": float(np.mean(counts == 1)),
    }


def greedy_cluster_metrics(
    mols: Sequence[Chem.Mol],
    thresholds: Sequence[float] = CASF_GEOMETRIC_CLUSTER_THRESHOLDS,
) -> dict[str, float]:
    out: dict[str, float] = {"conformer_count": float(len(mols))}
    if not mols:
        for threshold in thresholds:
            tag = threshold_tag(threshold)
            out[f"greedy_clusters_{tag}"] = math.nan
            out[f"clusters_per_100_{tag}"] = math.nan
            out[f"cluster_entropy_{tag}"] = math.nan
            out[f"largest_cluster_fraction_{tag}"] = math.nan
            out[f"effective_clusters_{tag}"] = math.nan
            out[f"simpson_concentration_{tag}"] = math.nan
            out[f"singleton_fraction_{tag}"] = math.nan
        return out
    for threshold in thresholds:
        tag = threshold_tag(threshold)
        assignments = greedy_cluster_assignments(mols, threshold)
        dist = cluster_distribution_metrics(assignments)
        cluster_count = dist["cluster_count"]
        out[f"greedy_clusters_{tag}"] = cluster_count
        out[f"clusters_per_100_{tag}"] = 100.0 * cluster_count / len(mols)
        out[f"cluster_entropy_{tag}"] = dist["cluster_entropy"]
        out[f"largest_cluster_fraction_{tag}"] = dist["largest_cluster_fraction"]
        out[f"effective_clusters_{tag}"] = dist["effective_cluster_count"]
        out[f"simpson_concentration_{tag}"] = dist["simpson_concentration"]
        out[f"singleton_fraction_{tag}"] = dist["singleton_fraction"]
    return out


def pairwise_rmsd_stats(mols: Sequence[Chem.Mol], cap: int = PAIRWISE_SAMPLE_CAP) -> dict[str, float]:
    sampled = sample_mols(mols, cap=cap)
    values = []
    for mol_a, mol_b in itertools.combinations(sampled, 2):
        value = heavy_rmsd_same_topology(mol_a, mol_b)
        if math.isfinite(value):
            values.append(value)
    if not values:
        return {"pairwise_mean": math.nan, "pairwise_p90": math.nan}
    arr = np.asarray(values, dtype=float)
    return {
        "pairwise_mean": float(np.mean(arr)),
        "pairwise_p90": float(np.percentile(arr, 90)),
    }


def get_rotatable_torsions(mol: Chem.Mol) -> list[tuple[int, int, int, int]]:
    pattern = Chem.MolFromSmarts("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]")
    matches = mol.GetSubstructMatches(pattern)
    torsions: list[tuple[int, int, int, int]] = []
    seen: set[tuple[int, int]] = set()
    for j, k in matches:
        bond_key = tuple(sorted((j, k)))
        if bond_key in seen:
            continue
        seen.add(bond_key)
        atom_j = mol.GetAtomWithIdx(j)
        atom_k = mol.GetAtomWithIdx(k)
        neighbors_j = [
            atom.GetIdx()
            for atom in atom_j.GetNeighbors()
            if atom.GetIdx() != k and atom.GetAtomicNum() > 1
        ]
        neighbors_k = [
            atom.GetIdx()
            for atom in atom_k.GetNeighbors()
            if atom.GetIdx() != j and atom.GetAtomicNum() > 1
        ]
        if not neighbors_j:
            neighbors_j = [atom.GetIdx() for atom in atom_j.GetNeighbors() if atom.GetIdx() != k]
        if not neighbors_k:
            neighbors_k = [atom.GetIdx() for atom in atom_k.GetNeighbors() if atom.GetIdx() != j]
        if neighbors_j and neighbors_k:
            torsions.append((neighbors_j[0], j, k, neighbors_k[0]))
    return torsions


def mean_torsion_std_deg(mols: Sequence[Chem.Mol]) -> float:
    if not mols:
        return math.nan
    torsions = get_rotatable_torsions(mols[0])
    if not torsions:
        return math.nan
    profiles = []
    for mol in mols:
        conf = mol.GetConformer(0)
        profiles.append([rdMolTransforms.GetDihedralDeg(conf, *torsion) for torsion in torsions])
    return float(np.mean(np.std(np.asarray(profiles, dtype=float), axis=0)))


def forcefield_energy(mol: Chem.Mol, conf_id: int = 0, ff_variant: str = "MMFF94s") -> float:
    try:
        minimized = mol.GetProp("minimization_applied").strip().lower() == "true" if mol.HasProp("minimization_applied") else False
        if minimized and mol.HasProp("post_min_energy"):
            return float(mol.GetProp("post_min_energy"))
    except Exception:
        pass
    try:
        if not AllChem.MMFFHasAllMoleculeParams(mol):
            return math.nan
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=ff_variant)
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
        if ff is None:
            return math.nan
        return float(ff.CalcEnergy())
    except Exception:
        return math.nan


def energy_stats(mols: Sequence[Chem.Mol]) -> dict[str, float]:
    values = [forcefield_energy(mol) for mol in mols]
    finite = np.asarray([value for value in values if math.isfinite(value)], dtype=float)
    if finite.size == 0:
        return {
            "energy_count": 0.0,
            "energy_min": math.nan,
            "energy_max": math.nan,
            "energy_median": math.nan,
            "energy_std": math.nan,
        }
    return {
        "energy_count": float(finite.size),
        "energy_min": float(np.min(finite)),
        "energy_max": float(np.max(finite)),
        "energy_median": float(np.median(finite)),
        "energy_std": float(np.std(finite)),
    }


def safe_mean(values: Sequence[float]) -> float:
    finite = [float(v) for v in values if isinstance(v, (int, float)) and math.isfinite(v)]
    return float(statistics.mean(finite)) if finite else math.nan


def safe_median(values: Sequence[float]) -> float:
    finite = [float(v) for v in values if isinstance(v, (int, float)) and math.isfinite(v)]
    return float(statistics.median(finite)) if finite else math.nan


def safe_sum(values: Sequence[float]) -> float:
    finite = [float(v) for v in values if isinstance(v, (int, float)) and math.isfinite(v)]
    return float(sum(finite)) if finite else math.nan
