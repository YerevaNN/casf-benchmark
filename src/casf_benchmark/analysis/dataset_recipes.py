"""Fixed-budget conformer recipes; sampling probabilities and geometric summaries."""

from __future__ import annotations

import hashlib
import itertools
import math
from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln


@dataclass(frozen=True)
class Recipe:
    name: str
    weights: tuple[tuple[str, float], ...]


def recipe_grid(sources=("qwen", "loqi", "flowr")) -> list[Recipe]:
    sources = tuple(sources)
    if not sources or len(set(sources)) != len(sources) or not set(sources) <= {"qwen", "loqi", "flowr"}:
        raise ValueError("Recipe sources must be distinct supported learned generators")
    recipes = [Recipe(source, ((source, 1.0),)) for source in sources]
    for left, right in itertools.combinations(sources, 2):
        for percent in (75, 50, 25):
            recipes.append(Recipe(f"{left}{percent}_{right}{100-percent}",
                                  ((left, percent / 100), (right, 1 - percent / 100))))
    if len(sources) == 3:
        for favored in sources:
            recipes.append(Recipe("_".join(f"{source}{50 if source == favored else 25}" for source in sources),
                                  tuple((source, .5 if source == favored else .25) for source in sources)))
        recipes.append(Recipe("_".join(sources) + "_equal", tuple((source, 1 / 3) for source in sources)))
    recipes.extend(Recipe(source, ((source, 1.0),)) for source in ("rdkit_raw", "rdkit_minimized"))
    return recipes


def allocation_variants(recipe: Recipe, budget: int) -> list[dict[str, int]]:
    """Largest-remainder quotas; average all tied ways to assign a remainder."""
    if budget < 0 or not np.isclose(sum(weight for _, weight in recipe.weights), 1):
        raise ValueError("Invalid budget or weights")
    names = [name for name, _ in recipe.weights]
    exact = np.asarray([weight * budget for _, weight in recipe.weights])
    base = np.floor(exact + 1e-10).astype(int)
    remaining = int(budget - base.sum())
    if remaining == 0:
        return [dict(zip(names, map(int, base)))]
    remainders = exact - base
    boundary = np.sort(remainders)[-remaining]
    mandatory = set(np.flatnonzero(remainders > boundary + 1e-8))
    tied = np.flatnonzero(np.isclose(remainders, boundary, atol=1e-8, rtol=0))
    variants = []
    for extra in itertools.combinations(tied, remaining - len(mandatory)):
        counts = base.copy()
        counts[list(mandatory | set(extra))] += 1
        variants.append(dict(zip(names, map(int, counts))))
    return variants


def no_hit(population: int, hits: np.ndarray, draws: int) -> np.ndarray:
    """Vectorized hypergeometric no-hit probability, including empty pools."""
    hits = np.asarray(hits, dtype=int)
    if population < 0 or not 0 <= draws <= population or np.any((hits < 0) | (hits > population)):
        raise ValueError("Invalid population, hits, or draw count")
    result = np.ones(hits.shape, dtype=float)
    if draws == 0:
        return result
    impossible = population - hits < draws
    result[impossible] = 0
    selected = (hits > 0) & ~impossible
    misses = population - hits[selected]
    result[selected] = np.exp(gammaln(misses + 1) - gammaln(misses - draws + 1)
                              - gammaln(population + 1) + gammaln(population - draws + 1))
    return np.clip(result, 0, 1)


def recipe_statistics(recipe: Recipe, budget: int, populations: dict[str, int],
                      hits: dict[str, np.ndarray], valid_counts: dict[str, int]) -> dict:
    first = next(iter(hits.values()))
    coverage, support = [], []
    realized, valid = [], []
    allocations = []
    for requested in allocation_variants(recipe, budget):
        draws = {name: min(count, populations[name]) for name, count in requested.items()}
        miss = np.ones(np.asarray(first).shape)
        count = np.zeros_like(miss)
        valid_draws = 0.0
        for name, n in draws.items():
            population = populations[name]
            miss *= no_hit(population, hits[name], n)
            if population:
                count += n * np.asarray(hits[name]) / population
                valid_draws += n * valid_counts[name] / population
        coverage.append(1 - miss)
        support.append(count)
        realized.append(sum(draws.values()))
        valid.append(valid_draws)
        allocations.append(draws)
    return {"coverage": np.mean(coverage, axis=0), "support": np.mean(support, axis=0),
            "realized_draws": float(np.mean(realized)), "expected_valid": float(np.mean(valid)),
            "allocations": allocations}


def paired_gains_losses(left: Recipe, right: Recipe, budget: int,
                        populations: dict[str, int], hits: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Exact gains/losses with shared-source draws coupled as nested prefixes."""
    left_variants = allocation_variants(left, budget)
    right_variants = allocation_variants(right, budget)
    gains, losses = [], []
    count = math.lcm(len(left_variants), len(right_variants))
    for index in range(count):
        a = left_variants[index % len(left_variants)]
        b = right_variants[index % len(right_variants)]
        miss_a = np.ones_like(next(iter(hits.values())), dtype=float)
        miss_b = miss_a.copy()
        miss_both = miss_a.copy()
        for name in set(a) | set(b):
            na, nb = min(a.get(name, 0), populations[name]), min(b.get(name, 0), populations[name])
            miss_a *= no_hit(populations[name], hits[name], na)
            miss_b *= no_hit(populations[name], hits[name], nb)
            miss_both *= no_hit(populations[name], hits[name], max(na, nb))
        gains.append(np.maximum(0, miss_b - miss_both))
        losses.append(np.maximum(0, miss_a - miss_both))
    return np.mean(gains, axis=0), np.mean(losses, axis=0)


def stable_seed(seed: int, *parts: object) -> int:
    text = "|".join([str(seed), *map(str, parts)])
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "little")


def greedy_distance_clusters(distances: np.ndarray, indices: np.ndarray, radius: float) -> dict[str, float]:
    """The paper's nearest-fixed-representative greedy rule on a frozen matrix."""
    centers: list[int] = []
    counts: list[int] = []
    for index in indices:
        if centers:
            values = distances[index, centers]
            nearest = int(np.argmin(values))
            if values[nearest] < radius:
                counts[nearest] += 1
                continue
        centers.append(int(index))
        counts.append(1)
    if not counts:
        return {"clusters": 0., "largest_fraction": np.nan, "effective_clusters": 0.}
    probabilities = np.asarray(counts, dtype=float) / sum(counts)
    return {"clusters": float(len(counts)), "largest_fraction": float(probabilities.max()),
            "effective_clusters": float(np.exp(-np.sum(probabilities * np.log(probabilities))))}


def bootstrap_interval(values: np.ndarray, *, replicates: int, seed: int) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    means = values[rng.integers(0, len(values), size=(replicates, len(values)))].mean(axis=1)
    return tuple(map(float, np.quantile(means, [.025, .975])))
