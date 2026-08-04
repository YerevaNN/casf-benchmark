"""Shared CASF stratum bin definitions."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

import pandas as pd

ROTATABLE_BOND_EDGES = [0, 5, 10, 15, math.inf]
ROTATABLE_BOND_LABELS = ["0-4", "5-9", "10-14", "15+"]
HEAVY_ATOM_EDGES = [0, 20, 30, 40, math.inf]
HEAVY_ATOM_LABELS = ["<20", "20-29", "30-39", "40+"]


@dataclass(frozen=True)
class StratumDefinition:
    column: str
    edges: list[float]
    labels: list[str]


STRATA = {
    "rotatable_bonds": StratumDefinition("rotatable_bonds", ROTATABLE_BOND_EDGES, ROTATABLE_BOND_LABELS),
    "heavy_atoms": StratumDefinition("heavy_atoms", HEAVY_ATOM_EDGES, HEAVY_ATOM_LABELS),
}


def assign_stratum(values: Iterable[object], breakdown: str) -> pd.Series:
    if breakdown not in STRATA:
        raise ValueError(f"Unknown CASF breakdown {breakdown!r}")
    spec = STRATA[breakdown]
    if isinstance(values, pd.Series):
        numeric = pd.to_numeric(values, errors="coerce")
    else:
        numeric = pd.to_numeric(pd.Series(list(values)), errors="coerce")
    return pd.cut(numeric, bins=spec.edges, labels=spec.labels, right=False, include_lowest=True)
