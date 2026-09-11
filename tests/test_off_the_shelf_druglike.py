from __future__ import annotations

import importlib.util
import pickle
from collections import OrderedDict
from pathlib import Path
from types import ModuleType

import numpy as np

from casf_benchmark.paths import REPO_ROOT


def load_script(name: str) -> ModuleType:
    path = REPO_ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_merge_parts_preserves_druglike_order(tmp_path: Path, monkeypatch) -> None:
    inference = load_script("run_off_the_shelf_druglike")
    monkeypatch.setattr(inference, "normalize_conformers", lambda mols: mols)
    druglike = tmp_path / "druglike.pickle"
    with druglike.open("wb") as handle:
        pickle.dump(OrderedDict([("CC", {}), ("CO", {})]), handle)

    for index, (smiles, conformers) in enumerate((("CC", ["a"]), ("CO", ["b", "c"]))):
        path = inference.part_path(tmp_path / "run", index)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("wb") as handle:
            pickle.dump(OrderedDict([(smiles, conformers)]), handle)

    destination = inference.merge_parts(druglike, tmp_path / "run")
    with destination.open("rb") as handle:
        merged = pickle.load(handle)
    assert list(merged) == ["CC", "CO"]
    assert merged["CO"] == ["b", "c"]


def test_covmat_metrics_have_expected_recall_and_precision() -> None:
    covmat = load_script("eval_druglike_covmat")
    metrics = covmat.per_molecule_metrics(
        np.asarray([[0.0, 1.0], [2.0, 3.0]], dtype=float),
        threshold=0.75,
        dmax=3.0,
    )
    assert metrics["cov_r_075"] == 0.5
    assert metrics["cov_p_075"] == 0.5
    assert metrics["mat_r"] == 1.0
    assert metrics["mat_p"] == 0.5
    assert metrics["cmat_r"] == 1.0
    assert metrics["cmat_p"] == 0.5
