from __future__ import annotations

import itertools
import json
import math
import subprocess

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from casf_benchmark.analysis.pharmacophore import (
    PharmitBackend,
    distinct_centers,
    exact_feature_correspondence_match,
    exact_feature_correspondence_transform,
    exact_feature_membership,
    export_candidate_library,
    hit_probability,
    hypergeometric_no_hit_probability,
    maximum_triangle_area,
    mixture_hit_probability,
    molecule_coordinate_hash,
    proper_kabsch,
    transform_query,
)
from casf_benchmark.cli.analyze_pharmacophores import (
    _molecule_weighted_agreement,
    conda_package_inventory,
)


def test_conda_package_inventory_reads_metadata_without_environment_manager(tmp_path) -> None:
    metadata = tmp_path / "conda-meta"
    metadata.mkdir()
    (metadata / "zlib-1.3.1-0.json").write_text(
        json.dumps({"name": "zlib", "version": "1.3.1", "build": "0", "build_number": 0})
    )
    (metadata / "python-3.12-1.json").write_text(
        json.dumps(
            {
                "name": "python",
                "version": "3.12",
                "build": "1",
                "build_number": 1,
                "channel": "conda-forge",
            }
        )
    )

    packages = conda_package_inventory(tmp_path)

    assert [record["name"] for record in packages] == ["python", "zlib"]
    assert packages[0]["channel"] == "conda-forge"
    assert packages[1]["url"] is None


def test_pharmit_null_points_are_normalized(tmp_path, monkeypatch) -> None:
    backend = PharmitBackend(tmp_path / "pharmit")

    def fake_run(self, arguments, *, cwd=None):
        output = tmp_path / "query.json"
        output.write_text(json.dumps({"points": None}))
        return subprocess.CompletedProcess(arguments, 0, "", "")

    monkeypatch.setattr(PharmitBackend, "_run", fake_run)
    query = backend.extract_interaction_query(
        tmp_path / "ligand.sdf",
        tmp_path / "receptor.pdb",
        tmp_path / "query.json",
        0.5,
    )
    assert query["points"] == []


def test_pharmit_database_parent_is_created(tmp_path, monkeypatch) -> None:
    backend = PharmitBackend(tmp_path / "pharmit")
    database = tmp_path / "nested" / "method" / "database"

    def fake_run(self, arguments, *, cwd=None):
        assert database.parent.is_dir()
        return subprocess.CompletedProcess(arguments, 0, "created", "")

    monkeypatch.setattr(PharmitBackend, "_run", fake_run)
    assert backend.create_database(tmp_path / "library.sdf", database) == "created"


def test_pharmit_search_disables_lossy_prefilters(tmp_path, monkeypatch) -> None:
    backend = PharmitBackend(
        tmp_path / "pharmit",
        skip_triplet_fingerprints=True,
        skip_geometric_pruning=True,
    )
    captured = []

    def fake_run(self, arguments, *, cwd=None):
        captured.extend(map(str, arguments))
        return subprocess.CompletedProcess(arguments, 0, "searched", "")

    monkeypatch.setattr(PharmitBackend, "_run", fake_run)
    backend.search_many(
        database_path=tmp_path / "database",
        query_paths=[tmp_path / "query.json"],
        output_paths=[tmp_path / "matches.sdf"],
        max_hits=11,
        workers=2,
    )
    assert "--skip-fingers" in captured
    assert "--skip-geometric-pruning" in captured


def test_hypergeometric_probability_matches_enumeration() -> None:
    population, hits, draws = 8, 3, 4
    hit_indices = set(range(hits))
    subsets = list(itertools.combinations(range(population), draws))
    expected_no_hit = sum(not hit_indices.intersection(subset) for subset in subsets) / len(subsets)
    observed = hypergeometric_no_hit_probability(population, hits, draws)
    assert observed == pytest.approx(expected_no_hit)
    assert hit_probability(population, hits, draws) == pytest.approx(1 - expected_no_hit)


def test_all_hit_no_hit_and_complete_pool_endpoints() -> None:
    assert hit_probability(1000, 0, 1000) == 0.0
    assert hit_probability(1000, 1000, 1) == 1.0
    assert hit_probability(1000, 1, 1000) == 1.0
    probabilities = [hit_probability(1000, 17, budget) for budget in [10, 20, 50, 1000]]
    assert probabilities == sorted(probabilities)


def test_equal_mixture_probability_matches_independent_formula() -> None:
    observed = mixture_hit_probability(100, 7, 100, 11, 20)
    qwen_no_hit = hypergeometric_no_hit_probability(100, 7, 10)
    loqi_no_hit = hypergeometric_no_hit_probability(100, 11, 10)
    assert observed == pytest.approx(1 - qwen_no_hit * loqi_no_hit)
    with pytest.raises(ValueError, match="even"):
        mixture_hit_probability(100, 7, 100, 11, 21)


def test_equal_mixture_probability_matches_monte_carlo() -> None:
    rng = np.random.default_rng(1719)
    trials = 40_000
    qwen_hits = np.arange(5)
    loqi_hits = np.arange(8)
    observed_hits = 0
    for _ in range(trials):
        qwen_draw = rng.choice(40, size=10, replace=False)
        loqi_draw = rng.choice(40, size=10, replace=False)
        observed_hits += bool(
            np.intersect1d(qwen_draw, qwen_hits).size
            or np.intersect1d(loqi_draw, loqi_hits).size
        )
    exact = mixture_hit_probability(40, 5, 40, 8, 20)
    assert observed_hits / trials == pytest.approx(exact, abs=0.006)


def test_hit_probability_caps_requested_budget_at_supplied_pool() -> None:
    assert hit_probability(121, 1, 1000) == 1.0
    assert hit_probability(121, 0, 1000) == 0.0


def test_agreement_summary_weights_molecules_not_candidate_pairs() -> None:
    agreement = pd.DataFrame(
        [
            {
                "method": "m",
                "smiles": "A",
                "query_id": "a1",
                "candidate_count": 100,
                "both_hit_count": 100,
                "neither_count": 0,
            },
            {
                "method": "m",
                "smiles": "A",
                "query_id": "a2",
                "candidate_count": 100,
                "both_hit_count": 0,
                "neither_count": 100,
            },
            {
                "method": "m",
                "smiles": "B",
                "query_id": "b1",
                "candidate_count": 10,
                "both_hit_count": 0,
                "neither_count": 10,
            },
        ]
    )
    molecule, summary = _molecule_weighted_agreement(
        agreement,
        ["both_hit_count", "neither_count"],
        replicates=200,
        seed=7,
    )
    assert molecule.set_index("smiles").loc["A", "both_hit_count"] == 0.5
    observed = summary.set_index("category").loc["both_hit_count", "fraction"]
    assert observed == pytest.approx(0.25)


def test_exact_feature_correspondence_is_typed_permuted_and_monotone() -> None:
    query = np.asarray(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]
    )
    query_types = ["A", "B", "B", "A"]
    rotation = np.asarray([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    transformed = query @ rotation + np.asarray([3.0, 4.0, 5.0])
    permutation = [3, 2, 1, 0]
    candidate = transformed[permutation]
    candidate_types = [query_types[index] for index in permutation]
    assert exact_feature_correspondence_match(
        query, query_types, candidate, candidate_types, 0.01
    )
    transform = exact_feature_correspondence_transform(
        query, query_types, candidate, candidate_types, 0.01
    )
    assert transform is not None
    fitted = candidate[transform["candidate_indices"]] @ transform["rotation"]
    fitted += transform["translation"]
    assert np.allclose(fitted, query[transform["query_indices"]])
    perturbed = candidate.copy()
    perturbed[0, 0] += 0.2
    assert not exact_feature_correspondence_match(
        query, query_types, perturbed, candidate_types, 0.05
    )
    assert exact_feature_correspondence_match(
        query, query_types, perturbed, candidate_types, 0.2
    )


def test_exact_feature_membership_accepts_ragged_candidate_features() -> None:
    points = [
        {"name": "A", "x": 0.0, "y": 0.0, "z": 0.0},
        {"name": "B", "x": 1.0, "y": 0.0, "z": 0.0},
        {"name": "B", "x": 0.0, "y": 1.0, "z": 0.0},
        {"name": "A", "x": 0.0, "y": 0.0, "z": 1.0},
    ]
    first = np.asarray([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    second = np.vstack([first, [9, 9, 9]])
    coordinates = np.empty(2, dtype=object)
    coordinates[:] = [first, second]
    types = [("A", "B", "B", "A"), ("A", "B", "B", "A", "B")]
    assert exact_feature_membership(points, coordinates, types, 0.01).tolist() == [
        True,
        True,
    ]


def test_proper_kabsch_recovers_rotation_without_reflection() -> None:
    source = np.asarray(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    )
    angle = 0.41
    rotation = np.asarray(
        [
            [math.cos(angle), -math.sin(angle), 0.0],
            [math.sin(angle), math.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    target = source @ rotation + np.asarray([4.0, -2.0, 7.0])
    fitted_rotation, translation, rmsd = proper_kabsch(source, target)
    assert np.linalg.det(fitted_rotation) == pytest.approx(1.0)
    assert rmsd < 1e-10
    assert np.allclose(source @ fitted_rotation + translation, target)


def test_query_transform_rejects_reflection_and_preserves_distances() -> None:
    query = {
        "points": [
            {"name": "Aromatic", "x": 0.0, "y": 0.0, "z": 0.0},
            {"name": "Hydrophobic", "x": 2.0, "y": 0.0, "z": 0.0},
        ]
    }
    rotation = np.asarray([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    moved = transform_query(query, rotation, np.asarray([1.0, 2.0, 3.0]))
    before = np.linalg.norm(
        np.asarray([query["points"][0][axis] for axis in "xyz"])
        - np.asarray([query["points"][1][axis] for axis in "xyz"])
    )
    after = np.linalg.norm(
        np.asarray([moved["points"][0][axis] for axis in "xyz"])
        - np.asarray([moved["points"][1][axis] for axis in "xyz"])
    )
    assert after == pytest.approx(before)
    with pytest.raises(ValueError, match="proper rotation"):
        transform_query(query, np.diag([-1.0, 1.0, 1.0]), np.zeros(3))


def test_distinct_centers_and_noncollinearity() -> None:
    points = [
        {"x": 0.0, "y": 0.0, "z": 0.0},
        {"x": 0.05, "y": 0.0, "z": 0.0},
        {"x": 1.0, "y": 0.0, "z": 0.0},
        {"x": 0.0, "y": 1.0, "z": 0.0},
        {"x": 0.0, "y": 0.0, "z": 1.0},
    ]
    centers = distinct_centers(points, 0.1)
    assert len(centers) == 4
    assert maximum_triangle_area(centers) == pytest.approx(math.sqrt(3) / 2)


def test_candidate_export_preserves_coordinates_and_order(tmp_path) -> None:
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert AllChem.EmbedMolecule(mol, randomSeed=5) == 0
    source_hash = molecule_coordinate_hash(mol)
    records = export_candidate_library(
        method="test_method",
        smiles="CCO",
        mols=[mol, Chem.Mol(mol)],
        sdf_path=tmp_path / "library.sdf",
    )
    assert [record["candidate_index"] for record in records] == [0, 1]
    assert all(record["exported"] for record in records)
    assert molecule_coordinate_hash(mol) == source_hash
    supplier = [m for m in Chem.SDMolSupplier(str(tmp_path / "library.sdf"), removeHs=False) if m]
    assert [m.GetProp("_Name") for m in supplier] == [
        records[0]["candidate_id"],
        records[1]["candidate_id"],
    ]
