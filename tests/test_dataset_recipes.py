import itertools

import numpy as np
import pytest

from casf_benchmark.analysis.dataset_recipes import (
    Recipe, allocation_variants, greedy_distance_clusters, no_hit,
    paired_gains_losses, recipe_grid, recipe_statistics,
)


def test_grid_and_balanced_remainder():
    recipes = recipe_grid()
    assert len(recipes) == len({r.name for r in recipes}) == 18
    for recipe in recipes:
        for budget in (100, 500, 1000):
            assert all(sum(v.values()) == budget for v in allocation_variants(recipe, budget))
    balanced = next(r for r in recipes if r.name == 'qwen_loqi_flowr_equal')
    variants = allocation_variants(balanced, 1000)
    assert len(variants) == 3
    assert {name: sum(v[name] for v in variants) for name, _ in balanced.weights} == dict.fromkeys(['qwen','loqi','flowr'], 1000)


def test_two_source_grid_contains_only_the_seven_retained_arms():
    recipes = recipe_grid(['qwen', 'loqi'])
    assert [r.name for r in recipes] == [
        'qwen', 'loqi', 'qwen75_loqi25', 'qwen50_loqi50', 'qwen25_loqi75',
        'rdkit_raw', 'rdkit_minimized']
    for recipe in recipes:
        assert 'flowr' not in dict(recipe.weights)
        for budget in (100, 500, 1000):
            variants = allocation_variants(recipe, budget)
            assert len(variants) == 1
            assert sum(variants[0].values()) == budget
    for sources in [[], ['qwen', 'qwen'], ['qwen', 'unknown']]:
        with pytest.raises(ValueError, match='distinct supported'):
            recipe_grid(sources)


def test_active_configuration_excludes_flowr_without_changing_frozen_v1():
    import json
    from pathlib import Path
    import yaml
    from casf_benchmark.analysis.pharmacophore import sha256

    directory = Path(__file__).resolve().parents[1]/'src/casf_benchmark/config'
    config = yaml.safe_load((directory/'dataset_recipe_experiment.yaml').read_text())
    assert config['sources'] == ['qwen', 'loqi']
    assert len(recipe_grid(config['sources'])) == 7
    assert 'flowr' not in config
    for panel in config['panels'].values():
        assert set(panel['methods']) == {'qwen', 'loqi', 'rdkit_raw', 'rdkit_minimized'}
        assert 'flowr' not in panel.get('completed_pools', {})
    previous = Path(config['scope_revision']['supersedes'])/'protocol.json'
    if previous.is_file():
        assert sha256(directory/'dataset_recipe_experiment_flowr_v1.yaml') == json.loads(previous.read_text())['config_sha256']


def test_probabilities_against_exhaustive_tiny_libraries():
    # Enumerate every possible subset, with a different query for each source.
    matrices = {'qwen': np.array([[1,0,0,0],[0,0,0,0]], bool),
                'loqi': np.array([[0,0,0,0],[1,1,0,0]], bool),
                'flowr': np.array([[1,0,0,0],[0,1,0,0]], bool)}
    recipe = Recipe('test', (('qwen', .5), ('loqi', .25), ('flowr', .25)))
    populations = dict.fromkeys(matrices, 4)
    hits = {k: m.sum(axis=1) for k,m in matrices.items()}
    result = recipe_statistics(recipe, 4, populations, hits, populations)
    observations = []
    supports = []
    for a,b,c in itertools.product(itertools.combinations(range(4), 2), range(4), range(4)):
        count = matrices['qwen'][:,a].sum(axis=1) + matrices['loqi'][:,b] + matrices['flowr'][:,c]
        observations.append(count > 0)
        supports.append(count)
    np.testing.assert_allclose(result['coverage'], np.mean(observations, axis=0))
    np.testing.assert_allclose(result['support'], np.mean(supports, axis=0))


def test_pb_failures_consume_budget_and_empty_source_is_not_reallocated():
    recipe = Recipe('mix', (('qwen', .5), ('flowr', .5)))
    result = recipe_statistics(recipe, 4, {'qwen':4, 'flowr':0},
                               {'qwen':np.array([1]), 'flowr':np.array([0])},
                               {'qwen':1, 'flowr':0})
    assert result['coverage'][0] == pytest.approx(.5)
    assert result['expected_valid'] == .5
    assert result['realized_draws'] == 2
    np.testing.assert_array_equal(no_hit(0, np.array([0,0]), 0), [1,1])


def test_coupled_gains_minus_losses_is_coverage_difference():
    recipes = recipe_grid()[:16]
    populations = {'qwen':6, 'loqi':5, 'flowr':4}
    hits = {'qwen':np.array([0,1,5]), 'loqi':np.array([1,3,2]), 'flowr':np.array([2,1,0])}
    for left in recipes:
        for right in recipes:
            a = recipe_statistics(left, 5, populations, hits, populations)['coverage']
            b = recipe_statistics(right, 5, populations, hits, populations)['coverage']
            gained, lost = paired_gains_losses(left, right, 5, populations, hits)
            np.testing.assert_allclose(gained-lost, a-b, atol=1e-12)
    gain, loss = paired_gains_losses(recipes[0], recipes[0], 5, populations, hits)
    np.testing.assert_array_equal(gain, 0)
    np.testing.assert_array_equal(loss, 0)


def test_greedy_clusters_use_nearest_representative_and_strict_radius():
    distances = np.array([[0,.5,.9,.4],[.5,0,.4,.2],[.9,.4,0,.8],[.4,.2,.8,0]])
    result = greedy_distance_clusters(distances, np.arange(4), .5)
    assert result['clusters'] == 2
    assert result['largest_fraction'] == .75
    assert greedy_distance_clusters(distances, np.array([],int), .5)['clusters'] == 0


def test_distance_cache_is_symmetry_aware_and_does_not_move_sources(tmp_path):
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign
    from casf_benchmark.analysis.pharmacophore import molecule_coordinate_hash
    from casf_benchmark.cli.analyze_dataset_recipes import distance_matrix

    molecule = Chem.AddHs(Chem.MolFromSmiles('CCCO'))
    AllChem.EmbedMultipleConfs(molecule, numConfs=3, randomSeed=42)
    mols = []
    for conformer in molecule.GetConformers():
        mol = Chem.Mol(molecule)
        mol.RemoveAllConformers()
        mol.AddConformer(Chem.Conformer(conformer), assignId=True)
        mols.append(Chem.RemoveHs(mol))
    # Different atom order and a rigid translation must not create diversity.
    reordered = Chem.RenumberAtoms(mols[0], list(reversed(range(mols[0].GetNumAtoms()))))
    reordered.GetConformer().SetPositions(reordered.GetConformer().GetPositions() + [4, -2, 1])
    mols.append(reordered)
    hashes = [molecule_coordinate_hash(mol) for mol in mols]
    path = tmp_path/'distances.npz'
    actual = distance_matrix(mols, path, 'signature', 1)
    expected = np.zeros(actual.shape)
    for i, j in itertools.combinations(range(len(mols)), 2):
        expected[i,j] = expected[j,i] = rdMolAlign.GetBestRMS(Chem.Mol(mols[i]), Chem.Mol(mols[j]))
    np.testing.assert_allclose(actual, expected, atol=1e-6)
    assert actual[0,-1] < 1e-6
    assert [molecule_coordinate_hash(mol) for mol in mols] == hashes
    np.testing.assert_array_equal(distance_matrix(mols, path, 'signature', 1), actual)
    assert distance_matrix([], path, 'new_signature', 1).shape == (0,0)


def test_prepare_appends_only_new_candidates_and_rejects_changed_prefix(tmp_path, monkeypatch):
    import json
    import pickle
    from types import SimpleNamespace
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from casf_benchmark.cli import analyze_dataset_recipes as cli
    from casf_benchmark.generation import conformer_sets

    smiles = 'CCCO'
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    mol = Chem.RemoveHs(mol)
    pool = tmp_path/'pool.pickle'
    config = {'panels': {'test': {'methods': {'qwen':'original'}}}, 'radii':[.5,.75,1.], 'workers':1}
    old, root = tmp_path/'old', tmp_path/'new'
    old.mkdir()
    (old/'query_manifest.csv').write_text('fixture')
    monkeypatch.setattr(cli, 'pool_path', lambda *args: pool)
    calls = []
    def check(mols, reference, **kwargs):
        calls.append(len(mols))
        return SimpleNamespace(passes=[True]*len(mols), check_rows=[{'geometry':True} for _ in mols])
    monkeypatch.setattr(conformer_sets, 'posebusters_geometry_passes', check)
    group = pd.DataFrame(columns=['query_id'])
    for count in [2,3,3]:
        with pool.open('wb') as handle:
            pickle.dump({smiles:[Chem.Mol(mol) for _ in range(count)]}, handle)
        cli.prepare_source(config, 'test', 'qwen', smiles, group, {'methods':{}}, old, root)
    assert calls == [2,1]
    path = root/'prepared'/cli.stable_smiles_id(smiles)/'qwen.npz'
    with np.load(path, allow_pickle=True) as stored:
        assert stored['matches_r0p5'].shape == (0,3)
        assert stored['passes'].tolist() == [True]*3
    assert json.loads(path.with_suffix('.json').read_text())['source_sha256'] == cli.sha256(pool)
    changed = Chem.Mol(mol)
    changed.GetConformer().SetPositions(changed.GetConformer().GetPositions() + [1,0,0])
    with pool.open('wb') as handle:
        pickle.dump({smiles:[changed, mol, mol]}, handle)
    with pytest.raises(ValueError, match='prefix/query changed'):
        cli.prepare_source(config, 'test', 'qwen', smiles, group, {'methods':{}}, old, root)


def test_freeze_is_idempotent_and_refuses_changed_inputs_before_writing(tmp_path, monkeypatch):
    import pandas as pd
    from casf_benchmark.cli import analyze_dataset_recipes as cli

    config_path = tmp_path/'config.yaml'
    config_path.write_text('fixture')
    source_path = tmp_path/'source.yaml'
    source_path.write_text('source fixture')
    pool_path = tmp_path/'pool.pickle'
    pool_path.write_bytes(b'original pool')
    old, root = tmp_path/'old', tmp_path/'output'
    old.mkdir()
    (old/'query_manifest.csv').write_text('queries')
    (old/'query_freeze_manifest.json').write_text('frozen')
    source = {'_path':str(source_path), 'methods':{'qwen_id':{'path':str(pool_path)}}}
    manifest = pd.DataFrame({'expected_smiles':['CCO'], 'eligible':[True]})
    config = {'output_dir':str(root), 'panels':{'panel':{'methods':{'qwen':'qwen_id'}}},
              'sources':['qwen','loqi','flowr']}
    monkeypatch.setattr(cli, 'context', lambda *args: (source, old, manifest, root))
    cli.freeze(config, config_path)
    cli.freeze(config, config_path)
    protocol = (root/'protocol.json').read_bytes()
    inventory = (root/'immutable_source_inventory.csv').read_bytes()
    pool_path.write_bytes(b'changed source')
    with pytest.raises(ValueError, match='inventory changed'):
        cli.freeze(config, config_path)
    assert (root/'protocol.json').read_bytes() == protocol
    assert (root/'immutable_source_inventory.csv').read_bytes() == inventory
    config_path.write_text('changed config')
    with pytest.raises(ValueError, match='different configuration'):
        cli.freeze(config, config_path)


def test_reference_diagnostics_distinguish_unscreened_from_failed(tmp_path, monkeypatch):
    import pandas as pd
    from casf_benchmark.cli import analyze_dataset_recipes as cli

    smiles = 'CCO'
    manifest = pd.DataFrame({'expected_smiles':[smiles], 'name':['ethanol'], 'eligible':[True]})
    root = tmp_path/'panel'
    prepared = root/'prepared'/cli.stable_smiles_id(smiles)
    prepared.mkdir(parents=True)
    refs = pd.DataFrame({
        'query_id':['a','b','c'], 'reference_group':[0,1,2],
        'reference_geometry_status':['passed_archived_checks', 'not_reference_posebusters_screened', 'bond_lengths'],
        'minimum_ligand_occupancy':[1.,1.,1.]})
    refs.to_csv(prepared/'references.csv',index=False)
    refs[['query_id']].assign(smiles=smiles,recipe='qwen',coverage=.5,support=1.,feature_count=4).to_csv(
        root/'per_query_support.csv',index=False)
    coverage = pd.DataFrame([{'smiles':smiles,'recipe':'qwen','gate':'posebusters','budget':1000,
                             'scenario':'requested','weighting':'queries','radius':.5,'coverage':.5}])
    monkeypatch.setattr(cli,'context',lambda *args: ({},tmp_path/'old',manifest,root))
    cli.diagnostic_summaries({'bootstrap_replicates':10,'seed':1},'panel',coverage)
    diagnostic = pd.read_csv(root/'native_geometry_per_molecule.csv')
    assert set(diagnostic.native_geometry) == {
        'passed_archived_checks','not_screened','flagged_archived_checks'}
