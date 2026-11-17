import json

import pandas as pd
import pytest

from casf_benchmark.cli.report_dataset_recipes import paired_interval, read_completed, render_report


def test_paired_interval_correctly_reverses_a_comparison():
    rows = pd.DataFrame([{'left':'a','right':'b','difference':.2,'ci_low':-.1,'ci_high':.4}])
    assert paired_interval(rows,'a','b') == (.2,-.1,.4)
    assert paired_interval(rows,'b','a') == (-.2,-.4,.1)
    assert paired_interval(rows,'a','a') == (0,0,0)
    with pytest.raises(ValueError,match='Missing or ambiguous'):
        paired_interval(rows,'a','c')


def test_final_report_rejects_current_pool_or_partial_results(tmp_path):
    root = tmp_path/'plinder22'
    root.mkdir()
    config = {'output_dir':str(tmp_path),'panels':{'plinder22':{}}}
    for status in ['analysis_complete_on_current_pools','coverage_complete_diversity_partial']:
        (root/'run_manifest.json').write_text(json.dumps({'status':status}))
        with pytest.raises(RuntimeError,match='not ready'):
            read_completed(config)


def test_report_distinguishes_a_point_leader_from_a_joint_winner():
    import itertools
    from casf_benchmark.analysis.dataset_recipes import recipe_grid

    names = [r.name for r in recipe_grid()]
    table = pd.DataFrame([{'recipe':name,'coverage_0p5':.8-i*.01,'coverage_0p75':.9-i*.005,
                           'clusters_r0p5':100+i*10,'realized_draws':1000,
                           'expected_valid':950,'support':30} for i,name in enumerate(names)]).set_index('recipe')
    comparisons = pd.DataFrame([{'left':left,'right':right,'difference':.01,'ci_low':-.05,'ci_high':.06}
                                for left,right in itertools.combinations(names,2)])
    shapes = comparisons.assign(metric='clusters')
    sensitivity = pd.DataFrame([
        {'recipe':'qwen50_loqi25_flowr25','gate':'posebusters','radius':.5,'budget':1000,
         'scenario':scenario,'weighting':weighting,'coverage':.7}
        for scenario,weighting in [('equal_availability','queries'),('requested','reference_groups')]])
    row = {'table':table,'pairs':comparisons,'diversity_pairs':shapes,'coverage':sensitivity,
           'flowr':{'conformer_count':1000,'molecule_count':1,'minimum_count':1000,
                    'molecules':[{'attempts':1234,'shortfall':0,'elapsed_seconds':3600}]}}
    report = render_report({}, {'plinder22':row,'atc48':row}, {'accepted_additions':57,'attempts':57})
    assert 'No recipe dominates every alternative' in report
    assert 'not a demonstration' in report
    assert 'summed sampling time: 1.00 hours' in report
    assert 'not parallel elapsed wall time' in report
    assert '57 additions from 57 attempts' in report
    assert '0.75 A' in report
    assert 'exploratory paired 95% intervals' in report
    assert 600 < len(report.split()) < 1800


def test_seven_recipe_report_does_not_require_or_compare_flowr():
    import itertools
    from casf_benchmark.analysis.dataset_recipes import recipe_grid

    config = {'sources':['qwen', 'loqi'], 'scope_revision':{'excluded_sources':['flowr']}}
    names = [r.name for r in recipe_grid(config['sources'])]
    table = pd.DataFrame([{'recipe':name,'coverage_0p5':.8-i*.01,'coverage_0p75':.9-i*.005,
                           'clusters_r0p5':100+i*10,'realized_draws':1000,
                           'expected_valid':950,'support':30}
                          for i,name in enumerate(names)]).set_index('recipe')
    pairs = pd.DataFrame([{'left':left,'right':right,'difference':.01,'ci_low':-.05,'ci_high':.06}
                          for left,right in itertools.combinations(names,2)])
    sensitivity = pd.DataFrame([
        {'recipe':'qwen50_loqi50','gate':'posebusters','radius':.5,'budget':1000,
         'scenario':scenario,'weighting':weighting,'coverage':.7}
        for scenario,weighting in [('equal_availability','queries'),('requested','reference_groups')]])
    row = {'table':table, 'pairs':pairs, 'diversity_pairs':pairs.assign(metric='clusters'),
           'coverage':sensitivity}
    report = render_report(config, {'plinder22':row, 'atc48':row}, {'accepted_additions':57,'attempts':57})
    assert 'contains 7 recipes' in report
    assert 'All 7 rows' in report
    assert 'no further generation' in report
    assert 'not preregistered independent validation' in report
    assert 'qwen25_loqi75' in report
    assert 'final FlowR has' not in report
    assert 'qwen50_loqi25_flowr25' not in report


def test_completed_seven_recipe_report_reuses_external_qwen_completion(tmp_path):
    from casf_benchmark.analysis.dataset_recipes import recipe_grid
    from casf_benchmark.analysis.pharmacophore import sha256

    names = [r.name for r in recipe_grid(['qwen', 'loqi'])]
    original = tmp_path/'original_qwen.pickle'
    original.write_bytes(b'completed Qwen fixture')
    qwen_path = tmp_path/'v1'/'topup_manifest.json'
    qwen_path.parent.mkdir()
    qwen_path.write_text(json.dumps({'combined_cache_sha256':sha256(original)}))
    config = {'output_dir':str(tmp_path/'v2'), 'sources':['qwen','loqi'], '_sha256':'config',
              'qwen_completion_manifest':str(qwen_path), 'panels':{
                  p:{'methods':dict.fromkeys(['qwen','loqi','rdkit_raw','rdkit_minimized'])}
                  for p in ['plinder22','atc48']}}
    for panel in config['panels']:
        root = tmp_path/'v2'/panel
        root.mkdir(parents=True)
        table_path = root/'recipe_endpoint_table.csv'
        pd.DataFrame({'recipe':names}).to_csv(table_path,index=False)
        for name in ['coverage_summary','paired_recipe_differences','paired_diversity_differences']:
            pd.DataFrame({'fixture':[1]}).to_csv(root/f'{name}.csv',index=False)
        (root/'run_manifest.json').write_text(json.dumps({
            'status':'complete','config_sha256':'config','summary_sha256':sha256(table_path),
            'generation_pools':{name:{'path':str(original),'sha256':sha256(original)}
                                for name in config['panels'][panel]['methods']}}))
    results, qwen = read_completed(config)
    assert set(results) == {'plinder22','atc48'}
    assert all(row['flowr'] is None for row in results.values())
    assert qwen['combined_cache_sha256'] == sha256(original)
    manifest_path = tmp_path/'v2'/'plinder22'/'run_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    manifest['generation_pools']['flowr'] = {'path':str(original),'sha256':sha256(original)}
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(RuntimeError, match='configured active sources'):
        read_completed(config)


def test_two_source_submission_has_no_flowr_or_merge_dependency(tmp_path, monkeypatch, capsys):
    import importlib.util
    import sys
    from pathlib import Path
    import yaml
    from casf_benchmark.analysis.pharmacophore import sha256

    path = Path(__file__).resolve().parents[1]/'scripts/submit_dataset_recipe_analysis.py'
    spec = importlib.util.spec_from_file_location('recipe_submission', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    config_path = tmp_path/'config.yaml'
    config_path.write_text(yaml.safe_dump({'output_dir':str(tmp_path),'sources':['qwen','loqi'],
                                         'panels':{'plinder22':{},'atc48':{}}}))
    (tmp_path/'protocol.json').write_text(json.dumps({
        'config_sha256':sha256(config_path),
        'panels':{'plinder22':{'molecule_count':22},'atc48':{'molecule_count':48}}}))
    monkeypatch.setattr(sys, 'argv', [str(path),'--config',str(config_path)])
    module.main()
    output = capsys.readouterr().out
    assert 'generate_recipe_flowr' not in output
    assert '-merge' not in output
    assert '--dependency=afterok:' not in output.splitlines()[0]
    assert '--array=0-47%8' in output
    assert '--cpus-per-task=2 --gres=cpuonly:2' in output
    assert '--cpus-per-task=1 --gres=cpuonly:1' in output
    assert 'casf_benchmark.cli.report_dataset_recipes' in output
