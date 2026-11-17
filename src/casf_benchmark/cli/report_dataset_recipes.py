"""Produce a medium-length factual report after both recipe panels finish."""

from __future__ import annotations

import argparse
import json
import os
import io
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml

from casf_benchmark.analysis.dataset_recipes import recipe_grid
from casf_benchmark.analysis.pharmacophore import atomic_json_dump, sha256


def paired_interval(frame, left, right, metric=None):
    if left == right:
        return 0., 0., 0.
    if metric is not None:
        frame = frame[frame.metric == metric]
    rows = frame[((frame.left == left) & (frame.right == right)) |
                 ((frame.left == right) & (frame.right == left))]
    if len(rows) != 1:
        raise ValueError(f'Missing or ambiguous paired comparison: {left}, {right}, {metric}')
    row = rows.iloc[0]
    sign = 1 if row.left == left else -1
    low, high = sorted([sign * row.ci_low, sign * row.ci_high])
    return sign * row.difference, low, high


def read_completed(config):
    root = Path(config['output_dir'])
    results = {}
    expected_recipes = {r.name for r in recipe_grid(config.get('sources', ('qwen', 'loqi', 'flowr')))}
    for panel in config['panels']:
        path = root / panel
        manifest_path = path/'run_manifest.json'
        if not manifest_path.is_file():
            raise RuntimeError(f'{panel}: final analysis is not ready')
        manifest = json.loads(manifest_path.read_text())
        if (manifest['status'] != 'complete' or manifest.get('source_completions_pending') or
                manifest.get('diversity_molecules_missing')):
            raise RuntimeError(f'{panel}: completed-pool coverage and diversity are not ready')
        if manifest['config_sha256'] != config['_sha256']:
            raise RuntimeError(f'{panel}: configuration changed')
        table_path = path/'recipe_endpoint_table.csv'
        if sha256(table_path) != manifest['summary_sha256']:
            raise RuntimeError(f'{panel}: summary checksum changed')
        table = pd.read_csv(table_path).set_index('recipe')
        if not table.index.is_unique or set(table.index) != expected_recipes:
            raise RuntimeError(f'{panel}: the complete {len(expected_recipes)}-recipe table is required')
        active_sources = set(config['panels'][panel].get('methods', {}))
        if active_sources and set(manifest['generation_pools']) != active_sources:
            raise RuntimeError(f'{panel}: parent pools do not match the configured active sources')
        for pool in manifest['generation_pools'].values():
            if sha256(Path(pool['path'])) != pool['sha256']:
                raise RuntimeError(f'{panel}: parent pool changed since analysis')
        flowr = None
        if 'flowr' in manifest['generation_pools']:
            flowr = json.loads((root/'generation/flowr'/f'{panel}_manifest.json').read_text())
            if (flowr['status'] not in ['complete', 'complete_with_reported_shortfalls'] or
                    flowr['pool_sha256'] != manifest['generation_pools']['flowr']['sha256']):
                raise RuntimeError(f'{panel}: analysis does not use the final FlowR pool')
        results[panel] = {'manifest':manifest, 'table':table, 'flowr':flowr,
            'coverage':pd.read_csv(path/'coverage_summary.csv'),
            'pairs':pd.read_csv(path/'paired_recipe_differences.csv'),
            'diversity_pairs':pd.read_csv(path/'paired_diversity_differences.csv')}
    qwen_path = Path(config.get('qwen_completion_manifest', root/'generation/qwen_atc48/topup_manifest.json'))
    qwen = json.loads(qwen_path.read_text())
    if qwen['combined_cache_sha256'] != results['atc48']['manifest']['generation_pools']['qwen']['sha256']:
        raise RuntimeError('ATC48 analysis does not use the completed Qwen pool')
    runtime_path = qwen_path.parent/'runtime_sessions.json'
    qwen['runtime_sessions'] = json.loads(runtime_path.read_text()) if runtime_path.is_file() else []
    return results, qwen


def collect_job_runtime(root):
    ledger = json.loads((root/'analysis_jobs.json').read_text())
    jobs = set(ledger['generation_jobs'].values()) | set(ledger.get('generation_resume_jobs', {}).values())
    jobs.update(row['id'] for row in ledger['jobs'].values())
    try:
        result = subprocess.run(['sacct','-X','--parsable2','-j',','.join(sorted(jobs)),
            '--format=JobID,JobName,State,ExitCode,ElapsedRaw,AllocCPUS,Start,End'],
            check=True, capture_output=True, text=True, timeout=30)
        frame = pd.read_csv(io.StringIO(result.stdout), sep='|')
    except (OSError, subprocess.SubprocessError, pd.errors.EmptyDataError):
        return None
    if not len(frame):
        return None
    path = root/'job_runtime_audit.csv'
    temporary = path.with_name(f'.{path.name}.{os.getpid()}.tmp')
    frame.to_csv(temporary, index=False)
    os.replace(temporary, path)
    return path


def render_report(config, results, qwen, runtime_audit=None):
    sources = config.get('sources', ('qwen', 'loqi', 'flowr'))
    recipes = recipe_grid(sources)
    has_flowr = 'flowr' in sources
    grid_description = (
        'Qwen, LoQI and FlowR separately; every learned pair at 75/25, 50/50 and 25/75; '
        'three 50/25/25 permutations; equal thirds; and RDKit raw and minimized controls.' if has_flowr else
        'Qwen and LoQI separately; Qwen/LoQI at 75/25, 50/50 and 25/75; '
        'and RDKit raw and minimized controls.')
    sections = ['DATASET RECIPE EXPERIMENT: COMPLETED RESULTS', '',
        'Purpose and scope',
        'This experiment asks how to construct a protein-independent conformer training dataset that '
        'contains many different shapes while retaining experimentally observed pharmacophore geometries. '
        'It is a dataset-construction test, not a demonstration that a downstream model trained on the '
        'dataset will perform better. No model was retrained and the manuscript was not edited.', '',
        'What was compared',
        f'The frozen grid contains {len(recipes)} recipes: {grid_description} '
        'Torsional diffusion is excluded. The main allowance is 1,000 draws before PoseBusters '
        '(PB), with 100 and 500 as smaller-budget checks. Failed PB checks consume their draw rather than '
        'being replaced. These are stored-library comparisons, not equal raw-generation-attempt '
        'or equal-compute comparisons. No new geometry repair or stereo-based replenishment is '
        'introduced for Qwen or LoQI.', '',
        'Molecules, queries and matching',
        'PLINDER22 uses 932 frozen queries from 22 molecules. ATC48 uses 311 queries from 44 eligible '
        'molecules; its four molecules without eligible queries appear in supporting validity/diversity '
        'results, not the joint endpoint. ATC48 contains all eligible ATC/DrugBank-annotated identities '
        'in the historical PLINDER candidate table after excluding the original PLINDER-23 identities '
        'and unsupported elements. Up to eight distinct PDB references per molecule were selected '
        'by a seeded ordering without examining generator results; it is not a random sample of all '
        'PLINDER molecules. A query describes the full selected protein-contact feature '
        'pattern in an experimental ligand pose. Each candidate is tested against references of the '
        'same molecule using rigid fits, with no internal bond rotation during matching. Generation '
        'uses SMILES only; protein information is restricted to evaluation. Primary hits must pass '
        'ligand-only PB checks and a 0.5-A feature-center tolerance. Tolerances of 0.75 and 1.0 A are '
        'sensitivities; this is not docking, affinity prediction or cross-molecule virtual screening.', '',
        'How the numbers were obtained',
        'Coverage is the exact expected fraction of queries recovered by uniformly sampling the '
        'specified source quotas without replacement. Queries are averaged within each molecule, '
        'then molecules equally. Diversity is measured on actual combined sampled libraries, using '
        'symmetry-aware heavy-atom RMSD clusters at 0.5 and 1.0 A. Twenty reproducible draws per '
        'allocation are averaged; tied rounding allocations receive equal weight. Paired uncertainty '
        'resamples molecules 2,000 times, not individual queries or repeated libraries. Reference-group '
        'weighting, raw hits, equal source availability, flexibility and training overlap are diagnostics. '
        f'Intervals remain exploratory and unadjusted for searching {len(recipes)} recipes.', '',
        'What happened during execution',
        f'ATC48 Qwen initially lacked {qwen["accepted_additions"]} accepted conformers. The top-up retained '
        f'{qwen["accepted_additions"]} additions from {qwen["attempts"]} attempts, preserving every original '
        'coordinate and reaching 1,000 per molecule. This replaced missing graph/parsing-accepted '
        'records, not PB or coordinate-stereochemistry failures. The PLINDER Qwen pilot pool '
        'similarly contains 53 additions; its original paper benchmark pool remains unchanged. '
        'All other active source libraries are reused without new generation or coordinate optimization.']
    if config.get('scope_revision'):
        sections.append('Scope revision: FlowR was removed at the user\'s request after an audit '
                        'identified custom relaxation, mirror correction and stereo-rejection/refill '
                        'in its benchmark adapter. Its running generation and dependent jobs were '
                        'cancelled, while checkpoints and the original 18-recipe exploration were '
                        'preserved separately. The revised grid requires no further generation. '
                        'This revision followed inspection of preliminary results; it is not '
                        'preregistered independent validation. Molecules, queries, tolerances and '
                        'the no-PB-refill rule were not changed.')
    sessions = qwen.get('runtime_sessions', [])
    if sessions:
        seconds = sum(s['elapsed_seconds'] for s in sessions)
        sections.append(f'Qwen completion-session wall time: {seconds:.1f} seconds including model setup. '
                        'This excludes generation of the original pool.')
    sections.append('Per-job elapsed wall time, state, CPU allocation and start/end timestamps are '
                    'recorded in job_runtime_audit.csv. Parallel job durations must not be summed '
                    'and called overall elapsed wall time.' if runtime_audit else
                    'A per-job Slurm runtime audit was unavailable when this report ran; consult the job logs for wall time.')
    for panel, data in results.items():
        flowr = data.get('flowr')
        if flowr is None:
            continue
        short = sum(r['shortfall'] > 0 for r in flowr['molecules'])
        attempts = sum(r['attempts'] for r in flowr['molecules'])
        hours = sum(r.get('elapsed_seconds', 0) for r in flowr['molecules'])/3600
        sections.append(f'{panel}: final FlowR has {flowr["conformer_count"]:,} candidates for '
                        f'{flowr["molecule_count"]} molecules; minimum pool {flowr["minimum_count"]}; '
                        f'{short} molecules remain short after completion attempts. Recorded additional '
                        f'network attempts: {attempts:,}; summed sampling time: {hours:.2f} hours. '
                        'This time excludes model setup and is not parallel elapsed wall time. '
                        'Short source quotas are not reassigned; realized allowances and '
                        'equal-availability results accompany the main comparison.')
    sections += ['', 'Final comparisons (PB-valid; requested allowance 1,000)']
    common_frontier = set(r.name for r in recipes)
    joint_leaders = common_frontier.copy()
    for panel, data in results.items():
        table = data['table']
        ordered = table.sort_values('coverage_0p5', ascending=False)
        first, second = ordered.index[:2]
        delta, low, high = paired_interval(data['pairs'], first, second)
        sections.append(f'{panel}: the highest coverage point estimate is {first}, '
            f'{100*ordered.iloc[0].coverage_0p5:.1f}% at 0.5 A and '
            f'{100*ordered.iloc[0].coverage_0p75:.1f}% at 0.75 A, with '
            f'{ordered.iloc[0].clusters_r0p5:.1f} mean shape clusters. Versus the next coverage '
            f'point estimate ({second}), the paired difference is {100*delta:+.2f} percentage '
            f'points (95% interval {100*low:+.2f} to {100*high:+.2f}).')
        for single in sources:
            if single == first:
                continue
            delta, low, high = paired_interval(data['pairs'], first, single)
            shape, shape_low, shape_high = paired_interval(data['diversity_pairs'], first, single, 'clusters')
            sections.append(f'  {first} versus {single}: coverage {100*delta:+.2f} percentage points '
                            f'({100*low:+.2f}, {100*high:+.2f}); clusters {shape:+.1f} '
                            f'({shape_low:+.1f}, {shape_high:+.1f}). These are exploratory paired 95% intervals.')
        sections.append('  Recipe                              Coverage(0.5A)  Clusters(0.5A)  PB-valid  Match support  Draws')
        for name in [r.name for r in recipes]:
            row = table.loc[name]
            sections.append(f'  {name:35s} {100*row.coverage_0p5:7.1f}% {row.clusters_r0p5:16.1f} '
                            f'{row.expected_valid:9.1f} {row.support:14.1f} {row.realized_draws:6.0f}')
        sections.append('  PB-valid is the expected number passing PB per sampled library; match support '
                        'is the expected number of PB-valid matching candidates per query, averaged '
                        'within and then across molecules. All counts refer to the query-eligible cohort.')
        candidate = 'qwen50_loqi25_flowr25' if has_flowr else 'qwen50_loqi50'
        comparator = ('qwen50_loqi50' if has_flowr else
                      max(['qwen', 'loqi'], key=lambda name: table.loc[name, 'coverage_0p5']))
        comparison_label = ('Fixed three-way 50/25/25 versus Qwen/LoQI 50/50' if has_flowr else
                            f'Fixed Qwen/LoQI 50/50 versus the higher-coverage single ({comparator})')
        delta, low, high = paired_interval(data['pairs'], candidate, comparator)
        shape, shape_low, shape_high = paired_interval(data['diversity_pairs'], candidate, comparator, 'clusters')
        sections.append(f'  {comparison_label}: coverage '
                        f'{100*delta:+.2f} percentage points ({100*low:+.2f}, {100*high:+.2f}); '
                        f'clusters {shape:+.1f} ({shape_low:+.1f}, {shape_high:+.1f}). '
                        'Parentheses are exploratory paired 95% intervals.')
        coverage = data['coverage']
        sensitivity = coverage[(coverage.recipe == candidate) & (coverage.gate == 'posebusters') &
                               (coverage.radius == .5) & (coverage.budget == 1000)]
        equal = sensitivity[(sensitivity.scenario == 'equal_availability') & (sensitivity.weighting == 'queries')].iloc[0]
        groups = sensitivity[(sensitivity.scenario == 'requested') & (sensitivity.weighting == 'reference_groups')].iloc[0]
        sections.append(f'  {candidate} sensitivity coverage: {100*equal.coverage:.1f}% under equal '
                        f'availability; {100*groups.coverage:.1f}% with equal reference-group weighting.')
        frontier = set()
        for name,row in table.iterrows():
            dominated = ((table.coverage_0p5 >= row.coverage_0p5) &
                         (table.clusters_r0p5 >= row.clusters_r0p5) &
                         ((table.coverage_0p5 > row.coverage_0p5) | (table.clusters_r0p5 > row.clusters_r0p5))).any()
            if not dominated:
                frontier.add(name)
            if not ((row.coverage_0p5 >= table.coverage_0p5).all() and
                    (row.clusters_r0p5 >= table.clusters_r0p5).all()):
                joint_leaders.discard(name)
        common_frontier &= frontier
        raw, minimized = table.loc['rdkit_raw'], table.loc['rdkit_minimized']
        sections.append(f'  RDKit minimized versus raw: {minimized.clusters_r0p5-raw.clusters_r0p5:+.1f} '
                        f'clusters and {100*(minimized.coverage_0p5-raw.coverage_0p5):+.2f} coverage '
                        'percentage points. This describes these control arms, not every possible minimization strategy.')
    sections += ['', 'Is there a clear winner?',
        ('Descriptive joint leaders across both panels: '+', '.join(sorted(joint_leaders))+'. '
         'Check the paired intervals and sensitivities before treating that as a reliable advantage.' if joint_leaders else
         'No recipe dominates every alternative on both coverage and geometric diversity across both panels. '
         'The highest coverage point estimate is therefore not automatically the best training-library recipe.'),
        'Recipes on the descriptive coverage/diversity frontier in both panels: '+(', '.join(sorted(common_frontier)) or 'none')+'. '
        'These are candidates for choosing a practical compromise, not independently validated optima.', '',
        'What can be concluded for the paper',
        'Use the final numbers to distinguish coverage recovery from geometric breadth and to justify '
        'a transparent fixed-budget dataset recipe. A mixture can be motivated by complementary strengths '
        'if the gain survives the paired comparisons and availability/reference sensitivities; do not '
        'call it uniquely optimal merely because it ranks first in this exploratory grid. The RDKit '
        'control comparison quantifies the observed diversity/coverage trade-off of minimization. '
        'The defensible claim concerns coverage of experimental feature geometries and sampled shape '
        'diversity, not improved downstream pharmacophore-conditioned generation without a training test.', '',
        'Important limitations',
        'ATC48 has exact Qwen training-identity overlap for 20/48 molecules (17/44 query-eligible); '
        'no exact overlap does not establish scaffold separation or absence from other generators\' '
        'training data. The panels remain small, queries from the same molecule are correlated, '
        'and native-reference quality and source shortfalls can matter. Native references were not '
        'filtered by candidate PB validity: PLINDER22 retains 18 queries with archived geometry flags, '
        'while ATC48 native references were not PB-screened. Unscreened references must not be '
        'interpreted as failed reference checks. PB is a screening rule rather '
        'than a binding assay. Report source sizes, quality diagnostics and exploratory uncertainty '
        'alongside any preferred recipe. Completion sampling times are not hardware-matched full '
        f'generation speed comparisons. All {len(recipes)} rows and supporting strata remain in the panel CSV files.']
    return '\n\n'.join(sections)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=Path('src/casf_benchmark/config/dataset_recipe_experiment.yaml'))
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    config['_sha256'] = sha256(args.config)
    results, qwen = read_completed(config)
    root = Path(config['output_dir'])
    runtime_audit = collect_job_runtime(root)
    destination = root/'FINAL_REPORT.txt'
    temporary = destination.with_name(f'.{destination.name}.{os.getpid()}.tmp')
    temporary.write_text(render_report(config, results, qwen, runtime_audit))
    os.replace(temporary, destination)
    atomic_json_dump(root/'final_report_manifest.json', {
        'status':'complete', 'created_at_utc':datetime.now(timezone.utc).isoformat(),
        'report_sha256':sha256(destination), 'report_source_sha256':sha256(Path(__file__)),
        'job_runtime_audit_sha256':sha256(runtime_audit) if runtime_audit else None,
        'config_sha256':config['_sha256'], 'panels':{
            panel:sha256(root/panel/'run_manifest.json') for panel in config['panels']}})
    print(destination)


if __name__ == '__main__':
    main()
