"""Queue configured recipe analysis, with only the required source dependencies."""

from __future__ import annotations

import argparse
import json
import re
import shlex
import subprocess
from pathlib import Path

import yaml

from casf_benchmark.analysis.pharmacophore import atomic_json_dump, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=Path('src/casf_benchmark/config/dataset_recipe_experiment.yaml'))
    parser.add_argument('--flowr-atc-job')
    parser.add_argument('--flowr-plinder-job')
    parser.add_argument('--qwen-job', help='Optional dependency if ATC48 Qwen completion is still running')
    parser.add_argument('--max-parallel', type=int, default=8)
    parser.add_argument('--submit', action='store_true', help='Without this flag, print the plan only')
    args = parser.parse_args()
    for job in [args.flowr_atc_job, args.flowr_plinder_job, args.qwen_job]:
        if job is not None and not re.fullmatch(r'\d+', job):
            parser.error('Generation dependencies must be numeric Slurm job IDs')
    config = yaml.safe_load(args.config.read_text())
    root = Path(config['output_dir'])
    protocol = json.loads((root/'protocol.json').read_text())
    if protocol['config_sha256'] != sha256(args.config):
        raise ValueError('Freeze the current configuration before queueing analysis')
    if args.max_parallel <= 0:
        parser.error('--max-parallel must be positive')
    has_flowr = 'flowr' in config['sources']
    if has_flowr and not (args.flowr_atc_job and args.flowr_plinder_job):
        parser.error('Active FlowR sources require both generation job IDs')
    if not has_flowr and (args.flowr_atc_job or args.flowr_plinder_job):
        parser.error('FlowR dependencies cannot be supplied when FlowR is excluded')
    repository = Path(__file__).resolve().parents[1]
    ledger = root/'analysis_jobs.json'
    if args.submit and ledger.exists():
        raise RuntimeError(f'An analysis queue already exists; inspect {ledger} before resubmitting')
    jobs = {}
    generation_jobs = {
        'flowr_atc48': args.flowr_atc_job, 'flowr_plinder22': args.flowr_plinder_job,
        'qwen_atc48': args.qwen_job}
    record = {'config_sha256': sha256(args.config), 'generation_jobs': {
        name: job for name, job in generation_jobs.items() if job},
        'jobs': jobs, 'submission_complete': False}

    def queue(label, dependencies, options, command):
        dependency_args = [f'--dependency=afterok:{":".join(dependencies)}'] if dependencies else []
        invocation = ['sbatch', '--parsable', f'--chdir={repository}', *dependency_args,
                      f'--output={root}/{label}-%j.out', f'--job-name=recipe-{label}',
                      *options, *command]
        print(shlex.join(invocation), flush=True)
        job_id = subprocess.check_output(invocation, text=True).strip().split(';')[0] if args.submit else f'<{label}>'
        jobs[label] = {'id': job_id, 'dependencies': dependencies, 'command': invocation}
        if args.submit:
            atomic_json_dump(ledger, record)
        return job_id

    summaries = []
    for panel in config['panels']:
        dependencies = []
        if has_flowr:
            generation = args.flowr_plinder_job if panel == 'plinder22' else args.flowr_atc_job
            merge = queue(f'{panel}-merge', [generation],
                          ['--partition=research_cpu', '--cpus-per-task=1', '--gres=cpuonly:1', '--mem=8G', '--time=00:30:00'],
                          ['--wrap', shlex.join([str(repository/'.venv/bin/python'),
                           'scripts/generate_recipe_flowr.py', '--config', str(args.config.resolve()), '--panel', panel, '--merge'])])
            dependencies.append(merge)
        if panel == 'atc48' and args.qwen_job:
            dependencies.append(args.qwen_job)
        count = protocol['panels'][panel]['molecule_count']
        array = f'--array=0-{count-1}%{args.max_parallel}'
        prepare = queue(f'{panel}-prepare-final', dependencies, [array],
                        ['scripts/run_dataset_recipes.sbatch', 'prepare', '--config', str(args.config.resolve()), '--panel', panel])
        diversity = queue(f'{panel}-diversity-final', [prepare], [array, '--time=04:00:00'],
                          ['scripts/run_dataset_recipes.sbatch', 'diversity', '--config', str(args.config.resolve()), '--panel', panel])
        summaries.append(queue(f'{panel}-summary-final', [diversity],
              ['--cpus-per-task=2', '--gres=cpuonly:2', '--mem=16G', '--time=00:30:00'],
              ['scripts/run_dataset_recipes.sbatch', 'summarize', '--config', str(args.config.resolve()), '--panel', panel]))
    report = queue('final-report', summaries,
                   ['--partition=research_cpu', '--cpus-per-task=1', '--gres=cpuonly:1', '--mem=8G', '--time=00:30:00'],
                   ['--wrap', shlex.join([str(repository/'.venv/bin/python'), '-m',
                    'casf_benchmark.cli.report_dataset_recipes', '--config', str(args.config.resolve())])])
    record['submission_complete'] = True
    if args.submit:
        atomic_json_dump(ledger, record)
        atomic_json_dump(root/'final_report_job.json', {
            'job_id': report, 'dependencies': summaries, 'status': 'queued',
            'automatic_chat_notification_enabled': False,
            'report_path': str(root/'FINAL_REPORT.txt')})
    else:
        print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
