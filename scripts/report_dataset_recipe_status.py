"""Inspect this experiment's generation checkpoints, analysis results, and Slurm jobs."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import yaml
import pandas as pd

from casf_benchmark.analysis.pharmacophore import load_pickle
from casf_benchmark.cli.analyze_dataset_recipes import pool_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=Path('src/casf_benchmark/config/dataset_recipe_experiment.yaml'))
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    root = Path(config['output_dir'])
    protocol = json.loads((root/'protocol.json').read_text())
    for panel in config['panels']:
        source = yaml.safe_load(Path(config['panels'][panel]['source_config']).read_text())
        expected_smiles = set(pd.read_csv(Path(source['output_dir'])/'query_manifest.csv',
                                         usecols=['expected_smiles']).expected_smiles)
        inventory = {}
        for name in config['panels'][panel]['methods']:
            try:
                pool = load_pickle(pool_path(config, panel, name, source))
            except FileNotFoundError:
                pool = {}
            counts = {s: len(pool.get(s, [])) for s in expected_smiles}
            if name == 'flowr':
                for path in sorted((root/'generation/flowr'/panel).glob('*.pickle')):
                    record = load_pickle(path)
                    counts[record['smiles']] = len(record['mols'])
            inventory[name] = {'molecules_at_1000': sum(n >= 1000 for n in counts.values()),
                               'candidates': sum(counts.values()),
                               'minimum_pool_size': min(counts.values())}
        print(json.dumps({'panel': panel, 'expected_molecules': protocol['panels'][panel]['molecule_count'],
                          'active_sources': inventory}), flush=True)
        manifest = root/panel/'run_manifest.json'
        if manifest.is_file():
            result = json.loads(manifest.read_text())
            print(f'{panel} last completed analysis: {result["status"]}; source hashes in {manifest}', flush=True)
    ledger_path = root/'analysis_jobs.json'
    if ledger_path.is_file():
        ledger = json.loads(ledger_path.read_text())
        jobs = set(ledger['generation_jobs'].values())
        jobs.update(ledger.get('generation_resume_jobs', {}).values())
        jobs.update(row['id'] for row in ledger['jobs'].values())
        subprocess.run(['sacct', '-X', '-j', ','.join(sorted(jobs)),
                        '--format=JobID,JobName%30,State,ExitCode,Elapsed,AllocCPUS'], check=True)


if __name__ == '__main__':
    main()
