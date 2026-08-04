# CASF Benchmark

CASF–ChEMBL3D exact-intersection conformer benchmark: data prep → generation → analysis → Streamlit dashboard.

Bundled results under `data/results/` let you run the dashboard without regenerating anything.

## Quick start (dashboard)

```bash
git clone https://github.com/YerevaNN/casf-benchmark.git && cd casf-benchmark
conda env create -f environment-dashboard.yml && conda activate casf-benchmark-dashboard
streamlit run apps/dashboard/streamlit_app.py
```

Uses [`data/results/casf_analysis_dashboard.sqlite`](data/results/casf_analysis_dashboard.sqlite). Override with `CASF_DASHBOARD_DB`.

## Full pipeline (cluster)

```bash
conda env create -f environment-analysis.yml && conda activate casf-benchmark-analysis
export PYTHONPATH=src CASF_BENCHMARK_DATA_ROOT=/mnt/weka/mbedrosian
```

## Documentation

| Doc | Topic |
| --- | --- |
| [docs/data_preparation.md](docs/data_preparation.md) | CASF–ChEMBL3D intersection mapping |
| [docs/generation_methods.md](docs/generation_methods.md) | RDKit/torsion + external model output contract |
| [docs/analyzer.md](docs/analyzer.md) | Geometric analysis pipeline |
| [docs/materialization.md](docs/materialization.md) | External pool → tiered layout |
| [docs/dashboard.md](docs/dashboard.md) | SQLite build + Streamlit |
| [docs/extras.md](docs/extras.md) | Install, Weka paths, deployment |
| [docs/casf16_hypothesis.md](docs/casf16_hypothesis.md) | Hypothesis (H0/H1, metrics, predictions) |
| [docs/casf16-core-hypothesis-assessment.md](docs/casf16-core-hypothesis-assessment.md) | Core 94-ligand results vs hypothesis |

## Common commands

```bash
./scripts/submit_casf.sh generate core          # RDKit/torsion Slurm generation
./scripts/submit_casf.sh analyze core           # geometric analysis
python scripts/build_casf_analysis_dashboard_db.py
./scripts/push_with_mirror.sh                   # sync personal mirror for Streamlit Cloud
```

Public dashboard mirror: [MenuaB/casf-benchmark](https://github.com/MenuaB/casf-benchmark) (see [docs/extras.md](docs/extras.md#streamlit-cloud)).
