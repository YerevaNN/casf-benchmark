# Installation

Two environments are supported. Choose based on what you need.

## Dashboard only

**Use when:** exploring bundled precomputed results.

**Requirements:** Python 3.10, pandas, Streamlit.

```bash
cd casf-benchmark
conda env create -f environment-dashboard.yml
conda activate casf-benchmark-dashboard
streamlit run apps/dashboard/streamlit_app.py
```

Or with pip from an existing Python 3.10 environment:

```bash
pip install -e ".[dashboard]"
streamlit run apps/dashboard/streamlit_app.py
```

### Environment variables

| Variable | Default | Purpose |
| --- | --- | --- |
| `CASF_DASHBOARD_DB` | `data/results/casf_analysis_dashboard.sqlite` | Main dashboard SQLite |
| `CASF_EXTENDED_DB` | `data/results/extended_casf_analysis.sqlite` | Extended analysis SQLite |

## Full analysis

**Use when:** running data prep, classical generation, materializing external pools, geometric analysis, or rebuilding the dashboard.

**Requirements:** Python 3.10, RDKit, zarr, PoseBusters, NumPy, pandas, PyYAML, matplotlib.

RDKit and zarr are best installed via conda:

```bash
cd casf-benchmark
conda env create -f environment-analysis.yml
conda activate casf-benchmark-analysis
pip install -e ".[dev]"   # if not already installed by the env file
pytest tests/ -q
```

### PYTHONPATH

Scripts expect the package on `PYTHONPATH`. Slurm wrappers set:

```bash
export PYTHONPATH="${REPO_ROOT}/src"
```

For manual runs:

```bash
export PYTHONPATH="$(pwd)/src"
python scripts/analyze_casf_conformer_sets.py ...
```

### Cluster data

Production CASF16 MOL2 ligands, ChEMBL3D zarr/topologies, and **pre-generated conformer SDF pools** live on shared storage (not in this repo). On the analysis cluster:

```bash
export CASF_BENCHMARK_DATA_ROOT=/mnt/weka/mbedrosian
```

Use `config/casf_analysis_sources.weka.yaml` when rebuilding from Weka run roots instead of bundled `data/results/runs/`.

Full path reference (SDF run roots per generator, reference inputs, dashboard artifacts): [weka_data_paths.md](weka_data_paths.md).

**Cluster workflow (clone → ingest new Qwen gens → rebuild dashboard):** [cluster_quickstart.md](cluster_quickstart.md).

### Generator entry point

RDKit/torsion baseline generation is implemented in `src/casf_benchmark/generation/conformer_sets.py`:

```bash
export PYTHONPATH=src
python -m casf_benchmark.generation.conformer_sets --help
# or
python src/casf_benchmark/generation/conformer_sets.py --help
```

Slurm wrappers and `scripts/launch_casf_conformer_sets_parallel.py` invoke this script. Default output roots are under `/mnt/weka/mbedrosian/pharma_generation_analysis/` (see [weka_data_paths.md](weka_data_paths.md)).

## Git LFS

Large bundled artifacts may be tracked with Git LFS:

```bash
git lfs install
git lfs pull
```

## Troubleshooting

**Streamlit cannot find the database.** Confirm `data/results/casf_analysis_dashboard.sqlite` exists or set `CASF_DASHBOARD_DB`.

**RDKit import errors.** Use `environment-analysis.yml`; RDKit is not pip-installable on all platforms.

**ChEMBL3D zarr errors.** Install `zarr` and run with the analysis environment. Reference-mode analysis requires ChEMBL3D data on disk.

**PoseBusters missing.** `pip install posebusters` or use the analysis conda environment.
