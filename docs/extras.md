# Extras

Installation, cluster paths, deployment, and optional tooling not covered in the core method docs.

## Installation

| Use case | Command |
| --- | --- |
| Dashboard only | `conda env create -f environment-dashboard.yml` → `streamlit run apps/dashboard/streamlit_app.py` |
| Full pipeline | `conda env create -f environment-analysis.yml` → `pip install -e ".[dev]"` → `pytest tests/ -q` |

Requires Python 3.10+. Analysis env needs RDKit + zarr (via conda). Run `pip install -e ".[dev]"` once per environment; no `PYTHONPATH` needed.

Large bundled files: `git lfs pull` if needed.

## Weka paths

```bash
export CASF_BENCHMARK_DATA_ROOT=/mnt/weka/mbedrosian
```

| Role | Path |
| --- | --- |
| Core/ref mapping CSV | `$CASF_BENCHMARK_DATA_ROOT/data/casf16/casf16_*_exact_intersection.csv` (preferred default when present; bundled copy under `data/mapping/`) |
| Intersection MOL2 | `.../CASF16/core_chembl3d_exact_intersection_ligands/` (or ref) |
| ChEMBL3D topologies | `$CASF_BENCHMARK_DATA_ROOT/data/chembl3d/topologies/` |
| RDKit/torsion output | `.../pharma_generation_analysis/{core,ref}_pb_full_dynamic_chembl_count/` |
| External ML pools | `.../codex_dir/` (see [generation_methods.md](generation_methods.md)) |
| Production dashboard DB | `.../pharma_generation_analysis/casf_analysis_dashboard.sqlite` |

Defaults in [`src/casf_benchmark/paths.py`](../src/casf_benchmark/paths.py). Production analysis registry: [`src/casf_benchmark/config/casf_analysis_sources.weka.yaml`](../src/casf_benchmark/config/casf_analysis_sources.weka.yaml).

## Cluster ingest

After external inference on Weka:

```bash
export CASF_BENCHMARK_DATA_ROOT=/mnt/weka/mbedrosian
./scripts/ingest_external_generation.sh REF /path/to/run qwen_4b_bigdata
# or materialize + analyze separately — see materialization.md
```

Then copy `analysis/tables/geometric_per_ligand_long.csv` into [`data/results/runs/{run_id}/`](../data/results/runs/) if updating the bundled repo, and rebuild the dashboard.

Slurm helpers: [`scripts/submit_casf.sh`](../scripts/submit_casf.sh), [`scripts/run_pb_once_casf_analysis.sbatch`](../scripts/run_pb_once_casf_analysis.sbatch).

## Qwen checkpoint batches

Checkpoint labels and the inference output directory for each cohort live in one
manifest, [`src/casf_benchmark/config/qwen_generation_runs.yaml`](../src/casf_benchmark/config/qwen_generation_runs.yaml)
(read via `casf_benchmark.catalog.load_generation_runs`). Add a checkpoint there once
and every batch script picks it up.

Run the pickle → SDF/manifest → materialize → analyze pipeline for all checkpoints,
core cohort first:

```bash
export CASF_GENERATION_RESULTS_ROOT=/path/to/gen_results   # holds the manifest's dirs
export CASF_RUNS_ROOT=/path/to/casf_benchmark_runs         # materialized run outputs
./scripts/run_casf_batch_all.sh
```

Both scripts are idempotent — they skip steps whose output already exists, so a
partial batch can simply be re-run. Per-item logs and `batch_summary.log` land under
`$CASF_RUNS_ROOT/_logs/`. To run a single item:
`./scripts/run_casf_batch_item.sh LABEL core|ref PICKLE_DIR` (see its header for
`SHARD_PARALLEL`, `PB_WORKERS`, `ANALYZE_WORKERS`).

Raw pickles are converted by
[`scripts/convert_qwen_pickle_to_raw_generation.py`](../scripts/convert_qwen_pickle_to_raw_generation.py),
which maps generated SMILES back to CASF `ligand_id`s via the intersection mapping CSV.

## Druglike test set

A small hand-picked set of named drugs/candidates with no crystal conformers, so it
is scored on validity/diversity/energy only — no RMSD recovery. Three steps:

```bash
export CASF_GENERATION_RESULTS_ROOT=/path/to/gen_results
export CASF_DRUGLIKE_PICKLE=/path/to/druglike_smi.pickle
export CASF_EVAL_OUT_DIR=/path/to/druglike_eval

./scripts/run_druglike_eval_batch.sh                      # per-checkpoint CSVs
python scripts/build_druglike_eval_db.py --eval-dir "$CASF_EVAL_OUT_DIR"
python scripts/build_druglike_covmat.py \
  --generation-results-root "$CASF_GENERATION_RESULTS_ROOT" \
  --druglike-db "$CASF_EVAL_OUT_DIR/druglike_eval.sqlite"
```

The last step folds in COV/MAT results and publishes `extended_druglike_summary` /
`extended_druglike_per_molecule` into the extended-analysis sidecar DB, where the
dashboard shows them as "Druglike summary" and "Druglike per-molecule" tabs.

## Streamlit Cloud

Org repos cannot authorize Streamlit’s GitHub App. Deploy from personal mirror [`MenuaB/casf-benchmark`](https://github.com/MenuaB/casf-benchmark):

1. [share.streamlit.io](https://share.streamlit.io) → `MenuaB/casf-benchmark` → main file `apps/dashboard/streamlit_app.py` → Python 3.10  
2. Uses root [`requirements.txt`](../requirements.txt)

Sync mirror after YerevaNN pushes:

```bash
./scripts/push_with_mirror.sh
```

Optional auto-mirror: add `PERSONAL_MIRROR_PAT` secret → [`.github/workflows/mirror-personal.yml`](../.github/workflows/mirror-personal.yml).

## Extended analysis script (optional)

[`scripts/extended_casf_analysis.py`](../scripts/extended_casf_analysis.py) computes K-efficiency, paired deltas, frontier ranks from the dashboard DB → `extended_analysis/` + [`extended_casf_analysis.sqlite`](../data/results/extended_casf_analysis.sqlite). Results: [casf16_hypothesis.md](casf16_hypothesis.md) · [casf16-core-hypothesis-assessment.md](casf16-core-hypothesis-assessment.md).

```bash
casf-extended-analysis --skip-k-efficiency
casf-extended-analysis --recompute-k-rmsd --k-workers 8
```

K-efficiency and the energy windows read per-conformer SDFs under each run's
`generation/`, which the bundled run roots under `data/results/runs/` do not carry.
Point them at fuller copies without rebuilding the dashboard:

```bash
python scripts/build_generation_root_overrides.py      # probes the weka sources config
casf-extended-analysis \
  --generation-root-overrides src/casf_benchmark/config/casf_generation_root_overrides.yaml \
  --allow-k-efficiency-failures
```

The generated override map is machine-specific and gitignored. `--allow-k-efficiency-failures`
keeps the tasks that succeeded and reports the rest in `extended_k_efficiency_failures.csv`
instead of aborting the whole run.

## Repository layout

```
src/casf_benchmark/config/   YAML catalogs (packaged)
data/mapping/    intersection CSVs (bundled)
data/results/    dashboard SQLite, per-run analysis CSVs
docs/            this documentation set
scripts/         CLIs and Slurm
src/casf_benchmark/   Python package
apps/dashboard/  Streamlit
tests/
```
