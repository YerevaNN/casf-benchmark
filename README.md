# CASF Benchmark

Standalone pipeline for the CASF–ChEMBL3D exact-intersection conformer benchmark: data preparation, classical generation, external-generator ingestion, geometric analysis, and an interactive Streamlit dashboard.

Precomputed results are bundled under `data/` so you can explore the dashboard immediately after cloning.

## Quick start — view the dashboard (~5 minutes)

```bash
git clone <repo-url> casf-benchmark
cd casf-benchmark

conda env create -f environment-dashboard.yml
conda activate casf-benchmark-dashboard

streamlit run apps/dashboard/streamlit_app.py
```

The app opens `data/results/casf_analysis_dashboard.sqlite` by default. Override with:

```bash
export CASF_DASHBOARD_DB=/path/to/casf_analysis_dashboard.sqlite
```

Bundled artifacts:

| File | Description |
| --- | --- |
| `data/mapping/casf16_*_exact_intersection.csv` | Core/ref ligand mapping tables |
| `data/results/casf_analysis_dashboard.sqlite` | Dashboard database |
| `data/results/casf_analysis_master.csv` | Wide comparison summary |
| `data/results/casf_per_ligand_long.csv` | Global per-ligand long table |
| `data/results/extended_casf_analysis.sqlite` | Extended analysis tab |
| `data/results/runs/*/analysis/tables/geometric_per_ligand_long.csv` | Per-source analysis inputs |

Large files may use Git LFS — run `git lfs pull` after clone if needed.

## Public web dashboard (Streamlit Cloud)

The live app is deployed from a **personal mirror** ([menuab/casf-benchmark](https://github.com/menuab/casf-benchmark)) because YerevaNN org repos cannot authorize Streamlit’s GitHub App.

**Sync is not automatic** unless you use the mirror push script or enable the optional GitHub Action. After pushing to YerevaNN only, run:

```bash
./scripts/push_with_mirror.sh
```

Full setup: [docs/deployment.md](docs/deployment.md).

## Full analysis environment

For data prep, generation, materialization, analysis, and dashboard rebuild:

```bash
conda env create -f environment-analysis.yml
conda activate casf-benchmark-analysis
pytest tests/ -q
```

See [docs/installation.md](docs/installation.md) for details.

## Add your own generator (e.g. Qwen)

Inference runs **outside** this repo. Produce outputs matching the layout in [docs/adding_external_generators.md](docs/adding_external_generators.md), then:

```bash
python scripts/materialize_casf_generation_sets.py --root /path/to/your/root --method your_method
python scripts/analyze_casf_conformer_sets.py --run-root /path/to/your/root ...
# Update config/casf_analysis_sources.yaml, then:
python scripts/build_casf_analysis_master_csv.py
python scripts/build_casf_analysis_dashboard_db.py
streamlit run apps/dashboard/streamlit_app.py
```

## Method documentation

| Doc | Topic |
| --- | --- |
| [casf16_chembl3d_exact_match_method.md](docs/casf16_chembl3d_exact_match_method.md) | Data preparation |
| [casf16_chembl3d_conformer_generation_method.md](docs/casf16_chembl3d_conformer_generation_method.md) | RDKit/torsion generation |
| [casf16_materialize_generation_sets_method.md](docs/casf16_materialize_generation_sets_method.md) | External pool materialization |
| [casf16_generation_pipeline_and_models_method.md](docs/casf16_generation_pipeline_and_models_method.md) | Generator catalog + output contract |
| [casf16_benchmark_analysis_method.md](docs/casf16_benchmark_analysis_method.md) | Geometric analysis pipeline |
| [casf16_dashboard_method.md](docs/casf16_dashboard_method.md) | Dashboard architecture |

## Cluster / Weka production paths

Set `CASF_BENCHMARK_DATA_ROOT=/mnt/weka/mbedrosian` to use production data paths. A Weka-oriented analysis config is kept at `config/casf_analysis_sources.weka.yaml`.

## Repository layout

```
casf-benchmark/
├── config/           YAML catalogs and analysis source registry
├── data/             Bundled mapping CSVs and precomputed results
├── docs/             Method notes and guides
├── scripts/          CLI tools and Slurm wrappers
├── src/casf_benchmark/   Python package
├── apps/dashboard/   Streamlit UI
└── tests/
```
