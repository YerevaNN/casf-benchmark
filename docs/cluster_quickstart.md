# Cluster quickstart (Weka)

Self-contained guide for cloning this repo on the analysis cluster, ingesting **new external generator outputs** (e.g. Qwen ref inference), and updating the production dashboard — without custom help.

Prerequisites:

- Weka mount at `/mnt/weka/mbedrosian` (or set `CASF_BENCHMARK_DATA_ROOT`)
- Conda
- Git LFS (`git lfs install && git lfs pull` if you need bundled artifacts)

Path reference: [weka_data_paths.md](weka_data_paths.md).

---

## 1. Clone and install (once per machine)

```bash
git clone https://github.com/MenuaB/casf-benchmark.git
cd casf-benchmark

conda env create -f environment-analysis.yml
conda activate casf-benchmark-analysis

export CASF_BENCHMARK_DATA_ROOT=/mnt/weka/mbedrosian
export REPO_ROOT="$(pwd)"
export PYTHONPATH="${REPO_ROOT}/src"
```

Verify Weka inputs are visible:

```bash
ls /mnt/weka/mbedrosian/data/casf16/casf16_ref_chembl3d_exact_intersection.csv
ls /mnt/weka/mbedrosian/data/casf16/CASF16_REF/ref_chembl3d_exact_intersection_ligands | head
ls /mnt/weka/mbedrosian/data/chembl3d/topologies | head
```

---

## 2. What the repo expects from inference (outside this repo)

Inference is **not** implemented here. After running Qwen (or any model), outputs must look like:

```
{RUN_ROOT}/generation/manifest.tsv
{RUN_ROOT}/generation/{method}/{mol_id}.sdf
```

- `{method}` is **untiered** (e.g. `qwen_4b_bigdata`, not `qwen_4b_bigdata_fixed`)
- Up to 1000 conformers per ligand in multi-record SDFs
- Manifest: tab-separated, one row per `(generation_method, mol_id)`

See [adding_external_generators.md](adding_external_generators.md) for the full contract.

### Canonical Weka run roots

| Cohort | Generator | `{RUN_ROOT}` |
| --- | --- | --- |
| Core | Qwen (all checkpoints) | `/mnt/weka/mbedrosian/codex_dir/qwen_gens` |
| Ref | Qwen (all checkpoints) | `/mnt/weka/mbedrosian/codex_dir/qwen/generations/casf16_ref_qwen_1k` |
| Core | LOQI | `/mnt/weka/mbedrosian/codex_dir/loqi/generations/casf16_core_loqi_1k` |
| Ref | LOQI | `/mnt/weka/mbedrosian/codex_dir/loqi/generations/casf16_ref_loqi_1k` |

Constants also live in `src/casf_benchmark/paths.py` and `config/casf_generation_families.yaml`.

---

## 3. Example: new Qwen ref inference → dashboard

Assume inference wrote raw SDFs under the **ref** root for checkpoint `qwen_4b_bigdata`.

### Step A — Materialize + analyze (one command)

```bash
cd casf-benchmark
conda activate casf-benchmark-analysis
export PYTHONPATH=src

chmod +x scripts/ingest_external_generation.sh
./scripts/ingest_external_generation.sh REF \
  /mnt/weka/mbedrosian/codex_dir/qwen/generations/casf16_ref_qwen_1k \
  qwen_4b_bigdata
```

This runs:

1. `materialize_casf_generation_sets.py` → `_fixed`, `_dynamic`, `_chembl_count` tier dirs
2. `analyze_casf_conformer_sets.py` → `{RUN_ROOT}/analysis/tables/geometric_per_ligand_long.csv`

For **core** Qwen into the existing core pool:

```bash
./scripts/ingest_external_generation.sh CORE \
  /mnt/weka/mbedrosian/codex_dir/qwen_gens \
  qwen_4b_revisited
```

Multiple checkpoints in one invocation:

```bash
./scripts/ingest_external_generation.sh REF \
  /mnt/weka/mbedrosian/codex_dir/qwen/generations/casf16_ref_qwen_1k \
  qwen_4b_bigdata qwen_4b_revisited
```

### Step B — Rebuild production dashboard

`config/casf_analysis_sources.weka.yaml` already lists `qwen_ref` at the ref run root. After analysis completes:

```bash
./scripts/rebuild_dashboard_weka.sh
```

Outputs:

| Artifact | Path |
| --- | --- |
| Dashboard SQLite | `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite` |
| Master CSV | `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_master.csv` |
| Per-ligand long CSV | `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_per_ligand_long.csv` |

### Step C — Verify

```bash
sqlite3 /mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite \
  "SELECT run_id, method, ligands, casf_hit_0p75 FROM comparison_rows \
   WHERE run_id='qwen_ref' AND tier='fixed' ORDER BY method LIMIT 10;"
```

View in Streamlit:

```bash
CASF_DASHBOARD_DB=/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite \
  streamlit run apps/dashboard/streamlit_app.py
```

---

## 4. Adding a **new** checkpoint name

If inference uses a method string not in `config/casf_generation_families.yaml`:

1. Add a family block (copy an existing Qwen entry, change `id`, `method_prefix`, `display_label`, `variant`).
2. Set `core_root` / `ref_root` to the Weka paths above.
3. Run ingest (step 3A) with the new `--method` name.
4. Rebuild dashboard (step 3B).

No change to `casf_analysis_sources.weka.yaml` is needed unless you create a **new run root** (separate directory from existing Qwen pools).

---

## 5. Full multi-cohort refresh (all generators)

When all SDF pools on Weka are up to date and you want to re-analyze everything:

```bash
sbatch scripts/run_pb_once_casf_analysis.sbatch
```

That Slurm job:

- Analyzes reference baselines (core + ref)
- Analyzes RDKit/torsion, Qwen (core **and ref**), LOQI, DMT-L, Torsional Diffusion, MCF
- Rebuilds master CSV + dashboard DB on Weka using `casf_analysis_sources.weka.yaml`

---

## 6. Optional: bundle results for git / offline dashboard

To update the laptop clone (analysis tables only, no SDFs):

```bash
# Example: ship new Qwen ref analysis table into the repo
mkdir -p data/results/runs/qwen_ref/analysis/tables
cp /mnt/weka/mbedrosian/codex_dir/qwen/generations/casf16_ref_qwen_1k/analysis/tables/geometric_per_ligand_long.csv \
   data/results/runs/qwen_ref/analysis/tables/
```

Add to `config/casf_analysis_sources.yaml`:

```yaml
  - run_id: qwen_ref
    label: "Qwen (ref)"
    ligand_set: ref
    root: data/results/runs/qwen_ref
```

Rebuild bundled dashboard locally:

```bash
python scripts/build_casf_analysis_master_csv.py
python scripts/build_casf_analysis_dashboard_db.py
```

Commit CSV + updated SQLite (Git LFS for large files).

---

## 7. Troubleshooting

| Problem | Fix |
| --- | --- |
| `FileNotFoundError` for ligand MOL2 | Check `--casf-ligand-dir` / Weka mount; see [weka_data_paths.md](weka_data_paths.md) |
| Materialize finds no SDFs | Confirm `{RUN_ROOT}/generation/{method}/{mol_id}.sdf` and manifest rows |
| Analyzer missing methods | Run materialize first; tier dirs must exist |
| Dashboard missing new run | Use `rebuild_dashboard_weka.sh` (Weka config), not default bundled config |
| `qwen_ref` empty in dashboard | Ref analysis not finished, or rebuild ran before `geometric_per_ligand_long.csv` existed |

---

## Config files (which to edit when)

| File | Use on cluster | Use in git clone |
| --- | --- | --- |
| `config/casf_analysis_sources.weka.yaml` | Production rebuild (`rebuild_dashboard_weka.sh`) | — |
| `config/casf_analysis_sources.yaml` | — | Offline rebuild from `data/results/runs/` |
| `config/casf_generation_families.yaml` | Method catalog / display labels | Same |

Related: [adding_external_generators.md](adding_external_generators.md), [installation.md](installation.md), [casf16_dashboard_method.md](casf16_dashboard_method.md).
