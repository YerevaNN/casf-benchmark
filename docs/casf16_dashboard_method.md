# CASF Analysis Dashboard

## Purpose

This document describes the Streamlit dashboard at `apps/dashboard/streamlit_app.py`, the SQLite database it reads, how that database is built, and the commands required to refresh or launch the UI. The dashboard summarizes geometric analysis outputs for all configured generator cohorts and reference baselines on the CASF–ChEMBL3D intersection panel.

Analysis methodology is in [casf16_benchmark_analysis_method.md](casf16_benchmark_analysis_method.md). Data preparation is in [casf16_chembl3d_exact_match_method.md](casf16_chembl3d_exact_match_method.md).

---

## Architecture

```
Per-run analysis  ({run_root}/analysis/tables/geometric_per_ligand_long.csv)
        │
        ▼
build_casf_analysis_master_csv.py  →  casf_per_ligand_long.csv (global)
        │
        ▼
build_casf_analysis_dashboard_db.py  →  casf_analysis_dashboard.sqlite
        │
        ▼
streamlit run apps/dashboard/streamlit_app.py
        │
        optional ▼
extended_casf_analysis.py  →  extended_casf_analysis.sqlite
                              (or casf_analysis_dashboard_extended.sqlite)
```

The dashboard is **read-only**. It does not run generation, materialization, or geometric analysis.

---

## Data locations

### Primary dashboard database

| Path | Description |
| --- | --- |
| `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite` | Default SQLite database (~51 MB production build) |

Override at launch with environment variable `CASF_DASHBOARD_DB` or the sidebar path field.

### Global CSV (optional export)

| Path | Description |
| --- | --- |
| `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_per_ligand_long.csv` | Concatenated per-ligand long table (~58M rows file) |
| `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_master.csv` | Wide comparison summary (~91 rows) |

The dashboard builder can regenerate these while building the SQLite file.

### Extended analysis (optional second database)

| Path | Description |
| --- | --- |
| `/mnt/weka/mbedrosian/pharma_generation_analysis/extended_analysis/extended_casf_analysis.sqlite` | Sidecar DB with `extended_*` tables only |
| `/mnt/weka/mbedrosian/pharma_generation_analysis/extended_analysis/casf_analysis_dashboard_extended.sqlite` | Copy of dashboard DB with extended tables appended |

Override with sidebar **Extended DB** or `CASF_EXTENDED_ANALYSIS_DB`. If the primary dashboard DB already contains `extended_*` tables, the app uses it directly.

### Upstream per-run roots

Configured in `config/casf_analysis_sources.yaml`:

| run_id | ligand_set | root |
| --- | --- | --- |
| `rdkit_torsion_core` | core | `pharma_generation_analysis/core_pb_full_dynamic_chembl_count` |
| `rdkit_torsion_ref` | ref | `pharma_generation_analysis/ref_pb_full_dynamic_chembl_count` |
| `qwen_core` | core | `codex_dir/qwen_gens` |
| `loqi_core` / `loqi_ref` | core / ref | `codex_dir/loqi/generations/casf16_{core,ref}_loqi_1k` |
| `nextmol_dmt_l_*` | core / ref | `codex_dir/nextmol_dmt_l/generations/...` |
| `torsional_diffusion_*` | core / ref | `codex_dir/torsional_diffusion/generations/...` |
| `mcf_drugs_l_*` | core / ref | `codex_dir/mcf_drugs_l/generations/...` |
| `reference_core` / `reference_ref` | core / ref | `pharma_generation_analysis/reference_datasets/{core,ref}` |

Each run root must contain `analysis/tables/geometric_per_ligand_long.csv` from a completed `analyze_casf_conformer_sets.py` pass.

---

## Database schema (schema version 2)

Built by `scripts/build_casf_analysis_dashboard_db.py` (`PRAGMA user_version = 2`).

| Table | Rows (typical) | Contents |
| --- | --- | --- |
| `analysis_sources` | 13 | Configured run roots and labels |
| `catalog_families` | 19 | Generator families from `casf_generation_families.yaml` |
| `comparison_rows` | ~91 | Aggregated cohort metrics (generation + reference) |
| `comparison_strata` | ~750 | Same metrics split by rotatable-bond / heavy-atom bins |
| `per_ligand_long` | ~41k | Full per-ligand long-form metrics (all sources × ligands) |

### `comparison_rows` semantics

Each row is one **view** of a generator family or reference source:

- **Generation rows** (`row_type = generation`): keyed by `family`, `tier` (`fixed` / `dynamic` / `chembl_count`), and `method` (for example `loqi_raw_fixed`).
- **Reference rows** (`row_type = reference`): baselines `casf_crystal`, `casf_opt`, `chembl3d_sdf`, `chembl3d_gt`, `chembl3d_gt_pb` — tier is `reference`.

Key aggregated columns exposed in the UI:

| Tab | Columns |
| --- | --- |
| Overview | `total_confs`, `mean_confs_per_ligand`, `pairwise_mean`, `pairwise_p90`, `mean_torsion_std_deg` |
| Clustering | `mean_clusters_0p5` … `singleton_fraction_1p0`, `clusters_per_100_1p0`, … |
| Energy | `energy_min`, `energy_max`, `energy_median`, `energy_std` |
| CASF hits | `casf_best_rmsd`, `casf_median_rmsd`, `casf_hit_0p25` … `casf_hit_2p0` |
| CASF opt hits | Same metrics vs optimized poses when available |
| Funnel | Generation-only: `generated_ligands`, `selected_ligands`, `kept_ligands`, `pb_fail_rate_mean`, … |

### `comparison_strata`

Same structure as `comparison_rows` plus:

- `breakdown` — `rotatable_bonds` or `heavy_atoms`
- `stratum` — bin label (for example `0-5`, `6-10`)
- `view_tier` — tier used when aggregating generation rows

### `per_ligand_long`

One row per `(mol_id, method, ligand_set)`. Source of truth for rebuilding aggregates. Includes PoseBusters summaries, cluster counts, energies, CASF RMSD/hit flags, and manifest funnel fields for generation sources.

---

## Dashboard UI

### Sidebar filters

- **Dashboard DB** — path to SQLite file
- **Ligand set** — `core` or `ref`
- **Tier** — `fixed`, `dynamic`, or `chembl_count` (reference rows appear regardless when breakdown is Total)
- **Family** — catalog family id or All
- **Break down aggregates** — Total (whole cohort), Rotatable bonds, or Heavy atoms

### Main tabs

1. **Overview** — diversity and conformer counts  
2. **Clustering** — greedy cluster statistics at 0.5–3.0 Å  
3. **Energy** — MMFF94s pool statistics  
4. **CASF hits** — recovery vs crystal poses  
5. **CASF opt hits** — recovery vs optimized poses  
6. **Funnel** — generation filter funnel (PB pass rates, kept vs target)

Each tab renders a sortable table and a **Download CSV** button for the current view.

### Extended analysis section

When an extended SQLite file is found, additional tabs appear (paired deltas vs ChEMBL3D-PB and RDKit raw, frontier ranks, PB failure mechanisms, K-efficiency curves, rescue cases, sanity checks, etc.). Tables are defined in `EXTENDED_TABLES` inside `streamlit_app.py`.

---

## Commands

### Environment

```bash
cd /home/mbedrosian/code/casf-benchmark
export PYTHONPATH=src
PYTHON=/home/mbedrosian/.conda/envs/chembl3d/bin/python
```

Install dependencies including Streamlit (listed under `[project.optional-dependencies] dev` in `pyproject.toml`, or `uv sync --extra dev`).

### Launch the dashboard

```bash
streamlit run apps/dashboard/streamlit_app.py
```

With explicit database path:

```bash
CASF_DASHBOARD_DB=/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite \
CASF_EXTENDED_ANALYSIS_DB=/mnt/weka/mbedrosian/pharma_generation_analysis/extended_analysis/extended_casf_analysis.sqlite \
streamlit run apps/dashboard/streamlit_app.py
```

Default URL: `http://localhost:8501`.

On a remote machine, use SSH port forwarding or Streamlit’s `--server.address 0.0.0.0` as appropriate for your cluster policy.

### Rebuild dashboard database only

Requires per-run `geometric_per_ligand_long.csv` files under every root in `casf_analysis_sources.yaml`:

```bash
$PYTHON scripts/build_casf_analysis_master_csv.py
$PYTHON scripts/build_casf_analysis_dashboard_db.py
```

Custom paths:

```bash
$PYTHON scripts/build_casf_analysis_master_csv.py \
  --config config/casf_analysis_sources.yaml \
  --per-ligand-long-csv /mnt/weka/mbedrosian/pharma_generation_analysis/casf_per_ligand_long.csv \
  --output-csv /mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_master.csv

$PYTHON scripts/build_casf_analysis_dashboard_db.py \
  --config config/casf_analysis_sources.yaml \
  --output-db /mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite
```

The dashboard builder validates that stratum ligand counts sum to cohort totals; mismatches fail the build.

### Full analysis + dashboard rebuild (production)

Slurm orchestration for reference runs, all generation roots, master CSV, and dashboard DB:

```bash
sbatch scripts/run_pb_once_casf_analysis.sbatch
```

This job:

1. Reference-only analysis → `reference_datasets/{core,ref}`
2. Generation-only analysis for each configured root in `casf_analysis_sources.yaml`
3. `build_casf_analysis_master_csv.py`
4. `build_casf_analysis_dashboard_db.py`

Logs: `outputs/slurm_jobs/casf_pb_once/{jobid}.out`.

Per-cohort classical generation analysis (after RDKit/torsion Slurm merge):

```bash
./scripts/submit_casf.sh analyze core
# or
MERGE_JOB_ID=<id> ./scripts/submit_casf.sh analyze ref
```

### Extended analysis (optional, for Extended Analysis UI section)

```bash
$PYTHON scripts/extended_casf_analysis.py \
  --dashboard-db /mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite \
  --output-dir /mnt/weka/mbedrosian/pharma_generation_analysis/extended_analysis
```

K-efficiency curves (slow; reads SDF/MOL2 from disk):

```bash
$PYTHON scripts/extended_casf_analysis.py \
  --recompute-k-rmsd --k-workers 8
```

Skip K-efficiency for a faster pass:

```bash
$PYTHON scripts/extended_casf_analysis.py --skip-k-efficiency
```

Outputs under `extended_analysis/tables/`, `reports/`, `figures/`, plus `extended_casf_analysis.sqlite`.

---

## End-to-end refresh checklist

To reproduce the dashboard from scratch:

1. **Data prep** — intersection CSV + ligand directories ([casf16_chembl3d_exact_match_method.md](casf16_chembl3d_exact_match_method.md))
2. **Generate** — RDKit/torsion via `./scripts/submit_casf.sh generate {core|ref}`; external ML inference + [materialization](casf16_materialize_generation_sets_method.md) where applicable
3. **Analyze** — `run_pb_once_casf_analysis.sbatch` or per-root `analyze_casf_conformer_sets.py`
4. **Aggregate** — `build_casf_analysis_master_csv.py` + `build_casf_analysis_dashboard_db.py`
5. **Extended (optional)** — `extended_casf_analysis.py`
6. **View** — `streamlit run apps/dashboard/streamlit_app.py`

Adding a new generator family also requires an entry in `casf_generation_families.yaml`, a completed analysis run root, and a matching block in `casf_analysis_sources.yaml` before rebuilding the database.

---

## Troubleshooting

| Symptom | Likely cause |
| --- | --- |
| Dashboard DB not found | Analysis/build steps not run, or wrong path in sidebar |
| Missing family in strata breakdown | Rebuild DB after `casf_stratum_bins.assign_stratum` fix; older builds had incomplete strata |
| Extended section empty | Run `extended_casf_analysis.py` or point sidebar to `extended_casf_analysis.sqlite` |
| Row counts differ core vs ref | Expected: core ~94 ligands, ref ~1219; Qwen is core-only in current config |
| Stale data after re-analysis | Re-run both master CSV and dashboard DB builders; restart Streamlit (cached tables use DB mtime) |

---

## Related documentation

- [casf16_benchmark_analysis_method.md](casf16_benchmark_analysis_method.md) — geometric analyzer, PB-once, extended analysis details
- [casf_dashboard_database_audit.md](casf_dashboard_database_audit.md) — row-type semantics and schema audit notes
- [casf_extended_analysis_requirements.md](casf_extended_analysis_requirements.md) — extended table specifications
- [casf_pb_once_handoff.md](casf_pb_once_handoff.md) — implementation history and Slurm job record
