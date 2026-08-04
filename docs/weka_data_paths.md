# Weka production data paths

This document lists where CASF benchmark inputs and **pre-generated conformer SDF pools** live on the YerevaNN analysis cluster (Weka shared storage). These paths are **not** in the git clone; the repo ships precomputed **analysis tables** and dashboard SQLite under `data/results/` instead.

Set the root once:

```bash
export CASF_BENCHMARK_DATA_ROOT=/mnt/weka/mbedrosian
```

Canonical constants are defined in `src/casf_benchmark/paths.py`. Production run roots for dashboard rebuilds are in `config/casf_analysis_sources.weka.yaml`.

---

## Reference inputs (MOL2, ChEMBL3D, mapping CSVs)

| Role | Weka path |
| --- | --- |
| CASF-2016 core ligands | `/mnt/weka/mbedrosian/data/casf16/CASF16/ligands` |
| CASF-2016 ref ligands | `/mnt/weka/mbedrosian/data/casf16/CASF16_REF/ligands` |
| Core intersection MOL2 panel | `/mnt/weka/mbedrosian/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands` |
| Ref intersection MOL2 panel | `/mnt/weka/mbedrosian/data/casf16/CASF16_REF/ref_chembl3d_exact_intersection_ligands` |
| CASF optimized ligands | `/mnt/weka/mbedrosian/data/casf16/CASF16/ligands_opt` |
| Core intersection mapping CSV | `/mnt/weka/mbedrosian/data/casf16/casf16_core_chembl3d_exact_intersection.csv` |
| Ref intersection mapping CSV | `/mnt/weka/mbedrosian/data/casf16/casf16_ref_chembl3d_exact_intersection.csv` |
| ChEMBL3D topology SDF shards | `/mnt/weka/mbedrosian/data/chembl3d/topologies/{group}.sdf` |
| ChEMBL3D conformer zarr | `/mnt/weka/mbedrosian/data/chembl3d/zarr_database` |
| ChEMBL3D SMILES index | `/mnt/weka/mbedrosian/data/chembl3d_index/chembl3d_topology_smiles_index.csv` |

Bundled copies of the mapping CSVs (for offline work) are under `data/mapping/` in the repo.

---

## Pre-generated conformer SDF pools

Each run root follows:

```
{run_root}/
├── generation/
│   ├── manifest.tsv
│   └── {method}_{fixed|dynamic|chembl_count}/
│       └── {mol_id}.sdf
└── analysis/
    └── tables/
        └── geometric_per_ligand_long.csv
```

SDFs and manifests live under `generation/`. The git clone includes only `analysis/tables/` for dashboard rebuilds.

### RDKit random and torsion baselines (in-repo generator)

Produced by `src/casf_benchmark/generation/conformer_sets.py` (or `python -m casf_benchmark.generation.conformer_sets`).

| Cohort | Run root (SDFs + manifest) |
| --- | --- |
| Core (~94 ligands) | `/mnt/weka/mbedrosian/pharma_generation_analysis/core_pb_full_dynamic_chembl_count` |
| Ref (~1219 ligands) | `/mnt/weka/mbedrosian/pharma_generation_analysis/ref_pb_full_dynamic_chembl_count` |

Example SDF path:

```
/mnt/weka/mbedrosian/pharma_generation_analysis/core_pb_full_dynamic_chembl_count/generation/rdkit_random_raw_fixed/1bcu_1bcu_conf0.sdf
```

### Reference-only baselines (crystal, opt, ChEMBL3D)

| Cohort | Run root |
| --- | --- |
| Core | `/mnt/weka/mbedrosian/pharma_generation_analysis/reference_datasets/core` |
| Ref | `/mnt/weka/mbedrosian/pharma_generation_analysis/reference_datasets/ref` |

### External / learned generators (materialized tiers)

Raw inference outputs are materialized with `scripts/materialize_casf_generation_sets.py` into tiered `{method}_{fixed|dynamic|chembl_count}/` directories.

| Generator | Core run root | Ref run root |
| --- | --- | --- |
| Qwen (all checkpoints) | `/mnt/weka/mbedrosian/codex_dir/qwen_gens` | `/mnt/weka/mbedrosian/codex_dir/qwen/generations/casf16_ref_qwen_1k` |
| LOQI | `/mnt/weka/mbedrosian/codex_dir/loqi/generations/casf16_core_loqi_1k` | `/mnt/weka/mbedrosian/codex_dir/loqi/generations/casf16_ref_loqi_1k` |
| NExT-Mol DMT-L | `/mnt/weka/mbedrosian/codex_dir/nextmol_dmt_l/generations/casf16_core_nextmol_dmt_l_1k` | `/mnt/weka/mbedrosian/codex_dir/nextmol_dmt_l/generations/casf16_ref_nextmol_dmt_l_1k` |
| Torsional Diffusion | `/mnt/weka/mbedrosian/codex_dir/torsional_diffusion/generations/casf16_core_torsional_diffusion_1k` | `/mnt/weka/mbedrosian/codex_dir/torsional_diffusion/generations/casf16_ref_torsional_diffusion_1k` |
| MCF drugs-L | `/mnt/weka/mbedrosian/codex_dir/mcf_drugs_l/generations/casf16_core_mcf_drugs_l_1k` | `/mnt/weka/mbedrosian/codex_dir/mcf_drugs_l/generations/casf16_ref_mcf_drugs_l_1k` |

Example LOQI SDF path:

```
/mnt/weka/mbedrosian/codex_dir/loqi/generations/casf16_core_loqi_1k/generation/loqi_raw_fixed/1bcu_1bcu_conf0.sdf
```

Example Qwen ref SDF path:

```
/mnt/weka/mbedrosian/codex_dir/qwen/generations/casf16_ref_qwen_1k/generation/qwen_4b_bigdata_fixed/1ado_1ado_conf0.sdf
```

---

## Aggregated dashboard outputs (Weka production build)

| Artifact | Weka path |
| --- | --- |
| Dashboard SQLite | `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite` |
| Master CSV | `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_master.csv` |
| Global per-ligand long CSV | `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_per_ligand_long.csv` |
| Extended analysis outputs | `/mnt/weka/mbedrosian/pharma_generation_analysis/extended_analysis/` |

---

## Using Weka paths from the repo

**Rebuild dashboard from Weka run roots:**

```bash
export CASF_BENCHMARK_DATA_ROOT=/mnt/weka/mbedrosian
export PYTHONPATH=src
python scripts/build_casf_analysis_master_csv.py --config config/casf_analysis_sources.weka.yaml
python scripts/build_casf_analysis_dashboard_db.py --config config/casf_analysis_sources.weka.yaml
```

**Re-run geometric analysis** (needs SDFs + MOL2 on Weka):

```bash
python scripts/analyze_casf_conformer_sets.py \
  --run-root /mnt/weka/mbedrosian/codex_dir/loqi/generations/casf16_core_loqi_1k \
  --chembl-map-csv data/mapping/casf16_core_chembl3d_exact_intersection.csv \
  --ligand-dir /mnt/weka/mbedrosian/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands
```

**Extended analysis K-efficiency curves** read per-ligand SDFs from the run roots above. On a clone without SDFs, use the bundled `data/results/extended_casf_analysis.sqlite` for the dashboard Extended Analysis tab, or run on the cluster with `--skip-k-efficiency` when SDFs are unavailable.

See also: [installation.md](installation.md), [casf16_generation_pipeline_and_models_method.md](casf16_generation_pipeline_and_models_method.md).
