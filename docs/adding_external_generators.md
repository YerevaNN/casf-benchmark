# Adding External Generators

This guide explains how to ingest conformer outputs from your own model (Qwen, a custom checkpoint, etc.) into the CASF benchmark dashboard.

**This repository does not run model inference.** You generate conformers elsewhere, then use the scripts here to materialize tiers, analyze geometry, and rebuild the dashboard.

**On the analysis cluster:** start with [cluster_quickstart.md](cluster_quickstart.md) — it has copy-paste commands for Qwen core/ref on Weka.

Prerequisites: [installation.md](installation.md) (analysis environment), intersection mapping CSVs in `data/mapping/`.

---

## Step 1 — Run inference externally

Use ligands from:

| Cohort | Mapping CSV (repo) | Mapping CSV (Weka) | Ligand MOL2 dir (Weka) |
| --- | --- | --- | --- |
| Core (~94) | `data/mapping/casf16_core_chembl3d_exact_intersection.csv` | `/mnt/weka/mbedrosian/data/casf16/casf16_core_chembl3d_exact_intersection.csv` | `/mnt/weka/mbedrosian/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands` |
| Ref (~1219) | `data/mapping/casf16_ref_chembl3d_exact_intersection.csv` | `/mnt/weka/mbedrosian/data/casf16/casf16_ref_chembl3d_exact_intersection.csv` | `/mnt/weka/mbedrosian/data/casf16/CASF16_REF/ref_chembl3d_exact_intersection_ligands` |

Each row provides `mol_id`, `group`, ChEMBL3D identifiers, and `conformer_count` for tier sizing.

### Recommended Weka run roots

| Cohort | Qwen | Other learned models |
| --- | --- | --- |
| Core | `/mnt/weka/mbedrosian/codex_dir/qwen_gens` | e.g. `.../loqi/generations/casf16_core_loqi_1k` |
| Ref | `/mnt/weka/mbedrosian/codex_dir/qwen/generations/casf16_ref_qwen_1k` | e.g. `.../loqi/generations/casf16_ref_loqi_1k` |

Use **separate run roots** for core vs ref (do not mix cohort ligands in one manifest).

---

## Step 2 — Write outputs in the required layout

Place files under a single **generation root** (call it `{root}`):

```
{root}/
└── generation/
    ├── manifest.tsv
    └── {source_method}/
        ├── {mol_id_1}.sdf
        ├── {mol_id_2}.sdf
        └── ...
```

### `manifest.tsv`

Tab-separated, UTF-8. One row per `(generation_method, mol_id)` for **untiered** source methods only.

Required columns:

| Column | Description |
| --- | --- |
| `mol_id` | CASF ligand identifier (matches mapping CSV) |
| `generation_method` | Untiered method name, e.g. `qwen_4b_bigdata`, `my_model_raw` |
| `status` | `ok` when the SDF is present and valid |
| `set_tier` | Leave empty, or `raw` |

Optional columns: `num_target_confs` (default 1000), `group`, rotatable-bond counts.

**Do not** use tier suffixes (`_fixed`, `_dynamic`, `_chembl_count`) in `generation_method` before materialization.

Example row:

```
mol_id	generation_method	status	set_tier	num_target_confs
1a30	qwen_4b_bigdata	ok	raw	1000
```

### SDF files

- Path: `generation/{source_method}/{mol_id}.sdf`
- Multi-record SDF: one conformer per record, up to 1000 conformers per ligand
- Same heavy-atom topology as the ChEMBL3D reference for that ligand
- 3D coordinates required; hydrogens optional (analyzer strips as needed)

---

## Step 3 — Materialize tiers

Materialization applies PoseBusters-once validation and writes three tier directories:

```
generation/{source_method}_fixed/
generation/{source_method}_dynamic/
generation/{source_method}_chembl_count/
```

### Core example (Weka)

```bash
cd casf-benchmark
export PYTHONPATH=src

python scripts/materialize_casf_generation_sets.py \
  --root /mnt/weka/mbedrosian/codex_dir/qwen_gens \
  --chembl-map-csv /mnt/weka/mbedrosian/data/casf16/casf16_core_chembl3d_exact_intersection.csv \
  --ligand-dir /mnt/weka/mbedrosian/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands \
  --chembl-dataset-root /mnt/weka/mbedrosian/data/chembl3d \
  --method qwen_4b_bigdata
```

### Ref example (Weka)

```bash
python scripts/materialize_casf_generation_sets.py \
  --root /mnt/weka/mbedrosian/codex_dir/qwen/generations/casf16_ref_qwen_1k \
  --chembl-map-csv /mnt/weka/mbedrosian/data/casf16/casf16_ref_chembl3d_exact_intersection.csv \
  --ligand-dir /mnt/weka/mbedrosian/data/casf16/CASF16_REF/ref_chembl3d_exact_intersection_ligands \
  --chembl-dataset-root /mnt/weka/mbedrosian/data/chembl3d \
  --method qwen_4b_bigdata
```

Or use the wrapper script (materialize + analyze):

```bash
./scripts/ingest_external_generation.sh REF \
  /mnt/weka/mbedrosian/codex_dir/qwen/generations/casf16_ref_qwen_1k \
  qwen_4b_bigdata
```

See [casf16_materialize_generation_sets_method.md](casf16_materialize_generation_sets_method.md) for full options.

---

## Step 4 — Register in config

### Existing run root (new checkpoint only)

Add a family entry to `config/casf_generation_families.yaml` if the `method_prefix` is new.  
`config/casf_analysis_sources.weka.yaml` already lists standard run roots (`qwen_core`, `qwen_ref`, etc.) — no edit needed unless you use a **new directory**.

### New run root

Add a source block to `config/casf_analysis_sources.weka.yaml` (cluster) and optionally `config/casf_analysis_sources.yaml` (offline bundle):

```yaml
  - run_id: my_model_ref
    label: "My model (ref)"
    ligand_set: ref
    root: /mnt/weka/mbedrosian/codex_dir/my_model/generations/casf16_ref_my_model_1k
```

For bundled/offline rebuilds, copy `analysis/tables/geometric_per_ligand_long.csv` into `data/results/runs/my_model_ref/analysis/tables/` after analysis.

---

## Step 5 — Geometric analysis

```bash
export PYTHONPATH=src

# Ref cohort example
python scripts/analyze_casf_conformer_sets.py \
  --generation-only \
  --ligand-set ref \
  --output-dir /mnt/weka/mbedrosian/codex_dir/qwen/generations/casf16_ref_qwen_1k \
  --casf-ligand-dir /mnt/weka/mbedrosian/data/casf16/CASF16_REF/ref_chembl3d_exact_intersection_ligands \
  --chembl-map-csv /mnt/weka/mbedrosian/data/casf16/casf16_ref_chembl3d_exact_intersection.csv \
  --chembl-dataset-root /mnt/weka/mbedrosian/data/chembl3d \
  --workers 48 \
  --resume-parts \
  --quiet-rdkit-warnings
```

(`--output-dir` and `--run-root` are equivalent aliases.)

This writes `{root}/analysis/tables/geometric_per_ligand_long.csv`.

---

## Step 6 — Rebuild dashboard

**On Weka (production):**

```bash
./scripts/rebuild_dashboard_weka.sh
```

**Offline / bundled clone:**

```bash
python scripts/build_casf_analysis_master_csv.py
python scripts/build_casf_analysis_dashboard_db.py
streamlit run apps/dashboard/streamlit_app.py
```

---

## Checklist

- [ ] Untiered SDFs + `manifest.tsv` under `{root}/generation/`
- [ ] Method name registered in `casf_generation_families.yaml` (if new checkpoint)
- [ ] Materialization completed (three tier suffixes present)
- [ ] Geometric analysis completed (`analysis/tables/geometric_per_ligand_long.csv` exists)
- [ ] `casf_analysis_sources.weka.yaml` includes the run root (if new directory)
- [ ] `./scripts/rebuild_dashboard_weka.sh` completed
- [ ] (Optional) Analysis CSV copied to `data/results/runs/` for git bundle

---

## Reference method names (existing catalog)

| Model | Example `source_method` |
| --- | --- |
| LOQI | `loqi_raw` |
| NExT-Mol DMT-L | `nextmol_dmt_l_raw` |
| Torsional Diffusion | `torsional_diffusion_raw` |
| MCF drugs-L | `mcf_drugs_l_raw` |
| Qwen | `qwen_4b_bigdata`, `qwen_4b_revisited`, … |

Pick a unique `{source_method}` string for new checkpoints.
