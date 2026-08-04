# Adding External Generators

This guide explains how to ingest conformer outputs from your own model (Qwen, a custom checkpoint, etc.) into the CASF benchmark dashboard.

**This repository does not run model inference.** You generate conformers elsewhere, then use the scripts here to materialize tiers, analyze geometry, and rebuild the dashboard.

Prerequisites: [installation.md](installation.md) (analysis environment), intersection mapping CSVs in `data/mapping/`.

---

## Step 1 — Run inference externally

Use ligands from:

- `data/mapping/casf16_core_chembl3d_exact_intersection.csv` (core, ~94 ligands)
- `data/mapping/casf16_ref_chembl3d_exact_intersection.csv` (ref, ~1219 ligands)

Each row provides `mol_id`, `group`, ChEMBL3D identifiers, and `conformer_count` for tier sizing.

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

```bash
cd casf-benchmark
export PYTHONPATH=src

python scripts/materialize_casf_generation_sets.py \
  --root /path/to/your/root \
  --chembl-map-csv data/mapping/casf16_core_chembl3d_exact_intersection.csv \
  --ligand-dir /path/to/core_chembl3d_exact_intersection_ligands \
  --chembl-dataset-root /path/to/chembl3d \
  --method qwen_4b_bigdata
```

On a cluster with Weka data:

```bash
--ligand-dir /mnt/weka/mbedrosian/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands \
--chembl-dataset-root /mnt/weka/mbedrosian/data/chembl3d
```

See [casf16_materialize_generation_sets_method.md](casf16_materialize_generation_sets_method.md) for full options.

---

## Step 4 — Register in config

Add a family entry to `config/casf_generation_families.yaml` if needed, then add a source block to `config/casf_analysis_sources.yaml`:

```yaml
  - run_id: my_model_core
    label: "My model (core)"
    ligand_set: core
    root: /path/to/your/root
```

For bundled/offline rebuilds, copy `analysis/tables/geometric_per_ligand_long.csv` into `data/results/runs/my_model_core/analysis/tables/` after analysis.

---

## Step 5 — Geometric analysis

```bash
python scripts/analyze_casf_conformer_sets.py \
  --run-root /path/to/your/root \
  --chembl-map-csv data/mapping/casf16_core_chembl3d_exact_intersection.csv \
  --ligand-dir /path/to/ligands \
  --reference-mode pb_once
```

This writes `analysis/tables/geometric_per_ligand_long.csv` under the run root.

---

## Step 6 — Rebuild dashboard

```bash
python scripts/build_casf_analysis_master_csv.py
python scripts/build_casf_analysis_dashboard_db.py
streamlit run apps/dashboard/streamlit_app.py
```

---

## Checklist

- [ ] Untiered SDFs + `manifest.tsv` under `generation/`
- [ ] Method name registered in `casf_generation_families.yaml`
- [ ] Materialization completed (three tier suffixes present)
- [ ] Geometric analysis completed (`geometric_per_ligand_long.csv` exists)
- [ ] `casf_analysis_sources.yaml` updated
- [ ] Master CSV and dashboard SQLite rebuilt

---

## Reference method names (existing catalog)

| Model | Example `source_method` |
| --- | --- |
| LOQI | `loqi_raw` |
| NExT-Mol DMT-L | `nextmol_dmt_l_raw` |
| Torsional Diffusion | `torsional_diffusion_raw` |
| MCF drugs-L | `mcf_drugs_l_raw` |
| Qwen | `qwen_4b_bigdata`, `qwen_8b_bigdata`, … |

Pick a unique `{source_method}` string for new checkpoints.
