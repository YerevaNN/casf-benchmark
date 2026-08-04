# CASF Dashboard Database Audit

Dashboard DB: `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite`

## Scope

This audit checks whether dashboard tables mix different row meanings in a way that makes the UI hard to interpret. The main risk was mixing generation-method rows, such as `loqi_raw_dynamic`, with reference-baseline rows, such as `chembl3d_gt`, under the same visible table.

## Source Provenance

- `run_id` and `run_label` come from `config/casf_analysis_sources.yaml`.
- `source` comes from each analyzer output table. It is either a generation method (`*_fixed`, `*_dynamic`, `*_chembl_count`) or a reference baseline (`casf_crystal`, `casf_opt`, `chembl3d_sdf`, `chembl3d_gt`, `chembl3d_gt_pb`).
- `family`, `stage`, and `tier` are parsed from `source` by the analyzer/master builder.
- Reference rows are canonicalized to `reference_core` or `reference_ref`; they are not owned by LOQI, Qwen, DMT, or baseline generation runs.

## Table Audit

Pure generation tables:

- `generation_filter`: 63 rows, all generation rows.
- `generation`: 298 rows, all generation PoseBusters failure rows.

Pure reference tables:

- `reference_filter`: 32 rows, all reference rows.
- `reference`: 47 rows, all reference PoseBusters failure rows.

Mixed aggregate metric tables:

- `cluster`: 63 generation rows and 14 reference-baseline rows.
- `energy`: 63 generation rows and 14 reference-baseline rows.
- `casf_hit`: 63 generation rows and 14 reference-baseline rows.
- `casf_opt_hit`: core-only CASF optimized-ligand hit rows.
- `master_wide`: 63 generation rows and 32 reference rows by design. It is the full export table, not the main section UI.

After canonicalization:

- `master_wide`: 63 generation rows and 9 canonical reference rows.
- `reference_filter`: 9 canonical reference rows.
- `reference`: 14 canonical reference PB-failure rows.
- `cluster`, `energy`, and `casf_hit`: 63 generation rows and 4 canonical reference baseline rows each.
- `casf_opt_hit`: 45 core generation rows and 2 core reference baseline rows. Ref-set rows are excluded because CASF optimized ligand hits are core-only.

## Fixes Applied

The Streamlit dashboard keeps reference baselines visible next to generation rows in comparison tables, because visual side-by-side comparison is the point of those sections.

The mixed comparison tables now show a `row_type` column:

- `generation`: generated conformer-set rows such as `loqi_raw_fixed` or `qwen_4b_revisited_dynamic`.
- `baseline`: reference rows such as `chembl3d_gt` and `chembl3d_gt_pb`.

Reference baselines remain in the same tabs for:

- `Cluster Diversity`
- `Energy`
- `CASF Hits`
- `CASF Opt Hits`

The tier selector (`All`, `Fixed`, `Dynamic`, `ChEMBL count`) filters generation rows while keeping baseline rows visible in comparison tables. Reference-only tabs stay reference-only.

When a generation run is selected in the sidebar (for example `loqi_core`), the dashboard also keeps `reference_core` visible for comparison. This keeps reference baselines tied to the ligand set, not to the generation pipeline.

Display cleanup:

- Streamlit dataframe indexes are hidden.
- Display rows are sorted by `ligand_set`, `run_id`, tier order (`fixed`, `dynamic`, `chembl_count`, `reference`), then `source`.
- CSV downloads use the same sorted rows shown in the UI.

## Validation

The dashboard DB is aggregate-only. It does not include per-molecule `per_ligand_metrics`.

Current aggregate table counts:

- `analysis_sources`: 7
- `master_wide`: 72
- `generation_filter`: 63
- `generation`: 298
- `reference_filter`: 9
- `reference`: 14
- `cluster`: 67
- `energy`: 67
- `casf_hit`: 67
- `casf_opt_hit`: 47

## Second-Round Redundancy Audit

Checks performed:

- Duplicate row keys in every dashboard table.
- Reference rows attached to non-reference runs.
- Repeated `chembl3d_gt` / `chembl3d_gt_pb` baselines.
- Ref-set rows in CASF optimized-ligand hit tables.
- All-null and constant metric columns that would indicate misleading rows.

Findings:

- No duplicate row keys remain.
- No reference rows are attached to LOQI/Qwen/DMT/baseline generation runs.
- `chembl3d_gt` and `chembl3d_gt_pb` appear once per ligand set in `cluster`, `energy`, `casf_hit`, and `casf_opt_hit`.
- `casf_opt_hit` previously contained ref-set rows with all CASF-opt metric columns null; those rows were removed from both the master/dashboard builders and the rebuilt SQLite DB.


