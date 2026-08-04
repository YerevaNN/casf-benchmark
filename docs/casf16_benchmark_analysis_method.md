# CASF Benchmark Analysis Pipeline

## Purpose

This document describes the analysis pipeline that evaluates CASF–ChEMBL3D conformer generation outputs against crystal and optimized bound poses. The pipeline proceeds in three stages: per-run geometric analysis (`scripts/analyze_casf_conformer_sets.py`), cross-run aggregation into a dashboard database (`scripts/build_casf_analysis_master_csv.py` and `scripts/build_casf_analysis_dashboard_db.py`), and extended publication-oriented analyses (`scripts/extended_casf_analysis.py`).

Stage one computes per-ligand geometric recovery, diversity, clash, PoseBusters, and energy metrics from generation SDFs and reference structures. Stage two merges long-form metrics from many generator cohorts into a single SQLite database. Stage three performs deterministic paired comparisons, K-efficiency curves, Pareto frontier ranking, and failure-mechanism summaries without modifying the source database.

The geometric analyzer implements a PoseBusters-once strategy: generated conformers inherit PoseBusters labels from the generation manifest and are not re-validated at analysis time. Reference baselines are analyzed in separate reference-only runs where PoseBusters is executed against the CASF crystal ligand.

## Pipeline Overview

1. **Geometric analysis** — For each conformer-set run root, compute per-ligand metrics and cohort summary tables; write `geometric_per_ligand_long.csv` and `geometric_report.md` under `{run_root}/analysis/`.
2. **Dashboard build** — Aggregate long tables from all configured run roots into `casf_per_ligand_long.csv`, `casf_analysis_master.csv`, and `casf_analysis_dashboard.sqlite`.
3. **Extended analysis** — Read the dashboard database; emit paired deltas, K-efficiency curves, frontier summaries, and supplementary tables, reports, and figures under `extended_analysis/`.

Upstream inputs (intersection mapping, conformer generation, manifest production) are documented in [casf16_chembl3d_exact_match_method.md](casf16_chembl3d_exact_match_method.md) (data preparation), [casf16_chembl3d_conformer_generation_method.md](casf16_chembl3d_conformer_generation_method.md), [casf16_materialize_generation_sets_method.md](casf16_materialize_generation_sets_method.md), and [casf16_generation_pipeline_and_models_method.md](casf16_generation_pipeline_and_models_method.md) (generator models and reference datasets).

---

## Stage 1. Per-Run Geometric Analysis

Implementation: `scripts/analyze_casf_conformer_sets.py`.

### Input data

**Generation manifest and SDF outputs.** Generation-mode analysis requires the merged tab-separated manifest at `{output_dir}/generation/manifest.tsv`, or the union of parts under `generation/manifest_parts/` (concatenated on first load if no merged file exists). For each catalog method and eligible ligand, PoseBusters-passing conformers are read from `{generation_dir}/{method}/{mol_id}.sdf`. Empty SDF files are valid: they indicate a selected pool with zero post-PoseBusters survivors.

Default core generation root: `/mnt/weka/mbedrosian/pharma_generation_analysis/core_pb_full_dynamic_chembl_count`.

Pre-generated SDF pools for all dashboard generators (RDKit/torsion, LOQI, Qwen, etc.) are listed in [weka_data_paths.md](weka_data_paths.md). The git clone ships analysis CSVs only under `data/results/runs/`.

**Intersection mapping table.** Ligand eligibility follows the ChEMBL intersection CSV (`ligand_id`, `source_file`, ChEMBL3D identifiers, conformer counts). Chemical filters are not re-applied at analysis time.

Default path: `/mnt/weka/mbedrosian/data/casf16/casf16_core_chembl3d_exact_intersection.csv`.

**CASF ligand structures.** Crystal poses load from MOL2 files, searching manifest `source_input`, the configured ligand directory, the core intersection directory, and default CASF ligand paths. Optimized poses load from a separate directory when present; missing optimized structures omit optimized metrics without failing the run.

Defaults: crystal `/mnt/weka/mbedrosian/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands`; optimized `/mnt/weka/mbedrosian/data/casf16/CASF16/ligands_opt`.

**ChEMBL3D ground truth (reference mode).** Conformers load from topology SDF and zarr archive at `/mnt/weka/mbedrosian/data/chembl3d/topologies` and `/mnt/weka/mbedrosian/data/chembl3d/zarr_database`. Ligand `1tlp_1tlp_conf0` is deterministically subsampled to 2000 conformers before reference PoseBusters and metrics.

### Analysis modes

**Generation mode (default).** Discovers catalog generation methods in the manifest, evaluates SDF contents geometrically, and copies manifest `pb_*` fields without re-running PoseBusters on generation structures.

**Reference mode (`--reference-only`).** Produces five reference rows per mapped ligand: `casf_crystal`, `casf_opt`, `chembl3d_sdf`, `chembl3d_gt`, and `chembl3d_gt_pb`. PoseBusters runs at analysis time except for `chembl3d_gt_pb`, which filters the full ChEMBL3D PoseBusters pass mask. Requires `--ligand-set core` or `--ligand-set ref`.

`--reference-only` and `--generation-only` are mutually exclusive; `--generation-only` is an accepted alias for the default generation mode.

### Ligand selection and checks

Generation mode selects intersection CSV ligands present in the manifest with an SDF file for every discovered catalog method. Exclusions: `not_in_manifest`, `missing_generation_sdf`. Reference mode selects mapped ligands with a MOL2 in the ligand directory. Catalog methods are the ordered intersection of manifest methods and the generation catalog for the ligand set. Duplicate `ligand_id` rows in the mapping CSV are rejected. Post-computation validation ensures complete source coverage; reference mode requires successful ChEMBL3D loading for mapped ligands.

### Per-ligand metrics

Each ligand is processed independently; results cache as pickle files under `analysis/cache/geometric_generation_parts/` or `analysis/cache/geometric_reference_parts/` (reused with `--resume-parts`).

For every source and conformer list:

- **Steric clash** at DG cutoff 0.7 (bonds and valence angles excluded; hydrogens ignored).
- **PoseBusters summary** — copied from manifest for generation sources; run at analysis time for reference sources; filtered pass mask for `chembl3d_gt_pb`.
- **Conformational diversity** — mean torsion standard deviation, pairwise RMSD stats, greedy clustering at 0.5, 1.0, 2.0, and 3.0 Å with normalized cluster descriptors.
- **Energies** — per-conformer MMFF94s values in `geometric_energy_values.csv`.
- **CASF recovery** (comparison sources) — heavy-atom best-aligned RMSD to crystal and optimized references; best and median RMSD; hit indicators at 0.25, 0.5, 0.75, and 2.0 Å. Non-finite RMSD raises an error with atom-count diagnostics.

Generation workflow attaches manifest funnel fields per method. Reference workflow emits the five reference source rows described above; ChEMBL3D appears both unfiltered (`chembl3d_gt`) and PoseBusters-filtered (`chembl3d_gt_pb`) for fair comparison to filtered generation sets.

### Per-run outputs

Under `{run_root}/analysis/`:

- `tables/geometric_per_ligand_long.csv` — primary long-form metrics (alias: `geometric_per_ligand_metrics.csv`)
- `tables/geometric_generation_filter_summary.csv`, `geometric_generation_pb_failure_summary.csv`
- `tables/geometric_reference_filter_summary.csv`, `geometric_reference_pb_failure_summary.csv` (reference runs)
- `tables/geometric_cluster_summary.csv`, `geometric_energy_summary.csv`
- `tables/geometric_casf_hit_summary.csv`, `geometric_casf_opt_hit_summary.csv`, `geometric_energy_values.csv`
- `geometric_report.md` — consolidated Markdown report
- `cache/geometric_*_parts/` — per-ligand pickle caches

Cohort summary rates use per-ligand averaging before aggregation across ligands.

### Geometric analysis execution

Dependencies: RDKit, NumPy, pandas, PoseBusters (reference mode), `casf_benchmark.analysis.metrics`, `casf_benchmark.paths`, `casf_benchmark.catalog`. Set `PYTHONPATH` to include `src/`.

Key flags: `--output-dir`, `--casf-ligand-dir`, `--casf-opt-ligand-dir`, `--chembl-map-csv`, `--chembl-dataset-root`, `--ligand-set`, `--reference-only`, `--limit-molecules`, `--molecule-offset`, `--report-only`, `--resume-parts`, `--workers`, `--posebusters-workers`, `--posebusters-energy-threads`, `--quiet-rdkit-warnings`.

Production generation analysis: run with run root, ligand directory, and intersection CSV; add `--resume-parts`, `--workers 48`, and `--quiet-rdkit-warnings`. Reference analysis: add `--reference-only --ligand-set core` with `--workers 24`.

Slurm: per-run generation via `scripts/submit_casf.sh analyze {core|ref}` (`run_casf_ref_geometric_analysis.sbatch`, default 48 CPUs, 192 GB, 14 days). Full multi-cohort rebuild via `scripts/run_pb_once_casf_analysis.sbatch` (reference core/ref, then generation across all configured roots, then master CSV and dashboard DB builders).

---

## Stage 2. Dashboard Aggregation

Implementation: `scripts/build_casf_analysis_master_csv.py` and `scripts/build_casf_analysis_dashboard_db.py`.

These scripts are not invoked by the geometric or extended analyzers directly but form the required bridge between them. The master builder reads `geometric_per_ligand_long.csv` from each configured run root (defined in `config/casf_analysis_sources.yaml`), concatenates rows, and writes:

- `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_per_ligand_long.csv`
- `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_master.csv`

The dashboard builder loads the global long table and derives comparison rows and stratum breakdowns into SQLite:

- `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite`

Primary tables: `per_ligand_long`, `comparison_rows`, `comparison_strata`, `analysis_sources`, `catalog_families`.

The Streamlit dashboard is documented in [casf16_dashboard_method.md](casf16_dashboard_method.md). Implementation: `apps/dashboard/streamlit_app.py`.

---

## Stage 3. Extended Analysis

Implementation: `scripts/extended_casf_analysis.py`.

### Input data

**Primary database (read-only).** Loads `per_ligand_long`, `comparison_rows`, and `analysis_sources` from the dashboard SQLite file. Default: `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite`.

**Secondary filesystem inputs (K-efficiency only).** When K-efficiency curves are computed, the script may re-read generation SDFs, CASF MOL2 files, ChEMBL3D topology/zarr data, and per-run `geometric_energy_values.csv`. All other modules use database tables exclusively.

### Analysis modules

The script runs six modules in fixed order, then exports tables to CSV, Markdown, PNG, and SQLite.

**1. Sanity checks.** Validates required columns, duplicate keys, RMSD and hit value ranges, PoseBusters count accounting, and metadata warnings. Output: `extended_sanity_checks.csv`. Use `--fail-on-sanity-error` to abort on errors. Ligand-universe audit (Analysis 1) is skipped by default and not implemented here.

**2. K-efficiency curves.** Quantifies CASF crystal recovery as a function of sample size K from PoseBusters-passing pools. Eligible fixed-tier generation sources include RDKit, torsion, LOQI, torsional diffusion, NextMol, MCF-Drugs, and Qwen variants. Per-conformer CASF RMSD vectors cache under `cache/per_conformer_casf_rmsd_parts/`. For each K in the default grid (1, 2, 5, 10, 25, 50, 100, 250, 500, 1000), `--k-subsamples` (default 20) deterministic draws estimate mean best RMSD and hit rates at 0.5, 0.75, and 2.0 Å. Stratum aggregation uses rotatable-bond and heavy-atom bins. ChEMBL3D-PB K rows are explicitly skipped because the dashboard stores only per-ligand best RMSD, not per-conformer indices. Outputs include `extended_k_efficiency.csv`, stratum table, manifest, report, and PNG plots. Bypass with `--skip-k-efficiency`; force recompute with `--recompute-k-rmsd`.

**3. Paired delta comparisons.** Inner-joins each fixed-tier method and `chembl3d_gt` against baselines `chembl3d_gt_pb` and `rdkit_random_raw_fixed`. Reports common ligand counts, delta RMSD statistics, win/loss/tie rates, delta hit rates, and extreme rescued/worsened ligand lists. Stratified by rotatable bonds, heavy atoms, ChEMBL conformer count, and baseline difficulty bins.

**4. Diversity–validity–recovery frontier.** Computes clusters per 100 conformers, Hit@0.75 per 100 conformers, and valid hit score; assigns Pareto frontier rank maximizing Hit@0.75 and PB pass rate while minimizing conformer count.

**5. PoseBusters failure mechanisms.** Aggregates `pb_check_fail_counts_json` across ligands; identifies dominant failure tests and worst ligand tails. No manual chirality review.

**6. Bound-pose difficulty and rescue.** Stratifies by ChEMBL3D-PB baseline difficulty; reports rescue and catastrophic-loss rates; collects rescue cases for selected high-interest generators.

### Extended outputs

Default root: `/mnt/weka/mbedrosian/pharma_generation_analysis/extended_analysis`.

- `tables/` — all extended CSV tables prefixed `extended_*`
- `reports/` — Markdown summaries per module
- `figures/` — K-efficiency and frontier PNG plots
- `cache/per_conformer_casf_rmsd_parts/` — K-efficiency RMSD cache
- `extended_casf_analysis.sqlite` — sidecar database of extended tables only
- `casf_analysis_dashboard_extended.sqlite` — dashboard copy with extended tables appended
- `README.md` — reproduction command

### Extended analysis execution

Dependencies: pandas, NumPy, matplotlib, sqlite3; K-efficiency additionally requires RDKit and `analyze_casf_conformer_sets` helpers. Set `PYTHONPATH=src`.

Key flags: `--dashboard-db`, `--output-dir`, `--k-values`, `--k-subsamples`, `--k-seed`, `--k-workers`, `--recompute-k-rmsd`, `--skip-k-efficiency`, `--fail-on-sanity-error`, `--skip-analysis-1`.

Default invocation: `PYTHONPATH=src python scripts/extended_casf_analysis.py`. Rebuild K-efficiency: add `--recompute-k-rmsd --k-workers 8`.

Preconditions: populated dashboard database; for K-efficiency, reachable SDF and MOL2 paths from `analysis_sources`. Runtime is minutes without K-efficiency; full K-efficiency can require hours.

---

## Design Rationale

**Manifest-trusted PoseBusters for generation.** Avoids duplicating PoseBusters cost at analysis time while evaluating the same conformers written to SDF files.

**Separate reference runs.** Reference baselines compute once per ligand set and serve all generator comparisons with consistent ChEMBL3D loading and PoseBusters framing.

**Fair ChEMBL3D comparison via `chembl3d_gt_pb`.** PoseBusters-filtered ChEMBL3D baselines align validity assumptions with filtered generation sets.

**Per-ligand rate aggregation.** Cohort summaries average ligand-level rates before pooling, preventing large pools from dominating means.

**Read-only extended analysis.** Additive outputs preserve the canonical dashboard during exploratory and publication analyses.

**Ligand-paired extended comparisons.** Delta and rescue metrics require inner joins on `mol_id` within ligand set.

**Deterministic K-efficiency.** Fixed seeds and cached RMSD vectors yield reproducible curves without bootstrap resampling.

**Explicit ChEMBL3D-PB K skip.** The extended script documents rather than approximates K curves from aggregate best-RMSD rows.

**Fixed-tier focus in extended comparisons.** Paired deltas and frontier ranks use fixed-tier rows for comparable sample budgets.

---

## Related Documentation

- [casf16_chembl3d_exact_match_method.md](casf16_chembl3d_exact_match_method.md) — data preparation and intersection mapping
- [casf16_chembl3d_conformer_generation_method.md](casf16_chembl3d_conformer_generation_method.md) — conformer generation and manifests
- [casf16_materialize_generation_sets_method.md](casf16_materialize_generation_sets_method.md) — external ML tier materialization
- [casf16_dashboard_method.md](casf16_dashboard_method.md) — Streamlit dashboard and rebuild commands
- [casf_pb_once_handoff.md](casf_pb_once_handoff.md) — PB-once implementation and Slurm run record
- [casf_extended_analysis_requirements.md](casf_extended_analysis_requirements.md) — extended analysis requirements
- [casf_dashboard_database_audit.md](casf_dashboard_database_audit.md) — dashboard schema
- [casf_geometric_report_core_94_ligands.md](casf_geometric_report_core_94_ligands.md) — example geometric report (core cohort)
- [casf_diverse_conformer_hypothesis_analysis.md](casf_diverse_conformer_hypothesis_analysis.md) — narrative interpretation
