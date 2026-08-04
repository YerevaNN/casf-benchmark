# PB-once CASF implementation handoff

## Objective

Implement the PB-once CASF analysis strategy from the attached spec:

- Generation analysis trusts generation-time `pb_*` manifest columns and never runs PoseBusters on generation SDFs.
- Reference rows are computed once in separate reference-only runs.
- Per-ligand long CSVs are the source of truth for master/dashboard aggregation.
- The canonical catalog has 19 generation families and 3 tiers.

## Progress Log

- Started from a dirty worktree; existing changes included CASF scripts, dashboard files, tests, docs, `pyproject.toml`, and untracked app/config/test files. Treat these as user changes.
- Read the full attached spec from `/home/mbedrosian/.codex/attachments/8d162f60-9dc5-4bd6-85ee-8ddcbbd84e26/pasted-text.txt`.
- Inspected current analyzer, master builder, dashboard DB builder, Streamlit app, layout helper, and CASF analysis tests.
- Added canonical catalog files:
  - `config/casf_generation_families.yaml`
  - `src/casf_benchmark/casf_generation_catalog.py`
  - `src/casf_benchmark/casf_stratum_bins.py`
- Refactored `scripts/analyze_casf_conformer_sets.py`:
  - `discover_generation_methods()` now intersects manifest methods with the catalog and ignores untiered Qwen rows.
  - Generation work uses `analyze_generation_ligand()` and copies `pb_*` stats from the manifest.
  - Reference work uses `analyze_reference_ligand()` and keeps PoseBusters there.
  - Added `--reference-only`, `--generation-only`, and `--ligand-set`.
  - Removed deprecated attempted CSV writes and removed deprecated `--use-posebusters` / `--geometry-only` CLI flags.
  - Writes `geometric_per_ligand_long.csv` and keeps `geometric_per_ligand_metrics.csv` as a deprecated alias.
- Rewrote `config/casf_analysis_sources.yaml` with explicit non-baseline roots and `reference_sources`.
- Refactored `src/casf_benchmark/casf_analysis_dataset.py` so `build_master_frame()` reads per-root long tables and derives comparison rows by aggregation.
- Updated `scripts/build_casf_analysis_master_csv.py` to write/report the global long table and master CSV.
- Replaced `scripts/build_casf_analysis_dashboard_db.py` with the new SQLite schema:
  - `comparison_rows`
  - `comparison_strata`
  - `per_ligand_long`
  - `catalog_families`
  - `PRAGMA user_version = 2`
- Replaced `apps/dashboard/streamlit_app.py` with a unified tier/breakdown dashboard over `comparison_rows` and `comparison_strata`.
- Updated focused tests for catalog parsing, PB-once generation behavior, long-table master aggregation, dashboard strata counts, and the renamed generation cache parts directory.

## Validation Completed

- `PYTHONPATH=src:. uv run --with pytest pytest tests/pharmacophore/test_analyze_casf_conformer_sets.py tests/pharmacophore/test_casf_analysis_dataset.py tests/pharmacophore/test_conformer_sets_layout.py -q`
  - Result: 25 passed, 1 skipped.
- `git diff --check`
  - Result: passed.
- `python -m py_compile scripts/analyze_casf_conformer_sets.py scripts/build_casf_analysis_master_csv.py scripts/build_casf_analysis_dashboard_db.py apps/dashboard/streamlit_app.py src/casf_benchmark/casf_generation_catalog.py src/casf_benchmark/casf_analysis_dataset.py src/casf_benchmark/casf_stratum_bins.py src/casf_benchmark/paths.py`
  - Result: passed.
- `PYTHONPATH=src:. uv run --with pytest pytest -q`
  - Result: failed during test collection on unrelated environment/fixture issues:
    - `PermissionError` reading `/data/molgen/3DMolGen_data/pretokenized_prompts.json`.
    - Missing fixture `/home/mbedrosian/code/casf-benchmark/tests/evaluation/sample_by_smiles_df.csv`.
    - Missing module `molgen3D.pharmacophore.generate_chembl3d_etkdg_conformers`.

## Slurm Run Record

- Production analysis was moved to an independent Slurm job script:
  - `scripts/run_pb_once_casf_analysis.sbatch`
  - Submit command: `sbatch scripts/run_pb_once_casf_analysis.sbatch`
  - Python: `/home/mbedrosian/.conda/envs/chembl3d/bin/python`
  - Workers: `REFERENCE_WORKERS=24`, `GENERATION_WORKERS=48`
  - Slurm limits: `14-00:00:00`, `192G`, `48` CPUs (sized for heavy ref ligands like `1tlp` with 10951 ChEMBL3D conformers vs ref median 41 / mean 96)
- Reference-only `ref` uses the full 1236-ligand ref map. Ligand `1tlp_1tlp_conf0` (10951 ChEMBL3D conformers) is subsampled to 2000 deterministic random conformers before reference PoseBusters/metrics.
- Generation `ref` uses the full `casf16_ref_chembl3d_exact_intersection.csv` map (1236 ligands).
- Job history:
  - `206639`: first Slurm run, cancelled after stalling in `reference-ref`.
  - `206687`: second Slurm run, cancelled after stalling in `generation-rdkit_torsion_ref`.
  - `206814`: third Slurm run, cancelled after all 1224 selected rdkit/torsion ref parts were written but before summary completed; this exposed a generation resume bug.
  - `206819`: final Slurm run, completed successfully.
- Final Slurm status:
  - `sacct -j 206819 --format=JobID,JobName%32,State,ExitCode,Elapsed,MaxRSS -P`
  - Result: `206819|casf-pb-once|COMPLETED|0:0|02:25:01|`; batch max RSS was `24427984K`.
- Final output artifacts:
  - `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_per_ligand_long.csv`
    - 41,403 rows, 58M
  - `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_master.csv`
    - 91 rows, 104K
  - `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite`
    - 51M
- Final dashboard DB counts from Python `sqlite3`:
  - `analysis_sources`: 13
  - `catalog_families`: 19
  - `comparison_rows`: 91
  - `comparison_strata`: 178
  - `per_ligand_long`: 41403
- Final focused validation after the Slurm run:
  - `PYTHONPATH=src:. uv run --with pytest pytest tests/pharmacophore/test_analyze_casf_conformer_sets.py tests/pharmacophore/test_casf_analysis_dataset.py tests/pharmacophore/test_conformer_sets_layout.py -q`
  - Result: 25 passed, 1 skipped.
- Important implementation note from the Slurm run:
  - Generation resume checks must not require ChEMBL reference rows; generation part files contain generation sources only.
  - `scripts/analyze_casf_conformer_sets.py` was patched so `require_chembl` is only applied when `mode == "reference"`.
- Post-run dashboard fix:
  - The first populated dashboard DB had breakdown aggregates only for rdkit/torsion families.
  - Root cause: `assign_stratum()` rebuilt input values as a fresh `0..n` indexed pandas `Series`; assigning that result back to a filtered long-table slice aligned by index and dropped strata for later concatenated roots.
  - Fix: `src/casf_benchmark/casf_stratum_bins.py` now preserves an incoming `Series` index.
  - `scripts/build_casf_analysis_dashboard_db.py` validation now checks expected rows against strata rows, so missing family/tier strata fail fast.
  - Rebuilt `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite`; `comparison_strata` increased from 178 to 750 rows.
  - Verified strata coverage:
    - Core: 19 generation families for each of `fixed`, `dynamic`, `chembl_count`.
    - Ref: 8 generation families for each of `fixed`, `dynamic`, `chembl_count`.
    - No missing family/tier rows relative to `comparison_rows`.
  - Regression test added: `test_assign_stratum_preserves_series_index`.

## Current State

- The dashboard has been populated from the completed Slurm run.
- Focused CASF validation passes.
- Full `pytest` was attempted but did not reach test execution because collection failed on the three issues listed above.
- An untracked NFS artifact exists at `src/molgen3D/config/.nfs000000006423e5500000ffa8`; it contains the old `casf_analysis_sources.yaml` content. `rm` failed with `Device or resource busy`, so leave it until the filesystem releases it.

## Important Existing Behavior Before Refactor

- `scripts/analyze_casf_conformer_sets.py` currently has one `analyze_ligand()` that mixes generation and reference rows.
- It currently calls `posebusters_result()` on each generation method before calling `analyze_source()`.
- It writes deprecated attempted CSVs:
  - `geometric_generation_attempted_filter_summary.csv`
  - `geometric_generation_attempted_pb_failure_summary.csv`
- It writes `geometric_per_ligand_metrics.csv`, not `geometric_per_ligand_long.csv`.
- `src/casf_benchmark/casf_analysis_dataset.py` currently builds the master by merging summary CSVs, including reference de-duplication hacks.
- `apps/dashboard/streamlit_app.py` currently reads legacy detail tables and labels reference rows as `baseline`.

## Next Steps

1. Resolve the full-suite collection blockers if broad validation is required.
2. For any production rerun, prefer `scripts/run_pb_once_casf_analysis.sbatch`; it records the known stuck-ligand filters and uses `--resume-parts`.
3. Production rerun order remains:
   - Reference-only core/ref.
   - Generation-only for each configured root.
   - Master CSV.
   - Dashboard DB.
