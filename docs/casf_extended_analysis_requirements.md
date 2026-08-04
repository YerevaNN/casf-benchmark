# CASF Extended Analysis Implementation Requirements

## Scope

This document is an implementation brief for extending the CASF conformer-set analysis. The goal is to produce publication-grade evidence for the JCIM narrative:

> Quality-controlled geometric diversity improves recovery of CASF-like bound conformations, but useful diversity is constrained by chemical validity, sample budget, and ligand flexibility.

Do not turn this into a new unconstrained analysis project. Add narrowly scoped, reproducible analyses around the existing per-ligand tables and dashboard database.

Explicitly out of scope:

- Qwen ref generation/reruns.
- Bootstrap confidence intervals or multi-seed uncertainty analysis.

If uncertainty is needed later, implement it separately. The present deliverable should focus on deterministic extended analyses from the existing result files.

## Existing Inputs

Primary database:

- `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite`

Primary long-form table in the database:

- `per_ligand_long`

Relevant summary tables:

- `comparison_rows`
- `comparison_strata`
- `analysis_sources`

Per-run CSV source, when extra columns are needed:

- `<run_root>/analysis/tables/geometric_per_ligand_long.csv`
- `<run_root>/analysis/tables/geometric_generation_filter_summary.csv`
- `<run_root>/analysis/tables/geometric_generation_pb_failure_summary.csv`
- `<run_root>/generation/manifest.tsv`

The implementation should prefer the SQLite database for cross-run aggregation and only read per-run files for diagnostics not present in the database.

## Global Requirements

All analyses must be ligand-paired unless the output is explicitly labeled as unpaired.

Use `mol_id` as the pairing key. For each comparison, compute the common ligand universe first and report:

- `n_reference_ligands`
- `n_method_ligands`
- `n_common_ligands`
- `n_reference_only`
- `n_method_only`
- `excluded_reason_summary`, when available from manifest/selection audit

Never compare method means against ChEMBL3D-PB using different ligand sets without labeling the result as exploratory.

All output tables must include:

- `ligand_set`
- `run_id`
- `run_label`
- `source`
- `family`
- `tier`
- `method`
- `n_ligands`

For paired comparison tables, also include:

- `baseline_source`
- `n_common_ligands`
- `mean_delta_casf_best_rmsd`
- `median_delta_casf_best_rmsd`
- `win_rate_0p1A`
- `loss_rate_0p1A`
- `delta_hit_0p5`
- `delta_hit_0p75`
- `delta_hit_2p0`

Use CASF crystal RMSD as the main endpoint:

- `casf_best_rmsd`
- `casf_hit_0p5`
- `casf_hit_0p75`
- `casf_hit_2p0`

Use `casf_opt_*` only as supporting output.

## Required Analysis 1: Ligand Universe And Missingness Audit

Purpose:

Identify whether apparent missing ligands are caused by generation failure, analyzer prefiltering, missing SDFs, or stale analysis output.

Required output:

- `analysis/tables/extended_ligand_universe_audit.csv`
- `analysis/reports/extended_ligand_universe_audit.md`

Table columns:

- `ligand_set`
- `run_id`
- `run_label`
- `manifest_ligands`
- `analysis_generation_ligands`
- `reference_ligands`
- `missing_vs_reference`
- `extra_vs_reference`
- `zero_rotatable_bonds_excluded`
- `missing_generation_sdf_excluded`
- `manifest_ok_but_not_analyzed`
- `notes`

Implementation details:

1. Load reference ligand IDs from `reference_ref` and `reference_core` CASF crystal rows.
2. Load generation ligand IDs from `per_ligand_long where row_type='generation'`.
3. For each run root, load `generation/manifest.tsv` when present.
4. Reconstruct analyzer selection logic from `select_analysis_mol_ids()`:
   - start from the intersection CSV ligand IDs
   - exclude ligands missing from the generation manifest
   - exclude missing generation SDFs for discovered methods
   - do not re-filter on manifest `rotatable_bonds`; chemical eligibility is defined upstream in the CSV
5. Flag ligands that have manifest `status == ok`, nonzero `pb_pass_confs`, existing SDFs, but no row in `per_ligand_long`.

Acceptance checks:

- Core generation runs should report no missing ligands.
- Ref runs should report ligands excluded only when they are absent from the CSV, manifest, or generation SDFs.
- Ref generation analysis uses the full intersection CSV; large flexible ligands such as `6o0k`, `6o0m`, and `6o0p` are no longer skipped upstream.

## Required Analysis 2: K-Efficiency Curves

Purpose:

Test whether high-K wins are simply due to larger sample budget, and identify the K at which diverse methods overtake ChEMBL3D-PB.

Required output:

- `analysis/tables/extended_k_efficiency.csv`
- `analysis/tables/extended_k_efficiency_strata.csv`
- `analysis/figures/k_efficiency_hit_0p75_ref.png`
- `analysis/figures/k_efficiency_best_rmsd_ref.png`
- `analysis/reports/extended_k_efficiency.md`

Methods to include:

- ChEMBL3D-PB
- LOQI fixed
- torsional diffusion fixed
- RDKit raw fixed
- RDKit minimized fixed
- torsion raw fixed
- NextMol DMT-L fixed
- MCF drugs-L fixed, if present

K values:

- `1, 2, 5, 10, 25, 50, 100, 250, 500, 1000`

Required metrics at each K:

- `mean_best_rmsd_at_k`
- `median_best_rmsd_at_k`
- `hit_0p5_at_k`
- `hit_0p75_at_k`
- `hit_2p0_at_k`
- `mean_pb_pass_rate_at_k`, if available
- `mean_conformers_available`
- `n_ligands_with_at_least_k`

Sampling rules:

- Use deterministic subsampling with fixed seeds, but do not report seed variance.
- Use at least 20 deterministic subsamples per K for methods with more than K conformers per ligand.
- For each ligand and K, sample from that ligand's PB-passing conformer set.
- If per-conformer RMSDs are not currently saved, add a per-conformer RMSD export first. Do not approximate K curves from per-ligand best RMSD.

Required strata:

- rotatable bonds: `0-4`, `5-9`, `10-14`, `15+`
- heavy atoms: `<20`, `20-29`, `30-39`, `40+`

Interpretation requirements:

- Report crossover K versus ChEMBL3D-PB for Hit@0.75.
- Report whether crossover K shifts upward with rotatable-bond count.
- Avoid claiming a method is more sample-efficient unless it beats ChEMBL3D-PB at the same K.

## Required Analysis 3: Paired Delta And Win/Loss Tables

Purpose:

Convert aggregate rankings into ligand-level evidence. This should become one of the main manuscript tables or SI tables.

Required output:

- `analysis/tables/extended_paired_delta_vs_chembl3d_pb.csv`
- `analysis/tables/extended_paired_delta_vs_rdkit_raw.csv`
- `analysis/reports/extended_paired_delta_summary.md`

Baselines:

- ChEMBL3D-PB
- RDKit raw fixed

Methods:

- all fixed-tier generation methods present in the database
- ChEMBL3D all, as a reference sanity check

Required metrics:

- mean and median per-ligand delta in `casf_best_rmsd`
- win/loss/tie rates using thresholds:
  - `0.05 A`
  - `0.10 A`
  - `0.25 A`
- hit deltas at:
  - `0.5 A`
  - `0.75 A`
  - `2.0 A`
- count of catastrophic losses:
  - method worse than baseline by `>1.0 A`
- count of large rescues:
  - method better than baseline by `>1.0 A`

Required stratification:

- rotatable-bond bin
- heavy-atom bin
- ChEMBL3D conformer-count bin
- baseline CASF difficulty bin:
  - ChEMBL3D-PB best RMSD `<=0.5 A`
  - `0.5-0.75 A`
  - `0.75-2.0 A`
  - `>2.0 A`

Acceptance checks:

- Results must be computed on common ligand IDs only.
- The report must list the top 20 rescued ligands and top 20 worsened ligands for each main method.

## Required Analysis 4: Diversity-Validity-Recovery Frontier

Purpose:

Show that the result is not a leaderboard and not "diversity is always good." The useful target is a frontier balancing coverage, physical validity, and bound-pose recovery.

Required output:

- `analysis/tables/extended_frontier_summary.csv`
- `analysis/figures/frontier_hit_0p75_vs_clusters_ref.png`
- `analysis/figures/frontier_rmsd_vs_validity_ref.png`
- `analysis/reports/extended_frontier_summary.md`

Required x/y/color encodings:

- x-axis option 1: mean `greedy_clusters_1p0`
- x-axis option 2: mean `conformer_count`
- y-axis option 1: `casf_hit_0p75`
- y-axis option 2: `casf_best_rmsd`
- color: PB fail rate or PB pass rate
- point size: mean conformer count

Required derived metrics:

- `useful_clusters_per_100_confs`
- `hit_0p75_per_100_confs`
- `valid_hit_score = casf_hit_0p75 * mean_pb_pass_rate`
- `frontier_rank`, computed by non-dominated sorting over:
  - maximize `casf_hit_0p75`
  - maximize PB pass rate
  - minimize mean conformer count

Acceptance checks:

- RDKit raw should appear chemically clean but not necessarily best recovery.
- Raw torsion should show high diversity with weaker useful recovery.
- LOQI and torsional diffusion should separate as recovery/validity tradeoffs on the ref set.

## Required Analysis 5: PoseBusters Failure Mechanism Summary

Purpose:

Separate generic invalidity from specific chemical failure modes.

Required output:

- `analysis/tables/extended_pb_failure_mechanisms.csv`
- `analysis/tables/extended_pb_failure_tail_ligands.csv`
- `analysis/reports/extended_pb_failure_mechanisms.md`

Required metrics:

- total PB fail rate by method
- median per-ligand PB pass rate
- number of ligands with PB pass rate `<0.9`
- number of ligands with PB pass rate `<0.5`
- dominant PB failure test by method
- second dominant PB failure test by method
- fraction of failures attributable to:
  - `tetrahedral_chirality`
  - `internal_steric_clash`
  - `energy_ratio`
  - bond geometry tests

Required tail table:

For each method, list the 25 worst ligands by PB pass rate with:

- `mol_id`
- `rotatable_bonds`
- `heavy_atoms`
- `pb_pass_rate`
- `pb_input_confs`
- top three failed PB tests
- `casf_best_rmsd`
- `casf_hit_0p75`

Important:

This is not a manual chirality audit. Do not classify chirality failures into true/artifact categories here. Only summarize the automated PB labels.

## Required Analysis 6: Bound-Pose Difficulty And Rescue Analysis

Purpose:

Identify where diverse methods add value: easy ligands, hard ligands, flexible ligands, large ligands, or ligands under-covered by ChEMBL3D-PB.

Required output:

- `analysis/tables/extended_bound_pose_difficulty.csv`
- `analysis/tables/extended_rescue_cases.csv`
- `analysis/reports/extended_bound_pose_rescue.md`

Define baseline difficulty using ChEMBL3D-PB:

- easy: `casf_best_rmsd <= 0.5 A`
- moderate: `0.5 < casf_best_rmsd <= 0.75 A`
- hard: `0.75 < casf_best_rmsd <= 2.0 A`
- failed: `casf_best_rmsd > 2.0 A`

For each method and difficulty bin, report:

- `n_common_ligands`
- `mean_delta_casf_best_rmsd`
- `hit_0p75`
- `rescue_rate`
- `catastrophic_loss_rate`

Definitions:

- rescue: baseline misses Hit@0.75 and method hits Hit@0.75
- catastrophic loss: baseline hits Hit@0.75 and method misses Hit@2.0

Required rescue-case table:

- top 50 ChEMBL3D-PB misses rescued by LOQI fixed
- top 50 ChEMBL3D-PB misses rescued by torsional diffusion fixed
- top 50 ChEMBL3D-PB misses rescued by RDKit raw fixed

Include:

- `mol_id`
- `rotatable_bonds`
- `heavy_atoms`
- baseline best RMSD
- method best RMSD
- delta RMSD
- method PB pass rate
- method cluster count at 1.0 A

## Required Analysis 7: Metadata And Sanity Checks

Purpose:

Prevent manuscript figures from being built on stale or inconsistent data.

Required output:

- `analysis/tables/extended_sanity_checks.csv`
- `analysis/reports/extended_sanity_checks.md`

Checks:

- ligand counts by run, source, tier
- duplicate `mol_id/source/tier` rows
- missing required columns
- impossible metric values:
  - negative RMSD
  - hit values outside `[0, 1]`
  - PB pass/fail counts inconsistent with PB input
  - conformer count greater than PB input when it should not be
- Qwen rotatable-bond metadata warning:
  - flag but do not fix in this task
- stale analysis warning:
  - manifest has successful rows and SDFs, but per-ligand table lacks the ligand

Acceptance checks:

- The report must explicitly state which runs are manuscript-ready and which are exploratory.
- Any figure-generating command should fail if required sanity checks fail.

## Implementation Placement

Preferred implementation:

- Add a new script: `scripts/extended_casf_analysis.py`
- Keep it read-only with respect to existing analysis outputs.
- Write new outputs under:
  - `/mnt/weka/mbedrosian/pharma_generation_analysis/extended_analysis/tables`
  - `/mnt/weka/mbedrosian/pharma_generation_analysis/extended_analysis/figures`
  - `/mnt/weka/mbedrosian/pharma_generation_analysis/extended_analysis/reports`

Do not modify `scripts/analyze_casf_conformer_sets.py` unless per-conformer RMSD export is required for K-efficiency curves. If that export is required, make the smallest possible addition and gate it behind a CLI flag.

Suggested CLI:

```bash
PYTHONPATH=src python scripts/extended_casf_analysis.py \
  --dashboard-db /mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite \
  --output-dir /mnt/weka/mbedrosian/pharma_generation_analysis/extended_analysis
```

Suggested optional flags:

- `--analysis-sources-yaml config/casf_analysis_sources.yaml`
- `--k-values 1,2,5,10,25,50,100,250,500,1000`
- `--skip-k-efficiency`
- `--fail-on-sanity-error`

## Testing Requirements

Add focused tests rather than broad end-to-end tests.

Required tests:

- common-ligand pairing excludes non-overlapping ligands
- paired delta signs are correct
- win/loss threshold logic is correct
- rescue/catastrophic-loss definitions are correct
- non-dominated frontier ranking works on a small synthetic table
- PB failure JSON parsing handles empty and malformed JSON safely
- ligand missingness audit distinguishes:
  - zero-rotatable exclusion
  - missing SDF
  - manifest-ok-but-not-analyzed

Suggested test file:

- `tests/pharmacophore/test_extended_casf_analysis.py`

Run:

```bash
PYTHONPATH=src:. pytest tests/pharmacophore/test_extended_casf_analysis.py -q
```

If shared code is touched, also run:

```bash
PYTHONPATH=src:. pytest tests/pharmacophore/test_analyze_casf_conformer_sets.py tests/pharmacophore/test_casf_analysis_dataset.py -q
```

## Final Deliverables

The implementing agent should produce:

- one executable script
- one test file
- all required CSV tables
- all required Markdown reports
- key PNG figures for K-efficiency and frontier plots
- a short `README.md` inside the output directory describing how to regenerate the outputs

The final report should answer these questions directly:

1. Does each high-performing method still win on common ligands?
2. At what K does each diverse method beat ChEMBL3D-PB?
3. Are improvements concentrated in flexible or hard baseline ligands?
4. Which methods sit on the recovery-validity frontier?
5. Which failures are narrow and fixable versus broad chemical-validity problems?
6. Which runs should be used for manuscript figures versus treated as exploratory?
