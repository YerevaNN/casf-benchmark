# Qwen checkpoint batches and the druglike set: what was run, and where the code lives

Status report for the work on branch `qwen_runs_on_casf_and_druglike`. Covers the two
new evaluation sets added on top of the existing CASF16 pipeline — 11 new Qwen
checkpoints benchmarked on CASF16 core/ref, and a new **druglike** cohort of 23 marketed
drugs / clinical candidates — the code that runs them, what actually completed, and the
headline numbers.

Snapshot date: 2026-09-09. Commit range: [`a25bf7a`](https://github.com/YerevaNN/casf-benchmark/commit/a25bf7a) plus the batch drivers committed alongside this report.

---

## 1. What the two new sets are

**CASF16 checkpoint sweep.** 11 Qwen conformer-generation checkpoints, spanning three
model sizes (0.6B / 1.7B / 4B), two coordinate tokenizers (uniform *binned* vs *FSQ*),
and three training recipes (`bigdata` pretrain, `4e`/`6e` SFT, and SFT *from* the bigdata
pretrain). Each was run through the existing pipeline on both CASF16 cohorts:

| Cohort | Ligands | Source of truth |
| --- | --- | --- |
| `casf16_core` | 94 | CASF16 core ∩ ChEMBL3D exact-topology intersection |
| `casf16_ref` | 1236 | CASF16 ref ∩ ChEMBL3D exact-topology intersection |

~900 conformers per ligand, materialized into the usual three tiers
(`fixed` / `dynamic` / `chembl_count`) and scored with the standard geometric analyzer
(PoseBusters validity, diversity/clustering, MMFF94s energy, RMSD-to-crystal recovery).

**Druglike set.** A new, deliberately small and interpretable cohort: 23 named molecules
(Nicotine, Diclofenac, Imatinib, Staurosporine, Dexamethasone, Amprenavir, Colchicine, …),
11–37 heavy atoms, 1–11 rotatable bonds. Unlike CASF16 there is **no single crystal pose
per molecule**; instead each molecule carries an ensemble of experimentally observed
conformers pooled across every PDB entry it appears in — 23–716 conformers per molecule,
**2450 in total**. That makes it a *bioactive-ensemble recall* set rather than a
single-pose recovery set, and it is scored along two independent axes:

- **In-repo** ([`scripts/eval_druglike_conformers.py`](../scripts/eval_druglike_conformers.py)): PoseBusters
  validity, diversity (torsion std / pairwise RMSD / greedy clustering) and MMFF94s energy,
  using the SMILES topology as the PoseBusters reference. No RMSD metric here.
- **External** (molgen3D's COV/MAT eval, run outside this repo): COV-R / COV-P / MAT-R /
  MAT-P against the 2450-conformer ground truth at a 0.75 Å threshold, plus censored
  variants at dmax = 3 Å. This repo *parses and merges* those outputs; it does not compute them.

13 checkpoints are configured for druglike — the 11 CASF16 ones plus
`qwen_0p6b_bigdata_step111000` and `qwen_0p6b_paired_step22000`.

---

## 2. Code that was written

Everything is driven off a **single manifest**, so adding a checkpoint is a one-place edit:

| File | Role |
| --- | --- |
| [`src/casf_benchmark/config/qwen_generation_runs.yaml`](../src/casf_benchmark/config/qwen_generation_runs.yaml) | The manifest: checkpoint label → inference output dir, per cohort (`casf16_core`, `casf16_ref`, `druglike`). Checkpoints not run on a cohort simply omit it. |
| `casf_benchmark.catalog.load_generation_runs` ([catalog.py](../src/casf_benchmark/catalog.py)) | Reads the manifest, returns ordered `(label, dirname)` pairs for a cohort; rejects duplicate labels. Cached. |
| [`tests/test_generation_runs_catalog.py`](../tests/test_generation_runs_catalog.py) | 7 tests: bundled manifest covers all three cohorts, labels unique, CASF16 cohorts ⊆ druglike, config order preserved, error cases. |

**CASF16 path:**

| File | Role |
| --- | --- |
| [`scripts/convert_qwen_pickle_to_raw_generation.py`](../scripts/convert_qwen_pickle_to_raw_generation.py) | `generation_results.pickle` (SMILES → `list[Mol]`) → the raw external-generation contract (`generation/manifest.tsv` + per-ligand SDFs). Maps generated SMILES back to CASF `ligand_id`s through the intersection CSV; one topology fanning out to several ligand_ids (same compound, multiple PDB entries) is written under each. Hard-fails on any unmapped SMILES rather than silently dropping it. |
| [`scripts/run_casf_batch_item.sh`](../scripts/run_casf_batch_item.sh) | One (checkpoint, cohort) item end to end: convert → materialize → analyze. Shards materialization one molecule per process across the node (with retries for ChEMBL3D cache collisions), then merges manifest parts. **Idempotent** — every step skips if its output exists. |
| [`scripts/run_casf_batch_all.sh`](../scripts/run_casf_batch_all.sh) | Sweeps the manifest for both cohorts. Core first (94 ligands, cheap → early results), then ref (1236). Deliberately does *not* `set -e`: one failing item never aborts the batch. Logs per item plus a `batch_summary.log`. |

**Druglike path:**

| File | Role |
| --- | --- |
| [`scripts/eval_druglike_conformers.py`](../scripts/eval_druglike_conformers.py) | Scores one checkpoint's pool: PB validity + per-check fail counts, diversity (subsampled to 200 conformers, since greedy clustering is O(n × clusters) and high-rotatable-bond molecules produce hundreds of clusters), energy over the PB-passing subset only. Parallel across molecules. Writes `{label}_per_molecule.csv` + `{label}_summary.csv`. |
| [`scripts/run_druglike_eval_batch.sh`](../scripts/run_druglike_eval_batch.sh) | Runs the above for every druglike checkpoint in the manifest. Idempotent (skips labels with a summary CSV), and distinguishes `FAIL` from `MISSING` (no pickle produced by inference). |
| [`scripts/build_druglike_eval_db.py`](../scripts/build_druglike_eval_db.py) | Per-checkpoint CSVs → `druglike_eval.sqlite` (`summary`, `per_molecule`). Re-runnable. |
| [`scripts/build_druglike_covmat.py`](../scripts/build_druglike_covmat.py) | Parses the external `covmat_results.txt` / `rmsd_matrix.csv`, folds COV/MAT into both tables, derives `model_size` / `tokenizer` / `recipe` / `step` from the label, and publishes `extended_druglike_summary` / `extended_druglike_per_molecule` into the **extended-analysis sidecar** DB. Drops its own previously-added columns first, so re-runs never collide. |

**Dashboard + supporting:**

- [`apps/dashboard/streamlit_app.py`](../apps/dashboard/streamlit_app.py): two new Extended Analysis tabs, "Druglike summary" and "Druglike per-molecule".
  The sidecar (not the main DB) is the target on purpose — `extended_*` tables in the main DB
  would flip the dashboard's auto-detected default away from the sidecar and hide everything already there.
- [`src/casf_benchmark/config/casf_generation_families.yaml`](../src/casf_benchmark/config/casf_generation_families.yaml) and both `casf_analysis_sources*.yaml`: 11 new families and 21 new run entries wired in, with display labels like "Qwen 1.7B FSQ bigdata (step 47023)".
- [`scripts/build_generation_root_overrides.py`](../scripts/build_generation_root_overrides.py): K-efficiency and the energy windows need per-conformer SDFs that the bundled
  run roots don't carry. This probes the Weka roots and emits a machine-specific (gitignored)
  `run_id → root` override map, so those analyses read a fuller copy while the dashboard keeps its provenance.
- [`docs/extras.md`](extras.md): "Qwen checkpoint batches", "Druglike test set" and the override-map workflow, with copy-pasteable commands.

Reproduction commands are in [docs/extras.md](extras.md#qwen-checkpoint-batches); nothing here needs
arguments beyond three env vars.

---

## 3. What actually ran

**CASF16: 21 of 22 items PASS.**

| Item | Outcome |
| --- | --- |
| 11 core items | all PASS (282 rows each = 94 ligands × 3 tiers) |
| 10 of 11 ref items | PASS (3708 rows each = 1236 × 3) |
| `qwen_0p6b_4e_from_bigdata_step29600` / **ref** | **FAIL** in the analyze step |

The one failure is a single-ligand crash, not a systemic one: convert and materialize
both completed (1236/1236 SDFs), and the analyzer died on ligand `5csp_5csp_conf0` with
`ValueError: atoms k and l have identical 3D coordinates` — a degenerate-coordinate pair
that breaks torsion computation. It fails identically on re-run because the pipeline
recomputes only the one missing ligand and then aborts the whole run. Consequence:
`data/results/runs/` carries 11 core + 10 ref Qwen roots, and that checkpoint appears in
core-cohort dashboard views only.

**Druglike: 12 of 13 checkpoints evaluated.**

`qwen_0p6b_paired_step22000` is reported `MISSING` — its inference directory
(`20260901_144618_qw600_pre_binned_paired_step-22000-hf_druglike`) contains only
`logs.txt`; no `generation_results.pickle` was ever written. That is an upstream inference
gap, not an eval failure. The published tables therefore hold 12 summary rows and
276 per-molecule rows (12 × 23).

---

## 4. Results

### CASF16 core, fixed tier — new checkpoints vs. existing baselines

| run                                     |   confs/lig |   best RMSD |   median RMSD |   hit@0.5 |   hit@2.0 |   clusters@1.0 |   E median |
|:----------------------------------------|------------:|------------:|--------------:|----------:|----------:|---------------:|-----------:|
| qwen_0p6b_fsq_bigdata_step70534         |     928.574 |       0.286 |         1.321 |     0.822 |     1.000 |         76.222 |     51.051 |
| qwen_1p7b_fsq_bigdata_step47023         |     932.511 |       0.304 |         1.319 |     0.789 |     1.000 |         75.478 |     50.670 |
| qwen_4b_bigdata_step110000              |     873.766 |       0.314 |         1.312 |     0.769 |     1.000 |         72.648 |     65.350 |
| qwen_1p7b_bigdata_step74000             |     859.681 |       0.332 |         1.292 |     0.772 |     1.000 |         72.707 |     64.379 |
| loqi_core                               |     957.064 |       0.341 |         1.305 |     0.822 |     1.000 |         20.700 |     37.922 |
| torsional_diffusion_core                |     839.745 |       0.372 |         1.456 |     0.789 |     1.000 |         56.644 |     88.117 |
| mcf_drugs_l_core                        |     860.457 |       0.387 |         1.459 |     0.753 |     1.000 |         33.667 |     46.212 |
| nextmol_dmt_l_core                      |     860.915 |       0.395 |         1.454 |     0.753 |     0.989 |         31.032 |     41.271 |
| qwen_0p6b_fsq_6e_from_bigdata_step28320 |     908.926 |       0.416 |         1.396 |     0.763 |     0.978 |         43.527 |     49.560 |
| qwen_1p7b_fsq_6e_from_bigdata_step28320 |     919.319 |       0.429 |         1.468 |     0.750 |     0.978 |         42.120 |     49.158 |
| qwen_0p6b_4e_from_bigdata_step29600     |     834.915 |       0.430 |         1.484 |     0.725 |     0.989 |         49.681 |     67.060 |
| qwen_1p7b_4e_from_bigdata_step26400     |     841.415 |       0.433 |         1.462 |     0.707 |     0.989 |         41.239 |     65.343 |
| qwen_1p7b_4e_step29600                  |     796.543 |       0.441 |         1.479 |     0.728 |     0.989 |         42.293 |     66.214 |
| qwen_4b_4e_from_bigdata_step20000       |     863.745 |       0.443 |         1.475 |     0.703 |     0.978 |         43.330 |     64.605 |
| rdkit_torsion (best variant)            |     999.947 |       0.447 |         1.334 |     0.628 |     1.000 |         35.883 |     81.161 |
| qwen_4b_4e_step20000                    |     794.266 |       0.476 |         1.493 |     0.667 |     0.989 |         38.731 |     66.829 |

`best RMSD` / `median RMSD` are mean-over-ligands best and median RMSD to the crystal pose (Å);
`hit@x` is the fraction of ligands with a conformer within x Å; `clusters@1.0` is mean greedy
clusters at 1.0 Å; `E median` is median MMFF94s energy.

### CASF16 ref, fixed tier

| run                                     |   confs/lig |   best RMSD |   median RMSD |   hit@0.5 |   hit@2.0 |   clusters@1.0 |   E median |
|:----------------------------------------|------------:|------------:|--------------:|----------:|----------:|---------------:|-----------:|
| qwen_1p7b_fsq_bigdata_step47023         |     942.995 |       0.320 |         1.296 |     0.789 |     0.999 |         63.568 |     49.464 |
| qwen_0p6b_fsq_bigdata_step70534         |     938.831 |       0.325 |         1.306 |     0.782 |     1.000 |         65.873 |     50.084 |
| qwen_4b_bigdata_step110000              |     889.955 |       0.346 |         1.322 |     0.773 |     0.997 |         59.906 |     64.530 |
| qwen_1p7b_bigdata_step74000             |     878.364 |       0.351 |         1.322 |     0.766 |     0.997 |         59.176 |     64.090 |
| loqi_ref                                |     970.159 |       0.386 |         1.303 |     0.760 |     0.988 |         23.943 |     31.323 |
| torsional_diffusion_ref                 |     842.552 |       0.391 |         1.414 |     0.729 |     0.995 |         47.188 |     81.371 |
| mcf_drugs_l_ref                         |     869.816 |       0.454 |         1.436 |     0.709 |     0.970 |         27.539 |     39.875 |
| nextmol_dmt_l_ref                       |     874.909 |       0.458 |         1.401 |     0.699 |     0.973 |         27.805 |     33.538 |
| qwen_0p6b_fsq_6e_from_bigdata_step28320 |     912.398 |       0.471 |         1.435 |     0.699 |     0.974 |         34.336 |     48.007 |
| qwen_1p7b_fsq_6e_from_bigdata_step28320 |     912.874 |       0.473 |         1.443 |     0.695 |     0.974 |         33.256 |     47.617 |
| rdkit_torsion (best variant)            |     999.843 |       0.477 |         1.340 |     0.595 |     0.998 |         38.146 |     76.891 |
| qwen_4b_4e_from_bigdata_step20000       |     867.336 |       0.493 |         1.455 |     0.692 |     0.968 |         31.285 |     61.975 |
| qwen_1p7b_4e_from_bigdata_step26400     |     845.814 |       0.505 |         1.462 |     0.691 |     0.967 |         30.427 |     62.707 |
| qwen_4b_4e_step20000                    |     821.871 |       0.511 |         1.467 |     0.675 |     0.967 |         30.582 |     64.697 |
| qwen_1p7b_4e_step29600                  |     816.828 |       0.512 |         1.474 |     0.688 |     0.963 |         31.340 |     63.842 |

**What the CASF16 numbers say.** The two `fsq_bigdata` checkpoints are the best models in
the benchmark on both cohorts — 0.286 Å / 0.320 Å best-RMSD, ahead of LoQI (0.341 / 0.386),
torsional diffusion (0.372 / 0.391), NextMol-DMT-L, MCF-drugs-L and every RDKit-ETKDG
variant. The ordering is consistent core→ref, so it is not a 94-ligand artifact.

Three patterns hold across the sweep:

1. **The bigdata pretrain is what matters, and SFT on top of it costs accuracy.** Every
   `*_bigdata_step*` checkpoint beats every `*_4e*` / `*_6e*` SFT checkpoint. On ref, the
   gap is 0.320–0.351 Å (pretrain) vs 0.471–0.512 Å (SFT) — the SFT variants land at or
   below the RDKit baseline.
2. **FSQ tokenization > uniform binning**, holding recipe fixed: 0.320/0.325 Å (FSQ bigdata)
   vs 0.346/0.351 Å (binned bigdata) on ref, and the FSQ runs also sit ~15 kcal/mol lower
   in median energy (≈50 vs ≈64).
3. **Size buys nothing here.** 0.6B FSQ bigdata ties or beats 1.7B and 4B; the 4B SFT
   checkpoints are among the weakest. Data/recipe/tokenizer dominate parameter count at
   this scale.

The Qwen runs also produce far more distinct conformations than the diffusion/flow
baselines — 60–76 clusters at 1.0 Å vs 24–47 — i.e. the recovery advantage comes with
genuinely broader ensembles, not a narrow well-placed pool.

### Druglike (23 molecules, 2450 ground-truth conformers, 0.75 Å threshold)

| label                                   | model_size   | tokenizer   |   confs |   PB pass |   COV-R mean |   COV-R med |   MAT-R mean |   COV-P mean |   MAT-P mean |
|:----------------------------------------|:-------------|:------------|--------:|----------:|-------------:|------------:|-------------:|-------------:|-------------:|
| qwen_0p6b_fsq_bigdata_step70534         | 0.6B         | FSQ         |   22988 |     0.891 |        0.887 |       1.000 |        0.385 |        0.461 |        0.996 |
| qwen_1p7b_fsq_6e_from_bigdata_step28320 | 1.7B         | FSQ         |   22952 |     0.884 |        0.789 |       1.000 |        0.651 |        0.454 |        1.219 |
| qwen_0p6b_fsq_6e_from_bigdata_step28320 | 0.6B         | FSQ         |   22937 |     0.867 |        0.790 |       1.000 |        0.630 |        0.449 |        1.215 |
| qwen_1p7b_fsq_bigdata_step47023         | 1.7B         | FSQ         |   22982 |     0.837 |        0.910 |       1.000 |        0.375 |        0.468 |        1.000 |
| qwen_4b_bigdata_step110000              | 4B           | Binned      |   22949 |     0.753 |        0.863 |       1.000 |        0.435 |        0.414 |        1.027 |
| qwen_4b_4e_from_bigdata_step20000       | 4B           | Binned      |   22907 |     0.734 |        0.767 |       0.969 |        0.543 |        0.436 |        1.228 |
| qwen_0p6b_bigdata_step111000            | 0.6B         | Binned      |   22973 |     0.714 |        0.867 |       1.000 |        0.441 |        0.356 |        1.108 |
| qwen_1p7b_4e_from_bigdata_step26400     | 1.7B         | Binned      |   22903 |     0.689 |        0.750 |       1.000 |        0.653 |        0.413 |        1.238 |
| qwen_0p6b_4e_from_bigdata_step29600     | 0.6B         | Binned      |   22986 |     0.683 |        0.772 |       1.000 |        0.558 |        0.413 |        1.270 |
| qwen_1p7b_bigdata_step74000             | 1.7B         | Binned      |   22963 |     0.663 |        0.866 |       1.000 |        0.448 |        0.374 |        1.091 |
| qwen_4b_4e_step20000                    | 4B           | Binned      |   22958 |     0.639 |        0.770 |       1.000 |        0.566 |        0.369 |        1.293 |
| qwen_1p7b_4e_step29600                  | 1.7B         | Binned      |   22979 |     0.628 |        0.744 |       1.000 |        0.581 |        0.375 |        1.289 |

`PB pass` is the overall PoseBusters pass rate over all ~23k generated conformers (in-repo eval).
COV/MAT are means over molecules from the external eval; COV in [0,1], MAT in Å.

**What the druglike numbers say.** The two axes agree with CASF16 on the headline and
disagree on the details, which is the useful part:

- **FSQ dominates validity.** The four FSQ checkpoints take the top four PB pass rates
  (0.837–0.891); every binned checkpoint sits at 0.628–0.753. Same direction as CASF16,
  larger effect.
- **`fsq_bigdata` also wins on precision-side geometry**: MAT-R 0.375–0.385 Å and MAT-P
  ≈1.00, versus 0.54–0.65 Å / 1.22–1.29 Å for the SFT variants. Its COV-R (0.887, 0.910)
  is the highest in the table.
- **Coverage is near-saturated per molecule** — median COV-R is 1.000 for 11 of 12
  checkpoints, so at ~1000 conformers per molecule the models cover the observed bioactive
  ensemble for a typical molecule. The mean (0.74–0.91) is dragged down by a minority of
  hard molecules, which is where the per-molecule table is worth reading.
- **COV-P ≈ 0.36–0.47 everywhere**: roughly half to two-thirds of generated conformers
  are *not* near any experimentally observed pose. Expected for a 1000-conformer pool
  against a 23–716-conformer reference, but it bounds how much of the pool is useful.
- **Failure modes are concentrated.** Across all checkpoints: `tetrahedral_chirality`
  47,022 conformer-checks failed, `internal_steric_clash` 26,182, `bond_angles` 8,633,
  `bond_lengths` 6,399, `energy_ratio` 5,505. Stereochemistry and clashes are the two
  problems worth fixing; bonded geometry is comparatively fine.

---

## 5. Is the code on GitHub?

Repository: [YerevaNN/casf-benchmark](https://github.com/YerevaNN/casf-benchmark), branch
**`qwen_runs_on_casf_and_druglike`**.

**Yes, with one gap that this commit closes.** Commit `a25bf7a` ("Add 11 Qwen checkpoints
and druglike conformer evaluation") was pushed and carries the manifest, the catalog
loader, the converter, all four Python analysis/eval/DB scripts, the config wiring, the
dashboard tabs, the docs, and the per-run result tables. But three files were left
untracked locally and were therefore **not** on GitHub — including two that `docs/extras.md`
already links to:

- `scripts/run_casf_batch_all.sh`
- `scripts/run_casf_batch_item.sh`
- `scripts/run_druglike_eval_batch.sh`
- `tests/test_generation_runs_catalog.py`

They are committed together with this report, so the documented workflow is now
reproducible from a clean clone. The 7 catalog tests pass.

The branch has not been merged to `main`; it is 2 commits ahead of `main` and 0 behind.

**Deliberately not in git:**

| Artifact | Why | Where it lives |
| --- | --- | --- |
| `casf_analysis_dashboard.sqlite` (118 MB), `extended_casf_analysis.sqlite` (95 MB), `casf_per_ligand_long.csv` (129 MB) | Over GitHub's 100 MB blob limit; LFS is disabled for this repo (with LFS off, the filters silently commit pointer stubs instead of data) | GitHub Release assets; gitignored |
| `casf_generation_root_overrides.yaml` | Machine-specific absolute paths, generated | Regenerate with `build_generation_root_overrides.py` |
| Materialized run roots, SDF pools, raw pickles | ~4.6 GB / 76k files | Weka: `/mnt/weka/vtarasov/outputs/casf_benchmark_runs/`, `druglike_eval/` |

Only `analysis/tables/geometric_per_ligand_long.csv` per run is committed under
`data/results/runs/` — which is all the dashboard build needs.

---

## 6. Known gaps and next steps

1. **`qwen_0p6b_4e_from_bigdata_step29600` / ref is still missing.** One degenerate-coordinate
   ligand (`5csp_5csp_conf0`) aborts the entire analyze step. Worth making the analyzer
   tolerate per-ligand failures (an `--allow-ligand-failures` flag, mirroring
   `--allow-k-efficiency-failures` in the extended analysis) rather than special-casing
   this ligand — the same crash will recur on any future run containing it.
2. **`qwen_0p6b_paired_step22000` never produced a druglike pickle.** Re-run inference, then
   `run_druglike_eval_batch.sh` picks it up with no changes (it is already in the manifest).
3. **The COV/MAT eval is external.** `build_druglike_covmat.py` parses molgen3D's
   `covmat_results.txt` by regex, so a change to that report's format breaks the merge
   silently-ish (missing columns rather than an error). The parse is also the only place
   the druglike RMSD numbers come from — nothing in this repo can regenerate them.
4. **Inaccurate docstring in `eval_druglike_conformers.py`.** It states the druglike set has
   "no crystal conformer for any molecule". The pickle does carry `confs` — 2450 experimental
   conformers with `pdb_ids` / `system_ids`. The accurate statement is that the set has no
   *single* reference pose, so the in-repo eval skips RMSD and leaves recovery to the COV/MAT eval.
5. **Only the `fixed` tier is discussed above.** `dynamic` and `chembl_count` tiers were
   materialized and analyzed for every checkpoint and are in the dashboard, but no
   cross-tier comparison for the new checkpoints has been written up.
6. **No dashboard join between the two sets.** Every CASF16 checkpoint is also in druglike
   (enforced by `test_casf16_cohorts_are_a_subset_of_druglike`), so a per-checkpoint view
   pairing CASF16 recovery against druglike COV/MAT and PB pass rate is straightforward and
   not yet built.
