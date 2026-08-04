# CASF Conformer Generation Benchmark — Hypothesis

Scientific hypothesis for the CASF–ChEMBL3D conformer generation benchmark. Pipeline docs: [generation_methods.md](generation_methods.md) · [analyzer.md](analyzer.md) · [materialization.md](materialization.md).

---

## Problem

Existing large-scale conformer datasets (ChEMBL3D, OMol, GEOM-drugs) are limited and biased toward vacuum-minimized or vacuum-like ensembles. They do not scale cleanly with modern compute, and they under-represent ligand geometry as it appears in protein-bound complexes. To train a large model for pharmacophore-related tasks, we need a dataset that scales with compute and includes conformers representative of **bound-state geometry**, not vacuum geometry.

This benchmark asks which generative methods can produce such ensembles. **CASF crystal poses** are used as a proxy for bound-state geometry. The analysis primarily tests a **sampling hypothesis** (which generator to use for dataset construction). Training-data effects are interpreted when comparing learned models, but matched model training is out of scope here.

Throughout, **diverse** means **geometrically diverse** (spread in 3D structure). Energy spread across basins is a secondary concern (see Metrics).

---

## Hypotheses

**H1 (sampling — primary):** At matched conformer budget and PoseBusters filtering, methods that produce more geometrically diverse ensembles recover CASF-bound conformers more often.

**H0 (formal null):** Geometric diversity and CASF recovery are uncorrelated (or not positively correlated) after controlling for sample size and validity.

**Competing explanations**

| Explanation | Claim |
| --- | --- |
| **Single-basin bias** | Ensembles collapsed onto one or a few local minima (e.g. heavy minimization) give poor bound-state coverage, even if individual conformers are low in energy. |
| **Coverage bias (H1)** | Ensembles spanning multiple distinct geometric basins are more likely to include conformers close to CASF bound poses at fixed sample budget. |
| **Training-data bias** | Learned generators inherit the conformational distribution of their training corpus (ChEMBL3D, GEOM-drugs, OMol, etc.); differences between AI methods partly reflect training data, not sampling alone. |

---

## Methods

**Panel:** ligands in the exact intersection of CASF (core and ref cohorts) with ChEMBL3D.

**References (not generated):** `casf_crystal`, `casf_opt`, `chembl3d_gt`, `chembl3d_gt_pb`.

**In-house baselines:** RDKit random embedding (raw and minimized); torsion perturbation from a ChEMBL3D seed (raw and minimized) — a simple heuristic for geometric diversity.

**Physics-based baseline:** ChEMBL3D conformers on the intersection panel.

**AI generators:** conditional conformer generation methods (SMILES-conditioned 3D generators) covering most current approaches — mostly trained on GEOM-drugs; LOQI on ChEMBL3D; Qwen variants on GEOM-drugs or OMol-scale “bigdata”; “FSQ” in Qwen names indicates conformers encoded with the project’s custom conformer tokenizer.

---

## Generation budget and PoseBusters

**In-house baselines (RDKit, torsion).** Each method fills a pre-PoseBusters pool of up to **1000 conformers**, then applies PoseBusters. The **fixed** tier retains the full pool for validation; dynamic and chembl_count tiers subsample from it (PoseBusters run once on the fixed pool; see [generation_methods.md](generation_methods.md)).

**AI methods.** Each method **generates up to 1000 conformers**; **only PoseBusters-passing structures are kept**. The number of kept conformers per ligand can be below 1000. Raw pools are tier-normalized by [materialization.md](materialization.md).

**Evaluation.** Diversity and CASF recovery metrics are computed on **PB-passing conformers only**.

**PB pass rate as yield.** `pb_fail_rate` and `kept_vs_target_rate` are reported separately. For dataset construction at scale, generation and PB filtering are assumed cheap enough that a method with strong coverage and recovery but moderately lower pass rate can compensate by oversampling and filtering. PB remains important: it defines the valid training pool and may remove strained or clashed geometries that would otherwise affect diversity or recovery.

---

## Sampling tiers

All methods are evaluated under three tiers. Each tier answers a different question.

| Tier | Target size | Use |
| --- | --- | --- |
| **`fixed`** | 1000 conformers | Behavior at scale; allows generators to exceed ChEMBL3D volume. |
| **`chembl_count`** | ChEMBL3D `conformer_count` per ligand | Fair comparison to existing ChEMBL3D datasets; prevents winning by raw volume alone. |
| **`dynamic`** | `max(1, −20 + 22 × rotatable_bonds)` | Ligand-adaptive budget between fixed and matched count. |

---

## Metrics

### Bound-state recovery (primary dependent variables)

Heavy-atom aligned RMSD to CASF crystal. Cohort summaries use per-ligand averaging.

| Metric | Definition |
| --- | --- |
| `casf_best_rmsd` | Per ligand: minimum RMSD to CASF crystal; cohort mean over ligands. |
| `casf_median_rmsd` | Per ligand: median RMSD to CASF crystal. |
| `casf_hit_0p25`, `casf_hit_0p5`, `casf_hit_0p75`, `casf_hit_2p0` | Fraction of ligands whose best crystal RMSD ≤ threshold (Hit@X). |

CASF optimized ligand metrics (`casf_opt_*`) are computed in parallel. Hit rate is the main success metric; best and median RMSD capture partial near-misses.

### Geometric diversity (primary independent variables)

Greedy heavy-atom RMSD clustering in generation order at thresholds **0.5, 1.0, 2.0, 3.0 Å**.

| Metric | Definition |
| --- | --- |
| `greedy_clusters_{threshold}` | Cluster count per ligand at that RMSD threshold. |
| `clusters_per_100_{threshold}` | Cluster count normalized to 100 conformers (primary cross-method comparison). |
| `cluster_entropy_{threshold}` | Shannon entropy of cluster occupancies, normalized to [0, 1]; 1 = evenly spread, ~0 = one cluster dominates. |
| `largest_cluster_fraction_{threshold}` | Fraction of conformers in the largest cluster; high values indicate a collapsed ensemble. |
| `effective_clusters_{threshold}` | exp(entropy); effective number of populated clusters. |
| `singleton_fraction_{threshold}` | Fraction of clusters with exactly one conformer. |
| `mean_torsion_std_deg` | Mean standard deviation of rotatable torsion angles. |
| `pairwise_mean`, `pairwise_p90` | Mean / 90th percentile of pairwise heavy-atom RMSD within the ensemble. |

Primary diversity summaries: `clusters_per_100_1p0` and `cluster_entropy_1p0`, with sensitivity across thresholds.

### Energy spread (secondary)

MMFF94s energies on PB-passing conformers. Energy metrics describe **spread across basins**, not “lower is better.”

| Metric | Definition |
| --- | --- |
| `energy_std` | Per-ligand energy standard deviation; high values suggest multi-basin coverage. |
| `energy_min`, `energy_max`, `energy_median` | Per-ligand energy range and central tendency. |

Low `energy_std` together with low cluster count indicates single-basin collapse. Geometric diversity is the primary axis; energy diversity across low basins is secondary.

### Validity and yield

| Metric | Definition |
| --- | --- |
| `pb_fail_rate_mean` | Per-ligand PB failure rate, cohort-averaged. |
| `kept_vs_target_rate_mean` | Per-ligand yield (`kept / target_confs`), cohort-averaged. |

### Extended analyses

K-efficiency curves (recovery vs sample size K), Pareto frontier (`clusters_per_100_1p0` vs Hit@0.75 vs PB pass rate), and paired deltas vs `rdkit_random_raw_fixed` and `chembl3d_gt_pb` — see [analyzer.md](analyzer.md) and [extras.md](extras.md).

---

## Predictions (if H1 holds)

1. Torsion perturbation and wide-sampling methods show higher `clusters_per_100_*` and `cluster_entropy_*` than single-basin / heavily minimized pipelines at the same tier.
2. Single-basin ensembles (low cluster count, low `energy_std`) show lower Hit@X even when `energy_min` is low.
3. Under **`chembl_count`**, diversity advantages may shrink because volume is capped.
4. Under **`fixed`**, diversity should translate more clearly into Hit@X.
5. GEOM-drugs-style vacuum recall may not rank methods the same way as CASF Hit@X.

---

## Falsification and revision

| Outcome | Interpretation |
| --- | --- |
| No positive association between `clusters_per_100_1p0` (or `cluster_entropy_1p0`) and `casf_hit_0p75` / `casf_hit_2p0` at matched tier and validity | Reject H1. |
| Diversity helps only above a threshold, or only for certain method families | Refine H1 (e.g. diversity helps when combined with scale or binding-relevant training signal). |
| High Hit@X with low `greedy_clusters_*` | Refine the diversity metric or shift to K-efficiency and conditional-generation behavior. |

Core-set results assessment: [casf16-core-hypothesis-assessment.md](casf16-core-hypothesis-assessment.md).

---

## Scope

This benchmark selects generative methods for building a scalable, bound-state-oriented conformer training dataset. It does **not** by itself prove that models trained on diverse ensembles outperform models trained on minimized corpora — that requires matched training ablations with controlled splits, token budgets, and inference settings.
