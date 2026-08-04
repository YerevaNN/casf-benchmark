# CASF core results: hypothesis assessment

Assessment of [casf16_hypothesis.md](casf16_hypothesis.md) against the 94-ligand core intersection. Qwen ref generations do not exist; ref observations used only for non-Qwen replication.

**Data:** [`data/results/casf_analysis_dashboard.sqlite`](../data/results/casf_analysis_dashboard.sqlite)

## Verdict: partial support; revise the hypothesis

Binding-relevant, validity-filtered coverage improves best-case bound-pose recovery at high sampling budget. Maximizing geometric diversity or clusters per 100 is **not** sufficient.

---

## 1. Claim-by-claim decision

| Hypothesis component | Observed core result | Decision | Confidence |
| --- | --- | --- | --- |
| Wider methods generate more geometric diversity. | Torsion raw and wide Qwen variants produce roughly 74–77 mean 1 Å clusters, versus 20–25 after minimization. | Supported | High descriptively |
| Minimization can create single-basin bias that hurts bound-pose recovery. | Within both RDKit and torsion, minimization sharply lowers cluster count and worsens paired best RMSD at fixed high K. | Supported | Moderate-to-high |
| More diversity at matched budget produces more Hit@X. | Method-level fixed rows correlate positively, but count-matched Qwen does not beat ChEMBL3D-PB on the full 94-ligand universe. | Mixed | Moderate |
| Maximum diversity is the best strategy. | Raw torsion and best Qwen have similar cluster counts but 77.7% versus 91.5% Hit@0.75; LOQI succeeds with few clusters. | Contradicted | High |
| Clusters per 100 is the primary cross-method diversity predictor. | Within-ligand centered association with Hit@0.75 is −0.075. | Not supported | Moderate |
| Energy diversity across low basins is useful. | Energy std for several raw/AI methods dominated by extreme outliers. | Unresolved | Low |
| Moderate PB loss can be handled by cheap oversampling. | Mean pass rates ~1.04–1.25× for most AI methods, but systematic chirality failures violate the assumption. | Conditionally supported | Moderate |
| Training-data distribution affects learned generator behavior. | Big-data and FSQ-pretrained Qwen variants outperform revisited/plain variants. | Strong lead, not causal proof | Moderate |
| A better conformer source will improve a trained downstream model. | No matched training-data ablation has been run. | Not tested | None |

## 2. Primary core fixed-pool result

| Method | Evaluable | Mean valid yield | Mean confs | Mean 1 Å clusters | Mean best RMSD | Mean median RMSD | 94-lig Hit@0.75 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Qwen 1.7B FSQ + big-data pretrain | 91/94 | 93.4% | 934 | 75.7 | 0.301 | 1.303 | 91.5% |
| Qwen 1.7B big data | 93/94 | 85.9% | 856 | 74.5 | 0.343 | 1.301 | 88.3% |
| Qwen 4B big data | 91/94 | 87.4% | 872 | 76.0 | 0.332 | 1.308 | 86.2% |
| Qwen 0.6B FSQ + big-data pretrain | 91/94 | 92.7% | 926 | 76.5 | 0.334 | 1.328 | 86.2% |
| Torsional Diffusion | 90/94 | 84.0% | 840 | 56.6 | 0.372 | 1.456 | 85.1% |
| MCF drugs-L | 93/94 | 86.0% | 860 | 33.7 | 0.387 | 1.459 | 84.0% |
| NextMol DMT-L | 93/94 | 86.1% | 861 | 31.0 | 0.395 | 1.454 | 84.0% |
| LOQI | 90/94 | 95.7% | 957 | 20.7 | 0.341 | 1.305 | 81.9% |
| RDKit raw | 94/94 | 100.0% | 1,000 | 35.9 | 0.447 | 1.334 | 80.9% |
| ChEMBL3D-PB | 94/94 | PB-passing subset | 83 | 7.8 | 0.524 | 1.325 | 79.8% |
| Torsion perturbation raw | 94/94 | 89.6% | 896 | 74.1 | 0.497 | 1.346 | 77.7% |
| RDKit minimized | 94/94 | 100.0% | 1,000 | 20.3 | 0.528 | 1.348 | 76.6% |
| Torsion perturbation minimized | 94/94 | 100.0% | 1,000 | 25.0 | 0.554 | 1.299 | 75.5% |

Best Qwen: 86/94 ligands within 0.75 Å vs 75/94 ChEMBL3D-PB; mean best RMSD 0.524→0.301 Å (~934 vs 83 conformers/ligand).

### Paired raw vs minimized ablations

| Comparison | Raw−min clusters | Raw−min energy med | Raw−min best RMSD | Bootstrap 95% CI | Hit@0.75 change |
| --- | --- | --- | --- | --- | --- |
| RDKit raw vs minimized | +15.6 | +54.1 | −0.081 Å | −0.131 to −0.029 | +4 ligands / +4.3 pp |
| Torsion raw vs minimized | +49.1 | +47.5 | −0.057 Å | −0.109 to −0.005 | +2 ligands / +2.1 pp |

## 3. Matched-count tier

| Method | Evaluable | Mean confs | Mean 1 Å clusters | Clusters/100 | Mean best RMSD | 94-lig Hit@0.75 |
| --- | --- | --- | --- | --- | --- | --- |
| ChEMBL3D-PB | 94/94 | 83.3 | 7.8 | 28.5 | 0.524 | 79.8% |
| Qwen 1.7B FSQ + big-data pretrain | 90/94 | 73.0 | 23.8 | 39.9 | 0.514 | 76.6% |
| NextMol DMT-L | 92/94 | 65.1 | 14.0 | 36.6 | 0.552 | 74.5% |
| Qwen 1.7B big data | 92/94 | 64.9 | 22.7 | 43.5 | 0.556 | 74.5% |
| LOQI | 90/94 | 75.6 | 9.4 | 28.2 | 0.530 | 71.3% |
| Torsional Diffusion | 90/94 | 63.8 | 17.9 | 43.0 | 0.574 | 71.3% |
| MCF drugs-L | 93/94 | 66.5 | 14.1 | 36.9 | 0.576 | 70.2% |
| RDKit raw | 94/94 | 105.5 | 14.9 | 32.1 | 0.604 | 66.0% |

At ChEMBL-count scale, best Qwen: 72/94 vs ChEMBL 75/94. Matches hypothesis prediction that diversity gains shrink under count-matched tier.

### Qwen K-efficiency

| K | Hit@0.75 (91 evaluable) | Hit@0.75 (94 universe) |
| --- | --- | --- |
| 1 | 28.9% | 28.0% |
| 2 | 36.8% | 35.6% |
| 5 | 51.3% | 49.7% |
| 10 | 59.2% | 57.3% |
| 25 | 69.8% | 67.6% |
| 50 | 77.1% | 74.6% |
| 100 | 81.3% | 78.7% |
| 250 | 88.0% | 85.2% |
| 500 | 90.9% | 88.0% |
| 1,000 | 94.5% | 91.5% |

Crossover vs ChEMBL 79.8% between K=100 and K=250.

## 4. Diversity vs recovery correlations

| Scope | Metric | Corr w/ Hit@0.75 | Note |
| --- | --- | --- | --- |
| Across 19 fixed method rows | Absolute 1 Å clusters | +0.591 | Broad method-level support; Qwen variants are not independent replicates. |
| Across 19 fixed method rows | 1 Å clusters per 100 | +0.525 | Positive descriptive association, but weaker than absolute coverage. |
| Within ligand, method-centered | Absolute 1 Å clusters | +0.060 | Only a weak relation to crossing the 0.75 Å threshold. |
| Within ligand, method-centered | 1 Å clusters per 100 | −0.075 | The hypothesis document's primary normalized metric is not predictive here. |
| Within ligand, method-centered | Cluster entropy | +0.153 | Even occupancy is more informative than local cluster density. |
| Within ligand, method-centered | Largest-cluster fraction | −0.214 | Collapse into one dominant mode is associated with worse recovery. |
| Within ligand, method-centered | Pairwise mean RMSD | +0.136 | Broad spread has a weak positive signal. |
| Within ligand, method-centered | Conformer count | +0.142 | Best-of-set recovery still benefits from budget. |

Demote clusters/100 as primary predictor; keep absolute coverage, entropy, largest-cluster fraction.

## 5. vs ChEMBL3D-PB (statistics)

| Method | Net Hit Δ | Mean best-RMSD Δ | Bootstrap 95% CI | McNemar p |
| --- | --- | --- | --- | --- |
| Qwen 1.7B FSQ + big-data pretrain | +11 | −0.220 Å | −0.289 to −0.160 | 0.0127 |
| Qwen 1.7B big data | +8 | −0.182 Å | −0.252 to −0.113 | 0.0768 |
| Torsional Diffusion | +5 | −0.153 Å | −0.210 to −0.101 | 0.2668 |
| MCF drugs-L | +4 | −0.141 Å | −0.215 to −0.072 | 0.5034 |
| NextMol DMT-L | +4 | −0.133 Å | −0.201 to −0.071 | 0.4807 |
| LOQI | +2 | −0.184 Å | −0.238 to −0.129 | 0.8036 |
| RDKit raw | +1 | −0.077 Å | −0.152 to −0.004 | 1.0000 |

## 6–7. Tail coverage & energy

| Method | Energy median | Energy std | Reading |
| --- | --- | --- | --- |
| ChEMBL3D-PB | 36.8 | 2.6 | Compact, sample-efficient reference; PB subset only. |
| LOQI | 37.9 | 2.9 | Low spread and strong recovery: energy spread is not required for success. |
| Qwen 1.7B FSQ + big-data pretrain | 50.6 | 138.0 | Good recovery, but a substantial high-energy tail needs auditing. |
| Torsional Diffusion | 88.1 | 5,599.5 | Extreme outliers make raw energy standard deviation unsuitable as evidence. |
| RDKit raw | 81.2 | 10.2 | Higher-energy raw ensemble; better high-K recovery than minimized RDKit. |
| RDKit minimized | 27.1 | 2.8 | Low-energy collapse; worse best RMSD and Hit@0.75 at high K. |
| Torsion raw | 73.9 | 4,782.1 | Diversity includes severe energy outliers; not a clean dataset recipe. |
| Torsion minimized | 26.4 | 2.6 | Physicality improves, but the ensemble loses useful high-K tails. |

## 8. PoseBusters & oversampling

| Method | Pass/yield | Multiplier | Qualification |
| --- | --- | --- | --- |
| Qwen 1.7B FSQ + big-data pretrain | 93.4% | 1.07× | All 94 generation rows exist, but three ligands have zero PB-passing conformers because every conformer fails tetrahedral chirality. |
| LOQI | 95.7% | 1.04× | All 94 generation rows exist, but four ligands have zero PB-passing conformers; all failures are tetrahedral chirality. |
| Torsional Diffusion | 84.0% | 1.19× | Moderate oversampling can recover count if failures are stochastic. |
| MCF / NextMol | ≈86.0% | ≈1.16× | One ligand each has zero PB-passing conformers; chirality dominates reported failures. |
| Torsion perturbation raw | 89.6% after selection | 1.12× post-PB | Including preselection clashes, only 84,266 of 167,000 trials remain: about 1.98× raw trials per kept conformer. |

PB-empty ligands (all fail tetrahedral_chirality): Qwen — 1nc3, 3syr, 4f2w; LOQI adds 3g31.

## 9. Qwen checkpoint contrasts

| Contrast | Hit@0.75 | Mean best RMSD | Lead |
| --- | --- | --- | --- |
| 1.7B FSQ: big-data pretrain vs no pretrain | 91.5% vs 83.0% | 0.301 vs 0.441 Å | Large, coherent pretraining/representation signal. |
| 1.7B: big data vs revisited | 88.3% vs 79.8% | 0.343 vs 0.456 Å | Training corpus/checkpoint choice appears more important than size alone. |
| 4B: big data vs revisited | 86.2% vs 78.7% | 0.332 vs 0.460 Å | The same direction repeats at 4B. |
| Best 1.7B vs 4B big data | 91.5% vs 86.2% | 0.301 vs 0.332 Å | No monotonic parameter-scaling story. |

## 10. Data-quality limits

| Severity | Issue | Consequence |
| --- | --- | --- |
| Critical | Qwen rotatable-bond metadata equals 2 for every ligand. | Qwen flexibility strata and dynamic-tier interpretation are invalid. Some FSQ-pretrained dynamic counts also do not follow the same apparent materialization behavior. |
| High | Rotatable-bond counts disagree between source families. | The flexible-ligand story cannot be compared fairly until every row joins one canonical ligand metadata table. |
| High | PB-empty outcomes have no continuous RMSD and are excluded from ordinary means. | Qwen best is 94.5% on 91 PB-evaluable ligands but 91.5% on the full 94-ligand universe. Use the latter as the primary endpoint. |
| High | Qwen has a K curve, but ChEMBL3D-PB does not have a PB-filtered per-conformer K curve. | Qwen first exceeds ChEMBL's full-set 79.8% universe Hit between K=100 and K=250, but a fully matched sample-efficiency comparison is still unavailable. |
| High | PoseBusters truth handling differs by source. | Reported chirality failures may mix real stereochemical errors with topology/reference-mapping artifacts. |
| High | PB failures are concentrated by ligand and mechanism. | Oversampling works for random attrition, not for a method that systematically emits the wrong stereochemistry for a ligand. |
| High | Energy standard deviations contain extreme outliers. | Raw energy spread cannot currently support the 'diverse low-energy basins' part of the thesis. |
| Medium | The core set has only four canonical 10–14-rotor ligands and none above 14. | Even after metadata repair, core-only flexibility conclusions will be underpowered. |
| Medium | The stale-manifest and universe audit was explicitly skipped. | Numeric checks pass, but artifact completeness should be revalidated before manuscript use. |
| Medium | Twelve related Qwen variants were inspected and the best was highlighted. | The unadjusted p-value for the winner is post-selection evidence, not a preregistered single test. |

## 11. Manuscript-safe claim

High-budget learned generators recover bound conformations beyond ChEMBL3D-PB; raw-vs-minimized ablations show minimization removes useful high-K coverage; validity-filtered binding-relevant coverage beats maximum diversity alone.

**Avoid:** max diversity wins; Qwen sample-efficient vs ChEMBL without K curve; flexible-ligand strata from Qwen metadata; energy spread = useful basins; diverse training helps before matched ablation.

## 12. Priority experiments

| # | Experiment | Why |
| --- | --- | --- |
| 1 | Repair metadata and denominators | Join canonical rotatable bonds/heavy atoms by ligand; define all 94 core ligands as the endpoint universe; rerun stale-manifest and PB-empty-ligand audits. |
| 2 | Complete the matched K analysis | Qwen K curves now exist. Export PB-passing per-conformer RMSDs for ChEMBL3D-PB and compare both sources at identical K with paired deterministic subsamples. |
| 3 | Test useful rather than raw diversity | Measure Hit@K after selecting one medoid per 1 Å cluster, and compare with random subsets of exactly the same size. |
| 4 | Make energy chemically interpretable | Report geometric diversity inside per-ligand relative energy windows such as ΔE ≤ 5, 10, and 20 kcal/mol; audit extreme MMFF outliers and validate a subset with xTB. |
| 5 | Audit PoseBusters failures | Manually inspect chirality-dominated ligands; standardize mol_true/topology handling; report per-ligand yield tails in addition to the mean. |
| 6 | Run Qwen on the reference set | Test whether the best core observation survives on the larger cohort. Keep core as the primary discovery set and ref as validation. |
| 7 | Perform matched training ablations | Hold molecules, conformer counts, token budget, model, and inference K fixed; vary only minimized, raw-diverse, and PB/energy/cluster-balanced training targets. |
