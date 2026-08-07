# CASF16 conformer generation story report

Expert read of the current CASF16 results, structured around the original hypothesis in
[`casf16_hypothesis.md`](casf16_hypothesis.md). This report focuses on text and tables, not
figures. It is meant as an internal decision document for whether there is a JCIM-grade story
and which result blocks are worth carrying into a manuscript.

## Executive verdict

The results are publishable, but not as the simple story that "more diversity wins."

The strongest story supported by the current data is:

> ChEMBL3D-PB is sample-efficient but saturates at finite bound-state coverage. High-budget
> learned generators, especially Qwen 1.7B FSQ with big-data pretraining, scale beyond that
> coverage ceiling. The winning signal is not raw geometric diversity; it is PB-valid,
> energy-window-retained, binding-relevant coverage.

This is a stronger and safer version of the initial hypothesis. The original H1 was directionally
right that collapsed/minimized ensembles miss useful bound-state coverage, but it was too broad
in treating normalized geometric diversity as the primary explanatory variable.

## What changed since the initial read

Two new analyses materially improved the interpretation.

| New test | What it adds | Main consequence |
| --- | --- | --- |
| ChEMBL3D-PB K-efficiency | ChEMBL now has a matched K recovery curve, not just a full-pool endpoint. | Qwen is not clearly better at low K. Qwen wins by scaling beyond ChEMBL's finite pool. |
| Relative-energy-windowed diversity | Recovery and cluster metrics are recomputed inside per-ligand MMFF relative-energy windows. | Best Qwen keeps most of its advantage at dE <= 20 kcal/mol, so the result is not mainly high-energy junk. |

These tests turn the story from "diversity helps" into "useful scalable coverage matters."

## Data used

Primary result files:

| File | Use |
| --- | --- |
| `data/results/casf_analysis_master.csv` | Aggregate method/tier summaries. |
| `data/results/casf_per_ligand_long.csv` | Per-ligand paired comparisons, denominators, PB-empty behavior. |
| `data/results/tables/extended_k_efficiency.csv` | K-efficiency curves. |
| `data/results/tables/extended_k_efficiency_comparisons.csv` | Qwen versus ChEMBL K-matched comparisons. |
| `data/results/tables/extended_energy_window_summary.csv` | Energy-windowed recovery and diversity. |
| `data/results/tables/extended_energy_window_comparisons.csv` | Hit losses under energy filtering and high-energy-hit fractions. |
| `data/results/tables/extended_pb_failure_mechanisms.csv` | PB failure mechanism summaries. |

Core-set universe: 94 ligands.

Ref-set universe in these artifacts: 1236 ligands. Qwen ref generations are still not present, so
ref-set analysis only reflects ChEMBL3D-PB and non-Qwen methods.

## How to read the newer tests

### K-efficiency

K-efficiency asks:

> If each ligand is allowed only K conformers from a method, how often does the method recover
> the CASF crystal pose?

For a ligand and method, the analysis selects up to K PB-passing conformers and recomputes best
RMSD and Hit@X from that subset.

For generated methods, the subset is a deterministic/random-subsample style sample from the
PB-passing generated pool. For ChEMBL3D-PB, the new test uses per-conformer ChEMBL conformer RMSDs
and PB status, so ChEMBL can now be evaluated at the same requested K values as generated methods.

The primary ChEMBL summary here is `random_k` with `capped_at_available`:

| Term | Meaning |
| --- | --- |
| `random_k` | Randomly sample up to K PB-passing ChEMBL conformers per ligand, averaged deterministically across repeats in the output table. |
| `capped_at_available` | If a ligand has fewer than K ChEMBL3D-PB conformers, use all available conformers rather than dropping the ligand. |
| `Mean used` | Mean number of conformers actually evaluated per ligand. For ChEMBL this saturates around 83 on core. |

This matters because ChEMBL has a finite pool. At K=250, K=500, and K=1000, ChEMBL is effectively
at its full available set for nearly all core ligands.

### Relative-energy-windowed diversity

Energy-windowed diversity asks:

> Does recovery survive if we only keep conformers within a plausible relative MMFF energy range?

For each ligand and method, conformer energy is converted to relative energy:

```text
dE_i = MMFF_energy_i - min_j(MMFF_energy_j for the same ligand and method)
```

Then recovery and diversity are recomputed after filtering to:

| Window | Interpretation |
| --- | --- |
| `All` | All PB-passing conformers. |
| `dE <= 20` | Low-to-moderate relative-energy conformers. This is the main plausibility window. |
| `dE <= 10` | Stricter relative-energy conformers. This is conservative and often harsh for raw/generated pools. |

Important: these are MMFF relative-energy windows, not binding free energies and not definitive
solution-phase thermodynamics. The result should be described as "MMFF-windowed plausibility," not
as proof that conformers are energetically populated in a protein pocket.

## Story component 1: high-budget core-set recovery

At the fixed high-budget tier, the best learned generators improve CASF crystal recovery over
ChEMBL3D-PB. The conservative Hit@0.75 denominator below is the full 94-ligand core universe;
PB-empty generated ligands count as misses.

| Method | Eval | Mean confs | Mean 1A clusters | Clusters/100 | Best RMSD | Hit@0.75 (94) | PB pass |
|:-------------------------------|:-------|-------------:|-------------------:|---------------:|------------:|----------------:|----------:|
| Qwen 1.7B fsq+bigdata-pretrain | 91/94 | 933.809 | 75.670 | 9.164 | 0.301 | 0.915 | 0.934 |
| Qwen 1.7B bigdata | 93/94 | 856.170 | 74.516 | 11.637 | 0.343 | 0.883 | 0.859 |
| Qwen 4B bigdata | 91/94 | 871.819 | 75.967 | 9.757 | 0.332 | 0.862 | 0.874 |
| Qwen 0.6B fsq+bigdata-pretrain | 91/94 | 926.351 | 76.538 | 9.385 | 0.334 | 0.862 | 0.927 |
| Torsional Diffusion | 90/94 | 839.745 | 56.644 | 8.584 | 0.372 | 0.851 | 0.840 |
| Qwen 0.6B bigdata | 93/94 | 851.160 | 74.581 | 11.094 | 0.357 | 0.851 | 0.853 |
| MCF drugs-L | 93/94 | 860.457 | 33.667 | 6.407 | 0.387 | 0.840 | 0.860 |
| NextMol DMT-L | 93/94 | 860.915 | 31.032 | 5.706 | 0.395 | 0.840 | 0.861 |
| LOQI | 90/94 | 957.064 | 20.700 | 2.071 | 0.341 | 0.819 | 0.957 |
| RDKit random (raw) | 94/94 | 999.947 | 35.883 | 3.589 | 0.447 | 0.809 | 1.000 |
| ChEMBL3D ground truth PB | 94/94 | 83.287 | 7.798 | 28.502 | 0.524 | 0.798 | NA |
| Torsion perturb (raw) | 94/94 | 896.447 | 74.149 | 8.820 | 0.497 | 0.777 | 0.896 |
| RDKit random (minimized) | 94/94 | 999.968 | 20.277 | 2.028 | 0.528 | 0.766 | 1.000 |
| Torsion perturb (minimized) | 94/94 | 999.968 | 25.011 | 2.501 | 0.554 | 0.755 | 1.000 |

Main read:

| Observation | Interpretation |
| --- | --- |
| Best Qwen reaches 91.5% Hit@0.75 versus 79.8% for ChEMBL3D-PB. | Strong high-budget recovery result. |
| Best Qwen mean best RMSD is 0.301 A versus 0.524 A for ChEMBL3D-PB. | The result is not only threshold crossing; it improves continuous best RMSD. |
| Torsion raw has Qwen-like cluster count but worse than ChEMBL3D-PB Hit@0.75. | Raw diversity is not enough. |
| LOQI has low cluster count but strong RMSD and decent Hit@0.75. | A learned prior can be useful without maximal diversity. |

Tight-threshold behavior is also favorable for best Qwen:

| Method | Hit@0.25 | Hit@0.5 | Hit@0.75 | Hit@2.0 | Best RMSD |
|:-------------------------------|-----------:|----------:|-----------:|----------:|------------:|
| ChEMBL3D ground truth PB | 0.255 | 0.638 | 0.798 | 0.979 | 0.524 |
| Qwen 1.7B fsq+bigdata-pretrain | 0.489 | 0.798 | 0.915 | 0.968 | 0.301 |
| Qwen 1.7B bigdata | 0.521 | 0.723 | 0.883 | 0.989 | 0.343 |
| Qwen 4B bigdata | 0.521 | 0.723 | 0.862 | 0.968 | 0.332 |
| Qwen 0.6B fsq+bigdata-pretrain | 0.468 | 0.766 | 0.862 | 0.957 | 0.334 |
| LOQI | 0.532 | 0.787 | 0.819 | 0.957 | 0.341 |
| Torsional Diffusion | 0.415 | 0.755 | 0.851 | 0.957 | 0.372 |
| RDKit random (raw) | 0.309 | 0.628 | 0.809 | 1.000 | 0.447 |
| Torsion perturb (raw) | 0.309 | 0.564 | 0.777 | 0.989 | 0.497 |

The tight cutoff result is especially useful for a manuscript. Best Qwen improves Hit@0.25 and
Hit@0.5 substantially, suggesting a precision gain rather than only coarse basin recovery. Hit@2.0
is saturated for most methods and is less discriminating.

## Story component 2: ChEMBL is sample-efficient, Qwen scales beyond it

The K-efficiency result is one of the most important new additions. It weakens a simplistic
"Qwen is better at any budget" claim, but it creates a cleaner scaling story.

| K | ChEMBL3D-PB random-K | Qwen 1.7B FSQ+bigdata | LOQI | Torsional Diffusion | RDKit raw | Torsion raw |
|---------:|-----------------------:|------------------------:|-------:|----------------------:|------------:|--------------:|
| 10 | 0.641 | 0.573 | 0.622 | 0.563 | 0.549 | 0.480 |
| 25 | 0.730 | 0.676 | 0.696 | 0.654 | 0.617 | 0.553 |
| 50 | 0.771 | 0.746 | 0.734 | 0.701 | 0.665 | 0.609 |
| 100 | 0.790 | 0.787 | 0.762 | 0.754 | 0.707 | 0.658 |
| 250 | 0.798 | 0.852 | 0.793 | 0.807 | 0.756 | 0.724 |
| 500 | 0.798 | 0.880 | 0.805 | 0.840 | 0.787 | 0.760 |
| 1000 | 0.798 | 0.915 | 0.819 | 0.851 | 0.809 | 0.777 |

Main read:

| K region | What happens | Interpretation |
| --- | --- | --- |
| K <= 50 | ChEMBL3D-PB is best or near-best. | Existing ChEMBL conformers are highly sample-efficient. |
| K ~= 100 | Best Qwen and ChEMBL3D-PB are essentially tied. | This is the crossover region. |
| K >= 250 | Best Qwen clearly exceeds ChEMBL3D-PB. | Qwen's advantage is scalability beyond ChEMBL's finite pool. |

This is not a weakness if framed correctly. The result says ChEMBL3D-PB is a strong, compact
reference source, but it saturates because it has only about 83 PB-passing conformers per core
ligand. Learned generation can keep sampling and find additional bound-like geometries.

Manuscript-safe claim:

> ChEMBL3D-PB is more sample-efficient at low K, while learned generation provides better
> high-budget coverage after the ChEMBL pool saturates.

Avoid:

> Qwen is uniformly more sample-efficient than ChEMBL3D.

## Story component 3: energy-windowed recovery supports useful coverage

Energy-windowed analysis addresses a key reviewer concern:

> Are generated hits coming from unrealistic high-energy tails?

For the best Qwen variant, the answer is mostly no at dE <= 20 kcal/mol.

| Method | Window | Ligands w/conf | Mean confs in window | Mean clusters | Best RMSD | Hit@0.75 |
|:-------------------------------|:---------|:-----------------|-----------------------:|----------------:|------------:|-----------:|
| ChEMBL3D-PB | All | 94/94 | 83.287 | 7.798 | 0.524 | 0.798 |
| ChEMBL3D-PB | dE<=20 | 94/94 | 81.596 | 7.628 | 0.528 | 0.787 |
| ChEMBL3D-PB | dE<=10 | 94/94 | 68.202 | 6.638 | 0.554 | 0.766 |
| Qwen 1.7B fsq+bigdata-pretrain | All | 91/94 | 933.809 | 75.670 | 0.301 | 0.915 |
| Qwen 1.7B fsq+bigdata-pretrain | dE<=20 | 91/94 | 756.936 | 42.055 | 0.317 | 0.894 |
| Qwen 1.7B fsq+bigdata-pretrain | dE<=10 | 91/94 | 499.372 | 18.352 | 0.392 | 0.798 |
| Qwen 0.6B fsq+bigdata-pretrain | All | 91/94 | 926.351 | 76.538 | 0.334 | 0.862 |
| Qwen 0.6B fsq+bigdata-pretrain | dE<=20 | 91/94 | 744.149 | 40.714 | 0.350 | 0.840 |
| Qwen 0.6B fsq+bigdata-pretrain | dE<=10 | 91/94 | 481.287 | 18.198 | 0.412 | 0.819 |
| Qwen 1.7B bigdata | All | 93/94 | 856.170 | 74.516 | 0.343 | 0.883 |
| Qwen 1.7B bigdata | dE<=20 | 93/94 | 489.606 | 26.817 | 0.385 | 0.830 |
| Qwen 1.7B bigdata | dE<=10 | 93/94 | 185.670 | 9.892 | 0.487 | 0.734 |
| Qwen 4B bigdata | All | 91/94 | 871.819 | 75.967 | 0.332 | 0.862 |
| Qwen 4B bigdata | dE<=20 | 91/94 | 472.160 | 25.934 | 0.373 | 0.819 |
| Qwen 4B bigdata | dE<=10 | 91/94 | 173.543 | 10.066 | 0.504 | 0.734 |
| NextMol DMT-L | All | 93/94 | 860.915 | 31.032 | 0.395 | 0.840 |
| NextMol DMT-L | dE<=20 | 93/94 | 797.213 | 23.065 | 0.401 | 0.840 |
| NextMol DMT-L | dE<=10 | 93/94 | 652.213 | 14.548 | 0.448 | 0.819 |
| LOQI | All | 90/94 | 957.064 | 20.700 | 0.341 | 0.819 |
| LOQI | dE<=20 | 90/94 | 909.213 | 17.467 | 0.347 | 0.819 |
| LOQI | dE<=10 | 90/94 | 770.266 | 11.233 | 0.409 | 0.787 |
| Torsional Diffusion | All | 90/94 | 839.745 | 56.644 | 0.372 | 0.851 |
| Torsional Diffusion | dE<=20 | 90/94 | 432.894 | 11.144 | 0.469 | 0.787 |
| Torsional Diffusion | dE<=10 | 90/94 | 212.840 | 5.656 | 0.605 | 0.702 |
| RDKit random (raw) | All | 94/94 | 999.947 | 35.883 | 0.447 | 0.809 |
| RDKit random (raw) | dE<=20 | 94/94 | 527.330 | 11.702 | 0.489 | 0.745 |
| RDKit random (raw) | dE<=10 | 94/94 | 238.234 | 5.457 | 0.598 | 0.691 |
| Torsion perturb (raw) | All | 94/94 | 896.447 | 74.149 | 0.497 | 0.777 |
| Torsion perturb (raw) | dE<=20 | 94/94 | 357.021 | 12.457 | 0.594 | 0.723 |
| Torsion perturb (raw) | dE<=10 | 94/94 | 182.234 | 5.298 | 0.727 | 0.606 |
| RDKit random (minimized) | All | 94/94 | 999.968 | 20.277 | 0.528 | 0.766 |
| RDKit random (minimized) | dE<=20 | 94/94 | 975.532 | 17.564 | 0.534 | 0.766 |
| RDKit random (minimized) | dE<=10 | 94/94 | 902.043 | 14.011 | 0.575 | 0.713 |
| Torsion perturb (minimized) | All | 94/94 | 999.968 | 25.011 | 0.554 | 0.755 |
| Torsion perturb (minimized) | dE<=20 | 94/94 | 983.043 | 22.404 | 0.560 | 0.755 |
| Torsion perturb (minimized) | dE<=10 | 94/94 | 913.053 | 15.777 | 0.595 | 0.745 |

Key positives:

| Result | Why it matters |
| --- | --- |
| Best Qwen drops only from 91.5% to 89.4% Hit@0.75 at dE <= 20. | The high-budget gain mostly survives energy filtering. |
| Best Qwen dE <= 20 best RMSD is 0.317 A, still much better than ChEMBL3D-PB at 0.528 A. | Continuous RMSD supports the same conclusion. |
| Qwen 0.6B FSQ+bigdata and NextMol are strong at dE <= 10. | Strict energy windows reveal other robust learned methods, not just the top all-conformer Qwen. |
| Torsion raw collapses from 77.7% all to 60.6% at dE <= 10. | Naive diversity relies heavily on high-energy coverage. |

Key negatives and cautions:

| Result | Why it matters |
| --- | --- |
| Best Qwen drops from 89.4% at dE <= 20 to 79.8% at dE <= 10. | Very strict MMFF filtering erodes the headline advantage. |
| Qwen 1.7B bigdata and Qwen 4B bigdata are much weaker at dE <= 10. | Not all Qwen variants produce equally energy-window-robust hits. |
| NextMol and Qwen 0.6B FSQ+bigdata tie or beat best Qwen at dE <= 10 Hit@0.75. | The "best" method depends on whether the manuscript emphasizes high-budget all/PB coverage or strict energy plausibility. |

Energy-window hit losses and high-energy-hit fractions:

| Method | Comparison | Delta Hit@0.75 | High-energy hit fraction |
|:-------------------------------|:-----------------------|-----------------:|---------------------------:|
| Qwen 1.7B fsq+bigdata-pretrain | all_vs_deltaE_20 | -0.021 | 0.023 |
| Qwen 1.7B fsq+bigdata-pretrain | deltaE_20_vs_deltaE_10 | -0.096 | 0.023 |
| Qwen 1.7B bigdata | all_vs_deltaE_20 | -0.053 | 0.060 |
| Qwen 1.7B bigdata | deltaE_20_vs_deltaE_10 | -0.096 | 0.060 |
| Qwen 4B bigdata | all_vs_deltaE_20 | -0.043 | 0.049 |
| Qwen 4B bigdata | deltaE_20_vs_deltaE_10 | -0.085 | 0.049 |
| Qwen 0.6B fsq+bigdata-pretrain | all_vs_deltaE_20 | -0.021 | 0.025 |
| Qwen 0.6B fsq+bigdata-pretrain | deltaE_20_vs_deltaE_10 | -0.021 | 0.025 |
| LOQI | all_vs_deltaE_20 | 0.000 | 0.000 |
| LOQI | deltaE_20_vs_deltaE_10 | -0.032 | 0.000 |
| Torsional Diffusion | all_vs_deltaE_20 | -0.064 | 0.075 |
| Torsional Diffusion | deltaE_20_vs_deltaE_10 | -0.085 | 0.075 |
| RDKit random (raw) | all_vs_deltaE_20 | -0.064 | 0.079 |
| RDKit random (raw) | deltaE_20_vs_deltaE_10 | -0.053 | 0.079 |
| Torsion perturb (raw) | all_vs_deltaE_20 | -0.053 | 0.068 |
| Torsion perturb (raw) | deltaE_20_vs_deltaE_10 | -0.117 | 0.068 |
| ChEMBL3D-PB | all_vs_deltaE_20 | -0.011 | 0.013 |
| ChEMBL3D-PB | deltaE_20_vs_deltaE_10 | -0.021 | 0.013 |

The best Qwen high-energy-hit fraction is low. That supports a manuscript claim that the result is
not mainly driven by extreme energy artifacts. The dE <= 10 erosion should still be reported because
it is scientifically informative and will build trust.

## Story component 4: diversity helps descriptively, but the original metric is not causal enough

The initial hypothesis centered on geometric diversity, especially `clusters_per_100_1p0`.
The current results support a weaker version:

> Collapsed ensembles are bad, and broad coverage helps at high budget, but normalized cluster
> density is not a reliable predictor of bound-state recovery.

Across fixed method rows, diversity metrics correlate with Hit@0.75. But within a ligand, after
centering by ligand, the association is weak.

| Metric | Across-method corr with Hit@0.75 | Within-ligand centered corr |
|:-------------------------|-----------------------------------:|------------------------------:|
| Mean 1A clusters | 0.591 | 0.060 |
| Clusters/100 | 0.525 | -0.032 |
| Entropy | 0.271 | 0.128 |
| Largest cluster fraction | -0.449 | -0.197 |
| Pairwise mean RMSD | 0.445 | 0.136 |

Interpretation:

| Observation | Reading |
| --- | --- |
| Across-method cluster count is positively correlated with recovery. | Broad high-budget coverage tracks method performance descriptively. |
| Within-ligand `clusters_per_100` is near zero/negative. | The originally proposed normalized metric should not be the central mechanistic claim. |
| Largest-cluster fraction is negatively associated with recovery. | Collapse into one dominant mode is a more convincing negative signal than raw cluster density is a positive signal. |

Raw versus minimized ablations still support the single-basin/collapse concern:

| Comparison | Delta clusters | Delta best RMSD | Delta Hit@0.75 | Delta energy median |
|:------------------------|-----------------:|------------------:|-----------------:|----------------------:|
| RDKit raw - minimized | 15.606 | -0.081 | 0.043 | 54.068 |
| Torsion raw - minimized | 49.138 | -0.057 | 0.021 | 47.502 |

These ablations are worth including because they are easy to understand: minimization lowers energy
and reduces diversity, but also worsens best RMSD and Hit@0.75 at high K.

However, the torsion raw result is the counterexample that prevents overclaiming. Torsion raw has
very high cluster count but poor recovery, and it degrades badly under energy filtering. This is
evidence for "useful diversity," not "maximum diversity."

## Story component 5: checkpoint and representation effects are interesting but not causal proof

The Qwen variants suggest that training corpus/representation choices matter more than simple
parameter count.

| Contrast | Hit A | Hit B | Best RMSD A | Best RMSD B |
|:----------------------------------------|--------:|--------:|--------------:|--------------:|
| 1.7B FSQ + bigdata pretrain vs 1.7B FSQ | 0.915 | 0.830 | 0.301 | 0.441 |
| 1.7B bigdata vs 1.7B revisited | 0.883 | 0.798 | 0.343 | 0.456 |
| 4B bigdata vs 4B revisited | 0.862 | 0.787 | 0.332 | 0.460 |
| 1.7B FSQ+bigdata vs 4B bigdata | 0.915 | 0.862 | 0.301 | 0.332 |

Main read:

| Observation | Interpretation |
| --- | --- |
| FSQ + bigdata pretraining improves over plain FSQ. | Strong internal lead that representation/pretraining matters. |
| Bigdata variants beat revisited variants at both 1.7B and 4B. | Training data/checkpoint distribution likely matters. |
| Best 1.7B beats 4B bigdata. | No clean monotonic parameter scaling story. |

This can be included, but as an exploratory model-family result unless Qwen ref validation confirms it.

## Story component 6: validity and PB failure modes are important limitations

PoseBusters is not only a filter; it also exposes systematic method failures.

| Method | PB fail rate | Median ligand PB pass | Ligands PB pass <0.5 | Dominant failure | Frac chirality failures | Frac clash failures | Frac energy failures |
|:-------------------------------|---------------:|------------------------:|-----------------------:|:----------------------|--------------------------:|----------------------:|-----------------------:|
| Qwen 1.7B fsq+bigdata-pretrain | 0.066 | 0.996 | 4 | tetrahedral_chirality | 0.693 | 0.273 | 0.023 |
| Qwen 1.7B bigdata | 0.141 | 0.958 | 6 | tetrahedral_chirality | 0.346 | 0.297 | 0.085 |
| Qwen 4B bigdata | 0.126 | 0.949 | 6 | internal_steric_clash | 0.292 | 0.333 | 0.091 |
| LOQI | 0.043 | 1.000 | 4 | tetrahedral_chirality | 1.000 | 0.000 | 0.000 |
| Torsional Diffusion | 0.160 | 0.950 | 7 | internal_steric_clash | 0.215 | 0.475 | 0.260 |
| MCF drugs-L | 0.140 | 0.999 | 13 | tetrahedral_chirality | 0.740 | 0.206 | 0.022 |
| NextMol DMT-L | 0.139 | 1.000 | 13 | tetrahedral_chirality | 0.766 | 0.217 | 0.002 |
| Torsion perturb (raw) | 0.104 | 0.938 | 0 | energy_ratio | 0.000 | 0.000 | 1.000 |

Important negative result:

| Issue | Interpretation |
| --- | --- |
| Best Qwen still has PB-empty ligands. | Oversampling helps stochastic failure, but not systematic chirality failure. |
| Several learned methods have chirality-dominated failures. | This should be reported directly, not hidden. |
| Torsion raw failures are energy-ratio dominated. | This matches the energy-window story: raw diversity includes implausible conformers. |

For the manuscript, treat PB failures as part of the result: a scalable conformer source must
optimize both recovery and validity. This strengthens the benchmark framing.

## Ref-set reflection

Qwen ref generations are still absent, so the ref set cannot validate the main Qwen claim. It can,
however, test whether the non-Qwen baseline behavior seen on core is anomalous. It is not.

### Ref fixed-pool result

| Method | Eval | Mean confs | Mean 1A clusters | Best RMSD | Hit@0.75 | PB pass |
|:----------------------------|:----------|-------------:|-------------------:|------------:|-----------:|----------:|
| LOQI | 1206/1236 | 970.159 | 23.943 | 0.386 | 0.846 | 0.970 |
| Torsional Diffusion | 1203/1236 | 842.552 | 47.188 | 0.391 | 0.840 | 0.843 |
| MCF drugs-L | 1227/1236 | 869.816 | 27.539 | 0.454 | 0.835 | 0.870 |
| NextMol DMT-L | 1226/1236 | 874.909 | 27.805 | 0.458 | 0.833 | 0.875 |
| RDKit random (raw) | 1236/1236 | 999.843 | 38.146 | 0.477 | 0.815 | 1.000 |
| ChEMBL3D ground truth PB | 1236/1236 | 67.128 | 7.348 | 0.551 | 0.790 | NA |
| RDKit random (minimized) | 1235/1236 | 998.799 | 23.506 | 0.537 | 0.769 | 0.999 |
| Torsion perturb (raw) | 1236/1236 | 897.395 | 61.576 | 0.491 | 0.766 | 0.897 |
| Torsion perturb (minimized) | 1235/1236 | 999.170 | 23.977 | 0.549 | 0.761 | 0.999 |

Ref confirms several non-Qwen patterns:

| Pattern | Core | Ref |
| --- | --- | --- |
| LOQI is strong despite low cluster count. | Yes | Yes |
| Torsional Diffusion is strong at high K. | Yes | Yes |
| RDKit raw is respectable but not top. | Yes | Yes |
| Torsion raw is highly diverse but not strong. | Yes | Yes |
| Minimized baselines underperform raw RDKit at high K. | Yes | Yes |

### Ref K-efficiency

| K | ChEMBL3D-PB random-K | LOQI | Torsional Diffusion | MCF drugs-L | NextMol DMT-L | RDKit raw | Torsion raw |
|---------:|-----------------------:|-------:|----------------------:|--------------:|----------------:|------------:|--------------:|
| 25 | 0.722 | 0.682 | 0.658 | 0.666 | 0.668 | 0.634 | 0.567 |
| 50 | 0.764 | 0.726 | 0.713 | 0.721 | 0.723 | 0.687 | 0.620 |
| 100 | 0.785 | 0.761 | 0.756 | 0.760 | 0.763 | 0.727 | 0.667 |
| 250 | 0.790 | 0.798 | 0.802 | 0.800 | 0.800 | 0.769 | 0.719 |
| 500 | 0.790 | 0.823 | 0.827 | 0.821 | 0.821 | 0.794 | 0.747 |
| 1000 | 0.790 | 0.846 | 0.840 | 0.835 | 0.833 | 0.815 | 0.766 |

This mirrors core: ChEMBL is strong at low K and saturates around 79%; high-budget generated
methods exceed it once K is large enough. This is strong support for the "database-efficient but
finite, generators scalable" story even without Qwen ref.

### Ref energy-windowed result

| Method | Window | Hit@0.75 | Best RMSD | Mean clusters | Mean confs in window |
|:-------------------------|:---------|-----------:|------------:|----------------:|-----------------------:|
| ChEMBL3D ground truth PB | All | 0.790 | 0.551 | 7.348 | 67.128 |
| ChEMBL3D ground truth PB | dE<=20 | 0.788 | 0.559 | 7.113 | 65.054 |
| ChEMBL3D ground truth PB | dE<=10 | 0.761 | 0.587 | 6.300 | 54.528 |
| LOQI | All | 0.846 | 0.386 | 23.943 | 970.159 |
| LOQI | dE<=20 | 0.838 | 0.409 | 17.348 | 923.458 |
| LOQI | dE<=10 | 0.807 | 0.462 | 10.761 | 803.923 |
| Torsional Diffusion | All | 0.840 | 0.391 | 47.188 | 842.552 |
| Torsional Diffusion | dE<=20 | 0.790 | 0.491 | 11.332 | 441.918 |
| Torsional Diffusion | dE<=10 | 0.697 | 0.625 | 5.652 | 214.794 |
| MCF drugs-L | All | 0.835 | 0.454 | 27.539 | 869.816 |
| MCF drugs-L | dE<=20 | 0.824 | 0.478 | 19.406 | 784.281 |
| MCF drugs-L | dE<=10 | 0.790 | 0.530 | 12.028 | 578.767 |
| NextMol DMT-L | All | 0.833 | 0.458 | 27.805 | 874.909 |
| NextMol DMT-L | dE<=20 | 0.825 | 0.468 | 21.769 | 823.752 |
| NextMol DMT-L | dE<=10 | 0.794 | 0.519 | 13.569 | 674.634 |
| RDKit random (raw) | All | 0.815 | 0.477 | 38.146 | 999.843 |
| RDKit random (raw) | dE<=20 | 0.757 | 0.549 | 12.557 | 518.554 |
| RDKit random (raw) | dE<=10 | 0.673 | 0.660 | 5.694 | 221.486 |
| Torsion perturb (raw) | All | 0.766 | 0.491 | 61.576 | 897.395 |
| Torsion perturb (raw) | dE<=20 | 0.676 | 0.616 | 10.371 | 373.104 |
| Torsion perturb (raw) | dE<=10 | 0.593 | 0.735 | 4.845 | 200.662 |

Ref strengthens the energy-window story for non-Qwen methods:

| Method family | Ref energy-window behavior |
| --- | --- |
| LOQI | Strong retention under dE <= 20 and still good under dE <= 10. |
| MCF / NextMol | Strong retention under dE <= 20 and competitive under dE <= 10. |
| Torsional Diffusion | Strong all-conformer recovery but substantial degradation under energy filtering. |
| RDKit raw and torsion raw | Large degradation under energy filtering, especially torsion raw. |

This is useful even before Qwen ref results because it shows the new energy-window test is not
just overfitting the 94 core ligands.

## Revised hypothesis assessment

| Original hypothesis component | Current result | Decision |
| --- | --- | --- |
| More geometrically diverse ensembles recover CASF-bound conformers more often. | Broadly true across methods at high budget, but weak within-ligand and contradicted by torsion raw. | Partially supported; revise. |
| Single-basin/minimized ensembles hurt recovery. | RDKit and torsion raw-vs-minimized ablations support this. | Supported. |
| Clusters/100 is the primary cross-method diversity metric. | It correlates across methods but fails within-ligand centered analysis. | Not manuscript-safe as primary mechanism. |
| High-budget generation helps beyond ChEMBL3D. | Best Qwen substantially exceeds ChEMBL3D-PB at K >= 250 and fixed high budget. | Supported on core. |
| ChEMBL-count comparisons may shrink diversity advantages. | K-efficiency confirms ChEMBL is strong at low K and near K=100. | Supported. |
| Energy spread identifies useful basins. | Raw energy std is not reliable, but relative-energy-windowed analysis is useful. | Revise to energy-windowed plausibility. |
| Training data or representation affects learned generators. | Qwen contrasts strongly suggest this, but ref validation is not available. | Strong lead, not causal proof. |

Revised hypothesis:

> At high sampling budgets, PB-valid learned generators can exceed database-limited ChEMBL3D
> bound-state coverage, but only when their diversity is chemically plausible and binding-relevant.
> Minimization and naive perturbation fail in opposite ways: minimization collapses useful coverage,
> while naive raw diversity produces many conformers that do not survive energy-windowed scrutiny.

## What is actually worth including in a JCIM paper

### Strong inclusion candidates

| Result block | Include? | Why |
| --- | --- | --- |
| Core high-budget fixed-pool table | Yes | Main discovery result. |
| ChEMBL K-efficiency versus best Qwen and non-Qwen generators | Yes | Defines the central "efficient but finite vs scalable" story. |
| Energy-windowed recovery at All, dE <= 20, dE <= 10 | Yes | Addresses the strongest chemical plausibility concern. |
| Raw-vs-minimized ablations | Yes | Simple mechanistic support for the collapse argument. |
| Torsion raw counterexample | Yes | Prevents overclaiming and makes the paper more credible. |
| Ref non-Qwen replication | Yes | Shows baseline patterns generalize beyond core. |
| PB failure mechanisms | Yes, concise | Validity is central for dataset construction. |

### Include cautiously

| Result block | Caution |
| --- | --- |
| Qwen checkpoint contrasts | Good lead, but post-selection and no ref validation yet. |
| Core-only Qwen superiority | Strong but should be called core-set until ref Qwen finishes. |
| dE <= 10 method ranking | Useful but should be secondary; MMFF strict windows may over-penalize bound-like strained conformers. |
| Diversity correlations | Use to motivate revision, not as primary proof. |

### Probably leave out or move to supplement

| Result block | Reason |
| --- | --- |
| Full list of all Qwen variants in main text | Too much detail; use a compact checkpoint contrast table. |
| Hit@2.0 as main endpoint | Saturated and not very discriminating. |
| Raw `energy_std` tables | Outlier-dominated; relative-energy windows are better. |
| Dynamic-tier Qwen flexibility claims | Qwen rotatable-bond metadata warning remains. |

## Suggested manuscript framing

### Possible title

1. `Scalable Learned Conformer Generation Improves Bound-State Ligand Geometry Recovery Beyond ChEMBL3D Coverage`
2. `Useful Diversity, Not Maximal Diversity, Drives Bound-State Conformer Recovery in CASF Ligands`
3. `Benchmarking Scalable Conformer Generators for Protein-Bound Ligand Geometry Recovery`

The first is the safest if Qwen ref validates. The third is safest before Qwen ref is complete.

### Possible abstract-level claim

> On the CASF core intersection, ChEMBL3D-PB is highly sample-efficient and reaches 79.8%
> Hit@0.75 with a mean of 83 PB-passing conformers per ligand. High-budget learned generation
> surpasses this ceiling, with the best Qwen variant reaching 91.5% Hit@0.75 and 0.301 A mean
> best RMSD. K-efficiency analysis shows that ChEMBL3D-PB remains competitive at low K, whereas
> learned generation wins after the ChEMBL pool saturates. Relative MMFF energy-windowed analysis
> shows that most of the best Qwen gain persists within dE <= 20 kcal/mol, while naive perturbation
> loses much of its apparent diversity and recovery under energy filtering. These results suggest
> that scalable conformer dataset construction should optimize for PB-valid, chemically plausible,
> binding-relevant coverage rather than minimization or raw geometric diversity alone.

### Recommended main-text result order

1. Define the benchmark: CASF crystal as bound-state proxy, PB-passing conformer pools.
2. Show ChEMBL3D-PB is compact and strong, but finite.
3. Show high-budget learned generation exceeds ChEMBL3D-PB on core.
4. Show K-efficiency: ChEMBL wins/ties at low K, Qwen wins after saturation.
5. Show energy-windowed analysis: best Qwen gain mostly survives dE <= 20.
6. Show why the original simple diversity hypothesis is incomplete: torsion raw and LOQI.
7. Show raw-vs-minimized ablations for single-basin collapse.
8. Show ref-set non-Qwen replication.
9. Discuss limitations: Qwen ref pending, PB chirality failures, MMFF energy limitations, Qwen metadata warning.

## Final recommendation

There is a good JCIM story here if the team embraces the revised interpretation:

> The important contribution is not that diversity alone predicts bound-pose recovery. The important
> contribution is a benchmark showing that ChEMBL3D conformers are efficient but coverage-limited,
> while learned generators can scale PB-valid, energy-window-retained coverage beyond that limit.

The paper will be strongest if Qwen ref generations confirm the core-set finding. Even before that,
the non-Qwen ref results are useful because they validate the broader pattern: high-budget learned
or learned-like generators can exceed ChEMBL's finite endpoint, while naive raw diversity degrades
under energy filtering.

## Pre-submission checks still worth doing

| Check | Why |
| --- | --- |
| Complete Qwen ref generation and rerun this report's key tables. | Needed to validate the main Qwen claim beyond core. |
| Repair or exclude Qwen rotatable-bond metadata. | Current warning prevents fair Qwen flexibility-stratum analysis. |
| Run stale manifest/SDF audit before freezing numbers. | Current sanity table still warns this was skipped. |
| Manually inspect PB-empty/chirality-failure ligands. | Important for trust in PoseBusters failure interpretation. |
| Decide whether dE <= 20 or dE <= 10 is the main energy-window result. | dE <= 20 supports the scalable coverage story; dE <= 10 is stricter and changes method ranking. |

