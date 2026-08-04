# CASF16 Conformer Generator Results

## Summary

The current CASF16 analysis identifies **LOQI as the strongest validated generator for constructing a large conformer training corpus**. LOQI is not the most geometrically expansive method; its advantage comes from the best observed combination of:

- strict CASF-bound-pose recovery,
- K-efficiency among generated methods,
- high PoseBusters pass rate,
- strong paired improvement over ChEMBL3D-PB,
- measurable improvement over high-K RDKit raw,
- and the clearest gains in hard, flexible, and large ligands.

ChEMBL3D-PB remains a strong compact reference. Its PB filtering removes conformers rather than ligands, and its CASF recovery barely changes after filtering. LOQI's value is that it densifies conformational regions that compact ChEMBL3D-PB under-covers, while remaining cleaner than other learned high-coverage generators.

This result establishes generator selection for dataset construction. Downstream training superiority requires a separate model-training ablation.

## Final Method Ranking

Ref fixed-tier results:

Column notes: `confs/ligand` is the mean number of PB-passing conformers used in the evaluated set. `PB pass` is the fraction of generated/input conformers that passed PoseBusters. `1.0 A clusters` is the mean number of geometric clusters per ligand at a 1.0 A RMSD threshold. `best RMSD` is the mean ligand-level best RMSD to the CASF bound pose. `valid-hit score` is approximately `PB pass x Hit@0.75`, so it penalizes methods that recover bound-like poses but generate many invalid conformers.

| method | ligands | confs/ligand | PB pass | 1.0 A clusters | best RMSD | Hit@0.75 | Hit@2.0 | valid-hit score |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| LOQI | 1232 | 970 | 0.970 | 21.9 | 0.384 | 0.869 | 0.988 | 0.843 |
| torsional diffusion | 1224 | 843 | 0.843 | 47.2 | 0.390 | 0.864 | 0.996 | 0.728 |
| MCF drugs-L | 1224 | 871 | 0.871 | 27.7 | 0.450 | 0.842 | 0.972 | 0.733 |
| NextMol DMT-L | 1224 | 875 | 0.875 | 26.6 | 0.457 | 0.840 | 0.974 | 0.735 |
| RDKit raw | 1224 | 1000 | 1.000 | 36.2 | 0.476 | 0.815 | 0.999 | 0.815 |
| torsion raw | 1224 | 898 | 0.898 | 60.9 | 0.490 | 0.766 | 0.994 | 0.688 |
| RDKit minimized | 1224 | 999 | 0.999 | 21.6 | 0.537 | 0.770 | 0.997 | 0.769 |
| ChEMBL3D-PB | 1235 | 67 | 1.000 | 7.3 | 0.550 | 0.791 | 0.976 | 0.791 |
| torsion minimized | 1224 | 999 | 0.999 | 22.4 | 0.550 | 0.762 | 0.993 | 0.761 |

Result interpretation:

- LOQI has the best ref-set best RMSD and the highest Hit@0.75 among validated ref-set generators.
- Torsional diffusion nearly matches LOQI in raw recovery, but its PB pass rate is much lower.
- RDKit raw is the strongest classical control: almost perfectly PB-clean, high Hit@2.0, and competitive strict recovery.
- Raw torsion has the highest cluster count but weak strict recovery, showing that geometric spread alone is not sufficient.
- Minimized ensembles are chemically clean but weaker in strict bound-pose recovery.

## K-Efficiency

Generated-method K-efficiency is now computed from cached per-conformer CASF RMSDs. ChEMBL3D-PB K rows remain unavailable because the dashboard currently exposes per-ligand PB-filtered best RMSD for ChEMBL3D-PB rather than per-conformer PB-passing RMSDs.

Ref set Hit@0.75 by requested K:

Column notes: each K column reports the mean Hit@0.75 after deterministically subsampling up to K PB-passing conformers per ligand and taking the ligand-level best RMSD within that subset. Higher values mean a method places useful conformers earlier in the sample stream.

| method | K=1 | K=5 | K=10 | K=50 | K=100 | K=250 | K=1000 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| LOQI | 0.340 | 0.547 | 0.626 | 0.746 | 0.782 | 0.820 | 0.869 |
| torsional diffusion | 0.290 | 0.495 | 0.581 | 0.732 | 0.777 | 0.824 | 0.864 |
| NextMol DMT-L | 0.297 | 0.497 | 0.582 | 0.728 | 0.770 | 0.808 | 0.840 |
| MCF drugs-L | 0.291 | 0.491 | 0.578 | 0.726 | 0.766 | 0.807 | 0.842 |
| RDKit raw | 0.285 | 0.469 | 0.544 | 0.687 | 0.726 | 0.769 | 0.815 |
| RDKit minimized | 0.312 | 0.499 | 0.570 | 0.686 | 0.715 | 0.744 | 0.770 |
| torsion raw | 0.259 | 0.419 | 0.489 | 0.619 | 0.666 | 0.719 | 0.766 |

Threshold crossing:

Column notes: each cell gives the smallest sampled K at which the method reaches the specified Hit@0.75 threshold. `not reached` means the method did not reach that threshold at any evaluated K.

| method | first K with Hit@0.75 >= 0.75 | first K with Hit@0.75 >= 0.80 | first K with Hit@0.75 >= 0.85 |
| --- | ---: | ---: | ---: |
| LOQI | 100 | 250 | 1000 |
| torsional diffusion | 100 | 250 | 500 |
| NextMol DMT-L | 100 | 250 | not reached |
| MCF drugs-L | 100 | 250 | not reached |
| RDKit raw | 250 | 1000 | not reached |
| RDKit minimized | 500 | not reached | not reached |
| torsion raw | 1000 | not reached | not reached |

Validity-weighted Hit@0.75 by K:

Column notes: values are `Hit@0.75 at K x mean PB pass rate`. This is not a replacement for reporting raw recovery and PB pass separately; it is a compact way to compare useful recovery after chemical-validity attrition.

| method | K=100 | K=250 | K=1000 |
| --- | ---: | ---: | ---: |
| LOQI | 0.777 | 0.816 | 0.865 |
| RDKit raw | 0.726 | 0.769 | 0.815 |
| RDKit minimized | 0.715 | 0.743 | 0.770 |
| NextMol DMT-L | 0.679 | 0.712 | 0.741 |
| torsional diffusion | 0.673 | 0.714 | 0.749 |
| MCF drugs-L | 0.672 | 0.708 | 0.739 |
| torsion raw | 0.598 | 0.645 | 0.688 |

Result interpretation:

- LOQI leads all generated methods from K=1 through K=100.
- Torsional diffusion slightly exceeds LOQI in raw Hit@0.75 at K=250 and K=500, but loses strongly after validity weighting.
- LOQI reaches Hit@0.75 >= 0.80 by K=250; RDKit raw reaches that threshold only at K=1000.
- LOQI is the best validity-weighted method at K=100, K=250, and K=1000.
- At K=1000, torsional diffusion has only 233 ref ligands with at least 1000 PB-passing conformers, compared with 1034 for LOQI and 1200 for RDKit raw. Its high-K curve therefore reflects stronger attrition.

## Paired Improvement Versus ChEMBL3D-PB

Ref fixed-tier paired comparison against ChEMBL3D-PB:

Column notes: this table uses only ligands present in both the method and ChEMBL3D-PB. `mean/median RMSD delta` is `method best RMSD - ChEMBL3D-PB best RMSD`; negative values mean the method is closer to the CASF bound pose. `wins >0.1 A` is the fraction of common ligands where the method improves RMSD by more than 0.1 A. `losses >0.1 A` is the fraction where it is worse by more than 0.1 A. `Hit@0.75 delta` is method Hit@0.75 minus ChEMBL3D-PB Hit@0.75 on the same ligand set. `large rescues >1 A` counts ligands improved by more than 1 A; `catastrophic losses >1 A` counts ligands worsened by more than 1 A.

| method | common ligands | mean RMSD delta | median RMSD delta | wins >0.1 A | losses >0.1 A | Hit@0.75 delta | large rescues >1 A | catastrophic losses >1 A |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| LOQI | 1202 | -0.166 | -0.095 | 0.477 | 0.064 | +0.078 | 33 | 6 |
| torsional diffusion | 1191 | -0.164 | -0.089 | 0.476 | 0.107 | +0.074 | 36 | 0 |
| MCF drugs-L | 1214 | -0.104 | -0.062 | 0.436 | 0.141 | +0.052 | 33 | 17 |
| NextMol DMT-L | 1213 | -0.097 | -0.041 | 0.407 | 0.153 | +0.049 | 30 | 13 |
| RDKit raw | 1223 | -0.075 | -0.011 | 0.316 | 0.253 | +0.025 | 37 | 0 |
| torsion raw | 1223 | -0.061 | -0.035 | 0.377 | 0.293 | -0.025 | 33 | 10 |
| RDKit minimized | 1222 | -0.013 | +0.011 | 0.237 | 0.326 | -0.021 | 35 | 2 |
| torsion minimized | 1222 | -0.001 | +0.015 | 0.232 | 0.331 | -0.029 | 30 | 12 |

Result interpretation:

- LOQI and torsional diffusion have almost identical mean improvements over ChEMBL3D-PB.
- LOQI has fewer routine losses than torsional diffusion and much higher PB pass rate.
- RDKit raw has the largest number of large rescues but a smaller aggregate Hit@0.75 gain because it also has many routine losses.
- Minimized sets do not improve strict recovery relative to ChEMBL3D-PB.

## Paired Improvement Versus RDKit Raw

Ref fixed-tier paired comparison against RDKit raw:

Column notes: this table uses only ligands present in both the method and RDKit raw. `mean/median RMSD delta vs RDKit raw` is `method best RMSD - RDKit raw best RMSD`; negative values mean the method is closer to the CASF bound pose. `wins >0.1 A`, `losses >0.1 A`, and `Hit@0.75 delta` are computed on that same paired ligand set.

| method | common ligands | mean RMSD delta vs RDKit raw | median RMSD delta | wins >0.1 A | losses >0.1 A | Hit@0.75 delta |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| LOQI | 1194 | -0.085 | -0.088 | 0.481 | 0.130 | +0.054 |
| torsional diffusion | 1191 | -0.082 | -0.059 | 0.404 | 0.137 | +0.050 |
| MCF drugs-L | 1215 | -0.026 | -0.054 | 0.402 | 0.172 | +0.027 |
| NextMol DMT-L | 1214 | -0.018 | -0.045 | 0.376 | 0.187 | +0.026 |
| torsion raw | 1224 | +0.014 | -0.003 | 0.322 | 0.311 | -0.049 |
| RDKit minimized | 1223 | +0.061 | +0.025 | 0.207 | 0.357 | -0.046 |
| torsion minimized | 1223 | +0.074 | +0.037 | 0.212 | 0.383 | -0.054 |

Result interpretation:

- LOQI and torsional diffusion improve over high-K RDKit raw on paired ligands.
- LOQI has the strongest median improvement over RDKit raw while retaining much better validity than torsional diffusion.
- RDKit raw remains a strong baseline because it is chemically clean and has high Hit@2.0, but it is less efficient at strict Hit@0.75.

## Bound-Pose Difficulty

Difficulty is defined by ChEMBL3D-PB best RMSD: easy <=0.5 A, moderate 0.5-0.75 A, hard 0.75-2.0 A, failed >2.0 A.

Ref set, paired against ChEMBL3D-PB:

Column notes: ligands are first binned by ChEMBL3D-PB best RMSD, then each method is compared on common ligands within that bin. `mean RMSD delta` is `method best RMSD - ChEMBL3D-PB best RMSD`; negative values indicate improvement. `rescue rate` is the fraction of baseline misses that become Hit@0.75 under the method. `catastrophic loss rate` is the fraction of ligands worsened by more than 1 A.

| baseline difficulty | method | ligands | mean RMSD delta | Hit@0.75 | rescue rate | catastrophic loss rate |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| easy | LOQI | 717 | -0.066 | 0.986 | 0.000 | 0.000 |
| easy | torsional diffusion | 708 | -0.039 | 0.983 | 0.000 | 0.000 |
| easy | RDKit raw | 727 | +0.078 | 0.945 | 0.000 | 0.000 |
| moderate | LOQI | 234 | -0.142 | 0.897 | 0.000 | 0.012 |
| moderate | torsional diffusion | 233 | -0.155 | 0.884 | 0.000 | 0.000 |
| moderate | RDKit raw | 241 | -0.066 | 0.809 | 0.000 | 0.000 |
| hard | LOQI | 221 | -0.422 | 0.570 | 0.560 | 0.000 |
| hard | torsional diffusion | 220 | -0.421 | 0.559 | 0.547 | 0.000 |
| hard | RDKit raw | 225 | -0.388 | 0.489 | 0.489 | 0.000 |
| failed | LOQI | 30 | -0.858 | 0.067 | 0.067 | 0.000 |
| failed | torsional diffusion | 30 | -1.290 | 0.133 | 0.133 | 0.000 |
| failed | RDKit raw | 30 | -1.525 | 0.200 | 0.200 | 0.000 |

Result interpretation:

- Easy ligands are already well covered by ChEMBL3D-PB.
- Hard ligands are the strongest recovery signal: LOQI and torsional diffusion rescue about 55% of ChEMBL3D-PB hard misses at Hit@0.75.
- Failed ligands behave differently. RDKit raw has the highest rescue rate in the failed bin, but even there strict recovery remains low.
- The largest aggregate value of learned diverse generators comes from ligands where compact ChEMBL3D-PB under-covers the bound pose.

## Flexible And Large-Ligand Strata

LOQI fixed versus ChEMBL3D-PB, paired by rotatable-bond bin:

Column notes: bins are computed from ligand properties, then LOQI and ChEMBL3D-PB are compared only on common ligands in each bin. Delta columns use `LOQI - ChEMBL3D-PB`, so negative RMSD deltas and positive Hit@0.75 deltas favor LOQI.

| rotatable-bond bin | common ligands | mean RMSD delta | median RMSD delta | wins >0.1 A | losses >0.1 A | Hit@0.75 delta |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0-4 | 509 | -0.133 | -0.070 | 0.418 | 0.018 | +0.029 |
| 5-9 | 563 | -0.150 | -0.096 | 0.478 | 0.101 | +0.082 |
| 10-14 | 121 | -0.361 | -0.222 | 0.686 | 0.083 | +0.256 |
| 15+ | 9 | -0.408 | -0.594 | 0.889 | 0.111 | +0.222 |

LOQI fixed and RDKit raw fixed versus ChEMBL3D-PB, paired by heavy-atom bin:

Column notes: delta columns are computed relative to ChEMBL3D-PB within each heavy-atom bin. Negative mean RMSD deltas indicate lower RMSD than ChEMBL3D-PB; positive Hit@0.75 deltas indicate higher strict recovery.

| heavy-atom bin | common ligands | LOQI mean delta | LOQI Hit@0.75 delta | RDKit raw mean delta | RDKit raw Hit@0.75 delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| <20 | 534 | -0.138 | +0.041 | -0.027 | +0.017 |
| 20-29 | 498 | -0.153 | +0.088 | -0.080 | +0.016 |
| 30-39 | 149 | -0.267 | +0.154 | -0.192 | +0.086 |
| 40+ | 21 | -0.473 | +0.238 | -0.389 | +0.000 |

K-resolved ref results for high-rotor ligands:

Column notes: these are K-efficiency summaries restricted to the indicated rotatable-bond strata. `K=1000 mean RMSD` is the mean best RMSD after using up to 1000 PB-passing conformers.

| stratum | method | K=10 Hit@0.75 | K=100 Hit@0.75 | K=1000 Hit@0.75 | K=1000 mean RMSD |
| --- | --- | ---: | ---: | ---: | ---: |
| 10-14 rotors | LOQI | 0.212 | 0.410 | 0.645 | 0.798 |
| 10-14 rotors | RDKit raw | 0.008 | 0.076 | 0.216 | 0.999 |
| 10-14 rotors | torsional diffusion | 0.006 | 0.053 | 0.139 | 1.249 |
| 15+ rotors | LOQI | 0.089 | 0.339 | 0.444 | 0.967 |
| 15+ rotors | RDKit raw | 0.000 | 0.000 | 0.000 | 1.591 |
| 15+ rotors | torsional diffusion | 0.000 | 0.000 | 0.000 | 1.845 |

K-resolved ref results for high-heavy-atom ligands:

Column notes: these are K-efficiency summaries restricted to the indicated heavy-atom strata. Low Hit@0.75 values in the 40+ bin indicate that strict bound-pose recovery remains hard even at high K.

| stratum | method | K=10 Hit@0.75 | K=100 Hit@0.75 | K=1000 Hit@0.75 | K=1000 mean RMSD |
| --- | --- | ---: | ---: | ---: | ---: |
| 30-39 heavy atoms | LOQI | 0.276 | 0.504 | 0.658 | 0.660 |
| 30-39 heavy atoms | RDKit raw | 0.109 | 0.349 | 0.586 | 0.745 |
| 30-39 heavy atoms | torsional diffusion | 0.132 | 0.397 | 0.611 | 0.735 |
| 40+ heavy atoms | LOQI | 0.024 | 0.152 | 0.333 | 1.130 |
| 40+ heavy atoms | RDKit raw | 0.002 | 0.010 | 0.095 | 1.214 |
| 40+ heavy atoms | torsional diffusion | 0.005 | 0.033 | 0.048 | 1.428 |

Result interpretation:

- LOQI's advantage grows with rotatable-bond count and heavy-atom count.
- The gain is not just an RMSD shift. In the 40+ heavy-atom bin, RDKit raw improves mean RMSD but does not improve Hit@0.75, while LOQI improves both.
- High-rotor and high-heavy-atom K curves show that LOQI's advantage appears across the sampled K range, not only at K=1000.

## ChEMBL3D Conformer-Count Bins

LOQI fixed versus ChEMBL3D-PB, paired by ChEMBL3D conformer-count bin:

Column notes: ligands are binned by how many PB-passing ChEMBL3D conformers they have. Delta columns compare LOQI to ChEMBL3D-PB within each bin, so this tests whether LOQI's gain persists even when ChEMBL3D-PB already has many conformers.

| ChEMBL3D conformer-count bin | common ligands | mean RMSD delta | Hit@0.75 delta | losses >0.1 A |
| --- | ---: | ---: | ---: | ---: |
| 1-9 | 378 | -0.166 | +0.042 | 0.000 |
| 10-49 | 295 | -0.198 | +0.081 | 0.054 |
| 50-99 | 199 | -0.161 | +0.106 | 0.095 |
| 100+ | 330 | -0.141 | +0.100 | 0.127 |

Result interpretation:

- LOQI improves every ChEMBL3D conformer-count bin.
- The improvement persists even when ChEMBL3D-PB has 100+ conformers.
- The benefit is therefore not explained only by ChEMBL3D-PB having too few conformers. The conformer distribution also matters.

## Rescue Cases

Largest ref rescue examples:

Column notes: `ChEMBL3D-PB RMSD` is the baseline best RMSD for that ligand. `method RMSD` is the method's best RMSD. `delta` is `method RMSD - ChEMBL3D-PB RMSD`; more negative values indicate larger rescues. `PB pass` and `clusters` are ligand-level method summaries.

| method | ligand | rot bonds | heavy atoms | ChEMBL3D-PB RMSD | method RMSD | delta | PB pass | clusters |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| RDKit raw | 6fo5 | 7 | 37 | 3.688 | 0.510 | -3.179 | 1.000 | 26 |
| RDKit raw | 2iuz | 3 | 28 | 2.420 | 0.198 | -2.222 | 1.000 | 8 |
| torsional diffusion | 2iuz | 3 | 28 | 2.420 | 0.300 | -2.120 | 0.950 | 19 |
| torsional diffusion | 6g2o | 4 | 27 | 2.149 | 0.144 | -2.005 | 0.973 | 5 |
| LOQI | 6g2o | 7 | 27 | 2.149 | 0.224 | -1.925 | 1.000 | 4 |
| LOQI | 6eya | 7 | 31 | 2.051 | 0.357 | -1.694 | 1.000 | 80 |
| LOQI | 6c7q | 9 | 34 | 1.832 | 0.222 | -1.610 | 1.000 | 18 |

Result interpretation:

- RDKit raw produces some of the largest individual rescues.
- LOQI and torsional diffusion produce strong rescues while also improving median paired behavior.
- Rescue behavior likely mixes generic high-K coverage with learned distributional bias.

## PoseBusters Failure Mechanisms

Ref fixed-tier PB failures:

Column notes: `PB fail rate` is the conformer-level failure fraction. `median pass` is the median ligand-level pass rate. `ligands <0.9 pass` and `<0.5 pass` count ligands with broad or severe PB attrition. `key fractions` summarize which PoseBusters checks dominate failed conformers; categories can overlap.

| method | PB fail rate | median pass | ligands <0.9 pass | ligands <0.5 pass | dominant failure | key fractions |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| RDKit raw | 0.0002 | 1.000 | 0 | 0 | internal clash | tiny absolute count |
| RDKit minimized | 0.0012 | 1.000 | 1 | 1 | bond angles / double bond stereo | 67.3% bond geometry |
| torsion minimized | 0.0008 | 1.000 | 1 | 1 | bond angles | 97.5% bond geometry |
| LOQI | 0.0298 | 1.000 | 43 | 34 | tetrahedral chirality | 95.0% chirality |
| torsion raw | 0.1025 | 0.952 | 469 | 18 | energy ratio | 100.0% energy_ratio |
| NextMol DMT-L | 0.1254 | 1.000 | 267 | 141 | tetrahedral chirality | 55.8% chirality, 29.8% clash |
| MCF drugs-L | 0.1291 | 0.999 | 270 | 146 | tetrahedral chirality | 51.8% chirality, 30.0% clash |
| torsional diffusion | 0.1570 | 0.969 | 453 | 133 | internal clash | 51.4% clash, 27.2% energy_ratio |

Result interpretation:

- LOQI is the cleanest learned generator in the ref-set analysis.
- LOQI's main failure mode is chirality rather than general strain/clash.
- Torsional diffusion has strong recovery but substantial clash and energy-ratio failures.
- NextMol and MCF have good recovery but larger severe-failure tails.
- Raw torsion has broad energy-ratio failure, consistent with high-strain perturbation.
- RDKit raw/minimized are almost entirely PB-clean.

## ChEMBL3D-PB Audit

The apparent `67` for ChEMBL3D-PB on the ref set is mean PB-passing conformers per ligand, not ligand count.

Column notes: `all ChEMBL3D mean confs` is the pre-PB conformer count per ligand. `ChEMBL3D-PB mean confs` is the post-PB conformer count per ligand. `zero PB ligands` counts ligands for which no conformer survived PB filtering.

| set | ChEMBL3D ligands | all ChEMBL3D mean confs | ChEMBL3D-PB ligands | ChEMBL3D-PB mean confs | zero PB ligands |
| --- | ---: | ---: | ---: | ---: | ---: |
| core | 94 | 115.2 | 94 | 83.3 | 0 |
| ref | 1235 | 87.5 | 1235 | 67.2 | 0 |

For ChEMBL3D, PB filtering barely changes CASF recovery:

Column notes: `PB fail rate` is the fraction of ChEMBL3D conformers removed by PoseBusters. `best RMSD` and `Hit@0.75` are computed after taking the best conformer per ligand in the corresponding set.

| set | source | mean confs | PB fail rate | best RMSD | Hit@0.75 |
| --- | --- | ---: | ---: | ---: | ---: |
| core | ChEMBL3D all | 115.2 | 0.277 | 0.524 | 0.798 |
| core | ChEMBL3D-PB | 83.3 | 0.000 | 0.524 | 0.798 |
| ref | ChEMBL3D all | 87.5 | 0.232 | 0.549 | 0.793 |
| ref | ChEMBL3D-PB | 67.2 | 0.000 | 0.550 | 0.791 |

ChEMBL3D PB failures:

Column notes: `PB input` is the number of conformers evaluated by PoseBusters. `PB failures` is the number that failed at least one PB check. `dominant failure` reports the most frequent failing check and its count.

| set | PB input | PB failures | failure rate | dominant failure |
| --- | ---: | ---: | ---: | --- |
| core | 10,828 | 2,999 | 0.277 | tetrahedral_chirality: 2,999 |
| ref | 108,073 | 25,110 | 0.232 | tetrahedral_chirality: 24,955 |

Result interpretation:

- PB filtering removes ChEMBL3D conformers, not ligands.
- ChEMBL3D-PB remains a strong compact reference after filtering.
- ChEMBL3D PB failures are overwhelmingly chirality-related.

## Qwen Core Result

Core fixed-tier excerpt:

Column notes: these are core-set results only. `PB pass`, `1.0 A clusters`, `best RMSD`, and `Hit@0.75` have the same definitions as in the ref-set ranking table.

| method | PB pass | 1.0 A clusters | best RMSD | Hit@0.75 |
| --- | ---: | ---: | ---: | ---: |
| Qwen 1.7B FSQ + bigdata pretrain | 0.934 | 75.7 | 0.301 | 0.945 |
| Qwen 0.6B FSQ + bigdata pretrain | 0.927 | 76.5 | 0.334 | 0.890 |
| LOQI | 0.957 | 20.7 | 0.341 | 0.856 |
| torsional diffusion | 0.840 | 56.6 | 0.372 | 0.889 |
| RDKit raw | 1.000 | 35.9 | 0.447 | 0.809 |
| ChEMBL3D-PB | 1.000 | 7.8 | 0.524 | 0.798 |

Result interpretation:

- Qwen 1.7B FSQ + bigdata pretrain is the strongest core-set learned generator by recovery.
- Qwen is currently core-only in the dashboard.
- Qwen stratified interpretation is limited by rotatable-bond metadata problems in at least one variant.
- LOQI remains the strongest validated ref-set generator in the available ref-set results.

## Current Result Limitations

| issue | current status |
| --- | --- |
| ligand universe | ref aggregate rows differ slightly: ChEMBL3D-PB 1235, LOQI 1232, several generated methods 1224; paired tables use common ligand IDs |
| ChEMBL3D-PB K curve | unavailable because per-conformer PB-passing RMSDs are not exported for ChEMBL3D-PB |
| Qwen ref validation | not available in the current dashboard |
| Qwen rotatable-bond metadata | at least one core Qwen variant has incorrect rotatable-bond metadata |
| chirality failures | automated PB labels; chemical interpretation depends on stereochemistry/atom-mapping audit |
| seed variance | not represented for close rankings such as LOQI versus torsional diffusion |
| compute cost | not included in the current result tables |

## Final Result Statement

LOQI is the best validated generator for building a large conformer training dataset in the available results. It improves paired CASF16 bound-pose recovery over ChEMBL3D-PB, improves over high-K RDKit raw, reaches high Hit@0.75 with fewer generated conformers than RDKit raw, and maintains the strongest validity-weighted performance among generated methods.

The current result establishes generator selection, not downstream training superiority.
