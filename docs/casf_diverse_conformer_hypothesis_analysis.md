# CASF16 Diverse-Conformer Hypothesis Analysis

## Executive Claim

The current results support a publishable benchmark story, but the strongest defensible claim is narrower than "geometrically diverse conformers are better training data than energy-minimized conformers."

The data currently support:

> Quality-controlled conformational diversity improves recovery of CASF16-like bound conformations, especially in the flexible-ligand tail. Energy-minimized ChEMBL3D-like ensembles are very sample-efficient, but they under-cover some protein-bound geometries that high-K diverse ensembles recover.

The data do not yet prove:

> A generative model trained on diverse conformer sets will be better than one trained on energy-minimized conformers.

That training-data claim needs matched model training where molecule split, model, token budget, inference K, filtering, and conformer count are controlled.

## Data Reviewed

Dashboard database:

- `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite`

Tables used:

- `master_wide`
- `cluster`
- `generation_filter`
- `generation`
- per-run `analysis/tables/geometric_per_ligand_metrics.csv`

Runs present in the dashboard:

- baseline RDKit/torsion core
- baseline RDKit/torsion ref
- Qwen core
- LOQI core/ref
- NextMol-DMT-L core/ref
- ChEMBL3D and CASF reference baselines

Important bookkeeping caveat: older prose refers to a 1219-ligand reference intersection, while the current dashboard/per-ligand files show 1227 generated/reference-baseline molecules for most ref runs and 1235 molecules in the LOQI ref/casf-crystal per-ligand file. Resolve this before publication; otherwise reviewers will catch the mismatch.

## Main Result

The experiment is not simply a diversity ranking. It shows a three-way tradeoff:

1. Low-energy/minimized conformers are physically clean and sample-efficient.
2. Raw or model-generated diverse conformers recover more CASF-bound geometries at large K.
3. Naive diversity introduces failure modes: high-energy torsion artifacts, stereochemistry errors, internal clashes, and local geometry failures.

This is a good paper direction because conformer-generation papers often conflate these objectives:

- local chemical validity
- low-energy plausibility
- geometric ensemble diversity
- recovery of protein-bound conformations

Your data show these are separable.

## Aggregate Findings

### Core 94-Ligand Set

Best fixed-tier generated methods by CASF best RMSD:

| method | mean confs | PB fail mean | 1.0 A clusters | best RMSD | Hit@0.5 | Hit@0.75 | Hit@2.0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Qwen 4B bigdata fixed | 872 | 0.126 | 74.2 | 0.323 | 0.747 | 0.901 | 1.000 |
| Qwen 1.7B bigdata fixed | 856 | 0.141 | 72.8 | 0.333 | 0.742 | 0.914 | 1.000 |
| LOQI fixed | 957 | 0.043 | 20.7 | 0.341 | 0.822 | 0.856 | 1.000 |
| Qwen 0.6B bigdata fixed | 851 | 0.147 | 71.8 | 0.357 | 0.720 | 0.871 | 1.000 |
| NextMol-DMT-L fixed | 861 | 0.139 | 31.0 | 0.395 | 0.753 | 0.849 | 0.989 |
| RDKit raw fixed | 1000 | 0.000 | 35.9 | 0.447 | 0.628 | 0.809 | 1.000 |
| ChEMBL3D-PB | 83 | n/a | 7.8 | 0.524 | 0.638 | 0.798 | 0.979 |

The AI generators, especially Qwen bigdata fixed, are substantially better than the classical baseline on this small core set. But this is currently only a core-set observation because Qwen ref results are not in the dashboard.

Per-ligand fixed-tier deltas versus ChEMBL3D-PB on core:

| method | mean RMSD delta | wins >0.1 A | losses >0.1 A | Hit@0.75 delta | PB pass |
| --- | ---: | ---: | ---: | ---: | ---: |
| Qwen 4B bigdata fixed | -0.204 | 0.532 | 0.064 | +0.103 | 0.997 |
| Qwen 1.7B bigdata fixed | -0.192 | 0.521 | 0.053 | +0.116 | 0.993 |
| LOQI fixed | -0.184 | 0.500 | 0.074 | +0.058 | 1.000 |
| NextMol-DMT-L fixed | -0.133 | 0.447 | 0.149 | +0.052 | 0.922 |
| RDKit raw fixed | -0.077 | 0.330 | 0.202 | +0.011 | 0.957 |
| torsion raw fixed | -0.027 | 0.372 | 0.309 | -0.021 | 0.936 |
| torsion minimized fixed | +0.030 | 0.266 | 0.426 | -0.043 | 0.957 |

On the core set, Qwen bigdata and LOQI are the strongest story. Qwen has broad geometric reach and high CASF recovery; LOQI has cleaner outputs and surprisingly strong recovery with much lower cluster count.

### Reference Set

Best fixed-tier generated methods by CASF best RMSD:

| method | mean confs | PB fail mean | 1.0 A clusters | best RMSD | Hit@0.5 | Hit@0.75 | Hit@2.0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| LOQI fixed | 970 | 0.030 | 24.0 | 0.386 | 0.760 | 0.867 | 0.988 |
| NextMol-DMT-L fixed | 874 | 0.125 | 28.0 | 0.460 | 0.697 | 0.838 | 0.973 |
| RDKit raw fixed | 1000 | 0.000 | 38.4 | 0.480 | 0.592 | 0.813 | 0.998 |
| torsion raw fixed | 897 | 0.103 | 62.0 | 0.494 | 0.577 | 0.764 | 0.993 |
| RDKit minimized fixed | 999 | 0.001 | 23.7 | 0.540 | 0.549 | 0.768 | 0.997 |
| torsion minimized fixed | 999 | 0.001 | 24.1 | 0.552 | 0.529 | 0.760 | 0.993 |
| ChEMBL3D-PB | 68 | n/a | 7.4 | 0.554 | 0.593 | 0.789 | 0.975 |

Per-ligand fixed-tier deltas versus ChEMBL3D-PB on ref:

| method | mean RMSD delta | wins >0.1 A | losses >0.1 A | Hit@0.75 delta | PB pass |
| --- | ---: | ---: | ---: | ---: | ---: |
| LOQI fixed | -0.166 | 0.465 | 0.063 | +0.077 | 0.998 |
| NextMol-DMT-L fixed | -0.095 | 0.403 | 0.153 | +0.049 | 0.945 |
| RDKit raw fixed | -0.074 | 0.315 | 0.256 | +0.024 | 0.973 |
| torsion raw fixed | -0.060 | 0.377 | 0.294 | -0.024 | 0.960 |
| torsion minimized fixed | -0.001 | 0.233 | 0.331 | -0.029 | 0.973 |

This validates the core-set conclusion that learned/conformer-generative methods can beat ChEMBL3D-PB for bound-conformer recovery. LOQI is the cleanest large-set winner in the current dashboard.

## Flexible-Ligand Tail

The rotatable-bond stratification is the most publishable evidence for your original intuition.

For ref ligands with 9+ rotatable bonds, fixed-tier deltas versus ChEMBL3D-PB:

| method | n | mean RMSD delta | wins >0.1 A | losses >0.1 A | Hit@0.75 | ChEMBL3D-PB Hit@0.75 | 1.0 A clusters |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| RDKit raw fixed | 91 | -0.403 | 0.604 | 0.231 | 0.275 | 0.253 | 253.0 |
| LOQI fixed | 91 | -0.375 | 0.648 | 0.099 | 0.472 | 0.253 | 161.2 |
| torsion minimized fixed | 91 | -0.301 | 0.582 | 0.231 | 0.363 | 0.253 | 164.2 |
| torsion raw fixed | 91 | -0.271 | 0.538 | 0.341 | 0.209 | 0.253 | 334.3 |
| NextMol-DMT-L fixed | 91 | +0.006 | 0.451 | 0.363 | 0.319 | 0.253 | 135.9 |

This says:

- ChEMBL3D-PB is strong overall, but it struggles on flexible molecules.
- Diversity helps the flexible tail, but only if quality remains acceptable.
- Raw torsion has enormous diversity, but too much of it is not useful at strict thresholds.
- LOQI has the best balance on the large set: lower diversity than raw torsion, better useful recovery.

For rigid ref ligands with 0-3 rotatable bonds:

| method | n | mean RMSD delta | wins >0.1 A | losses >0.1 A | Hit@0.75 | ChEMBL3D-PB Hit@0.75 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| LOQI fixed | 456 | -0.121 | 0.368 | 0.018 | 0.982 | 0.958 |
| NextMol-DMT-L fixed | 456 | -0.087 | 0.336 | 0.055 | 0.976 | 0.958 |
| torsion raw fixed | 456 | -0.048 | 0.320 | 0.213 | 0.947 | 0.958 |
| RDKit raw fixed | 456 | -0.010 | 0.195 | 0.193 | 0.965 | 0.958 |
| torsion minimized fixed | 456 | +0.044 | 0.116 | 0.274 | 0.945 | 0.958 |

For rigid molecules, the performance gap is smaller and ChEMBL3D-PB is already near saturated. The clean narrative is not "diversity always helps"; it is "diversity matters when the conformational search problem is hard."

## Diversity Is Not Sufficient

Across generated aggregate rows, 1.0 A cluster count correlates with better CASF recovery, but confounders matter.

Aggregate correlations:

| set | metric | corr with best RMSD | corr with Hit@0.75 |
| --- | --- | ---: | ---: |
| core | mean 1.0 A clusters | -0.712 | +0.700 |
| core | mean conformers | -0.752 | +0.738 |
| core | clusters per 100 conformers | +0.746 | -0.721 |
| ref | mean 1.0 A clusters | -0.519 | +0.387 |
| ref | mean conformers | -0.769 | +0.725 |
| ref | clusters per 100 conformers | +0.861 | -0.818 |

Interpretation:

- Absolute coverage matters: more total conformers and more total clusters generally improve best-case recovery.
- Normalized cluster density can look worse because high-density rows are usually smaller-K dynamic or ChEMBL-count selections.
- This argues for conformer-count efficiency curves rather than single aggregate tables.

## Quality Failure Modes

Classical methods and AI methods fail differently.

Fixed-tier dominant PoseBusters failures:

| method family | dominant failure mode | interpretation |
| --- | --- | --- |
| RDKit raw/minimized | near-zero failures | strong chemistry prior, but not always best bound recovery |
| raw torsion | energy_ratio | perturbations explore strained conformers; diversity is partly high-energy artifact |
| LOQI | tetrahedral_chirality | otherwise clean; likely fixable with stereochemistry constraints or post-filtering |
| NextMol-DMT-L | tetrahedral_chirality plus internal_steric_clash | stronger validity problem than LOQI |
| Qwen bigdata | internal_steric_clash, tetrahedral_chirality, bond angles/lengths, energy_ratio | excellent bound recovery on core but chemical validity needs targeted repair |

This is a useful story: neural generators are not just reproducing RDKit. They find bound-like geometries, but their chemical error profile is different and must be constrained.

## Strongest Paper Storylines

### Storyline 1: Bound-Pose Recovery Requires Quality-Filtered Diversity

Core claim:

> Standard energy-minimized conformer datasets are compact and sample-efficient, but they underrepresent protein-bound conformational tails. Diverse generators recover these tails, especially for flexible ligands, provided outputs are filtered or constrained for chemical validity.

Key figures:

- Hit@0.75 and best RMSD versus rotatable-bond bin.
- 1.0 A cluster count versus best RMSD, colored by PB fail rate.
- Per-ligand delta versus ChEMBL3D-PB, sorted by rotatable bonds.

Why it is publishable:

- It avoids claiming minimized conformers are bad.
- It reframes conformer generation around bioactive-pose coverage, not only low-energy reproduction.

### Storyline 2: AI Conformer Generators Learn Bioactive Geometry But Need Chemistry Constraints

Core claim:

> Learned generators can outperform RDKit/ChEMBL3D on CASF-bound conformer recovery, but current unconstrained pipelines trade recovery for stereochemical and local-geometry failures.

Evidence:

- Qwen bigdata fixed is best on core: best RMSD 0.323-0.333 and Hit@0.75 0.901-0.914.
- LOQI fixed validates on ref: best RMSD 0.386 and Hit@0.75 0.867.
- Failure profiles are concentrated and actionable: stereochemistry, clashes, bond geometry.

Needed before publication:

- Qwen ref results.
- Same input molecule set for LOQI/NextMol/baseline.
- A post-filtered or constrained-decoding version of Qwen to show recovery can be retained while failures drop.

### Storyline 3: Diversity Has a Useful Frontier, Not a Monotonic Optimum

Core claim:

> Maximum diversity is not optimal. Useful conformer sets lie on a Pareto frontier of diversity, validity, and bound-pose recovery.

Evidence:

- Raw torsion has the largest cluster counts but worse Hit@0.75 than LOQI and sometimes ChEMBL3D-PB.
- Minimization removes failures but collapses diversity and can hurt bound-pose recovery.
- LOQI appears near the current frontier: moderate diversity, low PB failures, strong CASF recovery.

Key figure:

- Pareto plot: x = mean 1.0 A clusters, y = Hit@0.75 or best RMSD, point size = conformer count, color = PB fail rate.

### Storyline 4: ChEMBL3D Is a Strong But Biased Reference

Core claim:

> ChEMBL3D-like energy-minimized references are not weak baselines; they are strong sample-efficient references that nevertheless under-cover flexible bioactive conformers.

Evidence:

- Ref ChEMBL3D-PB: only 68 conformers/ligand, best RMSD 0.554, Hit@0.75 0.789.
- LOQI fixed improves to best RMSD 0.386, Hit@0.75 0.867, but uses about 970 conformers.
- For 9+ rotatable-bond ligands, ChEMBL3D-PB Hit@0.75 is 0.253; LOQI fixed reaches 0.472.

This is reviewer-safe because it treats ChEMBL3D as a strong baseline instead of attacking it.

## Recommended Next Experiments

### 1. Matched Training-Data Ablation

Train identical generative models where only conformer target distribution changes.

Candidate training sets:

- ChEMBL3D-PB/minimized baseline.
- RDKit minimized, ChEMBL-count matched.
- RDKit raw PB-filtered, ChEMBL-count matched.
- LOQI-generated PB-filtered, ChEMBL-count matched.
- mixed minimized + diverse.
- cluster-medoids selected from diverse ensembles.
- torsion perturbation with short relaxation and energy-window filtering.

Controls:

- same molecules
- same conformers per molecule
- same total token budget
- same split by standardized molecule identity
- same inference K
- same output filters

This is the experiment that can actually prove or refute the training-data hypothesis.

### 2. Conformer-Count Efficiency Curves

Evaluate K = 1, 2, 5, 10, 25, 50, 100, 250, 500, 1000.

Report:

- Hit@0.5/0.75/1.0/2.0 versus K
- best RMSD versus K
- PB-passing yield versus K
- stratification by rotatable-bond bin

Expected result:

- ChEMBL3D/minimized references win early.
- diverse methods win late, especially for flexible ligands.
- the crossover point is a clean quantitative result.

### 3. Cluster-Representative Utility

For each ligand:

1. Cluster at 1.0 A.
2. Select one medoid per cluster.
3. Compare CASF recovery using medoids only against random conformer subsampling with the same count.

This answers whether clusters are chemically useful modes or just noisy geometric dispersion.

### 4. Energy-Windowed Diversity

Report cluster counts and hit rates within energy windows:

- within 5 kcal/mol of the best observed conformer
- within 10 kcal/mol
- within 20 kcal/mol
- all PB-passing

This will prevent reviewers from saying "you are only rewarding unrealistic conformers."

### 5. Bound-Pose Strain Analysis

For each CASF crystal ligand:

1. Minimize the CASF bound ligand in isolation.
2. Measure displacement from the crystal pose.
3. Estimate strain energy if feasible.
4. Ask whether diverse methods recover high-strain bound poses better than minimized datasets.

This directly tests the mechanism behind the hypothesis: protein-bound conformations may live outside the low-energy isolated-ligand ensemble.

### 6. Stereochemistry Audit

Many ChEMBL3D, LOQI, NextMol, and Qwen failures are tetrahedral-chirality failures.

Audit a stratified sample:

- true stereochemical inversion
- atom-ordering/mapping issue
- unspecified stereocenter issue
- PoseBusters/reference-standardization artifact

Do not publish a strong conclusion about stereochemistry until this is done.

### 7. Qwen Ref Validation

Qwen bigdata fixed is the strongest core result, but the dashboard only has Qwen core. Run the same Qwen pipelines on the ref intersection before using Qwen as a central paper claim.

Priority:

1. Qwen 1.7B bigdata fixed
2. Qwen 4B bigdata fixed
3. Qwen 1.7B FSQ fixed
4. dynamic and ChEMBL-count tiers for sample-efficiency curves

## Proposed Paper Framing

Suggested title direction:

> Quality-Filtered Geometric Diversity Improves Bioactive Conformer Coverage Beyond Energy-Minimized Ensembles

Suggested abstract-level narrative:

1. Conformer-generation datasets are often built from low-energy minimized ensembles.
2. Protein-bound ligand conformations are not always represented by compact low-energy ensembles, particularly for flexible ligands.
3. On CASF16 intersections with ChEMBL3D, diverse generation methods improve bound-conformer recovery over ChEMBL3D-PB and minimized baselines.
4. The gain is strongest in flexible molecules and at high K.
5. Unfiltered diversity is not sufficient: raw torsion and neural generators expose different physical-validity failure modes.
6. The useful target is quality-filtered diversity, ideally selected by cluster medoids under energy and stereochemical constraints.

## Bottom Line

The current data give a strong benchmark paper and a strong motivation for a training-data paper. The benchmark paper can be written now around CASF-bound conformer coverage and the diversity/quality frontier. The training-data paper requires matched model training to show that these improved conformer sets actually transfer into better generative models.
