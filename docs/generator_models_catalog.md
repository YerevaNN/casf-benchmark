# Generator models catalog

Detailed write-ups for every generator and reference dataset in the dashboard. Operational generation steps: [generation_methods.md](generation_methods.md). Data loading: [casf16_chembl3d_data_loader_method.md](casf16_chembl3d_data_loader_method.md).

Catalog definitions: [`config/casf_generation_families.yaml`](../config/casf_generation_families.yaml).

---

## In-House Classical Baselines

These four families are implemented in `src/casf_benchmark/generation/conformer_sets.py`. Procedural detail: [generation_methods.md](generation_methods.md).

### RDKit random (raw and minimized)

**Method.** Conformers are embedded with RDKit ETKDGv3 from a topology-only copy of the ChEMBL3D reference graph (`useRandomCoords=True`, enforced chirality, no pruning). The raw family accumulates embeddings that pass a finite-coordinate check. The minimized family applies MMFF94s minimization (default 500 iterations) before pool accumulation.

**Training data.** None; rule-based cheminformatics.

**Compute.** CPU-only; parallel MMFF workers default to eight processes per ligand on Slurm array jobs.

**Role in the benchmark.** Matched-budget stochastic baselines on the same 1000-conformer fixed pool and tier subsampling scheme.

### Torsion perturbation (raw and minimized)

**Method.** Starting from the ChEMBL3D topology reference, all rotatable torsions (SMARTS `[!$(*#*)&!D1]-!@[!$(*#*)&!D1]`) are perturbed independently by uniform draws in ±120°. The raw family retains candidates passing a distance-geometry steric clash filter at cutoff 0.7. The minimized family adds MMFF94s minimization and a post-minimization clash filter at the same cutoff.

**Training data.** None; deterministic torsion sampling around the ChEMBL3D topology seed.

**Compute.** CPU-only, with optional parallel torsion-batch workers.

---

## External Learned Generators

Frozen checkpoints run on intersection ligands, producing up to 1000 raw conformers per molecule, then tier-normalized for dashboard ingestion. See [materialization.md](materialization.md).

### LOQI

**Reference.** Nikitin et al., *Scalable Low-Energy Molecular Conformer Generation with Quantum Mechanical Accuracy*, ChemRxiv 2025, [doi:10.26434/chemrxiv-2025-k4h7v](https://doi.org/10.26434/chemrxiv-2025-k4h7v). Code: [isayevlab/LoQI](https://github.com/isayevlab/LoQI).

**Method.** Sterochemistry-aware equivariant diffusion (Megalodon + EGNN). Generates 3D coordinates from a 2D/3D input graph; 25 DDPM steps reported sufficient. Targets low-energy conformers, not full thermodynamic ensembles.

**Training data.** **ChEMBL3D**: ~1.8M molecules, ~250M AIMNet2-optimized conformers (Omega Classic initial sampling → AIMNet2 opt).

**Compute.** No total GPU-hours in paper; repo configs use 1–2 GPUs, 800 epochs, batch 150. Inference ~0.1 s/conformer on RTX 3090 (batched).

**Primary achievement.** Joint ChEMBL3D + LoQI release; single low-energy conformers with explicit stereochemistry including macrocycles.

**Reported performance (ChEMBL3D hold-out, median ΔE kcal/mol):** LoQI **0.33** vs DMT-L 1.72, Torsional Diffusion 2.62, OMEGA 1.42. R/S accuracy 0.956, E/Z 0.992.

**CASF deployment.** `codex_dir/loqi/generations/casf16_{core,ref}_loqi_1k`.

### NExT-Mol DMT-L

**Reference.** Qiao et al., ICLR 2025, [arXiv:2502.12638](https://arxiv.org/abs/2502.12638).

**Method.** MoLlama (960M, Llama-2) + DMT diffusion head (150M for conformer prediction). Relational attention over atom/pair features with adaptive layer norm on diffusion time.

**Training data.** MoLlama: ~1.8B ZINC-15 molecules. DMT: **GEOM-DRUGS** (243k train), one conformer/molecule/epoch.

**Compute.** MoLlama: 555k steps, 4× A100-40G, ~2 weeks. DMT: 3000 epochs, batch 256.

**Primary achievement.** Billion-scale 1D LM + 3D diffusion transformer for conformer prediction.

**Reported performance (GEOM-DRUGS, δ=0.75 Å):** COV-R mean **85.8%**, med 92.3%; AMR-R mean **0.375 Å**.

**CASF deployment.** `codex_dir/nextmol_dmt_l/generations/casf16_{core,ref}_nextmol_dmt_l_1k`.

### Torsional Diffusion

**Reference.** Jing et al., NeurIPS 2022. Code: [gcorso/torsional-diffusion](https://github.com/gcorso/torsional-diffusion).

**Method.** Diffusion on rotatable-bond torsion angles; RDKit fixes bond lengths/angles. **20 denoising steps** at inference (vs 5000 for GeoDiff).

**Training data.** **GEOM-DRUGS** + GEOM-QM9 (GeoMol splits, CREST/GFN2-xTB ensembles).

**Compute.** ~**2000 GPU-hours** (NeurIPS checklist).

**Primary achievement.** First ML generator to consistently beat OMEGA on GEOM-DRUGS COV/AMR.

**Reported performance (GEOM-DRUGS test, δ=0.75 Å):** COV-R mean/med **72.7 / 80.0**; AMR-R **0.582 / 0.565**.

**CASF deployment.** `codex_dir/torsional_diffusion/generations/casf16_{core,ref}_torsional_diffusion_1k`.

### MCF drugs-L

**Reference.** Wang et al., ICML 2024, [arXiv:2311.17932](https://arxiv.org/abs/2311.17932). Code: [apple/ml-mcf](https://github.com/apple/ml-mcf).

**Method.** PerceiverIO coordinate diffusion on graph Laplacian eigenvector fields. Dashboard uses **L** checkpoint (242M params).

**Training data.** **GEOM-DRUGS** (GeoMol splits), 300 epochs, batch 128.

**Compute.** **8× A100 GPUs**.

**Primary achievement.** Scaling domain-agnostic coordinate diffusion beats torsion-specialized models on GEOM-DRUGS.

**Reported performance (GEOM-DRUGS, δ=0.75 Å):** COV-R **81.6 / 89.2**; AMR-R **0.468 / 0.438**.

**CASF deployment.** `codex_dir/mcf_drugs_l/generations/casf16_{core,ref}_mcf_drugs_l_1k`.

---

## In-House Qwen Conformer Generators

Continued pretraining from **Qwen3** base checkpoints (`0p6b` → 0.6B, `1p7b` → 1.7B, `4b` → 4B). **Enriched SMILES** serialization (or **FSQ tokenizer** when `fsq` in name). Autoregressive generation conditioned on 2D SMILES; decode via `decode_cartesian_v2`.

### Training data suffixes

| Suffix | Data |
| --- | --- |
| *(default)* | GEOM-DRUGS |
| `revisited` | GEOM-DRUGS-Revisited ([doi:10.1039/D5DD00206K](https://doi.org/10.1039/D5DD00206K)) |
| `bigdata` | Filtered **OMol** (test-set overlap removed) |
| `bigdata_to_revisited` | OMol pretrain → Revisited finetune |
| `fsq` | GEOM-DRUGS + FSQ coordinate tokenizer |
| `fsq_bigdata_pretrain` | OMol pretrain → GEOM-DRUGS + FSQ |

### Dashboard Qwen variants

| Family ID | Base | Training | Tokenization |
| --- | --- | --- | --- |
| `qwen_0p6b_bigdata` | Qwen3-0.6B | OMol (bigdata) | Enriched SMILES |
| `qwen_0p6b_bigdata_to_revisited` | Qwen3-0.6B | OMol → Revisited | Enriched SMILES |
| `qwen_0p6b_fsq` | Qwen3-0.6B | GEOM-DRUGS | FSQ |
| `qwen_0p6b_fsq_bigdata_pretrain` | Qwen3-0.6B | OMol → GEOM-DRUGS | FSQ |
| `qwen_1p7b_bigdata` | Qwen3-1.7B | OMol (bigdata) | Enriched SMILES |
| `qwen_1p7b_bigdata_to_revisited` | Qwen3-1.7B | OMol → Revisited | Enriched SMILES |
| `qwen_1p7b_fsq` | Qwen3-1.7B | GEOM-DRUGS | FSQ |
| `qwen_1p7b_fsq_bigdata_pretrain` | Qwen3-1.7B | OMol → GEOM-DRUGS | FSQ |
| `qwen_1p7b_revisited` | Qwen3-1.7B | Revisited only | Enriched SMILES |
| `qwen_4b_bigdata` | Qwen3-4B | OMol (bigdata) | Enriched SMILES |
| `qwen_4b_revisited` | Qwen3-4B | Revisited | Enriched SMILES |

**CASF inference:** intersection SMILES from mapping CSV → up to 1000 conformers. Core: `codex_dir/qwen_gens`; ref: `codex_dir/qwen/generations/casf16_ref_qwen_1k`.

---

## Reference Datasets

Not generators; geometric baselines for comparison.

### CASF-2016

**Reference.** Su et al., *J. Chem. Inf. Model.* 2019, [doi:10.1021/acs.jcim.8b00545](https://doi.org/10.1021/acs.jcim.8b00545).

**How gathered.** 285 ligands from PDBbind v.2016 refined set (4,057 complexes). Cluster at ≥90% protein sequence similarity → **57 clusters** → 5 representatives each spanning binding affinity (≥100-fold span within cluster). Excludes duplicate compounds/stereoisomers.

**Benchmark use:** `casf_crystal` (bound MOL2, primary RMSD target), `casf_opt` (re-optimized poses).

### ChEMBL3D

**Reference.** Nikitin et al., ChemRxiv 2025, [doi:10.26434/chemrxiv-2025-k4h7v](https://doi.org/10.26434/chemrxiv-2025-k4h7v). Dataset: [doi:10.1184/R1/31428449](https://doi.org/10.1184/R1/31428449).

**How gathered.** ChEMBL v34 → FixpKa → Omega Classic → Flipper stereoisomers → AIMNet2 optimization (CPCM solvation) → energy-window filtering. **~1.8M molecules, ~250M conformers.**

**Benchmark use:** `chembl3d_sdf` (topology SDF), `chembl3d_gt` (full zarr ensemble), `chembl3d_gt_pb` (PB-filtered ensemble vs CASF crystal). Loader details: [casf16_chembl3d_data_loader_method.md](casf16_chembl3d_data_loader_method.md).

### GEOM-DRUGS

**Reference.** Axelrod & Gómez-Bombarelli, *Scientific Data* 2022, [doi:10.1038/s41597-022-01288-4](https://doi.org/10.1038/s41597-022-01288-4).

**How gathered.** ~304k drug-like molecules; CREST metadynamics + GFN2-xTB energies; ~37M conformers total. Standard training corpus for Torsional Diffusion, MCF, DMT, and default Qwen path.

See also: [data_preparation.md](data_preparation.md) · [analyzer.md](analyzer.md)
