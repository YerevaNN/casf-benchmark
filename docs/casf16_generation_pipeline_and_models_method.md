# CASF Conformer Generation Pipeline and Generator Models

## Purpose

This document describes the conformer generators and reference datasets that appear in the CASF benchmark dashboard. It complements the operational method notes for [exact CASF–ChEMBL3D matching](casf16_chembl3d_exact_match_method.md), [RDKit/torsion generation](casf16_chembl3d_conformer_generation_method.md), and [downstream analysis](casf16_benchmark_analysis_method.md).

The dashboard compares methods on ligands drawn from the CASF–ChEMBL3D exact-intersection panel. Each generator produces multi-conformer ensembles that are subsampled into three tiers (`fixed`, `dynamic`, `chembl_count`), filtered with PoseBusters, and evaluated against CASF crystal poses and ChEMBL3D reference conformers. Reference-only rows (`casf_crystal`, `casf_opt`, `chembl3d_sdf`, `chembl3d_gt`, `chembl3d_gt_pb`) are analyzed alongside generated sets but are not produced by the generation pipeline itself.

Catalog definitions live in `config/casf_generation_families.yaml` and `config/casf_analysis_sources.yaml`.

## Repository and execution environment

### Clone location

The benchmark code lives in the **casf-benchmark** standalone repository:

```
/home/mbedrosian/code/casf-benchmark
```

Bundled mapping CSVs and precomputed dashboard results ship under `data/`. Slurm batch scripts default to this path via `REPO_ROOT`. If you clone elsewhere, export `REPO_ROOT` before calling `scripts/submit_casf.sh`.

### Bundled results (clone-and-run dashboard)

| Role | Repo path |
| --- | --- |
| Intersection CSVs | `data/mapping/casf16_*_exact_intersection.csv` |
| Dashboard SQLite | `data/results/casf_analysis_dashboard.sqlite` |
| Master / long CSVs | `data/results/casf_analysis_master.csv`, `casf_per_ligand_long.csv` |
| Per-run analysis | `data/results/runs/{run_id}/analysis/tables/` |

See [adding_external_generators.md](adding_external_generators.md) for ingesting new model outputs.

## External generator output contract

Learned models (Qwen, LOQI, DMT-L, etc.) run **outside** this repo. To add outputs to the dashboard, write the following layout before calling `scripts/materialize_casf_generation_sets.py`:

```
{root}/generation/manifest.tsv
{root}/generation/{source_method}/{mol_id}.sdf
```

- **`manifest.tsv`**: tab-separated; one row per untiered `(generation_method, mol_id)`; method names must not include `_fixed`, `_dynamic`, or `_chembl_count` yet.
- **`{mol_id}.sdf`**: multi-record SDF with up to 1000 conformers per ligand, same topology as the ChEMBL3D reference.

After materialization the analyzer expects tiered methods `{source_method}_fixed`, `_dynamic`, `_chembl_count`. Full step-by-step: [adding_external_generators.md](adding_external_generators.md).

### Shared cluster data (optional regeneration)

For full pipeline reruns on the analysis cluster, external data lives on Weka (not in the git clone). **Pre-generated SDF pools, manifests, and MOL2 inputs** are documented in [weka_data_paths.md](weka_data_paths.md).

| Role | Path |
| --- | --- |
| CASF + intersection CSVs | `/mnt/weka/mbedrosian/data/casf16/` |
| ChEMBL3D topologies / zarr | `/mnt/weka/mbedrosian/data/chembl3d/` |
| RDKit/torsion generation output (SDFs) | `/mnt/weka/mbedrosian/pharma_generation_analysis/{core,ref}_pb_full_dynamic_chembl_count/` |
| External ML raw + materialized pools (SDFs) | `/mnt/weka/mbedrosian/codex_dir/` |
| Dashboard SQLite | `/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite` |

Complete data preparation (mapping CSV, intersection ligand directories) is in [casf16_chembl3d_exact_match_method.md](casf16_chembl3d_exact_match_method.md).

---

## Running conformer generation

Generation is split by generator type. All paths below assume the repository and environment above.

### 1. RDKit random and torsion baselines (in-repo)

**Production (Slurm, one array task per ligand).** From the repo root:

```bash
cd /home/mbedrosian/code/casf-benchmark

# Core cohort (~94 ligands) → 12 methods × 3 tiers per ligand
./scripts/submit_casf.sh generate core

# Reference cohort (~1219 ligands)
./scripts/submit_casf.sh generate ref
```

The wrapper submits Slurm array jobs (`scripts/run_casf_ref_conformer_molecule.sbatch`) plus a dependent merge job (`run_casf_ref_conformer_merge.sbatch`) that runs `--merge_manifest`. Default outputs:

- Core: `/mnt/weka/mbedrosian/pharma_generation_analysis/core_pb_full_dynamic_chembl_count`
- Ref: `/mnt/weka/mbedrosian/pharma_generation_analysis/ref_pb_full_dynamic_chembl_count`

Optional environment overrides: `CHEMBL_MAP_CSV`, `CASF_LIGAND_DIR`, `OUTPUT_DIR`, `CHEMBL3D_TOPOLOGY_ROOT`, `NUM_THREADS` (default 8), `MINIMIZE_WORKERS`, `PARTITION` (default `research`).

After the merge job finishes, note `merge_job_id` and optionally start analysis:

```bash
MERGE_JOB_ID=<id> ./scripts/submit_casf.sh analyze core
```

Full procedural detail: [casf16_chembl3d_conformer_generation_method.md](casf16_chembl3d_conformer_generation_method.md).

**Single-ligand smoke test (local, no Slurm):**

```bash
cd /home/mbedrosian/code/casf-benchmark
export PYTHONPATH=src
$PYTHON src/casf_benchmark/generation/conformer_sets.py \
  --chembl_map_csv /mnt/weka/mbedrosian/data/casf16/casf16_core_chembl3d_exact_intersection.csv \
  --ligand_dir /mnt/weka/mbedrosian/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands \
  --chembl3d_topology_root /mnt/weka/mbedrosian/data/chembl3d/topologies \
  --output_dir /mnt/weka/mbedrosian/pharma_generation_analysis/core_pb_full_dynamic_chembl_count \
  --molecule_offset 0 \
  --limit_molecules 1 \
  --num_threads 8 \
  --minimize_workers 8
```

**Local parallel launcher (alternative to Slurm on a multi-core node):**

```bash
cd /home/mbedrosian/code/casf-benchmark
export PYTHONPATH=src
$PYTHON scripts/launch_casf_conformer_sets_parallel.py \
  --chembl_map_csv /mnt/weka/mbedrosian/data/casf16/casf16_core_chembl3d_exact_intersection.csv \
  --ligand_dir /mnt/weka/mbedrosian/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands \
  --output_dir /mnt/weka/mbedrosian/pharma_generation_analysis/core_pb_full_dynamic_chembl_count \
  --workers 8 \
  --num_threads 1 \
  --minimize_workers 4
```

This processes pending ligands in parallel and merges `manifest_parts/` when done.

### 2. External learned generators (LOQI, DMT-L, Torsional Diffusion, MCF)

Inference is **not** implemented in this repository. Checkpoints are run separately; raw 1000-conformer pools land under `/mnt/weka/mbedrosian/codex_dir/` (see deployment paths in each model section below).

After raw SDFs and an untiered `generation/manifest.tsv` exist, **materialize** tiered outputs with the in-repo script:

```bash
cd /home/mbedrosian/code/casf-benchmark
export PYTHONPATH=src

# Example: LOQI core
$PYTHON scripts/materialize_casf_generation_sets.py \
  --root /mnt/weka/mbedrosian/codex_dir/loqi/generations/casf16_core_loqi_1k \
  --chembl-map-csv /mnt/weka/mbedrosian/data/casf16/casf16_core_chembl3d_exact_intersection.csv \
  --ligand-dir /mnt/weka/mbedrosian/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands \
  --method loqi_raw

# Example: LOQI ref (use ref CSV and ligand dir)
$PYTHON scripts/materialize_casf_generation_sets.py \
  --root /mnt/weka/mbedrosian/codex_dir/loqi/generations/casf16_ref_loqi_1k \
  --chembl-map-csv /mnt/weka/mbedrosian/data/casf16/casf16_ref_chembl3d_exact_intersection.csv \
  --ligand-dir /mnt/weka/mbedrosian/data/casf16/CASF16_REF/ref_chembl3d_exact_intersection_ligands \
  --method loqi_raw
```

Repeat with `--root` and `--method` for:

| Model | `--method` | Core `--root` | Ref `--root` |
| --- | --- | --- | --- |
| LOQI | `loqi_raw` | `codex_dir/loqi/generations/casf16_core_loqi_1k` | `.../casf16_ref_loqi_1k` |
| NExT-Mol DMT-L | `nextmol_dmt_l_raw` | `codex_dir/nextmol_dmt_l/generations/casf16_core_nextmol_dmt_l_1k` | `.../casf16_ref_nextmol_dmt_l_1k` |
| Torsional Diffusion | `torsional_diffusion_raw` | `codex_dir/torsional_diffusion/generations/casf16_core_torsional_diffusion_1k` | `.../casf16_ref_torsional_diffusion_1k` |
| MCF drugs-L | `mcf_drugs_l_raw` | `codex_dir/mcf_drugs_l/generations/casf16_core_mcf_drugs_l_1k` | `.../casf16_ref_mcf_drugs_l_1k` |

Full materialization reference: [casf16_materialize_generation_sets_method.md](casf16_materialize_generation_sets_method.md).

### 3. Qwen variants (in-house)

Qwen CASF inference is documented separately when added.

| Cohort | Raw pool root (Weka) |
| --- | --- |
| Core | `/mnt/weka/mbedrosian/codex_dir/qwen_gens` |
| Ref | `/mnt/weka/mbedrosian/codex_dir/qwen/generations/casf16_ref_qwen_1k` |

**Recommended cluster workflow:** [cluster_quickstart.md](cluster_quickstart.md) (`ingest_external_generation.sh` + `rebuild_dashboard_weka.sh`).

Core materialization example:

```bash
cd casf-benchmark
export PYTHONPATH=src
$PYTHON scripts/materialize_casf_generation_sets.py \
  --root /mnt/weka/mbedrosian/codex_dir/qwen_gens \
  --chembl-map-csv /mnt/weka/mbedrosian/data/casf16/casf16_core_chembl3d_exact_intersection.csv \
  --ligand-dir /mnt/weka/mbedrosian/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands
```

Ref: same command with ref paths (see [weka_data_paths.md](weka_data_paths.md)). Restrict to one checkpoint with `--method qwen_4b_bigdata` (repeat per variant as needed).

### 4. After generation

Run geometric analysis and rebuild the dashboard ([casf16_dashboard_method.md](casf16_dashboard_method.md)):

```bash
cd /home/mbedrosian/code/casf-benchmark
export PYTHONPATH=src

# Full multi-cohort rebuild (Slurm)
sbatch scripts/run_pb_once_casf_analysis.sbatch

# Or dashboard DB only (if per-run analysis/ tables already exist)
$PYTHON scripts/build_casf_analysis_master_csv.py
$PYTHON scripts/build_casf_analysis_dashboard_db.py

# View results
streamlit run apps/dashboard/streamlit_app.py
```

---

## Benchmark Context

### Ligand cohorts

Two intersection cohorts are used.

**Core cohort.** Ligands from the CASF-2016 core set that pass exact heavy-atom isomeric SMILES matching to ChEMBL3D and that contain at least one rotatable bond in the ChEMBL3D topology. Default mapping table: `/mnt/weka/mbedrosian/data/casf16/casf16_core_chembl3d_exact_intersection.csv`. Default ligand directory: `/mnt/weka/mbedrosian/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands`.

**Reference (ref) cohort.** The same matching protocol applied to the broader CASF-2016 reference ligand panel. Default mapping table: `/mnt/weka/mbedrosian/data/casf16/casf16_ref_chembl3d_exact_intersection.csv`. Default ligand directory: `/mnt/weka/mbedrosian/data/casf16/CASF16_REF/ref_chembl3d_exact_intersection_ligands`.

Both cohorts require a resolvable ChEMBL3D topology SDF entry for torsion-seeding and PoseBusters reference loading. CASF MOL2 coordinates remain the geometric ground truth for bound-pose recovery metrics.

### Sampling tiers

All dashboard generators expose three tiers per family:

- **Fixed** — a pool of up to 1000 conformers (default `fixed_set_size`).
- **Dynamic** — a random subset of size `max(1, −20 + 22 × rotatable_bonds)`, capped by the fixed pool.
- **ChEMBL-count** — a random subset of size equal to the ChEMBL3D conformer count for that molecule, capped by the fixed pool.

For RDKit/torsion baselines, tier subsampling and PoseBusters-once validation are applied inside `src/casf_benchmark/generation/conformer_sets.py`. For external and Qwen generators, raw 1000-conformer pools are normalized into the three tiers by `scripts/materialize_casf_generation_sets.py`; see [casf16_materialize_generation_sets_method.md](casf16_materialize_generation_sets_method.md).

### Dashboard generator families

| Family | Dashboard label | Core | Ref |
| --- | --- | --- | --- |
| `rdkit_random_raw` | RDKit random (raw) | yes | yes |
| `rdkit_random_minimized` | RDKit random (minimized) | yes | yes |
| `torsion_raw` | Torsion perturb (raw) | yes | yes |
| `torsion_minimized` | Torsion perturb (minimized) | yes | yes |
| `loqi_raw` | LOQI | yes | yes |
| `nextmol_dmt_l_raw` | NextMol DMT-L | yes | yes |
| `torsional_diffusion_raw` | Torsional Diffusion | yes | yes |
| `mcf_drugs_l_raw` | MCF drugs-L | yes | yes |
| `qwen_*` (10 variants) | Qwen (size/data/tokenizer suffix) | yes | yes (ref root configured; ingest when inference completes) |

Qwen **core** results are in the bundled dashboard. Qwen **ref** uses the same checkpoint names under `/mnt/weka/mbedrosian/codex_dir/qwen/generations/casf16_ref_qwen_1k` — see [cluster_quickstart.md](cluster_quickstart.md).

---

## In-House Classical Baselines

These four families are implemented in `src/casf_benchmark/generation/conformer_sets.py`. Full procedural detail is in [casf16_chembl3d_conformer_generation_method.md](casf16_chembl3d_conformer_generation_method.md).

### RDKit random (raw and minimized)

**Method.** Conformers are embedded with RDKit ETKDGv3 from a topology-only copy of the ChEMBL3D reference graph (`useRandomCoords=True`, enforced chirality, no pruning). The raw family accumulates embeddings that pass a finite-coordinate check. The minimized family applies MMFF94s minimization (default 500 iterations) before pool accumulation.

**Training data.** None; rule-based cheminformatics.

**Compute.** CPU-only; parallel MMFF workers default to eight processes per ligand on Slurm array jobs.

**Role in the benchmark.** These families provide matched-budget stochastic baselines against which learned generators are compared on the same 1000-conformer fixed pool and tier subsampling scheme.

### Torsion perturbation (raw and minimized)

**Method.** Starting from the ChEMBL3D topology reference, all rotatable torsions (SMARTS `[!$(*#*)&!D1]-!@[!$(*#*)&!D1]`) are perturbed independently by uniform draws in ±120°. The raw family retains candidates passing a distance-geometry steric clash filter at cutoff 0.7. The minimized family adds MMFF94s minimization and a post-minimization clash filter at the same cutoff.

**Training data.** None; deterministic torsion sampling around the ChEMBL3D topology seed.

**Compute.** CPU-only, with optional parallel torsion-batch workers.

**Design rationale.** Torsion perturbation tests whether a fixed ChEMBL3D seed plus local exploration can cover bioactive conformations without full coordinate re-embedding. It shares the ChEMBL3D reference frame with external structure-based generators while remaining fully classical.

---

## External Learned Generators

The following models were trained by their respective authors on public conformer datasets. For the CASF benchmark they are run as frozen checkpoints on intersection ligands, producing up to 1000 raw conformers per molecule that are then tier-normalized for dashboard ingestion.

### LOQI

**Reference.** Nikitin et al., *Scalable Low-Energy Molecular Conformer Generation with Quantum Mechanical Accuracy*, ChemRxiv 2025, [doi:10.26434/chemrxiv-2025-k4h7v](https://doi.org/10.26434/chemrxiv-2025-k4h7v). Code: [isayevlab/LoQI](https://github.com/isayevlab/LoQI).

**Method.** LoQI (Low-energy QM Informed conformer generative model) is a stereochemistry-aware equivariant diffusion model built on the Megalodon architecture: a Diffusion Transformer over invariant atom/bond features coupled to an EGNN coordinate update. It generates 3D coordinates in a single forward pass from a 2D/3D input graph. The authors report that 25 DDPM steps suffice for peak conformer quality. The model is trained to emit low-energy conformers rather than full thermodynamic ensembles.

**Training data.** End-to-end training on **ChEMBL3D** (see Reference Datasets below): ~1.8M unique drug-like molecules and ~250M AIMNet2-optimized conformers. Training conformers for individual molecules were initially sampled with OpenEye Omega Classic, then geometry-optimized with AIMNet2 under implicit solvation.

**Compute.** The paper does not report total GPU-hours. The public repository documents CUDA training with configurable GPU count (example configs use one or two GPUs; default training recipe lists 800 epochs and batch size 150). Inference on a typical ChEMBL molecule is ~0.1 s per conformer on an RTX 3090 when batched.

**Primary achievement.** Joint release of ChEMBL3D and a generative model that targets **single low-energy conformers** with explicit stereochemistry handling, including macrocycles, without fragment templates or predefined rotatable-bond lists.

**Reported benchmark performance (ChEMBL3D hold-out).** On the authors' low-energy conformer benchmark (median relative energy ΔE, topology preservation, post-generation relaxation energy):

| Method | Median ΔE (kcal/mol) | % within 0.1 kcal/mol of dataset min | Median ΔE_relax (kcal/mol) |
| --- | --- | --- | --- |
| RDKit ETKDG | 3.39 | 10.4 | 64.50 |
| OMEGA Classic | 1.42 | 15.7 | 13.19 |
| Torsional Diffusion | 2.62 | 15.3 | 63.09 |
| NExT-Mol DMT-L | 1.72 | 20.8 | 86.43 |
| **LoQI** | **0.33** | **38.5** | **3.48** |

LoQI also reports strong stereochemistry fidelity (R/S accuracy 0.956, E/Z accuracy 0.992 on held-out stereoisomers) and competitive RMSD on the Platinum Diverse protein-bound ligand set (median 0.45 Å at N = 250 vs 0.52 Å for RDKit).

**CASF deployment.** Raw outputs under `/mnt/weka/mbedrosian/codex_dir/loqi/generations/casf16_{core,ref}_loqi_1k`.

### NExT-Mol DMT-L

**Reference.** Qiao et al., *NExT-Mol: 3D Diffusion Meets 1D Language Modeling for 3D Molecule Generation*, ICLR 2025, [arXiv:2502.12638](https://arxiv.org/abs/2502.12638).

**Method.** NExT-Mol couples a pretrained 1D molecular language model (MoLlama, 960M parameters, Llama-2 decoder architecture) with a 3D diffusion head called Diffusion Molecular Transformer (DMT). DMT uses relational multi-head self-attention over atom and pair features with adaptive layer norm conditioning on diffusion time. For conformer **prediction**, a three-stage recipe is used: (1) train standalone DMT on GEOM-DRUGS, (2) warm up a cross-modal projector and LoRA adapters on MoLlama while DMT is frozen, (3) fine-tune the integrated model. **DMT-L** is the large 150M-parameter diffusion backbone used for conformer prediction (DMT-B, 55M, is used for de novo generation tasks).

**Training data.**

- MoLlama pretraining: ~1.8B molecules from ZINC-15 (90B SELFIES tokens after filtering; MW ≤ 500, logP ≤ 5).
- DMT conformer prediction: **GEOM-DRUGS** train split (243,473 / 30,433 / 1,000 molecules), one randomly sampled conformer per molecule per epoch, lowest-energy conformer selected for training following Huang et al. 2024.
- Transfer learning: MoLlama representations fused into pretrained DMT-B/L.

**Compute.**

- MoLlama pretraining: 555k global steps on 4× NVIDIA A100-40G, ~two weeks.
- DMT-B alone: ~52 s/epoch on one A100; DMT-B + MoLlama ~210 s/epoch.
- GEOM-DRUGS DMT training: 3000 epochs, batch size 256 (both DMT-B and DMT-L).

**Primary achievement.** First foundation-style pipeline to combine billion-scale 1D LM pretraining with a scalable 3D diffusion transformer, improving both de novo 3D generation distributional metrics and GEOM-DRUGS conformer-prediction COV/AMR over prior diffusion and MCF baselines.

**Reported benchmark performance (GEOM-DRUGS conformer prediction, δ = 0.75 Å, 2K generated vs K reference conformers).**

| Model | Params | COV-R mean | COV-R med | AMR-R mean | AMR-P mean |
| --- | --- | --- | --- | --- | --- |
| GeoDiff | 1.6M | 42.1 | 37.8 | 0.835 | 1.136 |
| Torsional Diffusion | 1.6M | 72.7 | 80.0 | 0.582 | 0.778 |
| TD w/ particle guidance | 1.6M | 77.0 | 82.6 | 0.543 | 0.656 |
| MCF-L | 242M | 81.6 | 89.2 | 0.468 | 0.705 |
| **DMT-L (NExT-Mol)** | **150M** | **85.8** | **92.3** | **0.375** | **0.598** |

On GEOM-QM9, DMT-L reports COV-R mean 96.0% and AMR-R mean 0.138 Å.

**CASF deployment.** Raw outputs under `/mnt/weka/mbedrosian/codex_dir/nextmol_dmt_l/generations/casf16_{core,ref}_nextmol_dmt_l_1k`.

### Torsional Diffusion

**Reference.** Jing et al., *Torsional Diffusion for Molecular Conformer Generation*, NeurIPS 2022, [proceedings link](https://proceedings.neurips.cc/paper_files/paper/2022/hash/994545b2308bbbbc97e3e687ea9e464f-Abstract-Conference.html). Code: [gcorso/torsional-diffusion](https://github.com/gcorso/torsional-diffusion).

**Method.** Diffusion on the hypertorus of rotatable-bond torsion angles rather than full 3D coordinates. RDKit provides initial local 3D structures (fixed bond lengths and angles); the model diffuses and denoises torsion angles with an SE(3)-aware score network (extrinsic-to-intrinsic formulation). Exact likelihoods enable Boltzmann-generator variants. Default GEOM-DRUGS inference uses **20 denoising steps** (vs 5000 for Euclidean GeoDiff).

**Training data.** **GEOM-DRUGS** and GEOM-QM9 using the GeoMol train/val/test splits (304k drug-like molecules, CREST/GFN2-xTB conformer ensembles). Featurization caches rotatable-bond graphs and local RDKit embeddings.

**Compute.** ~**2000 GPU-hours** on an internal cluster (reported in the NeurIPS paper checklist). CPU featurization ~2 h on first training run.

**Primary achievement.** First ML conformer generator to **consistently outperform OpenEye OMEGA** on GEOM-DRUGS ensemble RMSD metrics, while requiring two orders of magnitude fewer diffusion steps than coordinate-space GeoDiff.

**Reported benchmark performance (GEOM-DRUGS test set, δ = 0.75 Å, Table 1 mean/median).**

| Method | COV-R mean / med | AMR-R mean / med | COV-P mean / med | AMR-P mean / med |
| --- | --- | --- | --- | --- |
| RDKit ETKDG | 38.4 / 28.6 | 1.058 / 1.002 | 40.9 / 30.8 | 0.995 / 0.895 |
| OMEGA | 53.4 / 54.6 | 0.841 / 0.762 | 40.5 / 33.3 | 0.946 / 0.854 |
| GeoMol | 44.6 / 41.4 | 0.875 / 0.834 | 43.0 / 36.4 | 0.928 / 0.841 |
| GeoDiff | 42.1 / 37.8 | 0.835 / 0.809 | 24.9 / 14.5 | 1.136 / 1.090 |
| **Torsional Diffusion (20 steps)** | **72.7 / 80.0** | **0.582 / 0.565** | **55.2 / 56.9** | **0.778 / 0.729** |

At 20 steps, median runtime ~4.9 core-seconds per conformer vs 305 for GeoDiff (5000 steps) on CPU evaluation.

**CASF deployment.** Raw outputs under `/mnt/weka/mbedrosian/codex_dir/torsional_diffusion/generations/casf16_{core,ref}_torsional_diffusion_1k`.

### MCF drugs-L (Molecular Conformer Fields, Large)

**Reference.** Wang et al., *Swallowing the Bitter Pill: Simplified Scalable Conformer Generation* (MCF), ICML 2024, [arXiv:2311.17932](https://arxiv.org/abs/2311.17932). Code: [apple/ml-mcf](https://github.com/apple/ml-mcf).

**Method.** MCF treats each conformer as a **field** mapping graph Laplacian eigenvector positions to 3D coordinates. A PerceiverIO-based diffusion model denoises context/query atom-coordinate pairs without explicit torsion-angle or SE(3)-equivariant inductive biases. GEOM-DRUGS checkpoints are released in three sizes: **S** (13M), **B** (64M), and **L** (242M parameters). The dashboard uses the **L** checkpoint (`mcf_drugs_l`).

**Training data.** **GEOM-DRUGS** (and GEOM-QM9 for smaller models) with GeoMol splits. Training uses 300 epochs on DRUGS (Table 5 hyperparameters: batch size 128, model dim 1024, 6 Perceiver blocks).

**Compute.** **8× NVIDIA A100 GPUs**, 500 epochs reported in the appendix compute section (main DRUGS table lists 300 epochs; repository configs should be treated as authoritative for the released checkpoint).

**Primary achievement.** Shows that **scaling a domain-agnostic coordinate diffusion model** beats torsion-specialized architectures on GEOM-DRUGS COV/AMR without hand-built geometric assumptions; MCF-L improves precision coverage by ~20% over MCF-S.

**Reported benchmark performance (GEOM-DRUGS, δ = 0.75 Å, MCF-L / 1000 DDPM steps vs Torsional Diffusion / 20 steps).**

| Method | COV-R | COV-P | AMR-R | AMR-P |
| --- | --- | --- | --- | --- |
| Torsional Diffusion | 72.7 / 80.0 | 55.2 / 56.9 | 0.582 / 0.565 | 0.778 / 0.729 |
| **MCF-L (DDPM 1000 steps)** | **81.6 / 89.2** | **61.6 / 62.5** | **0.468 / 0.438** | **0.705 / 0.650** |

MCF also reports superior median absolute errors vs ground-truth ensemble energies and dipole moments relative to OMEGA and GeoDiff on GEOM-DRUGS.

**CASF deployment.** Raw outputs under `/mnt/weka/mbedrosian/codex_dir/mcf_drugs_l/generations/casf16_{core,ref}_mcf_drugs_l_1k`.

---

## In-House Qwen Conformer Generators

All Qwen variants are trained in-house by continued pretraining (chemical domain adaptation) from public **Qwen3** base checkpoints. Model size is encoded in the family name: `0p6b` → Qwen3-0.6B, `1p7b` → Qwen3-1.7B, `4b` → Qwen3-4B. Training uses the TorchTitan stack documented in [pretraining_runbook.md](pretraining_runbook.md); molecular instances are serialized in the enriched SMILES representation ([enriched_smiles.md](enriched_smiles.md)), where canonical isomeric SMILES topology is copied verbatim and 3D coordinates are appended as `<x,y,z>` blocks per atom.

### Training data variants

| Name suffix | Data source | Description |
| --- | --- | --- |
| *(no suffix beyond size)* | GEOM-DRUGS | Default conformer pretraining corpus: CREST-sampled, GFN2-xTB-annotated drug-like ensembles from the GEOM dataset (Axelrod & Gómez-Bombarelli, *Scientific Data* 2022, [doi:10.1038/s41597-022-01288-4](https://doi.org/10.1038/s41597-022-01288-4)). |
| `revisited` | GEOM-DRUGS-Revisited | Refinement of the GEOM-DRUGS splits and preprocessing from Nikitin et al., *GEOM-drugs revisited: toward more chemically accurate benchmarks for 3D molecule generation*, *Digital Discovery* 2025, [doi:10.1039/D5DD00206K](https://doi.org/10.1039/D5DD00206K). Excludes molecules fractured by GFN2-xTB replay and applies corrected valency/stability handling. |
| `bigdata` | Filtered OMol | Large-scale **OMol** conformer corpus with filters to remove overlap with held-out test sets and other exclusion rules (exact filter manifest is run-specific). Provides broader chemical coverage than GEOM-DRUGS alone. |
| `bigdata_to_revisited` | OMol → revisited | Two-stage schedule: pretrain on filtered OMol (`bigdata`), then continue on the GEOM-DRUGS-Revisited split. |
| `fsq` | GEOM-DRUGS (+ FSQ tokenizer) | Same GEOM-DRUGS data but uses the in-house **FSQ (finite scalar quantization) tokenizer** for coordinate-aware tokenization instead of plain enriched-SMILES decimal coordinates. |
| `fsq_bigdata_pretrain` | OMol → GEOM-DRUGS (+ FSQ) | FSQ tokenizer with OMol pretraining followed by GEOM-DRUGS fine-tuning. |

### Dashboard Qwen catalog

| Family ID | Display label |
| --- | --- |
| `qwen_0p6b_bigdata` | Qwen 0.6B bigdata |
| `qwen_0p6b_bigdata_to_revisited` | Qwen 0.6B bigdata→revisited |
| `qwen_0p6b_fsq` | Qwen 0.6B fsq |
| `qwen_0p6b_fsq_bigdata_pretrain` | Qwen 0.6B fsq+bigdata-pretrain |
| `qwen_1p7b_bigdata` | Qwen 1.7B bigdata |
| `qwen_1p7b_bigdata_to_revisited` | Qwen 1.7B bigdata→revisited |
| `qwen_1p7b_fsq` | Qwen 1.7B fsq |
| `qwen_1p7b_fsq_bigdata_pretrain` | Qwen 1.7B fsq+bigdata-pretrain |
| `qwen_1p7b_revisited` | Qwen 1.7B revisited |
| `qwen_4b_bigdata` | Qwen 4B bigdata |
| `qwen_4b_revisited` | Qwen 4B revisited |

Same `method_prefix` values are used for core and ref; cohort is determined by run root and `--ligand-set`.

**Method (inference).** Qwen models autoregressively generate enriched SMILES strings conditioned on a 2D SMILES prompt; coordinates are decoded back to RDKit molecules with `decode_cartesian_v2`. For CASF, models receive intersection ligand SMILES (from the mapping table) and produce conformer pools that are tier-normalized like other external generators.

**Compute.** Training runs on the project's H100/A100 TorchTitan fleet; exact GPU-hours vary by size and data variant. Generation for CASF is CPU/GPU mixed depending on checkpoint serving setup.

**Role in the benchmark.** Qwen variants test whether LM-style 3D coordinate prediction trained on GEOM-scale (or OMol-scale) data can recover **CASF crystal poses**. Core results are bundled in the git clone; ref results are ingested from Weka when available.

**CASF deployment.** Core: `/mnt/weka/mbedrosian/codex_dir/qwen_gens`. Ref: `/mnt/weka/mbedrosian/codex_dir/qwen/generations/casf16_ref_qwen_1k`. Tier materialization via `scripts/materialize_casf_generation_sets.py` or `scripts/ingest_external_generation.sh`.

---

## Reference Datasets

Reference rows in the dashboard are not generators. They supply geometric baselines and ground-truth ensembles for comparison.

### CASF-2016

**Reference.** Su et al., *Comparative Assessment of Scoring Functions: The CASF-2016 Update*, *Journal of Chemical Information and Modeling* 2019, [doi:10.1021/acs.jcim.8b00545](https://doi.org/10.1021/acs.jcim.8b00545).

**How the dataset was gathered.** CASF-2016 is a scoring-function benchmark decoupled from docking. The ligand test set is carved from the **PDBbind v.2016 refined set** (4,057 complexes). Protein–ligand complexes are clustered at ≥90% sequence similarity; from each of **57 clusters**, five representatives are chosen spanning the binding-affinity range (lowest, highest, and three evenly spaced intermediates, with ≥100-fold affinity span within a cluster and ≥1-fold between consecutive picks). Ligands are inspected to exclude identical compounds and stereoisomers. The result is **285 high-quality crystal structures** with reliable binding constants, distributed as MOL2/PDB files via PDBbind-CN.

**Use in this benchmark.**

- `casf_crystal` — experimental bound ligand coordinates from the intersection ligand directory (primary RMSD reference).
- `casf_opt` — optionally re-optimized ligand poses when present under `ligands_opt`.

CASF provides **protein-context-bound** geometries. It is complementary to vacuum conformer benchmarks such as GEOM-DRUGS COV/AMR.

### ChEMBL3D

**Reference.** Nikitin et al., *Scalable Low-Energy Molecular Conformer Generation with Quantum Mechanical Accuracy* (ChEMBL3D + LoQI), ChemRxiv 2025, [doi:10.26434/chemrxiv-2025-k4h7v](https://doi.org/10.26434/chemrxiv-2025-k4h7v). Dataset DOI: [10.1184/R1/31428449](https://doi.org/10.1184/R1/31428449).

**How the dataset was gathered.**

1. Start from **ChEMBL v34** drug-like structures.
2. Expand protonation/tautomer states with **OpenEye FixpKa**.
3. Generate initial 3D ensembles with **OpenEye Omega Classic**.
4. Enumerate stereoisomers with **OpenEye Flipper**.
5. Geometry-optimize every conformer with **AIMNet2** (MLIP trained to reproduce DFT energies with CPCM implicit solvent).
6. Filter broken topologies and duplicates; retain conformers within defined energy windows above the per-molecule minimum.

The public release contains **~1.8M unique molecules** and **~250M optimized conformers** (~180M within 6 kcal/mol of the minimum; subsets at 2.5 and 1 kcal/mol are also quoted). This is roughly two orders of magnitude larger than GEOM-DRUGS and emphasizes **low-energy single-molecule conformers** rather than full CREST thermodynamic ensembles.

**Local layout for the benchmark.**

- Topology SDF shards: `/mnt/weka/mbedrosian/data/chembl3d/topologies/{group}.sdf`
- Conformer zarr archive: `/mnt/weka/mbedrosian/data/chembl3d/zarr_database`
- SMILES index: `/mnt/weka/mbedrosian/data/chembl3d_index/chembl3d_topology_smiles_index.csv`

**Use in this benchmark.**

- `chembl3d_sdf` — topology reference loaded from the SDF shard (single entry per mapped molecule).
- `chembl3d_gt` — full stored conformer ensemble from zarr for that molecule (unfiltered).
- `chembl3d_gt_pb` — subset of `chembl3d_gt` passing PoseBusters against the CASF crystal reference (fair comparison to PoseBusters-filtered generated sets).

ChEMBL3D conformers anchor **torsion-seeded** classical pipelines and provide a low-energy vacuum ensemble baseline. They are not guaranteed to reproduce protein-bound CASF poses.

### GEOM-DRUGS (context for training corpora)

**Reference.** Axelrod & Gómez-Bombarelli, *GEOM, energy-annotated molecular conformations for property prediction and molecular generation*, *Scientific Data* 2022, [doi:10.1038/s41597-022-01288-4](https://doi.org/10.1038/s41597-022-01288-4).

**How the dataset was gathered.** ~304k drug-like molecules from medicinal-chemistry and MoleculeNet sources. Conformers were generated with the **CREST** metadynamics sampler using **GFN2-xTB** energies; selected subsets received higher-accuracy DFT labels. The release includes ~37M conformers across QM9, DRUGS, and related subsets with energy annotations and statistical weights.

GEOM-DRUGS is the standard training benchmark for Torsional Diffusion, MCF, NExT-Mol DMT, and the default Qwen pretraining path. **GEOM-DRUGS-Revisited** (see Qwen `revisited` suffix) corrects evaluation splits and stability metrics but remains rooted in the same underlying GEOM structures.

---

## End-to-End Pipeline Summary

```
CASF MOL2 ligands + ChEMBL3D index
        │
        ▼
Exact SMILES intersection mapping  ──►  core / ref CSV maps
        │
        ├─► RDKit random + torsion families (in-repo generator, 12 methods × 3 tiers)
        │
        ├─► LOQI / DMT-L / Torsional Diffusion / MCF-L  (external checkpoints, 1k pools)
        │
        └─► Qwen variants  (in-house checkpoints, core only)
                │
                ▼
        Tier normalization + PoseBusters  (materialize script for ML outputs)
                │
                ▼
        Geometric + CASF recovery analysis  ──►  dashboard SQLite
```

**Related documents.**

- [casf16_chembl3d_exact_match_method.md](casf16_chembl3d_exact_match_method.md) — data preparation and intersection mapping
- [casf16_materialize_generation_sets_method.md](casf16_materialize_generation_sets_method.md) — external ML tier materialization
- [casf16_chembl3d_conformer_generation_method.md](casf16_chembl3d_conformer_generation_method.md) — RDKit/torsion generation and validation details
- [casf16_benchmark_analysis_method.md](casf16_benchmark_analysis_method.md) — analysis, dashboard build, and extended statistics
- [casf16_dashboard_method.md](casf16_dashboard_method.md) — Streamlit UI and rebuild commands
- [pretraining_runbook.md](pretraining_runbook.md) — Qwen TorchTitan training configuration
- [enriched_smiles.md](enriched_smiles.md) — molecular serialization for Qwen models
