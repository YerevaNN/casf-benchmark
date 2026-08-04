# CASF–ChEMBL3D Intersection Conformer Generation

## Purpose

This document describes the procedure implemented in `src/casf_benchmark/generation/conformer_sets.py`. The script generates multi-conformer ensembles for ligands that appear in the CASF–ChEMBL3D exact-intersection mapping table. For each ligand it produces twelve conformer sets spanning four generation families and three sampling tiers, applies geometric and chemical validity filters, and records generation statistics in tab-separated manifest files.

The method is designed for benchmark construction rather than production docking. CASF crystal ligand coordinates serve as the PoseBusters reference geometry, while ChEMBL3D topology structures supply the molecular graph and torsion-seed coordinates for structure-based generators. Random-coordinate and torsion-perturbation pipelines are evaluated under matched sample budgets so that downstream geometric analysis can compare methods on a common ligand panel.

For the full generator catalog (learned models, Qwen variants, and reference dataset provenance), see [casf16_generation_pipeline_and_models_method.md](casf16_generation_pipeline_and_models_method.md).

## Production output locations (Weka)

Completed runs write SDFs and manifests under `{run_root}/generation/` on Weka shared storage. Default roots:

| Cohort | Run root |
| --- | --- |
| Core | `/mnt/weka/mbedrosian/pharma_generation_analysis/core_pb_full_dynamic_chembl_count` |
| Ref | `/mnt/weka/mbedrosian/pharma_generation_analysis/ref_pb_full_dynamic_chembl_count` |

Example: `.../core_pb_full_dynamic_chembl_count/generation/rdkit_random_raw_fixed/{mol_id}.sdf`.

Full table (external generators, reference baselines, inputs): [weka_data_paths.md](weka_data_paths.md).

Run locally:

```bash
export PYTHONPATH=src
python -m casf_benchmark.generation.conformer_sets --help
```

### Intersection mapping table

Each run begins from a comma-separated mapping file produced by `scripts/match_casf16_chembl3d_exact.py`. See [casf16_chembl3d_exact_match_method.md](casf16_chembl3d_exact_match_method.md) for full data preparation (CASF/ChEMBL3D layout, mapping, and intersection ligand curation). Rows link one CASF ligand identifier to one ChEMBL3D molecule. Required fields are `ligand_id`, `source_file`, `chembl3d_group`, `chembl3d_mol_id`, and `conformer_count`. SMILES are read from `chembl3d_isomeric_smiles` when present, otherwise from `casf_heavy_isomeric_smiles`.

Default path for the core cohort: `/mnt/weka/mbedrosian/data/casf16/casf16_core_chembl3d_exact_intersection.csv`.

A parallel reference cohort uses `/mnt/weka/mbedrosian/data/casf16/casf16_ref_chembl3d_exact_intersection.csv` with ligands drawn from the CASF-2016 reference set.

### CASF ligand structures

Three-dimensional coordinates for each mapped ligand are read from MOL2 files in a ligand directory. The script resolves paths as `{ligand_dir}/{source_file}` using the `source_file` column from the mapping table. These files provide fallback coordinates when ChEMBL3D topology loading fails and serve as the PoseBusters reference molecule when a conformer-bearing reference is required.

Default path for the core cohort: `/mnt/weka/mbedrosian/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands`.

### ChEMBL3D topology structures

For each mapped ligand the script loads a ChEMBL3D topology molecule from the on-disk SDF archive at `{topology_root}/{group}.sdf`, matching on `mol_id` or `_Name`. The loaded structure is prepared for torsion operations by ensuring explicit hydrogens and at least one conformer. If the topology SDF entry is unavailable, the script attempts to load the CASF MOL2 file as a fallback reference.

Default topology root: `/mnt/weka/mbedrosian/data/chembl3d/topologies`.

## Per-Ligand Workflow

Processing proceeds independently for each row in the mapping table. The following steps apply to every ligand before any conformer generation begins.

**Mapping validation.** If `chembl3d_group` or `chembl3d_mol_id` is empty, all twelve method outputs are marked with status `missing_chembl3d_mapping` and no conformers are generated.

**Reference loading.** The torsion reference is loaded via `load_torsion_ref`, which prefers the ChEMBL3D topology SDF and falls back to the CASF MOL2 path. Failure yields status `chembl3d_ref_load_failed` for all methods.

**Base topology for embedding.** A copy of the torsion reference with all conformers removed defines the topology-only molecule (`base_mol`) used as the RDKit embedding template. Rotatable-bond count is computed with RDKit `Descriptors.NumRotatableBonds`.

**Target set sizes.** Three sampling targets are computed per ligand. The fixed tier uses `--fixed_set_size` (default 1000). The dynamic tier uses `max(1, -20 + 22 × rotatable_bonds)`. The ChEMBL-count tier uses `max(1, conformer_count)` from the mapping table. When either dynamic or ChEMBL-count targets exceed the fixed pool size, subsets are capped at the fixed pool and a diagnostic message is printed; this cap is reflected in manifest status fields `dynamic_target_capped` or `chembl_count_target_capped`.

**Rotatable-torsion gate for torsion methods.** Torsion-based pipelines require at least one rotatable torsion as detected by the SMARTS pattern `[!$(*#*)&!D1]-!@[!$(*#*)&!D1]`. Ligands without rotatable torsions still run the RDKit random pipelines but receive empty torsion-method outputs with status `no_rotatable_bonds`.

**Output path assignment.** For each ligand the script prepares twelve output paths under `{output_dir}/{method}/{mol_id}.sdf`, one per generation method.

## Generation Families

Four independent generation families are executed per ligand. Each family accumulates a fixed-size pre-PoseBusters pool, then applies tier-specific subsampling and PoseBusters filtering as described in the Sampling and Validation section.

### RDKit random raw (`rdkit_random_raw`)

Conformers are embedded with RDKit ETKDGv3 using random initial coordinates (`useRandomCoords=True`), enforced chirality, no pruning (`pruneRmsThresh=-1.0`), and a per-batch seed derived deterministically from the global seed and ligand identifier. Embedding runs in batches of `--generation_batch_size` (default 1000) until the fixed pool reaches `--fixed_set_size`.

Within each batch, conformers with non-finite coordinates are rejected. No distance-geometry clash filter is applied in this family. The only geometric gate before pool accumulation is the finite-coordinate check.

### RDKit random minimized (`rdkit_random_minimized`)

This family follows the same ETKDGv3 embedding protocol as the raw family but minimizes each accepted embed before pool accumulation. Minimization uses MMFF94s (configurable via `--ff_variant`) for up to `--max_minimize_iters` iterations (default 500). Minimization outcomes with status outside `{0, 1}` are rejected. No pre- or post-minimization DG clash filter is applied. Minimization runs in parallel across worker processes controlled by `--minimize_workers` (defaulting to `--num_threads`).

Pool filling proceeds in tranches: each iteration embeds up to `max(remaining, 2500)` candidates (capped by batch size), minimizes them, and appends survivors until the fixed target is met or embedding fails.

### Torsion raw (`torsion_raw`)

Conformers are generated by copying the ChEMBL3D torsion reference and perturbing all detected rotatable torsions. Each torsion receives an independent uniform perturbation in `[-max_torsion_delta_deg, +max_torsion_delta_deg]` (default ±120°). The perturbation fraction defaults to 1.0, meaning every torsion is perturbed on every trial.

Candidates are filtered with a distance-geometry steric clash test at cutoff 0.7 before entering the pool. Clash detection uses RDKit `GetMoleculeBoundsMatrix` with set15 bounds, scaled van der Waals radii, triangle smoothing, and macrocycle-14 configuration disabled. Direct bonds and valence-angle (1,3) atom pairs are excluded; hydrogens are ignored. A pair is rejected when interatomic distance falls below 0.7 times the lower bound from the bounds matrix.

Torsion batch generation is parallelized across worker processes when `--minimize_workers` exceeds one.

### Torsion minimized (`torsion_minimized`)

This family combines torsion perturbation, pre-minimization clash filtering, MMFF94s minimization, and post-minimization clash filtering. The workflow operates in tranches:

1. Generate torsion-perturbed candidates and retain only those passing the pre-minimization DG clash filter at cutoff 0.7.
2. The first tranche targets `--torsion_min_pre_pool_size` pre-minimize passers (default 1500); subsequent refill tranches target 500 passers.
3. Minimize survivors with MMFF94s under the same acceptance criteria as the RDKit minimized family.
4. Reject minimized structures that fail the post-minimization DG clash filter at cutoff 0.7.
5. Append post-minimize passers to the final pool until the fixed target is reached.

## Sampling Tiers

Each generation family produces three output sets distinguished by the `set_tier` field in the manifest.

**Fixed tier.** The full pre-PoseBusters pool (default 1000 conformers) is retained for validation and output. Sidecar files record indices `0 … N-1` mapping back to the fixed pool.

**Dynamic tier.** A random subset of size `min(dynamic_target, fixed_pool_size)` is drawn without replacement from the fixed pool using a deterministic seed. The subset is chosen before PoseBusters so that failed conformers remain in the sampling frame. Sidecar files record the selected fixed-pool indices.

**ChEMBL-count tier.** A random subset of size `min(conformer_count, fixed_pool_size)` is drawn under the same protocol as the dynamic tier, using the ChEMBL3D conformer count from the mapping table as the target size.

## Validation and Filtering

Validation is applied in a fixed order that differs by family. The active filter chain for each method can be reconstructed from manifest fields and the generation funnel printed after single-molecule runs.

### Finite-coordinate check (RDKit families only)

All embedded conformers are inspected for finite x, y, and z coordinates on every atom. Non-finite coordinates are counted in `finite_rejected`.

### Distance-geometry clash filter (torsion families)

Applied at cutoff 0.7 for torsion raw (pre-pool only) and torsion minimized (both pre- and post-minimization). Rejections increment `clash_rejected` during pool accumulation and `post_min_clash_rejected` after minimization.

### MMFF94s minimization (minimized families only)

Each candidate is minimized with MMFF94s. Rejections are categorized as minimization failure (status not in `{0, 1}`), minimization error (exception during force-field setup or minimization), or post-minimization clash. Counts are recorded in `minimization_failed`, `minimization_error`, and `post_min_clash_rejected`. The manifest also reports minimization input size and derived rates.

### PoseBusters conformer validity

PoseBusters is executed once on the full fixed-tier pool. Pass/fail labels for dynamic and ChEMBL-count tiers are derived by indexing into the fixed-tier results, ensuring that PoseBusters is not re-run on subsets. This design preserves failed conformers in the sampling frame while avoiding redundant validation work.

The reference molecule passed to PoseBusters is the loaded torsion reference (ChEMBL3D topology or CASF MOL2 fallback). If the reference lacks coordinates, coordinates from the first predicted conformer are copied onto the reference for identity checks.

PoseBusters is mandatory. If the package is not importable, the script raises an error rather than silently marking conformers as valid. On batch failures the script falls back to a divide-and-conquer isolation strategy that assigns failure to individual conformers without discarding the entire pool.

The inline PoseBusters configuration activates the following binary checks:

- File loading for predicted and reference structures
- RDKit sanitisation
- InChI convertibility
- All atoms connected
- No radicals
- Identity checks against the reference: molecular formula, bond connectivity, tetrahedral chirality, and double-bond stereochemistry (InChI options `w`)
- Distance-geometry geometry: bond lengths (threshold 0.25), bond angles (threshold 0.25), internal steric clash (threshold 0.3, hydrogens ignored, sanitize enabled)
- Planar aromatic rings (5- and 6-membered aromatic SMARTS, flatness threshold 0.25 Å)
- Planar double bonds (trigonal carbon–carbon double-bond SMARTS, flatness threshold 0.25 Å)
- Energy ratio (threshold 100.0, ensemble size 50, non-strict InChI)

Conformers passing all checks are written to the output SDF. Per-check failure counts and rates are serialized as JSON in the manifest.

## Output

### Directory layout

The script resolves the generation directory via `resolve_generation_dir`. If the output path already contains method subdirectories, files are written directly there; otherwise a `generation/` subdirectory is created.

Default root for the core cohort: `/mnt/weka/mbedrosian/pharma_generation_analysis/core_pb_full_dynamic_chembl_count`.

Under the generation directory, outputs are organized as follows:

- `{method}/{mol_id}.sdf` — multi-record SDF containing PoseBusters-passing conformers for that method and ligand
- `{method}/{mol_id}.indices.tsv` — sidecar mapping output conformers back to fixed-pool indices (for dynamic and chembl_count tiers)
- `manifest_parts/{mol_id}.tsv` — per-ligand manifest written during distributed runs
- `manifest.tsv` — merged manifest produced after all parts complete

Each SDF record carries molecule- and conformer-level properties including `mol_id`, `input_smiles`, `source_input`, `generation_method`, `set_tier`, `conf_id`, `num_rotatable_bonds`, `num_target_confs`, `minimization_applied`, and minimization metadata when applicable.

Empty outputs (zero conformers after filtering or early failure) are written as empty files. Manifest rows still record the failure status.

### Manifest fields

The tab-separated manifest records one row per ligand per method (twelve rows per successfully processed ligand). Key fields include:

- Identification: `mol_id`, `input_smiles`, `source_input`, `generation_method`, `set_tier`
- Targets: `num_target_confs`, `rotatable_bonds`
- Generation funnel: `generated_candidates`, `finite_rejected`, `clash_rejected`, `pre_clash_passed`, `generation_batches`
- Minimization: `minimization_applied`, `minimization_input_confs`, `minimization_failed`, `minimization_error`, `post_min_clash_rejected`, and derived rates
- PoseBusters: `pb_input_confs`, `pb_pass_confs`, `pb_fail_confs`, `pb_pass_rate`, `pb_check_fail_counts_json`, `pb_check_fail_rates_json`
- Outcome: `kept_confs`, `selected_confs`, `waste_ratio`, `status`, `walltime_seconds`, `sdf_path`

Status values include `ok`, `failed_to_fill_pool`, `embedding_failed`, `all_candidates_invalid`, `empty_dynamic_subset`, `empty_chembl_count_subset`, `dynamic_target_capped`, `chembl_count_target_capped`, `missing_chembl3d_mapping`, `chembl3d_ref_load_failed`, and `no_rotatable_bonds`.

### Console diagnostics

After processing, the script prints aggregate counts per method: kept conformers after PoseBusters, PoseBusters pass fraction, and for fixed-tier raw methods the median number of embedding trials required to fill the pool. Single-molecule runs additionally print a generation funnel table and the path to the manifest part file.

## Execution

### Repository location

All generation commands assume a checkout at:

```
/home/mbedrosian/code/casf-benchmark
```

Slurm wrappers set `REPO_ROOT` to this path and `PYTHONPATH=${REPO_ROOT}/src`. See [casf16_generation_pipeline_and_models_method.md](casf16_generation_pipeline_and_models_method.md) for the full command reference across classical, external, and Qwen generators.

### Environment and dependencies

The script requires RDKit, PoseBusters, and the `molgen3D.pharmacophore` package from this repository. Run with the repository root on `PYTHONPATH` (the Slurm wrappers set `PYTHONPATH={REPO_ROOT}/src`).

A conda environment with RDKit and PoseBusters installed is required. The codebase references the chembl3d environment (`/home/mbedrosian/.conda/envs/chembl3d/bin/python`) as a known working configuration. Parallel minimization and torsion batch workers set `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and `NUMEXPR_NUM_THREADS` to 1 inside worker processes to avoid oversubscription.

No GPU is used. Network access is not required at runtime beyond what the environment already provides.

### Local invocation

Single-ligand test run with defaults: invoke the generator with `--molecule_offset 0`, `--limit_molecules 1`, `--num_threads 8`, and `--minimize_workers 8`. Full local runs over all mapped ligands use the same thread settings without offset or limit flags; this mode is not recommended for large panels due to runtime.

### Command-line parameters

The following arguments control inputs and behavior:

- `--chembl_map_csv` — intersection mapping table (default: core CSV path above)
- `--ligand_dir` — directory containing CASF MOL2 files (default: core intersection ligand directory)
- `--chembl3d_topology_root` — ChEMBL3D topology SDF root
- `--output_dir` — conformer-set root; generation artifacts are written under `generation/`
- `--limit_molecules` — maximum ligands to process, or `all` (default)
- `--molecule_offset` — skip this many CSV rows before applying the limit (used for array scheduling)
- `--merge_manifest` — merge `manifest_parts/*.tsv` into `manifest.tsv` and exit
- `--seed` — global random seed (default 42); combined with ligand-specific keys for reproducibility
- `--num_threads` — RDKit embedding threads and PoseBusters energy-ratio threads
- `--minimize_workers` — parallel MMFF minimization and torsion-batch workers (defaults to `--num_threads`)
- `--fixed_set_size` — fixed-tier pool size before PoseBusters (default 1000)
- `--generation_batch_size` — embed/perturb trials per accumulation batch (default 1000)
- `--torsion_min_pre_pool_size` — initial pre-minimize pool for torsion minimized (default 1500)
- `--max_torsion_delta_deg` — maximum torsion perturbation magnitude in degrees (default 120)
- `--perturb_fraction` — fraction of torsions perturbed per trial (default 1.0)
- `--ff_variant` — force field for minimization, `MMFF94` or `MMFF94s` (default `MMFF94s`)
- `--max_minimize_iters` — MMFF minimization iteration limit (default 500)

The clash cutoff arguments `--pre_clash_cutoff` and `--post_clash_cutoff` exist but are overridden internally by each pipeline family; they should not be relied upon to change production behavior.

### Slurm batch execution

From the repository root (`/home/mbedrosian/code/casf-benchmark`):

```bash
cd /home/mbedrosian/code/casf-benchmark
./scripts/submit_casf.sh generate core   # or: generate ref
```

Production runs use `scripts/submit_casf.sh`, which submits one Slurm array task per ligand and a dependent merge job. The core cohort is launched with `./scripts/submit_casf.sh generate core`; the reference cohort with `./scripts/submit_casf.sh generate ref`.

The wrapper validates that the mapping CSV, ligand directory, topology root, and Python interpreter exist. It counts rows in the mapping CSV, submits array jobs in chunks of up to 1000 tasks (limited to 80 concurrent tasks by default), and schedules a merge job that runs `--merge_manifest` after all array tasks complete.

Each array task invokes the generator with `--molecule_offset {base + task_id}`, `--limit_molecules 1`, and the thread settings from environment variables `NUM_THREADS` and `MINIMIZE_WORKERS` (default 8). Tasks skip processing if the corresponding manifest part already exists, supporting restart after partial completion.

Environment variables accepted by the submission wrapper include `CHEMBL_MAP_CSV`, `CASF_LIGAND_DIR`, `OUTPUT_DIR`, `CHEMBL3D_TOPOLOGY_ROOT`, `PYTHON`, `PARTITION`, `NUM_THREADS`, and `MINIMIZE_WORKERS`.

Default output directories are `/mnt/weka/mbedrosian/pharma_generation_analysis/core_pb_full_dynamic_chembl_count` for the core cohort and `/mnt/weka/mbedrosian/pharma_generation_analysis/ref_pb_full_dynamic_chembl_count` for the reference cohort.

After generation completes, geometric analysis can be submitted with `./scripts/submit_casf.sh analyze {core|ref}`, optionally depending on the merge job identifier.

### Preconditions for an independent run

1. A current intersection mapping CSV produced by the exact-match script.
2. CASF MOL2 files for every mapped ligand in the ligand directory.
3. ChEMBL3D topology SDF files for every mapped `(group, mol_id)` pair.
4. Write access to the output directory on shared storage.
5. An RDKit- and PoseBusters-enabled Python environment with this repository available on `PYTHONPATH`.

Per-ligand runtime depends strongly on rotatable-bond count and method family; torsion minimized is the most expensive due to clash filtering and iterative minimization tranches. Distributed array execution is the intended production mode.

## Downstream use

Merged manifest files and method SDFs feed the geometric analysis pipeline (`scripts/analyze_casf_conformer_sets.py`, submitted via `scripts/submit_casf.sh analyze`). Analysis compares generated ensembles against CASF crystal coordinates, CASF optimized ligands, and ChEMBL3D ground-truth conformers under matched sample tiers.

Example analysis artifacts for the core 94-ligand intersection are documented in `docs/casf_geometric_report_core_94_ligands.md`.

## Design rationale

**ChEMBL3D topology as the embedding template.** Using the matched ChEMBL3D topology (conformers removed) as `base_mol` ensures that RDKit random embedding operates on the same molecular graph as the ChEMBL3D reference, independent of CASF MOL2 atom ordering or hydrogen treatment.

**Single PoseBusters pass on the fixed pool.** Running PoseBusters once on the full fixed-tier pool and projecting results onto dynamic and ChEMBL-count subsets avoids repeated validation cost while keeping failed conformers visible in subset statistics.

**Deterministic subsampling.** Dynamic and ChEMBL-count subsets use seeded random sampling keyed by ligand and method, so repeated runs with the same global seed yield identical index selections.

**Conservative torsion eligibility.** Torsion methods are skipped when no rotatable torsions are detected, matching the eligibility filter applied during intersection mapping and preventing undefined perturbation behavior.

**Manifest-part architecture.** Writing one manifest part per ligand enables embarrassingly parallel Slurm execution with idempotent restart: existing parts are not regenerated.

## Related documentation

- `docs/casf16_chembl3d_exact_match_method.md` — construction of the intersection mapping table consumed by this script
- `src/casf_benchmark/paths.py` — canonical default paths and directory layout constants
- `docs/generate_casf_smiles_conformer_sets_report.md` — compact pipeline summary used during development
- `docs/casf16_benchmark_analysis_method.md` — geometric analysis, dashboard aggregation, and extended analysis pipeline
