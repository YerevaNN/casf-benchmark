# External Generation Set Materialization

## Purpose

This document describes `scripts/materialize_casf_generation_sets.py` and the library implementation in `src/casf_benchmark/casf_generation_normalizer.py`. These tools convert **raw fixed-pool conformer outputs** from external learned generators (LOQI, NExT-Mol DMT-L, Torsional Diffusion, MCF, Qwen, etc.) into the same **three-tier, PoseBusters-validated layout** that the in-house RDKit generator produces directly.

RDKit/torsion baselines do not use this step; they emit `{method}_{fixed|dynamic|chembl_count}` SDFs and manifest rows from `generate_casf_smiles_conformer_sets.py`. External checkpoints typically ship one SDF per ligand under a single method name (for example `loqi_raw` or `qwen_4b_bigdata`) representing a pre-PoseBusters pool of up to 1000 conformers.

Data preparation prerequisites are in [casf16_chembl3d_exact_match_method.md](casf16_chembl3d_exact_match_method.md).

---

## Role in the pipeline

```
External inference (out of repo)
        │
        ▼
Raw fixed pool:  generation/{source_method}/{mol_id}.sdf
                 generation/manifest.tsv  (one row per source method × ligand)
        │
        ▼
materialize_casf_generation_sets.py
        │
        ▼
Tiered outputs:  generation/{source_method}_fixed/{mol_id}.sdf
                 generation/{source_method}_dynamic/{mol_id}.sdf
                 generation/{source_method}_chembl_count/{mol_id}.sdf
                 + indices sidecars for dynamic/chembl_count tiers
                 + updated manifest.tsv (three rows per source method × ligand)
        │
        ▼
analyze_casf_conformer_sets.py  (PB-once: copies pb_* from manifest)
```

---

## Input contract

### Directory layout

Materialization operates on a **generation root** (default `{root}/generation` unless method subdirectories already sit directly under `--root`).

Required before materialization:

| Artifact | Description |
| --- | --- |
| `generation/manifest.tsv` | Tab-separated manifest with **untiered** source methods only |
| `generation/{source_method}/{mol_id}.sdf` | Multi-record SDF, one file per ligand, up to `fixed_set_size` conformers |

The manifest must list source methods **without** `_fixed`, `_dynamic`, or `_chembl_count` suffixes. Rows whose `generation_method` already ends with a tier suffix or whose `set_tier` is `fixed`/`dynamic`/`chembl_count` are skipped as already materialized.

Typical source method names: `loqi_raw`, `nextmol_dmt_l_raw`, `torsional_diffusion_raw`, `mcf_drugs_l_raw`, `qwen_4b_bigdata`.

### Intersection inputs

The script reads the same mapping table and ligand directory as classical generation:

- `--chembl-map-csv` — default core intersection CSV
- `--ligand-dir` — default core intersection MOL2 directory
- `--chembl-dataset-root` — ChEMBL3D root (`topologies/` used for PoseBusters reference loading)

Ligands are loaded via `load_intersection_molecules` with optional `--molecule-offset` and `--limit-molecules` for array-style partial runs.

---

## Method

For each discovered source pool and each intersection ligand, `materialize_source_pool` performs the following.

### 1. Load raw fixed pool

The multi-record SDF at `generation/{source_method}/{mol_id}.sdf` is loaded. Empty files or parse failures yield manifest status `fixed_pool_load_failed:*` for all three tier methods.

### 2. Load PoseBusters reference

ChEMBL3D topology is loaded via `load_torsion_ref(group, mol_id, topologies/, casf_mol2_fallback)`. Failure yields `reference_load_failed`.

### 3. Compute tier targets

Same formulas as classical generation:

- **Fixed target** — `num_target_confs` from manifest row, or `--fixed-set-size` (default 1000)
- **Dynamic target** — `max(1, -20 + 22 × rotatable_bonds)`, capped by pool size
- **ChEMBL-count target** — `conformer_count` from mapping CSV, capped by pool size

Rotatable bonds come from the manifest row or are recomputed from the reference molecule.

### 4. Subsample indices

Deterministic index draws (seed `--seed`, default 1729) select dynamic and chembl_count subsets **from the fixed pre-PoseBusters pool** before validation, matching the RDKit generator’s tier logic.

### 5. PoseBusters-once validation

`finalize_pipeline_pair` from `generate_casf_smiles_conformer_sets.py` runs PoseBusters **once** on the full fixed pool, then derives pass/fail labels for dynamic and chembl_count tiers by indexing into those results. This mirrors the PB-once strategy documented in [casf16_benchmark_analysis_method.md](casf16_benchmark_analysis_method.md).

### 6. Write outputs

For source stem `loqi_raw`, output methods are:

- `loqi_raw_fixed`
- `loqi_raw_dynamic`
- `loqi_raw_chembl_count`

Each tier writes `{generation_dir}/{method}/{mol_id}.sdf` and, for dynamic/chembl_count, an indices sidecar mapping output conformers back to fixed-pool indices.

Manifest rows include the same funnel fields as classical generation (`pb_input_confs`, `pb_pass_confs`, `kept_confs`, `set_tier`, etc.) so the analyzer can copy `pb_*` without re-running PoseBusters.

---

## Output

### Primary products

- Tiered SDF trees under `generation/`
- Updated `generation/manifest.tsv` (or parts under `generation/manifest_parts_normalized/` when using `--write-manifest-part`)
- Optional validation pass via `--validate-only`

### Manifest merge

`--merge-manifest-parts` concatenates `manifest_parts_normalized/*.tsv` into `manifest.tsv` (same pattern as classical generation merge jobs).

---

## Execution

### Environment

```bash
cd /home/mbedrosian/code/casf-benchmark
export PYTHONPATH=src
PYTHON=/home/mbedrosian/.conda/envs/chembl3d/bin/python
```

Requires RDKit, PoseBusters, pandas, and `molgen3D.pharmacophore` modules. CPU-only unless PoseBusters energy-ratio checks spawn threads (`--posebusters-energy-threads`).

### Full materialization (single process)

Example for LOQI core outputs:

```bash
$PYTHON scripts/materialize_casf_generation_sets.py \
  --root /mnt/weka/mbedrosian/codex_dir/loqi/generations/casf16_core_loqi_1k \
  --chembl-map-csv /mnt/weka/mbedrosian/data/casf16/casf16_core_chembl3d_exact_intersection.csv \
  --ligand-dir /mnt/weka/mbedrosian/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands \
  --method loqi_raw \
  --fixed-set-size 1000 \
  --seed 1729 \
  --posebusters-workers 1
```

Example for Qwen root (all untiered methods in manifest):

```bash
$PYTHON scripts/materialize_casf_generation_sets.py \
  --root /mnt/weka/mbedrosian/codex_dir/qwen_gens \
  --chembl-map-csv /mnt/weka/mbedrosian/data/casf16/casf16_core_chembl3d_exact_intersection.csv
```

Restrict methods with repeated `--method` flags or `--methods-csv`.

### Array-style partial runs

For large panels, mirror classical generation array semantics:

```bash
$PYTHON scripts/materialize_casf_generation_sets.py \
  --root /mnt/weka/mbedrosian/codex_dir/qwen_gens \
  --molecule-offset 42 \
  --limit-molecules 1 \
  --write-manifest-part
```

Then merge:

```bash
$PYTHON scripts/materialize_casf_generation_sets.py \
  --root /mnt/weka/mbedrosian/codex_dir/qwen_gens \
  --merge-manifest-parts
```

### Validation

```bash
$PYTHON scripts/materialize_casf_generation_sets.py \
  --root /mnt/weka/mbedrosian/codex_dir/loqi/generations/casf16_core_loqi_1k \
  --validate-only
```

Checks: required manifest columns; every tier SDF exists; dynamic/chembl_count indices sidecars present.

### Qwen artifact recovery

If a prior run wrote incorrectly tiered Qwen directories, quarantine stale artifacts before rematerializing:

```bash
$PYTHON scripts/materialize_casf_generation_sets.py \
  --root /mnt/weka/mbedrosian/codex_dir/qwen_gens \
  --quarantine-wrong-qwen-artifacts
```

This moves tier-suffixed `qwen_*` method dirs, normalized manifest parts, and the current manifest into a timestamped quarantine folder.

Restore a manifest backup:

```bash
$PYTHON scripts/materialize_casf_generation_sets.py \
  --root /mnt/weka/mbedrosian/codex_dir/qwen_gens \
  --restore-manifest-backup manifest.tsv.bak
```

---

## Command-line reference

| Flag | Default | Purpose |
| --- | --- | --- |
| `--root` | `codex_dir/qwen_gens` | Run root; generation dir resolved as `{root}/generation` |
| `--generation-dir` | auto | Override generation directory |
| `--source-manifest` | `{generation}/manifest.tsv` | Input manifest with untiered methods |
| `--output-manifest` | same as source | Manifest write target |
| `--chembl-map-csv` | core intersection CSV | Ligand mapping |
| `--ligand-dir` | core intersection ligands | CASF MOL2 directory |
| `--chembl-dataset-root` | `/mnt/weka/mbedrosian/data/chembl3d` | Topology root parent |
| `--method` | all pools in manifest | Repeatable filter for source methods |
| `--methods-csv` | — | CSV column `generation_method` list |
| `--limit-molecules` | all | Cap ligands processed |
| `--molecule-offset` | 0 | Skip first N mapping rows |
| `--fixed-set-size` | 1000 | Fixed-tier pool size |
| `--seed` | 1729 | Tier subsampling seed |
| `--posebusters-workers` | 1 | Parallel PoseBusters workers |
| `--write-manifest-part` | off | Write one ligand part to `manifest_parts_normalized/` |
| `--merge-manifest-parts` | off | Merge parts → manifest |
| `--validate-only` | off | Validate without writing |
| `--quarantine-wrong-qwen-artifacts` | off | Move bad Qwen tier artifacts aside |

---

## Preconditions

1. Intersection CSV and ligand directory prepared ([data preparation doc](casf16_chembl3d_exact_match_method.md)).
2. External inference completed: raw SDFs + untiered manifest rows for each source method.
3. ChEMBL3D topology shards reachable for every mapped `(group, mol_id)`.
4. PoseBusters importable in the active environment.

After materialization, run geometric analysis on the generation root before rebuilding the dashboard ([casf16_benchmark_analysis_method.md](casf16_benchmark_analysis_method.md)).

---

## Related documentation

- [casf16_chembl3d_conformer_generation_method.md](casf16_chembl3d_conformer_generation_method.md) — classical generator that materialization reuses for PB-once and tier logic
- [casf16_generation_pipeline_and_models_method.md](casf16_generation_pipeline_and_models_method.md) — external model catalog and deployment paths
- [casf16_dashboard_method.md](casf16_dashboard_method.md) — viewing materialized results in the Streamlit app
