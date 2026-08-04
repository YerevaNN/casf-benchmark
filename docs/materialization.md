# Materialization

Convert external raw conformer pools into the three-tier, PoseBusters-validated layout expected by the analyzer.

**Script:** [`scripts/materialize_casf_generation_sets.py`](../scripts/materialize_casf_generation_sets.py)  
**Library:** [`src/casf_benchmark/generation/normalizer.py`](../src/casf_benchmark/generation/normalizer.py)  
**Helper:** [`scripts/ingest_external_generation.sh`](../scripts/ingest_external_generation.sh) (materialize + analyze wrapper)

RDKit/torsion baselines skip this step — they write tiered outputs directly from [`conformer_sets.py`](../src/casf_benchmark/generation/conformer_sets.py).

## Flow

```
External inference  →  raw generation/{source_method}/{mol_id}.sdf + manifest.tsv
        →  materialize  →  {source_method}_{fixed|dynamic|chembl_count}/
        →  analyze_casf_conformer_sets.py
```

## Input contract

See [generation_methods.md#external-learned-models](generation_methods.md#external-learned-models).

Manifest must list **untiered** methods only (no `_fixed/_dynamic/_chembl_count`).

## What materialization does

For each source pool and ligand ([`normalizer.py`](../src/casf_benchmark/generation/normalizer.py)):

1. Load raw fixed pool SDF  
2. Load PB reference (ChEMBL3D topology + CASF fallback)  
3. Compute tier targets (same formulas as classical generator)  
4. Deterministic subsample indices (`--seed`, default 1729)  
5. PoseBusters-once via `finalize_pipeline_pair` from conformer_sets  
6. Write tier SDFs + updated manifest rows  

## Commands

```bash
export PYTHONPATH=src

python scripts/materialize_casf_generation_sets.py \
  --root /path/to/run \
  --chembl-map-csv data/mapping/casf16_core_chembl3d_exact_intersection.csv \
  --ligand-dir /path/to/ligands \
  --method loqi_raw

# All methods in manifest
python scripts/materialize_casf_generation_sets.py --root /path/to/qwen_gens

# Array-style + merge
python scripts/materialize_casf_generation_sets.py --root ... --molecule-offset 0 --limit-molecules 1 --write-manifest-part
python scripts/materialize_casf_generation_sets.py --root ... --merge-manifest-parts

# Validate layout only
python scripts/materialize_casf_generation_sets.py --root ... --validate-only
```

Useful flags: `--fixed-set-size`, `--posebusters-workers`, `--quarantine-wrong-qwen-artifacts`, `--restore-manifest-backup`.

## Next steps

1. [`analyze_casf_conformer_sets.py`](../scripts/analyze_casf_conformer_sets.py) on `--root`  
2. Add run to [`config/casf_analysis_sources.yaml`](../config/casf_analysis_sources.yaml)  
3. Rebuild dashboard — [dashboard.md](dashboard.md)

See also: [generation_methods.md](generation_methods.md) · [extras.md](extras.md#cluster-ingest)
