# Data preparation

Build the CASF–ChEMBL3D exact-intersection ligand panel before generation or analysis.

**Code:** [`scripts/match_casf16_chembl3d_exact.py`](../scripts/match_casf16_chembl3d_exact.py) · [`src/casf_benchmark/paths.py`](../src/casf_benchmark/paths.py) · [`src/casf_benchmark/chembl3d/loader.py`](../src/casf_benchmark/chembl3d/loader.py)

## Pipeline

```
CASF MOL2 panel  +  ChEMBL3D index / topologies / zarr
        →  match script  →  intersection CSV
        →  copy/symlink MOL2s  →  intersection ligand directory
        →  generation / analysis
```

| Cohort | CASF ligands | Output CSV | Intersection ligands |
| --- | --- | --- | --- |
| Core | [`CASF16/ligands`](../data/mapping/) | `casf16_core_chembl3d_exact_intersection.csv` | `CASF16/core_chembl3d_exact_intersection_ligands/` |
| Ref | `CASF16_REF/ligands` | `casf16_ref_chembl3d_exact_intersection.csv` | `CASF16_REF/ref_chembl3d_exact_intersection_ligands/` |

Bundled CSVs: [`data/mapping/`](../data/mapping/). Cluster paths: [extras.md](extras.md#weka-paths).

## External inputs

**CASF-2016** — Tripos MOL2 per ligand. Crystal poses feed analysis (`casf_crystal`); optional optimized poses in `ligands_opt` (`casf_opt`).

**ChEMBL3D** (on Weka):

| Path | Role |
| --- | --- |
| `chembl3d/topologies/{group}.sdf` | Topology + seed conformer (three-digit shard `000`–`999`) |
| `chembl3d/zarr_database` | Full ensembles (reference analysis only) |
| `chembl3d_index/chembl3d_topology_smiles_index.csv` | SMILES lookup index (prebuilt upstream; not regenerated here) |

**Index columns used:** `group`, `mol_id`, `isomeric_canonical_smiles` (match key), `conformer_count` (feeds `chembl_count` tier). Parser accepts 7- or 9-column layouts; rows with non-three-digit `group` are dropped.

## Matching algorithm

Conservative gate: a ligand is written only if MOL2 parses, heavy isomeric SMILES matches the index, topology SDF resolves, and the topology has ≥1 rotatable torsion.

### Stage 1 — CASF SMILES regeneration

For each `{ligand_dir}/*.mol2`:

1. **Parse:** `Chem.MolFromMol2File(sanitize=True, removeHs=False)`; on failure retry `sanitize=False` then `Chem.SanitizeMol`. Fail → `mol2_parse_failed`.
2. **H normalization:** `RemoveHs` → `AddHs(addCoords=True)` with all added H coords `(0, 0, 0)` so explicit-H SMILES reflect topology, not bound-state H placement.
3. **Four canonical SMILES:**
   - `casf_explicit_isomeric_smiles`
   - `casf_explicit_nonisomeric_smiles`
   - `casf_heavy_isomeric_smiles` — **intersection match key**
   - `casf_heavy_nonisomeric_smiles`

### Stage 2 — Index lookup

Match distinct CASF heavy isomeric SMILES to `isomeric_canonical_smiles`. Duplicate index rows for the same `(smiles, group, mol_id)` collapse to the row with largest `conformer_count`.

### Stage 3 — Topology resolution

Load `topologies/{group}.sdf`; resolve entry by `mol_id` or `_Name`. Only required shards are read.

### Stage 4 — Eligibility

For each matched SMILES, take the first hit with resolvable topology. Run `prepare_torsion_ref_mol` and rotatable-torsion SMARTS `[!$(*#*)&!D1]-!@[!$(*#*)&!D1]`. No torsions → exclude as `no_rotatable_bonds`.

**All exclusion reasons:** `mol2_parse_failed`, `no_smiles_match`, `missing_topology_sdf`, `no_rotatable_bonds`.

**Output columns:** `ligand_id`, `source_file`, four CASF SMILES variants, `chembl3d_group`, `chembl3d_mol_id`, `chembl3d_isomeric_smiles`, `conformer_count`.

## Commands

```bash
cd casf-benchmark && export PYTHONPATH=src

# Core (~94 ligands)
python scripts/match_casf16_chembl3d_exact.py \
  --ligand-dir $CASF_BENCHMARK_DATA_ROOT/data/casf16/CASF16/ligands \
  --output-csv $CASF_BENCHMARK_DATA_ROOT/data/casf16/casf16_core_chembl3d_exact_intersection.csv

# Ref (~1219 ligands) — must pass --output-csv (script default is core path)
python scripts/match_casf16_chembl3d_exact.py \
  --ligand-dir $CASF_BENCHMARK_DATA_ROOT/data/casf16/CASF16_REF/ligands \
  --output-csv $CASF_BENCHMARK_DATA_ROOT/data/casf16/casf16_ref_chembl3d_exact_intersection.csv
```

Overrides: `--chembl-index`, `--topology-root`. CPU only; requires RDKit + `casf_benchmark.generation.conformer_sets.get_rotatable_torsions`.

**Retain console diagnostics:** `casf_ligands_processed`, `exact_matched_ligands`, `excluded_*` counts and ligand lists.

## Intersection ligand directory

Matching writes the CSV only. Populate one MOL2 per row:

```bash
CSV=$CASF_BENCHMARK_DATA_ROOT/data/casf16/casf16_core_chembl3d_exact_intersection.csv
SRC=$CASF_BENCHMARK_DATA_ROOT/data/casf16/CASF16/ligands
DST=$CASF_BENCHMARK_DATA_ROOT/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands
mkdir -p "$DST"
tail -n +2 "$CSV" | while IFS=, read -r ligand_id source_file _; do
  ln -sf "$SRC/$source_file" "$DST/$source_file"
done
```

Use `cp` instead of `ln -sf` for portable copies. **Validation:** `*.mol2` count == `exact_matched_ligands`; every `ligand_id` resolves; file stem == `ligand_id`.

## Design rationale

- Regenerated SMILES enforce identical RDKit normalization on both sides.
- Heavy isomeric key avoids explicit-H mismatches while keeping stereochemistry.
- Topology SDF check guarantees loaders can instantiate a 3D reference.
- Rotatable-bond filter aligns with torsion-generator eligibility.

## Preconditions

1. CASF MOL2 panel on disk  
2. ChEMBL3D index + topologies (+ zarr for reference analysis)  
3. Intersection CSV + populated ligand directory  
4. RDKit; `PYTHONPATH=src`

Optional panel characterization: [`scripts/analyze_casf16_ligands.py`](../scripts/analyze_casf16_ligands.py) → `ligand_descriptors.csv` (not required).

See also: [generation_methods.md](generation_methods.md) · [analyzer.md](analyzer.md)
