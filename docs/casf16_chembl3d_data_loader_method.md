# ChEMBL3D Data Loader

## Purpose

This document describes how the CASF benchmark reads ChEMBL3D structures from on-disk storage. Implementation lives in `src/casf_benchmark/chembl3d/loader.py`. The loader is used at three distinct stages of the pipeline:

1. **Intersection mapping** — verify that a ChEMBL3D topology SDF entry exists for each CASF ligand candidate (`scripts/match_casf16_chembl3d_exact.py`).
2. **Conformer generation** — load a torsion-seed reference molecule for RDKit embedding and torsion perturbation (`load_torsion_ref`).
3. **Geometric analysis** — load topology SDF entries and full zarr conformer ensembles for reference baselines (`chembl3d_sdf`, `chembl3d_gt`, `chembl3d_gt_pb`).

Dataset provenance and how ChEMBL3D was assembled are in [generator_models_catalog.md](generator_models_catalog.md#chembl3d). Intersection mapping that produces the `(group, mol_id)` keys consumed here is in [data_preparation.md](data_preparation.md).

---

## On-disk layout

ChEMBL3D is stored as a sharded topology SDF tree plus a parallel zarr coordinate archive. Both are keyed by a three-digit **group** shard (000–999) and a string **mol_id** within that shard.

| Path | Format | Contents |
| --- | --- | --- |
| `{topology_root}/{group}.sdf` | Multi-record SDF | One topology molecule per `mol_id`; includes at least one 3D conformer used as the torsion-seed template |
| `{zarr_root}/{group}/mol_id` | zarr array (bytes) | Encoded molecule identifiers; rows index into coordinate arrays |
| `{zarr_root}/{group}/coord` | zarr array (float32, N×A×3) | Per-conformer atom coordinates aligned to the topology SDF atom order |
| `{zarr_root}/{group}/numbers` | zarr array (int, N×A) | Per-conformer atomic numbers; validated against topology on load |
| `{index_csv}` | CSV | SMILES lookup table: `group`, `mol_id`, `isomeric_canonical_smiles`, `conformer_count`, … |

Default mount paths (override with `CASF_BENCHMARK_DATA_ROOT` or CLI flags):

```
/mnt/weka/mbedrosian/data/chembl3d/topologies/{000..999}.sdf
/mnt/weka/mbedrosian/data/chembl3d/zarr_database/{000..999}/mol_id|coord|numbers
/mnt/weka/mbedrosian/data/chembl3d_index/chembl3d_topology_smiles_index.csv
```

The intersection mapping CSV (`casf16_*_exact_intersection.csv`) stores the `(chembl3d_group, chembl3d_mol_id, conformer_count)` triple that downstream code uses to locate each molecule.

---

## Dependencies and environment

The loader requires:

- **RDKit** — SDF parsing, conformer construction, mol2 fallback
- **NumPy** — coordinate array handling
- **zarr** — reading the conformer archive

Call `require_dependencies()` before zarr access; it raises a clear error if NumPy or zarr is missing. The project documents the `chembl3d` conda environment as a known working configuration:

```bash
/home/mbedrosian/.conda/envs/chembl3d/bin/python
```

No GPU is required. All I/O is local filesystem reads.

---

## API reference (step-by-step behavior)

### `load_topology_mol(group, mol_id, topology_root)`

**Step 1.** Resolve the SDF path as `{topology_root}/{int(group):03d}.sdf`. Return `None` if the file does not exist.

**Step 2.** Iterate records in the SDF supplier (`removeHs=False`, `sanitize=False`).

**Step 3.** Match on `mol_id` property or `_Name` field equal to the requested `mol_id`.

**Step 4.** Sanitize the matched molecule. Raise `ValueError` on sanitization failure (topology corruption).

**Step 5.** Return a deep copy of the matched `Chem.Mol` with conformers intact.

This function returns a single topology entry — typically the lowest-energy or first stored conformer in the SDF shard. It does **not** read the zarr archive.

### `prepare_torsion_ref_mol(mol)`

**Step 1.** Return `None` if the input is `None` or has zero conformers.

**Step 2.** Deep-copy the molecule.

**Step 3.** If no explicit hydrogens are present, call `Chem.AddHs(prepared, addCoords=True)`.

**Step 4.** Return the prepared molecule.

Explicit hydrogens are required for MMFF minimization and torsion SMARTS detection in the generation pipeline.

### `load_torsion_ref_from_chembl3d(group, mol_id, topology_root)`

Calls `load_topology_mol` then `prepare_torsion_ref_mol`. This is the preferred torsion-seed source for generation.

### `load_torsion_ref_from_mol2(mol2_path)`

Fallback when no ChEMBL3D topology entry resolves:

**Step 1.** Parse MOL2 with `sanitize=True`; retry with `sanitize=False` + manual sanitize on failure.

**Step 2.** Require at least one conformer.

**Step 3.** Run through `prepare_torsion_ref_mol`.

Used when a mapped ligand's ChEMBL3D SDF entry is missing but the CASF crystal MOL2 has usable 3D coordinates.

### `load_torsion_ref(group, mol_id, topology_root, mol2_path=None)`

**Step 1.** Try `load_torsion_ref_from_chembl3d`.

**Step 2.** If that returns `None` and `mol2_path` is provided, try `load_torsion_ref_from_mol2`.

**Step 3.** Return `(mol, source_tag)` where `source_tag` is one of:
- `"chembl3d_topology_sdf"` — ChEMBL3D topology used
- `"casf_mol2_fallback"` — CASF MOL2 used
- `"unavailable"` — no reference loaded

The generation script (`conformer_sets.py`) treats `"unavailable"` as a hard failure for all twelve method outputs on that ligand.

### `find_mol_id_indices(mol_id_array, mol_id)`

**Step 1.** Open the zarr `mol_id` array for the group (read-only).

**Step 2.** Encode the requested `mol_id` as UTF-8 bytes.

**Step 3.** Return all row indices where the stored bytes match exactly (handles duplicate rows if present).

### `load_chembl3d_conformers(group, mol_id, topology_root, zarr_root, limit=None, row_indices=None)`

This is the main ensemble loader used in reference-mode analysis.

**Step 1. Load topology template.**

Call `load_topology_mol`. Raise `FileNotFoundError` if the topology SDF entry is missing — the zarr coordinates cannot be interpreted without a matching atom graph.

**Step 2. Open zarr arrays.**

For group `{group:03d}`, open:
- `mol_id` — byte-encoded identifiers
- `coord` — shape `(N, n_atoms, 3)`
- `numbers` — shape `(N, n_atoms)` atomic numbers

Raise `FileNotFoundError` listing any missing array paths.

**Step 3. Resolve row indices.**

If `row_indices` is not provided:
- Call `find_mol_id_indices` on the `mol_id` array.
- Return an empty list if no rows match.
- Apply `limit` (first N rows) if specified.

If `row_indices` is provided explicitly, use those rows directly.

**Step 4. Validate and materialize conformers.**

For each row index:
1. Read the `numbers` row and compare to the topology's atomic-number sequence. Raise `ValueError` on mismatch (data integrity guard).
2. Read the `coord` row (length must equal topology atom count).
3. Clone the topology template, remove all conformers, attach a new conformer with the zarr coordinates.
4. Set properties: `_Name`, `chembl3d_group`, `chembl3d_mol_id`.

**Step 5.** Return the list of RDKit molecules (one per zarr row).

Each returned molecule shares the same bond topology and atom ordering as the SDF template; only coordinates differ.

---

## Usage in the benchmark pipeline

### During intersection mapping

`match_casf16_chembl3d_exact.py` calls `prepare_torsion_ref_mol` on the matched topology to verify rotatable bonds exist. Ligands with zero rotatable torsions are excluded with reason `no_rotatable_bonds`.

### During RDKit/torsion generation

For each intersection CSV row:

```
load_torsion_ref(group, mol_id, topology_root, casf_mol2_path)
    → torsion reference (ChEMBL3D preferred)
    → base_mol = copy with conformers removed (ETKDG embedding template)
    → PoseBusters reference molecule
```

The ChEMBL3D topology defines the molecular graph for embedding. CASF MOL2 coordinates are **not** used as the embedding template unless ChEMBL3D loading fails entirely.

### During external pool materialization

`materialize_casf_generation_sets.py` uses the same `load_torsion_ref` call to obtain the PoseBusters reference for tier validation of LOQI, DMT, MCF, Torsional Diffusion, and Qwen outputs.

### During reference-mode analysis

For each mapped ligand, `analyze_reference_ligand` in `scripts/analyze_casf_conformer_sets.py`:

| Reference source | Loader call | Result |
| --- | --- | --- |
| `chembl3d_sdf` | `load_topology_mol` | Single topology SDF entry |
| `chembl3d_gt` | `load_chembl3d_conformers` (full ensemble) | All zarr rows for `(group, mol_id)` |
| `chembl3d_gt_pb` | Same load, then PoseBusters filter | Subset passing validity vs CASF crystal |

**Special case:** ligand `1tlp_1tlp_conf0` is deterministically subsampled to 2000 conformers before PoseBusters and metrics (seed 42). This cap prevents excessive runtime on molecules with very large ChEMBL3D ensembles.

---

## Data flow diagram

```
chembl3d_topology_smiles_index.csv
        │  (SMILES lookup during mapping)
        ▼
intersection CSV  ──►  (group, mol_id, conformer_count)
        │
        ├─► topologies/{group}.sdf  ──► load_topology_mol / load_torsion_ref
        │                                      │
        │                                      ▼
        │                              generation seed + PB reference
        │
        └─► zarr_database/{group}/     ──► load_chembl3d_conformers
                 mol_id / coord / numbers          │
                                                   ▼
                                         chembl3d_gt / chembl3d_gt_pb
                                         (reference-mode analysis)
```

---

## Error handling and validation

| Condition | Behavior |
| --- | --- |
| Missing SDF shard | `load_topology_mol` returns `None` |
| Missing zarr group directory | `load_chembl3d_conformers` raises `FileNotFoundError` |
| Atomic number mismatch | `ValueError` with row index and observed vs expected numbers |
| Coordinate length ≠ atom count | `ValueError` with lengths |
| No zarr rows for mol_id | Returns empty list (reference analysis raises `RuntimeError`) |
| RDKit not installed | `RuntimeError` on any loader call requiring RDKit |

The loader does not cache open zarr arrays or SDF suppliers across calls. Each invocation opens files fresh. For large batch jobs (reference analysis on 1200+ ligands), this is acceptable because per-ligand work dominates I/O.

---

## Relationship to ChEMBL3D dataset construction

ChEMBL3D was built by Nikitin et al. (ChemRxiv 2025, [doi:10.26434/chemrxiv-2025-k4h7v](https://doi.org/10.26434/chemrxiv-2025-k4h7v)) from ChEMBL v34 drug-like structures:

1. Protomer enumeration (OpenEye FixpKa)
2. Initial 3D sampling (OpenEye Omega Classic)
3. Stereoisomer expansion (OpenEye Flipper)
4. AIMNet2 geometry optimization under implicit solvation
5. Filtering of broken topologies and high-energy outliers

The public release contains ~1.8M unique molecules and ~250M optimized conformers. The on-disk layout in this benchmark (SDF shards + zarr coordinates + SMILES index) is a project-specific materialization of that release for fast random access during Slurm array jobs.

The SMILES index CSV is a precomputed artifact on the shared mount. There is no index-build script in this repository; regenerating it requires the upstream ChEMBL3D release workflow.

---

## Preconditions checklist

Before calling the loader in production:

1. ChEMBL3D topology SDF tree mounted at `topology_root`.
2. Zarr archive mounted at `zarr_root` with matching group shards.
3. Intersection CSV with valid `(chembl3d_group, chembl3d_mol_id)` for each ligand.
4. Python environment with RDKit, NumPy, and zarr.
5. For generation: CASF MOL2 paths available as fallback when topology SDF entries are missing.

---

## Related documentation

- [data_preparation.md](data_preparation.md) — produces the mapping CSV and validates topology availability
- [generation_methods.md](generation_methods.md) — uses `load_torsion_ref` during generation
- [analyzer.md](analyzer.md) — uses ensemble loading in reference mode
- [generator_models_catalog.md](generator_models_catalog.md) — ChEMBL3D dataset provenance and LOQI training context
- [extras.md](extras.md#weka-paths) — production mount paths
