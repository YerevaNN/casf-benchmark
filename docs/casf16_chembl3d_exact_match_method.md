# CASF Benchmark Data Preparation

## Purpose

This document describes how benchmark inputs are assembled before conformer generation and analysis. It covers external dataset provenance, the ChEMBL3D on-disk layout, exact CASF–ChEMBL3D intersection mapping (`scripts/match_casf16_chembl3d_exact.py`), and curation of per-cohort ligand directories. Downstream generation is documented in [casf16_chembl3d_conformer_generation_method.md](casf16_chembl3d_conformer_generation_method.md); generator models in [casf16_generation_pipeline_and_models_method.md](casf16_generation_pipeline_and_models_method.md).

Canonical default paths live in `src/casf_benchmark/paths.py`.

---

## Overview

```
CASF-2016 MOL2 (core or ref panel)
        +
ChEMBL3D topology index + SDF shards + zarr conformers
        │
        ▼
Exact SMILES intersection mapping  →  core/ref CSV
        │
        ▼
Intersection ligand MOL2 directory  (copy or symlink)
        │
        ▼
Generation / analysis pipelines
```

Two ligand cohorts are prepared independently:

| Cohort | CASF ligand source | Intersection CSV | Intersection ligand directory |
| --- | --- | --- | --- |
| **Core** | `CASF16/ligands` (285-complex benchmark panel subset) | `casf16_core_chembl3d_exact_intersection.csv` | `CASF16/core_chembl3d_exact_intersection_ligands` |
| **Ref** | `CASF16_REF/ligands` (broader reference panel) | `casf16_ref_chembl3d_exact_intersection.csv` | `CASF16_REF/ref_chembl3d_exact_intersection_ligands` |

---

## External datasets

### CASF-2016

**Reference.** Su et al., *Comparative Assessment of Scoring Functions: The CASF-2016 Update*, *J. Chem. Inf. Model.* 2019, [doi:10.1021/acs.jcim.8b00545](https://doi.org/10.1021/acs.jcim.8b00545).

**Provenance.** CASF-2016 provides 285 protein–ligand complexes selected from the PDBbind v.2016 refined set: complexes are clustered at ≥90% protein sequence similarity and five ligands per cluster are chosen to span binding-affinity range. Ligand coordinates are distributed as Tripos MOL2 files together with the benchmark package (PDBbind-CN / CASF download).

**Local layout.**

| Path | Role |
| --- | --- |
| `/mnt/weka/mbedrosian/data/casf16/CASF16/ligands` | Core-panel MOL2 files (`*.mol2`) |
| `/mnt/weka/mbedrosian/data/casf16/CASF16/ligands_opt` | Optional re-optimized ligand poses (used as `casf_opt` in analysis) |
| `/mnt/weka/mbedrosian/data/casf16/CASF16_REF/ligands` | Reference-panel MOL2 files |

CASF structures are bound or prepared ligand geometries. They serve as crystal (or optimized) ground truth in recovery metrics; they are not modified by the preparation scripts except when copied into intersection directories.

### ChEMBL3D

**Reference.** Nikitin et al., ChemRxiv 2025 / ChEMBL3D dataset [doi:10.1184/R1/31428449](https://doi.org/10.1184/R1/31428449).

**Provenance.** Starting from ChEMBL v34 drug-like molecules: protomer states (OpenEye FixpKa), initial 3D ensembles (Omega Classic), stereoisomer enumeration (Flipper), then AIMNet2 geometry optimization under implicit solvation. The public release contains ~1.8M unique molecules and ~250M optimized conformers.

**Local layout.**

| Path | Role |
| --- | --- |
| `/mnt/weka/mbedrosian/data/chembl3d/topologies/{group}.sdf` | One multi-record SDF per three-digit shard (000–999); topology + at least one conformer for torsion seeding |
| `/mnt/weka/mbedrosian/data/chembl3d/zarr_database` | Full conformer ensembles for reference analysis (`chembl3d_gt`) |
| `/mnt/weka/mbedrosian/data/chembl3d_index/chembl3d_topology_smiles_index.csv` | SMILES index for fast intersection lookup |

The topology SMILES index is a precomputed artifact on the shared mount. There is no index-build script in this repository; regenerating it requires the upstream ChEMBL3D release workflow.

**Index CSV fields used by matching.**

- `group` — three-digit shard id
- `mol_id` — molecule identifier within the shard
- `isomeric_canonical_smiles` — match key against regenerated CASF SMILES
- `conformer_count` — stored ensemble size (feeds `chembl_count` sampling tier)

The parser accepts a seven-field layout (group, mol_id, original SMILES, isomeric SMILES, …, conformer count) or an extended nine-column layout; rows with non-three-digit groups are discarded.

---

## Exact CASF–ChEMBL3D intersection mapping

Implementation: `scripts/match_casf16_chembl3d_exact.py`.

The script identifies CASF ligands whose regenerated canonical heavy-atom isomeric SMILES exactly match a ChEMBL3D index entry, with additional structural eligibility checks. The output CSV is the authoritative map between CASF `ligand_id` and ChEMBL3D `(group, mol_id)`.

The method is intentionally conservative. A ligand is written only when (i) its MOL2 parses and normalizes, (ii) heavy isomeric SMILES matches the index, (iii) the ChEMBL3D topology SDF entry resolves on disk, and (iv) the topology has at least one rotatable torsion for downstream torsion-based generators.

### Stage 1. CASF SMILES regeneration

For each CASF MOL2 file, four canonical SMILES are computed after a fixed normalization protocol.

**Parsing.** `Chem.MolFromMol2File` with `sanitize=True`, `removeHs=False`; on failure, retry with `sanitize=False` then `Chem.SanitizeMol`. Failures are excluded as `mol2_parse_failed`.

**Hydrogen normalization.** `RemoveHs` → `AddHs(addCoords=True)` with all added hydrogen coordinates set to `(0, 0, 0)` so explicit-H SMILES reflect topology rather than bound-state H placement.

**SMILES recorded.**

1. `casf_explicit_isomeric_smiles`
2. `casf_explicit_nonisomeric_smiles`
3. `casf_heavy_isomeric_smiles` — **intersection match key**
4. `casf_heavy_nonisomeric_smiles`

### Stage 2. ChEMBL3D index lookup

Distinct CASF heavy isomeric SMILES are matched against `isomeric_canonical_smiles` in the index. Duplicate index rows for the same `(smiles, group, mol_id)` collapse to the row with largest `conformer_count`.

### Stage 3. Topology indexing

For each hit, the script loads `topologies/{group}.sdf` and resolves `mol_id` (or `_Name`). Only required shards are read.

### Stage 4. Eligibility filtering

For each matched SMILES, the first hit with a resolvable topology is selected. `prepare_torsion_ref_mol` and rotatable-torsion SMARTS `[!$(*#*)&!D1]-!@[!$(*#*)&!D1]` gate eligibility.

**Exclusion reasons:** `no_smiles_match`, `missing_topology_sdf`, `no_rotatable_bonds`.

### Mapping output

**Core default:** `/mnt/weka/mbedrosian/data/casf16/casf16_core_chembl3d_exact_intersection.csv`

**Ref default:** `/mnt/weka/mbedrosian/data/casf16/casf16_ref_chembl3d_exact_intersection.csv` (pass `--output-csv` explicitly; the script default is the core path).

Each row includes: `ligand_id`, `source_file`, four CASF SMILES variants, `chembl3d_group`, `chembl3d_mol_id`, `chembl3d_isomeric_smiles`, `conformer_count`.

### Mapping execution

```bash
cd /home/mbedrosian/code/casf-benchmark
export PYTHONPATH=src

# Core cohort (~94 ligands after filters)
/home/mbedrosian/.conda/envs/chembl3d/bin/python scripts/match_casf16_chembl3d_exact.py \
  --ligand-dir /mnt/weka/mbedrosian/data/casf16/CASF16/ligands \
  --output-csv /mnt/weka/mbedrosian/data/casf16/casf16_core_chembl3d_exact_intersection.csv

# Reference cohort (~1219 ligands)
/home/mbedrosian/.conda/envs/chembl3d/bin/python scripts/match_casf16_chembl3d_exact.py \
  --ligand-dir /mnt/weka/mbedrosian/data/casf16/CASF16_REF/ligands \
  --output-csv /mnt/weka/mbedrosian/data/casf16/casf16_ref_chembl3d_exact_intersection.csv
```

Optional overrides: `--chembl-index`, `--topology-root`.

**Dependencies.** RDKit; `casf_benchmark.chembl3d.loader` and `casf_benchmark.generation.conformer_sets` (rotatable-torsion helper). CPU only.

**Console diagnostics to retain:** `casf_ligands_processed`, `exact_matched_ligands`, `excluded_*` counts and ligand lists.

### Design rationale

- **Regenerated SMILES** ensure both sides use the same RDKit normalization.
- **Heavy isomeric key** avoids explicit-H representation mismatches while preserving stereochemistry.
- **Topology SDF verification** guarantees loaders can instantiate a 3D reference.
- **Rotatable-bond filter** excludes ligands that cannot use torsion-perturbation generators.

---

## Intersection ligand directory curation

Matching writes the CSV only; it does not copy MOL2 files. Generation and analysis expect a dedicated directory containing one MOL2 per mapped ligand so paths stay stable and panels are self-contained.

**Target directories.**

- Core: `/mnt/weka/mbedrosian/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands`
- Ref: `/mnt/weka/mbedrosian/data/casf16/CASF16_REF/ref_chembl3d_exact_intersection_ligands`

**Procedure.** For each row in the intersection CSV, copy or symlink `{source_casf_dir}/{source_file}` into the intersection directory. The `source_file` column names the original MOL2 under the cohort ligand directory (`CASF16/ligands` or `CASF16_REF/ligands`). The file stem must equal `ligand_id`.

Example (core):

```bash
CSV=/mnt/weka/mbedrosian/data/casf16/casf16_core_chembl3d_exact_intersection.csv
SRC=/mnt/weka/mbedrosian/data/casf16/CASF16/ligands
DST=/mnt/weka/mbedrosian/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands
mkdir -p "$DST"
tail -n +2 "$CSV" | while IFS=, read -r ligand_id source_file _; do
  ln -sf "$SRC/$source_file" "$DST/$source_file"
done
```

Use `cp` instead of `ln -sf` if the intersection directory must be portable off the Weka mount.

**Validation.** Count of `*.mol2` in the intersection directory should equal `exact_matched_ligands` from the mapping run. Every `ligand_id` in the CSV must resolve to an existing file.

---

## Optional exploratory step

`scripts/analyze_casf16_ligands.py` computes RDKit descriptors (rotatable bonds, heavy atoms, canonical SMILES) for all MOL2 files in a CASF ligand directory and writes `ligand_descriptors.csv` plus a summary markdown file. This script is not required for the benchmark pipeline; it supports panel characterization before matching.

---

## Preconditions checklist

Before running generation or analysis on a cohort:

1. CASF MOL2 panel present under the cohort ligand directory.
2. ChEMBL3D index CSV, topology SDF tree, and zarr archive readable on the mount.
3. Current intersection CSV produced by `match_casf16_chembl3d_exact.py`.
4. Intersection ligand directory populated from the CSV `source_file` column.
5. RDKit-enabled Python with repository `src/` on `PYTHONPATH`.

---

## Related documentation

- [casf16_chembl3d_conformer_generation_method.md](casf16_chembl3d_conformer_generation_method.md) — RDKit/torsion generation using the intersection CSV
- [casf16_materialize_generation_sets_method.md](casf16_materialize_generation_sets_method.md) — tier normalization for external ML pools
- [casf16_generation_pipeline_and_models_method.md](casf16_generation_pipeline_and_models_method.md) — generator catalog and reference dataset summaries
- [casf16_benchmark_analysis_method.md](casf16_benchmark_analysis_method.md) — analysis and dashboard build
