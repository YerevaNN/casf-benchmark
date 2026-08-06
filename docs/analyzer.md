# Geometric analyzer

Per-run evaluation of conformer sets against CASF crystal/optimized poses and ChEMBL3D references.

**Script:** [`scripts/analyze_casf_conformer_sets.py`](../scripts/analyze_casf_conformer_sets.py)  
**Metrics:** [`src/casf_benchmark/analysis/metrics.py`](../src/casf_benchmark/analysis/metrics.py)  
**Catalog:** [`src/casf_benchmark/catalog.py`](../src/casf_benchmark/catalog.py)  
**Slurm:** [`scripts/submit_casf.sh analyze`](../scripts/submit_casf.sh) · [`scripts/run_casf_ref_geometric_analysis.sbatch`](../scripts/run_casf_ref_geometric_analysis.sbatch) · [`scripts/run_pb_once_casf_analysis.sbatch`](../scripts/run_pb_once_casf_analysis.sbatch)

## PoseBusters-once policy

**Generation mode:** copies `pb_*` from the generation manifest — generated structures are **not** re-validated at analysis time.  
**Reference mode:** runs PoseBusters live against CASF crystal (except `chembl3d_gt_pb`, which filters the full ChEMBL3D PB pass mask).

## Inputs

| Input | Location / notes |
| --- | --- |
| Manifest | `{run_root}/generation/manifest.tsv` or concatenated `manifest_parts/` on first load |
| SDFs | `{generation_dir}/{method}/{mol_id}.sdf` — empty files valid (zero PB survivors) |
| Mapping CSV | `data/mapping/casf16_*_exact_intersection.csv` — chemical filters not re-applied |
| CASF crystal MOL2 | `--casf-ligand-dir`; also searches manifest `source_input`, core intersection dir, defaults |
| CASF optimized | `--casf-opt-ligand-dir` (default `CASF16/ligands_opt`); missing → omit opt metrics, no failure |
| ChEMBL3D | `chembl3d/topologies/`, `chembl3d/zarr_database` (reference mode) |

Default core generation root: `pharma_generation_analysis/core_pb_full_dynamic_chembl_count`.

## Modes

| Mode | Flags | Behavior |
| --- | --- | --- |
| Generation (default) | — or `--generation-only` | Catalog methods in manifest; geometric metrics on SDF contents; manifest `pb_*` copied |
| Reference only | `--reference-only --ligand-set {core\|ref}` | Five rows/ligand (see below) |

`--reference-only` and `--generation-only` are mutually exclusive.

### Reference sources (five rows per ligand)

| Source | Description |
| --- | --- |
| `casf_crystal` | CASF MOL2 crystal pose |
| `casf_opt` | CASF optimized ligand |
| `chembl3d_sdf` | ChEMBL3D topology SDF conformer |
| `chembl3d_gt` | Full ChEMBL3D zarr ensemble (unfiltered) |
| `chembl3d_gt_pb` | `chembl3d_gt` after PoseBusters pass mask |

Ligand `1tlp_1tlp_conf0`: ChEMBL3D ensemble subsampled to **2000** conformers (seed 42) before reference PB.

## Ligand selection

**Generation mode:** intersection CSV ligands that appear in manifest **and** have an SDF for every discovered catalog method. Exclusions: `not_in_manifest`, `missing_generation_sdf`. Methods = ordered intersection of manifest methods and generation catalog for the ligand set.

**Reference mode:** mapped ligands with MOL2 in ligand directory. Requires successful ChEMBL3D load for mapped ligands.

Duplicate `ligand_id` rows in mapping CSV → rejected. Post-run validation checks complete source coverage.

## Per-ligand metrics

Cached as pickles under `{run_root}/analysis/cache/geometric_{generation,reference}_parts/` (`--resume-parts`).

For each source and conformer list:

| Metric | Detail |
| --- | --- |
| Steric clash | DG cutoff **0.7**; bonds and valence angles excluded; H ignored |
| PoseBusters | Manifest copy (generation) or live run (reference); pass mask for `chembl3d_gt_pb` |
| Diversity | Mean torsion std; pairwise RMSD stats; greedy clusters at **0.5, 1.0, 2.0, 3.0 Å** |
| Energies | Per-conformer MMFF94s → `geometric_energy_values.csv` |
| CASF recovery | Heavy-atom best-aligned RMSD vs crystal and optimized; best/median RMSD; hits at **0.25, 0.5, 0.75, 2.0 Å** |

Non-finite RMSD → error with atom-count diagnostics. Generation workflow attaches manifest funnel fields per method.

Cohort summary rates: **average per-ligand rates first**, then aggregate across ligands.

## Outputs

```
{run_root}/analysis/
├── tables/geometric_per_ligand_long.csv      ← required for dashboard
├── tables/geometric_*_summary.csv
├── tables/geometric_generation_filter_summary.csv
├── tables/geometric_generation_pb_failure_summary.csv
├── tables/geometric_reference_*              # reference runs
├── tables/geometric_cluster_summary.csv
├── tables/geometric_energy_summary.csv
├── tables/geometric_casf_hit_summary.csv
├── tables/geometric_casf_opt_hit_summary.csv
├── tables/geometric_energy_values.csv
├── geometric_report.md                       ← per-run auto-generated tables (not a repo doc)
└── cache/
```

Bundled per-run CSVs: [`data/results/runs/{run_id}/`](../data/results/runs/). Results vs hypothesis: [casf16-core-hypothesis-assessment.md](casf16-core-hypothesis-assessment.md).

## Commands

```bash
# Single run root (generation)
casf-analyze-conformer-sets \
  --run-root /path/to/run \
  --chembl-map-csv data/mapping/casf16_core_chembl3d_exact_intersection.csv \
  --casf-ligand-dir /path/to/core_chembl3d_exact_intersection_ligands \
  --workers 48 --resume-parts --quiet-rdkit-warnings

# Reference baselines
casf-analyze-conformer-sets \
  --run-root /path/to/reference_datasets/core \
  --reference-only --ligand-set core --workers 24

# Slurm (classical cohort)
./scripts/submit_casf.sh analyze core
```

### CLI flags

| Flag | Role |
| --- | --- |
| `--output-dir` | Run root (alias: `--run-root`) |
| `--casf-ligand-dir` / `--casf-opt-ligand-dir` | Reference MOL2 dirs |
| `--chembl-map-csv` / `--chembl-dataset-root` | Mapping + ChEMBL3D paths |
| `--ligand-set` | `core` or `ref` (reference mode) |
| `--limit-molecules` / `--molecule-offset` | Partial runs |
| `--report-only` | Regenerate report from cached parts |
| `--workers` | Ligand-level parallelism |
| `--posebusters-workers` / `--posebusters-energy-threads` | PB threading |
| `--quiet-rdkit-warnings` | Suppress RDKit logs |

Slurm defaults ([`run_casf_ref_geometric_analysis.sbatch`](../scripts/run_casf_ref_geometric_analysis.sbatch)): 48 CPUs, 192 GB, 14 days. Full multi-cohort rebuild: [`run_pb_once_casf_analysis.sbatch`](../scripts/run_pb_once_casf_analysis.sbatch) (reference core/ref → all generation roots → master CSV + dashboard DB).

Worker threads limited via `OMP/MKL/OPENBLAS/NUMEXPR/RDKIT_NUM_THREADS=1`.

## Aggregation (downstream)

Per-run `geometric_per_ligand_long.csv` files merge via `casf-build-master-csv` and `casf-build-dashboard-db` using [`src/casf_benchmark/config/casf_analysis_sources.yaml`](../src/casf_benchmark/config/casf_analysis_sources.yaml). See [dashboard.md](dashboard.md).

Optional extended stats: `casf-extended-analysis` — see [extras.md](extras.md).

## Design rationale

- Manifest-trusted PB for generation avoids duplicate validation cost while analyzing the same SDF conformers.
- Separate reference runs give consistent ChEMBL3D loading and PB framing for all generator comparisons.
- `chembl3d_gt_pb` aligns validity assumptions with filtered generation sets.
- Per-ligand rate aggregation prevents large pools from dominating cohort means.

See also: [generation_methods.md](generation_methods.md) · [materialization.md](materialization.md)
