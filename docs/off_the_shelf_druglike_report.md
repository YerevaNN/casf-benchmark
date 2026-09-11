# Off-the-shelf conformer models on the druglike set

Status report for running the four cloned GitHub generators (LOQI, NExT-Mol DMT-L,
Torsional Diffusion, MCF drugs-L) on the 23-molecule druglike benchmark, with
PoseBusters validity and GEOM-style COV/MAT at threshold 0.75 Å.

Branch: `main` (pending push). Slurm job IDs and Weka paths below are on the
research cluster.

## Ground truth

Admin did provide reference structures. The evaluation uses the ensemble pickle
that aggregates them:

| Asset | Path | Notes |
| --- | --- | --- |
| Reference pickle | `/mnt/weka/vtarasov/druglike_smi.pickle` | 23 SMILES keys → `{geom_smiles, confs}`; **2450** reference conformers total |
| Raw per-molecule SDFs | `/mnt/weka/hrant/plinder_frequency/2024-06-v2/conformers/raw/` | 23 molecule directories; matches the pickle |
| Stale combined path | `/mnt/weka/hrant/plinder_frequency/2024-06-v2/conformers/combined/` | **Does not exist** on Weka; ignore |

COV/MAT and PoseBusters both read the pickle (via `eval_druglike_covmat.py` /
`eval_druglike_conformers.py`), not the raw SDF tree.

## What was run

- **Set:** 23 druglike molecules, **1000** conformers each (23000 generated per model).
- **Models:** LOQI, NExT-Mol DMT-L, Torsional Diffusion, MCF drugs-L — same upstream
  checkpoints documented in the off-the-shelf CASF manual under
  `/mnt/weka/mbedrosian/codex_dir/`.
- **Pipeline:** Slurm GPU array (`scripts/submit_off_the_shelf_druglike.sh`) → per-molecule
  `parts/*.pickle` → merge → PoseBusters + diversity/energy
  (`scripts/eval_druglike_conformers.py`) → in-repo COV/MAT
  (`scripts/eval_druglike_covmat.py`) → dashboard publish
  (`scripts/publish_druglike_dashboard.sh`).

### Slurm layout

| Stage | Script | Output root |
| --- | --- | --- |
| Inference | `scripts/run_off_the_shelf_druglike.sbatch` | `/mnt/weka/mbedrosian/codex_dir/druglike_off_the_shelf/gen_results/{model}_druglike/` |
| Finalize | `scripts/finalize_off_the_shelf_druglike.sbatch` | Merged pickle + `eval_covmat_075/` + CSVs under `.../evaluation/` |
| Logs | — | `/mnt/weka/mbedrosian/codex_dir/druglike_off_the_shelf/logs/` |

Array job **255879–255886** (inference + first eval wave). LOQI eval was re-run as
**255975** after adding `Chem.RemoveHs` normalization (explicit H in LOQI outputs
had driven PoseBusters to 0% pass).

## Code on GitHub

All new automation lives in **`YerevaNN/casf-benchmark`** (this repo):

| Path | Role |
| --- | --- |
| [`scripts/run_off_the_shelf_druglike.py`](../scripts/run_off_the_shelf_druglike.py) | Import-alias driver for the four codex adapters; array parts + merge |
| [`scripts/eval_druglike_covmat.py`](../scripts/eval_druglike_covmat.py) | In-repo COV-R/P, MAT-R/P at 0.75 Å vs the pickle GT |
| [`scripts/submit_off_the_shelf_druglike.sh`](../scripts/submit_off_the_shelf_druglike.sh) | Submit inference arrays + dependent finalize jobs |
| [`scripts/run_off_the_shelf_druglike.sbatch`](../scripts/run_off_the_shelf_druglike.sbatch) | One GPU task per molecule |
| [`scripts/finalize_off_the_shelf_druglike.sbatch`](../scripts/finalize_off_the_shelf_druglike.sbatch) | Merge, PB eval, COV/MAT |
| [`scripts/publish_druglike_dashboard.sh`](../scripts/publish_druglike_dashboard.sh) | Combine Qwen + OTS CSVs → `extended_druglike_*` tables |
| [`src/casf_benchmark/config/generation_runs.yaml`](../src/casf_benchmark/config/generation_runs.yaml) | Four new `*_druglike` catalog entries |
| [`tests/test_off_the_shelf_druglike.py`](../tests/test_off_the_shelf_druglike.py) | Merge order + COV/MAT unit checks |

Related prior work (Qwen checkpoints on the same set): [`docs/qwen_checkpoints_and_druglike_report.md`](qwen_checkpoints_and_druglike_report.md).

## Results (threshold 0.75 Å, 1000 confs/molecule)

### Headline metrics

| Model | PB pass (overall) | COV-R mean | COV-P mean | MAT-R mean | MAT-P mean | Status |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| LOQI | **90.9%** | **0.912** | **0.590** | **0.421** | **0.816** | Complete |
| NExT-Mol DMT-L | 62.2% | 0.827 | 0.370 | 0.464 | 1.172 | Complete |
| Torsional Diffusion | 83.4% | 0.786 | 0.426 | 0.455 | 1.189 | Complete |
| MCF drugs-L | — | — | — | — | — | **Inference in progress** (Slurm **255885**; ~1 h/molecule) |

LOQI leads on coverage-recall and coverage-precision among the finished off-the-shelf
runs; NExT-Mol is mid-pack on PB validity with competitive MAT-R.

### Per-model artifact paths

```
/mnt/weka/mbedrosian/codex_dir/druglike_off_the_shelf/
  gen_results/{loqi,nextmol_dmt_l,torsional_diffusion,mcf_drugs_l}_druglike/
    generation_results.pickle
    eval_covmat_075/{covmat_results.txt,rmsd_matrix.csv}
  evaluation/{label}_summary.csv
  evaluation/{label}_per_molecule.csv
  evaluation/druglike_eval.sqlite          # combined with Qwen after publish
```

## Dashboard

Extended Analysis tabs **`extended_druglike_summary`** and
**`extended_druglike_per_molecule`** in
`data/results/extended_casf_analysis.sqlite` include all **12 Qwen** druglike rows
plus the **three completed** off-the-shelf models. Re-publish after MCF finishes:

```bash
cd /home/mbedrosian/code/casf-benchmark
./scripts/publish_druglike_dashboard.sh
# then attach extended_casf_analysis.sqlite to a new dashboard-data-* release
```

Streamlit pin: set `CASF_DASHBOARD_RELEASE` to the release tag that carries the
updated sidecar (see [`src/casf_benchmark/release_data.py`](../src/casf_benchmark/release_data.py)).

## Fixes applied during the run

1. **LOQI cwd** — sampler expects `data/chembl3d_stereo/...` relative to the LOQI repo root.
2. **Merge Python** — finalize step uses each model's venv (torch required to unpickle parts).
3. **Explicit hydrogens** — `normalize_conformers()` strips H before PB and COV/MAT (LOQI).

## Short bullet summary (for Slack)

- Ran **4 cloned GitHub conformer models** on the **23-molecule druglike set** (1000 confs each) with ensemble GT from **`druglike_smi.pickle`** (2450 refs; raw SDFs under `plinder_frequency/.../conformers/raw/`).
- **LOQI, NExT-Mol, Torsional Diffusion** inference + PB + COV/MAT **done**; **MCF still generating** on Slurm.
- **Best finished COV-R:** LOQI **0.91**; best **COV-P:** LOQI **0.59**; PB leader among OTS: LOQI **91%**.
- New scripts + YAML entries are in **`casf-benchmark`**; dashboard publish via **`publish_druglike_dashboard.sh`** (15/16 rows live until MCF lands).
