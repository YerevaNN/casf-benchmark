#!/bin/bash -l
#SBATCH --job-name=druglike-ots-eval
#SBATCH --partition=research
#SBATCH --cpus-per-task=80
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/weka/mbedrosian/codex_dir/druglike_off_the_shelf/logs/%x_%j.out
#SBATCH --error=/mnt/weka/mbedrosian/codex_dir/druglike_off_the_shelf/logs/%x_%j.err

set -euo pipefail

REPO_ROOT="${CASF_REPO_ROOT:-/home/mbedrosian/code/casf-benchmark}"
MODEL="${MODEL:?MODEL must be set}"
LABEL="${LABEL:?LABEL must be set}"
DRUGLIKE_PICKLE="${DRUGLIKE_PICKLE:-/mnt/weka/vtarasov/druglike_smi.pickle}"
OUTPUT_DIR="${OUTPUT_DIR:?OUTPUT_DIR must be set}"
EVAL_DIR="${EVAL_DIR:?EVAL_DIR must be set}"
PYTHON="${PYTHON:-${REPO_ROOT}/.venv/bin/python}"

case "${MODEL}" in
  loqi)
    MERGE_PYTHON=/mnt/weka/mbedrosian/codex_dir/loqi/.venv/bin/python
    ;;
  nextmol_dmt_l)
    MERGE_PYTHON=/mnt/weka/mbedrosian/codex_dir/nextmol_dmt_l/.venv/bin/python
    ;;
  torsional_diffusion)
    MERGE_PYTHON=/mnt/weka/mbedrosian/codex_dir/torsional_diffusion/.venv_mnt/bin/python
    ;;
  mcf_drugs_l)
    MERGE_PYTHON=/mnt/weka/mbedrosian/codex_dir/mcf_drugs_l/.venv/bin/python
    ;;
  *)
    echo "Unknown MODEL=${MODEL}" >&2
    exit 2
    ;;
esac

export PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd "${REPO_ROOT}"
"${MERGE_PYTHON}" scripts/run_off_the_shelf_druglike.py \
  --model "${MODEL}" \
  --druglike-pickle "${DRUGLIKE_PICKLE}" \
  --output-dir "${OUTPUT_DIR}" \
  --merge-parts

"${PYTHON}" scripts/eval_druglike_conformers.py \
  --label "${LABEL}" \
  --druglike-pickle "${DRUGLIKE_PICKLE}" \
  --gen-pickle "${OUTPUT_DIR}/generation_results.pickle" \
  --out-dir "${EVAL_DIR}" \
  --molecule-workers 16 \
  --pb-workers-per-molecule 1 \
  --energy-num-threads 1

"${PYTHON}" scripts/eval_druglike_covmat.py \
  --ground-truth-pickle "${DRUGLIKE_PICKLE}" \
  --gen-pickle "${OUTPUT_DIR}/generation_results.pickle" \
  --out-dir "${OUTPUT_DIR}/eval_covmat_075" \
  --num-workers 80 \
  --threshold 0.75 \
  --dmax 3
