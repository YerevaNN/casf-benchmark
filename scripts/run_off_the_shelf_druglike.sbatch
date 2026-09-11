#!/bin/bash -l
#SBATCH --job-name=druglike-ots
#SBATCH --partition=research
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/weka/mbedrosian/codex_dir/druglike_off_the_shelf/logs/%x_%A_%a.out
#SBATCH --error=/mnt/weka/mbedrosian/codex_dir/druglike_off_the_shelf/logs/%x_%A_%a.err

set -euo pipefail

REPO_ROOT="${CASF_REPO_ROOT:-/home/mbedrosian/code/casf-benchmark}"
MODEL="${MODEL:?MODEL must be set}"
DRUGLIKE_PICKLE="${DRUGLIKE_PICKLE:-/mnt/weka/vtarasov/druglike_smi.pickle}"
OUTPUT_DIR="${OUTPUT_DIR:?OUTPUT_DIR must be set}"
NUM_CONFORMERS="${NUM_CONFORMERS:-1000}"
OFFSET="${SLURM_ARRAY_TASK_ID:-0}"

case "${MODEL}" in
  loqi)
    source /mnt/weka/mbedrosian/codex_dir/loqi/codex/env.sh
    BATCH_SIZE="${BATCH_SIZE:-512}"
    ;;
  nextmol_dmt_l)
    source /mnt/weka/mbedrosian/codex_dir/nextmol_dmt_l/codex/env.sh
    BATCH_SIZE="${BATCH_SIZE:-250}"
    ;;
  torsional_diffusion)
    source /mnt/weka/mbedrosian/codex_dir/torsional_diffusion/codex/env.sh
    BATCH_SIZE="${BATCH_SIZE:-128}"
    ;;
  mcf_drugs_l)
    source /mnt/weka/mbedrosian/codex_dir/mcf_drugs_l/codex/env.sh
    BATCH_SIZE="${BATCH_SIZE:-250}"
    ;;
  *)
    echo "Unknown MODEL=${MODEL}" >&2
    exit 2
    ;;
esac

export PYTHONPATH="${REPO_ROOT}/src:/home/mbedrosian/code/3DMolGen/src:${PYTHONPATH:-}"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd "${REPO_ROOT}"
python scripts/run_off_the_shelf_druglike.py \
  --model "${MODEL}" \
  --druglike-pickle "${DRUGLIKE_PICKLE}" \
  --output-dir "${OUTPUT_DIR}" \
  --num-conformers "${NUM_CONFORMERS}" \
  --molecule-offset "${OFFSET}" \
  --limit-molecules 1 \
  --batch-size "${BATCH_SIZE}"
