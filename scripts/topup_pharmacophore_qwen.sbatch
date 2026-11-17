#!/bin/bash -l
#SBATCH --partition=research
#SBATCH --gpus-per-node=h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=192G
#SBATCH --time=02:00:00
#SBATCH --job-name=pharm-qwen-topup

set -euo pipefail

REPO_ROOT="${CASF_REPO_ROOT:-/mnt/weka/mbedrosian/code/casf-benchmark}"
MOLGEN_ROOT="${MOLGEN_ROOT:-/mnt/weka/vtarasov/code/3DMolGen}"
MOLGEN_PYTHON="${MOLGEN_PYTHON:-/mnt/weka/mbedrosian/code/3DMolGen/.venv/bin/python}"

export PYTHONPATH="${MOLGEN_ROOT}/src:${REPO_ROOT}/src:${PYTHONPATH:-}"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd "${REPO_ROOT}"
exec "${MOLGEN_PYTHON}" scripts/topup_pharmacophore_qwen.py "$@"
