#!/bin/bash -l
#SBATCH --job-name=pharm-rdkit
#SBATCH --partition=research_cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=06:00:00

set -euo pipefail

REPO_ROOT="${CASF_REPO_ROOT:-/mnt/weka/mbedrosian/code/casf-benchmark}"
INPUT_PATH="${INPUT_PATH:?INPUT_PATH must be set}"
OUTPUT_DIR="${OUTPUT_DIR:?OUTPUT_DIR must be set}"
INDEX="${SLURM_ARRAY_TASK_ID:?Run this script as an array job}"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd "${REPO_ROOT}"
exec .venv/bin/python scripts/generate_pharmacophore_baselines.py generate \
  --input "${INPUT_PATH}" \
  --output "${OUTPUT_DIR}" \
  --index "${INDEX}" \
  --threads "${SLURM_CPUS_PER_TASK}"
