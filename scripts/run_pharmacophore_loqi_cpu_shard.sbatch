#!/bin/bash -l
#SBATCH --partition=research_cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --job-name=pharm-loqi-cpu

set -euo pipefail

REPO_ROOT="${CASF_REPO_ROOT:-/mnt/weka/mbedrosian/code/casf-benchmark}"
SHARD_COUNT="${LOQI_SHARD_COUNT:-4}"
TASK_INDEX="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
MOLECULE_INDEX=$((TASK_INDEX / SHARD_COUNT))
SHARD_INDEX=$((TASK_INDEX % SHARD_COUNT))

source /mnt/weka/mbedrosian/codex_dir/loqi/codex/env.sh
export PYTHONPATH="${REPO_ROOT}/src:/home/mbedrosian/code/3DMolGen/src:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

cd "${REPO_ROOT}"
python scripts/run_off_the_shelf_druglike.py \
  --model loqi \
  --druglike-pickle data/pharmacophore_atc48/generation_input.pickle \
  --output-dir data/results/pharmacophore_atc48/generation/loqi \
  --num-conformers 1000 \
  --molecule-offset "${MOLECULE_INDEX}" \
  --limit-molecules 1 \
  --batch-size 512 \
  --seed 20261117 \
  --conformer-shard-index "${SHARD_INDEX}" \
  --conformer-shard-count "${SHARD_COUNT}"
