#!/bin/bash -l
#SBATCH --partition=research_cpu
#SBATCH --cpus-per-task=8
#SBATCH --gres=cpuonly:8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --job-name=dataset-recipes

set -euo pipefail
RECIPE_ROOT=/mnt/weka/mbedrosian/code/casf-benchmark
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "${RECIPE_ROOT}"
source .venv/bin/activate
RECIPE_INDEX_ARGS=()
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    RECIPE_INDEX_ARGS=(--index "${SLURM_ARRAY_TASK_ID}")
fi
exec python -m casf_benchmark.cli.analyze_dataset_recipes "$@" "${RECIPE_INDEX_ARGS[@]}" --workers "${SLURM_CPUS_PER_TASK}"
