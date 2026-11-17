#!/bin/bash -l
#SBATCH --partition=research_cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --job-name=pharm-pb

set -euo pipefail

REPO_ROOT="${CASF_REPO_ROOT:-/mnt/weka/mbedrosian/code/casf-benchmark}"
CONFIG="${PHARMACOPHORE_CONFIG:?Set PHARMACOPHORE_CONFIG}"
METHODS_LIST="${PHARMACOPHORE_METHODS:?Set PHARMACOPHORE_METHODS}"
MOLECULE_COUNT="${PHARMACOPHORE_MOLECULE_COUNT:?Set PHARMACOPHORE_MOLECULE_COUNT}"

IFS=':' read -r -a METHODS <<< "${METHODS_LIST}"
TASK_INDEX="${SLURM_ARRAY_TASK_ID:?Run this launcher as an array}"
METHOD_INDEX=$((TASK_INDEX / MOLECULE_COUNT))
MOLECULE_INDEX=$((TASK_INDEX % MOLECULE_COUNT))

if (( METHOD_INDEX >= ${#METHODS[@]} )); then
    echo "Array task ${TASK_INDEX} does not map to a configured method" >&2
    exit 2
fi

cd "${REPO_ROOT}"
source .venv/bin/activate
exec python -m casf_benchmark.cli.analyze_pharmacophores \
    --config "${CONFIG}" \
    posebusters \
    --method "${METHODS[${METHOD_INDEX}]}" \
    --molecule-index "${MOLECULE_INDEX}" \
    --no-finalize
