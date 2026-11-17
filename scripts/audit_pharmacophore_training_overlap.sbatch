#!/bin/bash -l
#SBATCH --job-name=pharm-overlap
#SBATCH --partition=research_cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00

set -euo pipefail

REPO_ROOT="${CASF_REPO_ROOT:-/mnt/weka/mbedrosian/code/casf-benchmark}"
REFERENCE_AUDIT="${REFERENCE_AUDIT:?REFERENCE_AUDIT must be set}"
OUTPUT_DIR="${OUTPUT_DIR:?OUTPUT_DIR must be set}"

cd "${REPO_ROOT}"
exec .venv/bin/python scripts/audit_pharmacophore_training_overlap.py \
  --reference-audit "${REFERENCE_AUDIT}" \
  --output "${OUTPUT_DIR}"
