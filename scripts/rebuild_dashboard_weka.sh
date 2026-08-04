#!/usr/bin/env bash
# Rebuild production dashboard artifacts on Weka from configured run roots.
#
# Usage:
#   ./scripts/rebuild_dashboard_weka.sh
#
# Reads:  config/casf_analysis_sources.weka.yaml
# Writes: /mnt/weka/mbedrosian/pharma_generation_analysis/casf_{analysis_master,per_ligand_long,analysis_dashboard}.sqlite

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WEKA_ROOT="${CASF_BENCHMARK_DATA_ROOT:-/mnt/weka/mbedrosian}"
PYTHON="${PYTHON:-python}"

ANALYSIS_CONFIG="${REPO_ROOT}/config/casf_analysis_sources.weka.yaml"
DASHBOARD_DB="${WEKA_ROOT}/pharma_generation_analysis/casf_analysis_dashboard.sqlite"
MASTER_CSV="${WEKA_ROOT}/pharma_generation_analysis/casf_analysis_master.csv"
PER_LIGAND_LONG_CSV="${WEKA_ROOT}/pharma_generation_analysis/casf_per_ligand_long.csv"

export PYTHONPATH="${REPO_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"

echo "===== master csv $(date -Is) ====="
"${PYTHON}" "${REPO_ROOT}/scripts/build_casf_analysis_master_csv.py" \
  --config "${ANALYSIS_CONFIG}" \
  --output-csv "${MASTER_CSV}" \
  --per-ligand-long-csv "${PER_LIGAND_LONG_CSV}"

echo "===== dashboard db $(date -Is) ====="
"${PYTHON}" "${REPO_ROOT}/scripts/build_casf_analysis_dashboard_db.py" \
  --config "${ANALYSIS_CONFIG}" \
  --output-db "${DASHBOARD_DB}"

echo "===== done $(date -Is) ====="
echo "Dashboard DB: ${DASHBOARD_DB}"
echo "View locally: CASF_DASHBOARD_DB=${DASHBOARD_DB} streamlit run apps/dashboard/streamlit_app.py"
