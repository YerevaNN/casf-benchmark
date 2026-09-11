#!/usr/bin/env bash
# Merge Qwen + off-the-shelf druglike eval CSVs and publish extended_* tables.
#
# Usage:
#   ./scripts/publish_druglike_dashboard.sh
#
# Env overrides:
#   CASF_QWEN_DRUGLIKE_EVAL   Qwen PB/diversity CSV root
#   CASF_OTS_DRUGLIKE_EVAL    Off-the-shelf PB/diversity CSV root
#   CASF_DRUGLIKE_DB          Output SQLite (default: OTS eval dir)
#   CASF_GEN_RESULTS_ROOT     Qwen inference root (absolute OTS paths in YAML still resolve)
#   CASF_EXTENDED_DB          Sidecar DB the dashboard reads

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QWEN_EVAL="${CASF_QWEN_DRUGLIKE_EVAL:-/mnt/weka/vtarasov/outputs/casf_benchmark_runs/druglike_eval}"
OTS_EVAL="${CASF_OTS_DRUGLIKE_EVAL:-/mnt/weka/mbedrosian/codex_dir/druglike_off_the_shelf/evaluation}"
DRUGLIKE_DB="${CASF_DRUGLIKE_DB:-${OTS_EVAL}/druglike_eval.sqlite}"
GEN_ROOT="${CASF_GEN_RESULTS_ROOT:-/mnt/weka/vtarasov/outputs/outputs/gen_results}"
EXTENDED_DB="${CASF_EXTENDED_DB:-${REPO_ROOT}/data/results/extended_casf_analysis.sqlite}"

PYTHON="${PYTHON:-${REPO_ROOT}/.venv/bin/python}"
if [[ ! -x "${PYTHON}" ]]; then
  PYTHON="$(command -v python3)"
fi

echo "===== build druglike_eval.sqlite $(date -Is) ====="
"${PYTHON}" "${REPO_ROOT}/scripts/build_druglike_eval_db.py" \
  --eval-dir "${QWEN_EVAL}" \
  --eval-dir "${OTS_EVAL}" \
  --db-path "${DRUGLIKE_DB}"

echo "===== merge COV/MAT + publish extended tables $(date -Is) ====="
"${PYTHON}" "${REPO_ROOT}/scripts/build_druglike_covmat.py" \
  --generation-results-root "${GEN_ROOT}" \
  --druglike-db "${DRUGLIKE_DB}" \
  --extended-db "${EXTENDED_DB}"

echo "===== done $(date -Is) ====="
echo "Druglike DB: ${DRUGLIKE_DB}"
echo "Extended DB: ${EXTENDED_DB}"
echo "View: CASF_EXTENDED_DB=${EXTENDED_DB} streamlit run ${REPO_ROOT}/apps/dashboard/streamlit_app.py"
