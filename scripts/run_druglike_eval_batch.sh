#!/usr/bin/env bash
# Run scripts/eval_druglike_conformers.py for every druglike checkpoint listed in
# src/casf_benchmark/config/generation_runs.yaml.
#
# Idempotent: skips any label whose summary CSV is already present.
#
# Usage:
#   CASF_GENERATION_RESULTS_ROOT=... CASF_DRUGLIKE_PICKLE=... CASF_EVAL_OUT_DIR=... \
#     ./scripts/run_druglike_eval_batch.sh
#
# Required environment:
#   CASF_GENERATION_RESULTS_ROOT  root holding the per-checkpoint inference output dirs
#   CASF_DRUGLIKE_PICKLE          druglike test-set pickle (SMILES -> metadata)
#   CASF_EVAL_OUT_DIR             directory for the per-checkpoint eval CSVs
#
# Optional environment: PYTHON, MOLECULE_WORKERS, PB_WORKERS_PER_MOLECULE
#
# `set -e` is deliberately omitted: one failing checkpoint must not abort the batch.
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${PYTHON:-python}"
MOLECULE_WORKERS="${MOLECULE_WORKERS:-16}"
PB_WORKERS_PER_MOLECULE="${PB_WORKERS_PER_MOLECULE:-1}"

: "${CASF_GENERATION_RESULTS_ROOT:?set to the root holding the inference output dirs}"
: "${CASF_DRUGLIKE_PICKLE:?set to the druglike test-set pickle}"
: "${CASF_EVAL_OUT_DIR:?set to the directory for eval outputs}"

LOG_DIR="${CASF_EVAL_OUT_DIR}/_logs"
SUMMARY="${LOG_DIR}/batch_summary.log"
mkdir -p "${CASF_EVAL_OUT_DIR}" "${LOG_DIR}"

log() { echo "[$(date -Is)] $*" | tee -a "${SUMMARY}"; }

log "BATCH START"

while IFS=$'\t' read -r label dirname; do
  out_csv="${CASF_EVAL_OUT_DIR}/${label}_summary.csv"
  if [[ -s "${out_csv}" ]]; then
    log "SKIP ${label} (already done)"
    continue
  fi

  pickle_path="${CASF_GENERATION_RESULTS_ROOT}/${dirname}/generation_results.pickle"
  if [[ ! -f "${pickle_path}" ]]; then
    log "MISSING ${label} (no generation_results.pickle at ${pickle_path})"
    continue
  fi

  item_log="${LOG_DIR}/${label}.log"
  log "RUNNING ${label}"
  if "${PYTHON}" "${REPO_ROOT}/scripts/eval_druglike_conformers.py" \
      --label "${label}" \
      --druglike-pickle "${CASF_DRUGLIKE_PICKLE}" \
      --gen-pickle "${pickle_path}" \
      --out-dir "${CASF_EVAL_OUT_DIR}" \
      --molecule-workers "${MOLECULE_WORKERS}" \
      --pb-workers-per-molecule "${PB_WORKERS_PER_MOLECULE}" \
      >"${item_log}" 2>&1; then
    log "PASS ${label}"
  else
    log "FAIL ${label} (see ${item_log})"
  fi
done < <("${PYTHON}" -c '
from casf_benchmark.catalog import load_generation_runs
for label, dirname in load_generation_runs("druglike"):
    print(f"{label}\t{dirname}")
')

log "BATCH DONE"
