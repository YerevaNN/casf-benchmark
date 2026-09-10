#!/usr/bin/env bash
# Run scripts/run_casf_batch_item.sh for every (checkpoint, cohort) pair listed in
# src/casf_benchmark/config/generation_runs.yaml, sequentially. Each item
# internally shards across the node's CPUs.
#
# Core items run first (~94 ligands, cheap) so results land early; ref items
# (~1236 ligands each) run afterwards and take considerably longer.
#
# Usage:
#   CASF_GENERATION_RESULTS_ROOT=... CASF_RUNS_ROOT=... ./scripts/run_casf_batch_all.sh
#
# Required environment:
#   CASF_GENERATION_RESULTS_ROOT  root holding the per-checkpoint inference output dirs
#   CASF_RUNS_ROOT                root for the materialized/analyzed run outputs
#
# Optional environment: CASF_BENCHMARK_DATA_ROOT, PYTHON (see run_casf_batch_item.sh)
#
# `set -e` is deliberately omitted: one failing item must not abort the batch.
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${PYTHON:-python}"

: "${CASF_GENERATION_RESULTS_ROOT:?set to the root holding the inference output dirs}"
: "${CASF_RUNS_ROOT:?set to the root for materialized run outputs}"
export CASF_RUNS_ROOT  # run_casf_batch_item.sh reads it from the environment

LOG_DIR="${CASF_RUNS_ROOT}/_logs"
SUMMARY="${LOG_DIR}/batch_summary.log"
mkdir -p "${LOG_DIR}"

log() { echo "[$(date -Is)] $*" | tee -a "${SUMMARY}"; }

log "BATCH START"

for cohort in core ref; do
  while IFS=$'\t' read -r label dirname; do
    log "RUNNING ${label}/${cohort}"
    if bash "${REPO_ROOT}/scripts/run_casf_batch_item.sh" \
        "${label}" "${cohort}" "${CASF_GENERATION_RESULTS_ROOT}/${dirname}"; then
      log "PASS ${label}/${cohort}"
    else
      log "FAIL ${label}/${cohort}"
    fi
  done < <("${PYTHON}" -c '
import sys
from casf_benchmark.catalog import load_generation_runs
for label, dirname in load_generation_runs(f"casf16_{sys.argv[1]}"):
    print(f"{label}\t{dirname}")
' "${cohort}")
done

log "BATCH DONE"
