#!/usr/bin/env bash
# Run the full pickle -> SDF/manifest -> materialize -> analyze pipeline for one
# (checkpoint, cohort) item.
#
# Idempotent: skips any step whose expected output is already present, so it is safe
# to re-run after a partial failure.
#
# Usage:
#   CASF_RUNS_ROOT=... ./scripts/run_casf_batch_item.sh LABEL core|ref PICKLE_DIR
#
# Required environment:
#   CASF_RUNS_ROOT   root for the materialized/analyzed run outputs
#
# Optional environment:
#   CASF_BENCHMARK_DATA_ROOT  CASF/ChEMBL3D input root (default /mnt/weka/mbedrosian)
#   PYTHON                    python to run the converter with (default: python)
#   SHARD_PARALLEL            materialize shards in flight (default: 16 core / 24 ref)
#   PB_WORKERS                PoseBusters workers per shard (default: 2)
#   ANALYZE_WORKERS           geometric-analysis workers (default: SHARD_PARALLEL)
#
# `set -e` is deliberately omitted: step failures are reported and propagated by
# explicit exit-status checks below.
set -uo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $(basename "$0") LABEL core|ref PICKLE_DIR" >&2
  exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LABEL="$1"
COHORT="$(echo "$2" | tr '[:upper:]' '[:lower:]')"
PICKLE_DIR="$3"

case "${COHORT}" in
  core | ref) ;;
  *)
    echo "COHORT must be core or ref, got: $2" >&2
    exit 1
    ;;
esac

: "${CASF_RUNS_ROOT:?set to the root for materialized run outputs}"
WEKA_ROOT="${CASF_BENCHMARK_DATA_ROOT:-/mnt/weka/mbedrosian}"
PYTHON="${PYTHON:-python}"

if [[ "${COHORT}" == "core" ]]; then
  MAP_CSV="${WEKA_ROOT}/data/casf16/casf16_core_chembl3d_exact_intersection.csv"
  LIGAND_DIR="${WEKA_ROOT}/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands"
  SHARD_PARALLEL="${SHARD_PARALLEL:-16}"
else
  MAP_CSV="${WEKA_ROOT}/data/casf16/casf16_ref_chembl3d_exact_intersection.csv"
  LIGAND_DIR="${WEKA_ROOT}/data/casf16/CASF16_REF/ref_chembl3d_exact_intersection_ligands"
  SHARD_PARALLEL="${SHARD_PARALLEL:-24}"
fi

CHEMBL_ROOT="${WEKA_ROOT}/data/chembl3d"
PB_WORKERS="${PB_WORKERS:-2}"
ANALYZE_WORKERS="${ANALYZE_WORKERS:-${SHARD_PARALLEL}}"

METHOD="${LABEL}"
RUN_ROOT="${CASF_RUNS_ROOT}/${LABEL}_${COHORT}"
GEN_DIR="${RUN_ROOT}/generation"
LOG_DIR="${CASF_RUNS_ROOT}/_logs"
LOG="${LOG_DIR}/${LABEL}_${COHORT}.log"
N_LIGANDS=$(($(wc -l < "${MAP_CSV}") - 1))
mkdir -p "${LOG_DIR}"

log() { echo "[$(date -Is)] $*" | tee -a "${LOG}"; }

# One materialize shard: a single molecule, retried a few times because concurrent
# shards occasionally collide on the shared ChEMBL3D cache.
materialize_shard() {
  local offset="$1"
  local attempt
  for attempt in 1 2 3; do
    casf-materialize-generation-sets \
      --root "${RUN_ROOT}" \
      --chembl-map-csv "${MAP_CSV}" \
      --ligand-dir "${LIGAND_DIR}" \
      --chembl-dataset-root "${CHEMBL_ROOT}" \
      --method "${METHOD}" \
      --molecule-offset "${offset}" --limit-molecules 1 \
      --posebusters-workers "${PB_WORKERS}" --posebusters-energy-threads 1 \
      --write-manifest-part >>"${LOG}" 2>&1 && return 0
    echo "[retry] molecule-offset ${offset} attempt ${attempt} failed" >>"${LOG}"
    sleep $((attempt * 3))
  done
  echo "[gave up] molecule-offset ${offset} failed after 3 attempts" >>"${LOG}"
}
export -f materialize_shard
export RUN_ROOT MAP_CSV LIGAND_DIR CHEMBL_ROOT METHOD PB_WORKERS LOG

count_manifest_parts() {
  find "${GEN_DIR}/manifest_parts_normalized" -name '*.tsv' 2>/dev/null | wc -l
}

sweep_shards() {
  seq 0 $((N_LIGANDS - 1)) \
    | xargs -P "${SHARD_PARALLEL}" -I{} bash -c 'materialize_shard "$@"' _ {} \
    >>"${LOG}" 2>&1
}

log "=== START ${LABEL}/${COHORT} (n_ligands_in_mapping=${N_LIGANDS}) ==="

# --- Step 1: convert pickle -> raw SDFs + manifest ---
if [[ -s "${GEN_DIR}/manifest.tsv" ]]; then
  log "convert: manifest already present, skipping"
else
  log "convert: starting"
  if ! "${PYTHON}" "${REPO_ROOT}/scripts/convert_qwen_pickle_to_raw_generation.py" \
      --pickle "${PICKLE_DIR}/generation_results.pickle" \
      --chembl-map-csv "${MAP_CSV}" \
      --run-root "${RUN_ROOT}" \
      --method "${METHOD}" >>"${LOG}" 2>&1; then
    log "convert: FAILED"
    exit 1
  fi
  log "convert: done"
fi

# --- Step 2: materialize (sharded one molecule per process, then merge) ---
TIER_DIR_FIXED="${GEN_DIR}/${METHOD}_fixed"
FIXED_COUNT=$(find "${TIER_DIR_FIXED}" -name '*.sdf' 2>/dev/null | wc -l)

if [[ "${FIXED_COUNT}" -ge "${N_LIGANDS}" ]]; then
  log "materialize: tiers already complete (${FIXED_COUNT}/${N_LIGANDS}), skipping"
else
  log "materialize: sharding ${N_LIGANDS} ligands, parallel=${SHARD_PARALLEL} pb_workers=${PB_WORKERS}"
  sweep_shards
  PARTS_COUNT=$(count_manifest_parts)
  if [[ "${PARTS_COUNT}" -lt "${N_LIGANDS}" ]]; then
    log "materialize: ${PARTS_COUNT}/${N_LIGANDS} manifest parts after first pass, re-sweeping"
    sweep_shards
    log "materialize: $(count_manifest_parts)/${N_LIGANDS} manifest parts after re-sweep"
  fi

  log "materialize: shards done, merging manifest parts"
  casf-materialize-generation-sets \
    --root "${RUN_ROOT}" \
    --chembl-map-csv "${MAP_CSV}" \
    --ligand-dir "${LIGAND_DIR}" \
    --chembl-dataset-root "${CHEMBL_ROOT}" \
    --method "${METHOD}" \
    --merge-manifest-parts >>"${LOG}" 2>&1
  FIXED_COUNT=$(find "${TIER_DIR_FIXED}" -name '*.sdf' 2>/dev/null | wc -l)
  log "materialize: done (${FIXED_COUNT}/${N_LIGANDS} fixed-tier SDFs)"
fi

# --- Step 3: geometric analysis ---
LONG_CSV="${RUN_ROOT}/analysis/tables/geometric_per_ligand_long.csv"
if [[ -s "${LONG_CSV}" ]]; then
  log "analyze: output already present, skipping"
else
  log "analyze: starting (workers=${ANALYZE_WORKERS})"
  if ! casf-analyze-conformer-sets \
      --run-root "${RUN_ROOT}" \
      --chembl-map-csv "${MAP_CSV}" \
      --casf-ligand-dir "${LIGAND_DIR}" \
      --chembl-dataset-root "${CHEMBL_ROOT}" \
      --workers "${ANALYZE_WORKERS}" --resume-parts --quiet-rdkit-warnings \
      >>"${LOG}" 2>&1; then
    log "analyze: FAILED"
    exit 1
  fi
  log "analyze: done"
fi

log "=== DONE ${LABEL}/${COHORT} ==="
