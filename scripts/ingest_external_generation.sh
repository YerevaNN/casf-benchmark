#!/usr/bin/env bash
# Materialize tiers and run geometric analysis for one external generator run root.
#
# Usage:
#   ./scripts/ingest_external_generation.sh CORE|REF RUN_ROOT METHOD [METHOD ...]
#
# Examples:
#   # New Qwen ref checkpoint after inference wrote raw SDFs + manifest:
#   ./scripts/ingest_external_generation.sh REF \
#     /mnt/weka/mbedrosian/codex_dir/qwen/generations/casf16_ref_qwen_1k \
#     qwen_4b_bigdata
#
#   # New Qwen core checkpoint into the existing core root:
#   ./scripts/ingest_external_generation.sh CORE \
#     /mnt/weka/mbedrosian/codex_dir/qwen_gens \
#     qwen_4b_revisited
#
# Environment overrides: REPO_ROOT, PYTHON, WEKA_ROOT, WORKERS

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ingest_external_generation.sh CORE|REF RUN_ROOT METHOD [METHOD ...]

Runs, in order:
  1. scripts/materialize_casf_generation_sets.py  (tier materialization)
  2. scripts/analyze_casf_conformer_sets.py       (geometric analysis)

Prerequisites:
  - Analysis conda env (environment-analysis.yml)
  - Untiered SDFs + manifest.tsv under RUN_ROOT/generation/
  - Weka CASF/ChEMBL3D inputs mounted at WEKA_ROOT (default /mnt/weka/mbedrosian)

After success, rebuild the dashboard:
  ./scripts/rebuild_dashboard_weka.sh
EOF
}

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WEKA_ROOT="${CASF_BENCHMARK_DATA_ROOT:-/mnt/weka/mbedrosian}"
PYTHON="${PYTHON:-python}"
WORKERS="${WORKERS:-48}"

COHORT="${1:-}"
RUN_ROOT="${2:-}"
shift 2 || true
METHODS=("$@")

if [[ -z "${COHORT}" || -z "${RUN_ROOT}" || ${#METHODS[@]} -eq 0 ]]; then
  usage >&2
  exit 1
fi

case "${COHORT}" in
  CORE|core) COHORT=core ;;
  REF|ref) COHORT=ref ;;
  *)
    echo "COHORT must be CORE or REF, got: ${COHORT}" >&2
    exit 1
    ;;
esac

if [[ "${COHORT}" == "core" ]]; then
  CHEMBL_MAP="${WEKA_ROOT}/data/casf16/casf16_core_chembl3d_exact_intersection.csv"
  LIGAND_DIR="${WEKA_ROOT}/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands"
else
  CHEMBL_MAP="${WEKA_ROOT}/data/casf16/casf16_ref_chembl3d_exact_intersection.csv"
  LIGAND_DIR="${WEKA_ROOT}/data/casf16/CASF16_REF/ref_chembl3d_exact_intersection_ligands"
fi

CHEMBL_DATASET_ROOT="${WEKA_ROOT}/data/chembl3d"
export PYTHONPATH="${REPO_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"

MATERIALIZE=( "${PYTHON}" "${REPO_ROOT}/scripts/materialize_casf_generation_sets.py"
  --root "${RUN_ROOT}"
  --chembl-map-csv "${CHEMBL_MAP}"
  --ligand-dir "${LIGAND_DIR}"
  --chembl-dataset-root "${CHEMBL_DATASET_ROOT}"
)
for method in "${METHODS[@]}"; do
  MATERIALIZE+=( --method "${method}" )
done

echo "===== materialize $(date -Is) ====="
echo "cohort=${COHORT} root=${RUN_ROOT} methods=${METHODS[*]}"
"${MATERIALIZE[@]}"

echo "===== analyze $(date -Is) ====="
"${PYTHON}" "${REPO_ROOT}/scripts/analyze_casf_conformer_sets.py" \
  --generation-only \
  --ligand-set "${COHORT}" \
  --output-dir "${RUN_ROOT}" \
  --casf-ligand-dir "${LIGAND_DIR}" \
  --chembl-map-csv "${CHEMBL_MAP}" \
  --chembl-dataset-root "${CHEMBL_DATASET_ROOT}" \
  --workers "${WORKERS}" \
  --resume-parts \
  --quiet-rdkit-warnings

echo "===== done $(date -Is) ====="
echo "Analysis table: ${RUN_ROOT}/analysis/tables/geometric_per_ligand_long.csv"
echo "Next: ./scripts/rebuild_dashboard_weka.sh"
