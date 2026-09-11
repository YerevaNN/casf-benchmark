#!/usr/bin/env bash
# Submit off-the-shelf druglike inference arrays and dependent evaluations.

set -euo pipefail

REPO_ROOT="${CASF_REPO_ROOT:-/home/mbedrosian/code/casf-benchmark}"
ROOT="${OFF_THE_SHELF_ROOT:-/mnt/weka/mbedrosian/codex_dir/druglike_off_the_shelf}"
DRUGLIKE_PICKLE="${DRUGLIKE_PICKLE:-/mnt/weka/vtarasov/druglike_smi.pickle}"
NUM_CONFORMERS="${NUM_CONFORMERS:-1000}"
MAX_CONCURRENT="${MAX_CONCURRENT:-8}"
MODELS="${MODELS:-loqi nextmol_dmt_l torsional_diffusion mcf_drugs_l}"
EVAL_DIR="${ROOT}/evaluation"

declare -A LABELS=(
  [loqi]="loqi_druglike"
  [nextmol_dmt_l]="nextmol_dmt_l_druglike"
  [torsional_diffusion]="torsional_diffusion_druglike"
  [mcf_drugs_l]="mcf_drugs_l_druglike"
)

mkdir -p "${ROOT}/logs" "${EVAL_DIR}" "${ROOT}/gen_results"

for model in ${MODELS}; do
  label="${LABELS[${model}]:-}"
  if [[ -z "${label}" ]]; then
    echo "Unknown model: ${model}" >&2
    exit 2
  fi
  output_dir="${ROOT}/gen_results/${label}"
  mkdir -p "${output_dir}"

  array_job="$(
    sbatch --parsable \
      --array="0-22%${MAX_CONCURRENT}" \
      --export="ALL,CASF_REPO_ROOT=${REPO_ROOT},MODEL=${model},DRUGLIKE_PICKLE=${DRUGLIKE_PICKLE},OUTPUT_DIR=${output_dir},NUM_CONFORMERS=${NUM_CONFORMERS}" \
      "${REPO_ROOT}/scripts/run_off_the_shelf_druglike.sbatch"
  )"
  eval_job="$(
    sbatch --parsable \
      --dependency="afterok:${array_job}" \
      --export="ALL,CASF_REPO_ROOT=${REPO_ROOT},MODEL=${model},LABEL=${label},DRUGLIKE_PICKLE=${DRUGLIKE_PICKLE},OUTPUT_DIR=${output_dir},EVAL_DIR=${EVAL_DIR}" \
      "${REPO_ROOT}/scripts/finalize_off_the_shelf_druglike.sbatch"
  )"
  printf '%-24s inference=%s evaluation=%s output=%s\n' \
    "${model}" "${array_job}" "${eval_job}" "${output_dir}"
done
