#!/bin/bash
#SBATCH --job-name=casf-pb-once
#SBATCH --partition=research
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --gres=gpu:0
#SBATCH --gpus-per-node=0
#SBATCH --mem=192G
#SBATCH --time=14-00:00:00
#SBATCH --open-mode=append
#SBATCH --output=/home/mbedrosian/code/casf-benchmark/outputs/slurm_jobs/casf_pb_once/%j.out
#SBATCH --error=/home/mbedrosian/code/casf-benchmark/outputs/slurm_jobs/casf_pb_once/%j.err
#
# Full ref analysis includes heavy ligands such as 1tlp (10951 ChEMBL3D conformers)
# and 6o0k/6o0m/6o0p (large flexible generated ensembles). Prior cancelled runs were
# manual terminations, not Slurm timeouts; limits here are sized for completion.

set -euo pipefail

: "${REPO_ROOT:=/home/mbedrosian/code/casf-benchmark}"
: "${PYTHON:=/home/mbedrosian/.conda/envs/chembl3d/bin/python}"
: "${REFERENCE_WORKERS:=24}"
: "${GENERATION_WORKERS:=48}"

export PYTHONPATH="${REPO_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
export MOLGEN3D_QUIET_RDKIT_WARNINGS=1

ANALYZER="${REPO_ROOT}/scripts/analyze_casf_conformer_sets.py"
MASTER_BUILDER="${REPO_ROOT}/scripts/build_casf_analysis_master_csv.py"
DASHBOARD_BUILDER="${REPO_ROOT}/scripts/build_casf_analysis_dashboard_db.py"
ANALYSIS_CONFIG="${REPO_ROOT}/config/casf_analysis_sources.weka.yaml"
DASHBOARD_DB="/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_dashboard.sqlite"
MASTER_CSV="/mnt/weka/mbedrosian/pharma_generation_analysis/casf_analysis_master.csv"
PER_LIGAND_LONG_CSV="/mnt/weka/mbedrosian/pharma_generation_analysis/casf_per_ligand_long.csv"

CORE_MAP="/mnt/weka/mbedrosian/data/casf16/casf16_core_chembl3d_exact_intersection.csv"
REF_MAP="/mnt/weka/mbedrosian/data/casf16/casf16_ref_chembl3d_exact_intersection.csv"
CORE_LIGANDS="/mnt/weka/mbedrosian/data/casf16/CASF16/core_chembl3d_exact_intersection_ligands"
REF_LIGANDS="/mnt/weka/mbedrosian/data/casf16/CASF16_REF/ref_chembl3d_exact_intersection_ligands"

run_step() {
  local name="$1"
  shift
  echo "===== START ${name} $(date -Is) ====="
  "$@"
  echo "===== DONE ${name} $(date -Is) ====="
}

run_reference() {
  local ligand_set="$1"
  local output_dir="$2"
  local ligand_dir="$3"
  local chembl_map="$4"
  run_step "reference-${ligand_set}" \
    "${PYTHON}" -u "${ANALYZER}" \
      --reference-only \
      --ligand-set "${ligand_set}" \
      --output-dir "${output_dir}" \
      --casf-ligand-dir "${ligand_dir}" \
      --chembl-map-csv "${chembl_map}" \
      --workers "${REFERENCE_WORKERS}" \
      --resume-parts \
      --quiet-rdkit-warnings
}

run_generation() {
  local name="$1"
  local ligand_set="$2"
  local output_dir="$3"
  local ligand_dir="$4"
  local chembl_map="$5"
  run_step "generation-${name}" \
    "${PYTHON}" -u "${ANALYZER}" \
      --generation-only \
      --ligand-set "${ligand_set}" \
      --output-dir "${output_dir}" \
      --casf-ligand-dir "${ligand_dir}" \
      --chembl-map-csv "${chembl_map}" \
      --workers "${GENERATION_WORKERS}" \
      --resume-parts \
      --quiet-rdkit-warnings
}

echo "node=$(hostname)"
echo "python=${PYTHON}"
echo "repo=${REPO_ROOT}"
echo "reference_workers=${REFERENCE_WORKERS}"
echo "generation_workers=${GENERATION_WORKERS}"
echo "ref_map=${REF_MAP}"
echo "started=$(date -Is)"

mkdir -p "${REPO_ROOT}/outputs/slurm_jobs/casf_pb_once"

run_reference core \
  /mnt/weka/mbedrosian/pharma_generation_analysis/reference_datasets/core \
  "${CORE_LIGANDS}" \
  "${CORE_MAP}"

run_reference ref \
  /mnt/weka/mbedrosian/pharma_generation_analysis/reference_datasets/ref \
  "${REF_LIGANDS}" \
  "${REF_MAP}"

run_generation rdkit_torsion_core core \
  /mnt/weka/mbedrosian/pharma_generation_analysis/core_pb_full_dynamic_chembl_count \
  "${CORE_LIGANDS}" \
  "${CORE_MAP}"

run_generation rdkit_torsion_ref ref \
  /mnt/weka/mbedrosian/pharma_generation_analysis/ref_pb_full_dynamic_chembl_count \
  "${REF_LIGANDS}" \
  "${REF_MAP}"

run_generation qwen_core core \
  /mnt/weka/mbedrosian/codex_dir/qwen_gens \
  "${CORE_LIGANDS}" \
  "${CORE_MAP}"

run_generation qwen_ref ref \
  /mnt/weka/mbedrosian/codex_dir/qwen/generations/casf16_ref_qwen_1k \
  "${REF_LIGANDS}" \
  "${REF_MAP}"

run_generation loqi_core core \
  /mnt/weka/mbedrosian/codex_dir/loqi/generations/casf16_core_loqi_1k \
  "${CORE_LIGANDS}" \
  "${CORE_MAP}"

run_generation loqi_ref ref \
  /mnt/weka/mbedrosian/codex_dir/loqi/generations/casf16_ref_loqi_1k \
  "${REF_LIGANDS}" \
  "${REF_MAP}"

run_generation nextmol_dmt_l_core core \
  /mnt/weka/mbedrosian/codex_dir/nextmol_dmt_l/generations/casf16_core_nextmol_dmt_l_1k \
  "${CORE_LIGANDS}" \
  "${CORE_MAP}"

run_generation nextmol_dmt_l_ref ref \
  /mnt/weka/mbedrosian/codex_dir/nextmol_dmt_l/generations/casf16_ref_nextmol_dmt_l_1k \
  "${REF_LIGANDS}" \
  "${REF_MAP}"

run_generation torsional_diffusion_core core \
  /mnt/weka/mbedrosian/codex_dir/torsional_diffusion/generations/casf16_core_torsional_diffusion_1k \
  "${CORE_LIGANDS}" \
  "${CORE_MAP}"

run_generation torsional_diffusion_ref ref \
  /mnt/weka/mbedrosian/codex_dir/torsional_diffusion/generations/casf16_ref_torsional_diffusion_1k \
  "${REF_LIGANDS}" \
  "${REF_MAP}"

run_generation mcf_drugs_l_core core \
  /mnt/weka/mbedrosian/codex_dir/mcf_drugs_l/generations/casf16_core_mcf_drugs_l_1k \
  "${CORE_LIGANDS}" \
  "${CORE_MAP}"

run_generation mcf_drugs_l_ref ref \
  /mnt/weka/mbedrosian/codex_dir/mcf_drugs_l/generations/casf16_ref_mcf_drugs_l_1k \
  "${REF_LIGANDS}" \
  "${REF_MAP}"

run_step master-csv "${PYTHON}" -u "${MASTER_BUILDER}" \
  --config "${ANALYSIS_CONFIG}" \
  --output-csv "${MASTER_CSV}" \
  --per-ligand-long-csv "${PER_LIGAND_LONG_CSV}"
run_step dashboard-db "${PYTHON}" -u "${DASHBOARD_BUILDER}" \
  --config "${ANALYSIS_CONFIG}" \
  --output-db "${DASHBOARD_DB}"

echo "finished=$(date -Is)"
