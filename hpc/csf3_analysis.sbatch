#!/bin/bash --login
#SBATCH -p multicore
#SBATCH -n 4
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --job-name=dissertation-model
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail
module purge
module load apps/gcc/R/4.4.1

: "${ANALYSIS:?Set ANALYSIS before sbatch}"
: "${DISSERTATION_DATA_DIR:?Set DISSERTATION_DATA_DIR before sbatch}"
export R_LIBS_USER="${R_LIBS_USER:-${HOME}/R/csf3-4.4}"

exclude_contrastive="${EXCLUDE_CONTRASTIVE:-false}"
case "${exclude_contrastive}" in
  true|false) ;;
  *)
    echo "EXCLUDE_CONTRASTIVE must be true or false, got: ${exclude_contrastive}" >&2
    exit 2
    ;;
esac

repo_dir="${DISSERTATION_REPO_DIR:-${SLURM_SUBMIT_DIR}}"
output_dir="${DISSERTATION_OUTPUT_DIR:-${DISSERTATION_DATA_DIR}/results}"
mkdir -p "${output_dir}" "${DISSERTATION_DATA_DIR}/brm_cache"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Analysis: ${ANALYSIS}"
echo "Repository: ${repo_dir}"
echo "Data: ${DISSERTATION_DATA_DIR}"
echo "Output: ${output_dir}"
echo "Exclude S031/S032: ${exclude_contrastive}"
Rscript -e 'print(sessionInfo())'
sha256sum \
  "${DISSERTATION_DATA_DIR}/fixation_durations_word.csv" \
  "${DISSERTATION_DATA_DIR}/eye_measures_word.csv" \
  "${DISSERTATION_DATA_DIR}/monolingual_surprisal_word.csv" \
  "${DISSERTATION_DATA_DIR}/nmt_surprisal_soft_word.csv" \
  "${DISSERTATION_DATA_DIR}/nmt_alignment_mass_word.csv" \
  "${DISSERTATION_DATA_DIR}/attention_features_6_norm.csv" \
  "${DISSERTATION_DATA_DIR}/subtlex_us.csv"

run_r() {
  Rscript "$1" \
    "--data-dir=${DISSERTATION_DATA_DIR}" \
    "--output-dir=${output_dir}" \
    "--exclude-contrastive=${exclude_contrastive}"
}

case "${ANALYSIS}" in
  rq1_coef)
    run_r "${repo_dir}/RQ1/rq1_coef_maximal.R"
    ;;
  rq1_cv)
    run_r "${repo_dir}/RQ1/rq1_kfold_elpd.R"
    run_r "${repo_dir}/RQ1/rq1_direct_nmt_vs_mono.R"
    ;;
  rq1_robustness)
    run_r "${repo_dir}/RQ1/rq1_mass_stoplight_robustness.R"
    ;;
  rq1_locus)
    run_r "${repo_dir}/RQ1/rq_locus_kfold.R"
    ;;
  rq2_joint)
    run_r "${repo_dir}/RQ2/rq2_joint_maximal.R"
    ;;
  rq2_joint_stoplight)
    Rscript "${repo_dir}/RQ2/rq2_joint_maximal.R" \
      "--data-dir=${DISSERTATION_DATA_DIR}" \
      "--output-dir=${output_dir}" \
      "--exclude-contrastive=${exclude_contrastive}" \
      "--include-stoplight=true"
    ;;
  rq2_beyond)
    run_r "${repo_dir}/RQ2/rq2_beyond_kfold.R"
    ;;
  rq2_cv)
    run_r "${repo_dir}/RQ2/rq2_kfold_elpd.R"
    ;;
  rq2_reading_validation)
    run_r "${repo_dir}/RQ2/rq2_reading_cmono_validation.R"
    ;;
  rq3_cv)
    run_r "${repo_dir}/RQ3/rq3_kfold_elpd.R"
    ;;
  *)
    echo "Unknown ANALYSIS=${ANALYSIS}" >&2
    exit 2
    ;;
esac
