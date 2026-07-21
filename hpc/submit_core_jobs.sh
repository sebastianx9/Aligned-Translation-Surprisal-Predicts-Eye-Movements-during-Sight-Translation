#!/bin/bash

# Submit after csf3_check.sbatch has completed successfully.
# Usage: hpc/submit_core_jobs.sh /path/to/data [/path/to/output]
#
# The primary jobs retain S031/S032. A matched set of core sensitivity jobs
# refits the models after excluding the shared contrastive pair. The latter use
# a separate output directory, and the R scripts use distinct cache names.

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 DATA_DIR [OUTPUT_DIR]" >&2
  exit 2
fi

data_dir="$(cd "$1" && pwd)"
output_dir="${2:-${data_dir}/results}"
mkdir -p "${output_dir}"
sensitivity_output_dir="${output_dir}/exclude_contrastive"
mkdir -p "${sensitivity_output_dir}"
jobscript_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${jobscript_dir}/.." && pwd)"
jobscript="${jobscript_dir}/csf3_analysis.sbatch"

submit() {
  local analysis="$1"
  local dependency="${2:-}"
  local exclude_contrastive="${3:-false}"
  local job_name="${4:-${analysis}}"
  local analysis_output_dir="${5:-${output_dir}}"
  local job_export="ALL,DISSERTATION_REPO_DIR=${repo_dir},DISSERTATION_DATA_DIR=${data_dir},DISSERTATION_OUTPUT_DIR=${analysis_output_dir},ANALYSIS=${analysis},EXCLUDE_CONTRASTIVE=${exclude_contrastive}"
  local args=(--parsable --job-name="${job_name}" --export="${job_export}")
  if [[ -n "${dependency}" ]]; then
    args+=(--dependency="afterok:${dependency}")
  fi
  local submitted
  submitted="$(sbatch "${args[@]}" "${jobscript}")"
  printf '%s\n' "${submitted%%;*}"
}

rq1_coef_id="$(submit rq1_coef)"
rq1_cv_id="$(submit rq1_cv)"
rq1_joint_id="$(submit rq1_joint_predictive "${rq1_cv_id}")"
rq2_joint_id="$(submit rq2_joint)"
rq2_interaction_id="$(submit rq2_interaction_cv)"
rq3_cv_id="$(submit rq3_cv)"
rq1_locus_id="$(submit rq1_locus)"
rq2_reading_validation_id="$(submit rq2_reading_validation)"
rq1_robustness_id="$(submit rq1_robustness "${rq1_cv_id}")"
rq2_stoplight_id="$(submit rq2_joint_stoplight "${rq2_joint_id}")"

# Leave-pair-out sensitivity analyses. These are matched refits, not results
# obtained by subtracting S031/S032 from the primary pointwise ELPD values.
rq1_coef_excl_id="$(submit rq1_coef "${rq1_coef_id}" true rq1_coef_excl "${sensitivity_output_dir}")"
rq1_cv_excl_id="$(submit rq1_cv "${rq1_cv_id}" true rq1_cv_excl "${sensitivity_output_dir}")"
rq1_joint_excl_id="$(submit rq1_joint_predictive "${rq1_cv_excl_id}" true rq1_joint_excl "${sensitivity_output_dir}")"
rq2_joint_excl_id="$(submit rq2_joint "${rq2_joint_id}" true rq2_joint_excl "${sensitivity_output_dir}")"
rq2_interaction_excl_id="$(submit rq2_interaction_cv "${rq2_interaction_id}" true rq2_interaction_excl "${sensitivity_output_dir}")"
rq3_cv_excl_id="$(submit rq3_cv "${rq3_cv_id}" true rq3_cv_excl "${sensitivity_output_dir}")"
rq1_locus_excl_id="$(submit rq1_locus "${rq1_locus_id}" true rq1_locus_excl "${sensitivity_output_dir}")"
rq2_reading_validation_excl_id="$(submit rq2_reading_validation "${rq2_reading_validation_id}" true rq2_reading_excl "${sensitivity_output_dir}")"

printf '%-24s %s\n' \
  rq1_coef "${rq1_coef_id}" \
  rq1_cv "${rq1_cv_id}" \
  rq1_joint_predictive "${rq1_joint_id}" \
  rq1_robustness "${rq1_robustness_id}" \
  rq2_joint "${rq2_joint_id}" \
  rq2_joint_stoplight "${rq2_stoplight_id}" \
  rq2_interaction_cv "${rq2_interaction_id}" \
  rq3_cv "${rq3_cv_id}" \
  rq1_locus "${rq1_locus_id}" \
  rq2_reading_validation "${rq2_reading_validation_id}" \
  rq1_coef_exclude_pair "${rq1_coef_excl_id}" \
  rq1_cv_exclude_pair "${rq1_cv_excl_id}" \
  rq1_joint_exclude_pair "${rq1_joint_excl_id}" \
  rq2_joint_exclude_pair "${rq2_joint_excl_id}" \
  rq2_interaction_exclude_pair "${rq2_interaction_excl_id}" \
  rq3_cv_exclude_pair "${rq3_cv_excl_id}" \
  rq1_locus_exclude_pair "${rq1_locus_excl_id}" \
  rq2_reading_exclude_pair "${rq2_reading_validation_excl_id}"
