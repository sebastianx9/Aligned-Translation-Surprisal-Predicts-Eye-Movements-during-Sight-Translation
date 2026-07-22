#!/bin/bash

# Submit the small post-audit analysis set after the main CSF run:
#   * reading-stage c_nmt bridge (primary and leave-pair-out);
#   * targeted longer RQ3 coefficient refits (no k-fold rerun);
#   * leave-pair-out diagnostics omitted from the original dependency graph;
#   * posterior-predictive and residual checks after the primary RQ3 refits.
# Usage: hpc/submit_followup_jobs.sh DATA_DIR [OUTPUT_DIR]

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 DATA_DIR [OUTPUT_DIR]" >&2
  exit 2
fi

data_dir="$(cd "$1" && pwd)"
output_dir="${2:-${data_dir}/results}"
sensitivity_output_dir="${output_dir}/exclude_contrastive"
mkdir -p "${output_dir}" "${sensitivity_output_dir}" \
  "${output_dir}/diagnostics"

jobscript_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${jobscript_dir}/.." && pwd)"
analysis_jobscript="${jobscript_dir}/csf3_analysis.sbatch"
diagnostic_jobscript="${jobscript_dir}/csf3_diagnose_fit.sbatch"
assumption_jobscript="${jobscript_dir}/csf3_assumption_checks.sbatch"
repo_commit="$(git -C "${repo_dir}" rev-parse HEAD)"

if [[ -n "$(git -C "${repo_dir}" status --porcelain --untracked-files=no)" ]]; then
  echo "Commit or discard repository changes before submitting CSF jobs." >&2
  git -C "${repo_dir}" status --short >&2
  exit 2
fi
if [[ ! -f "${data_dir}/MANIFEST.sha256" ]]; then
  echo "Missing ${data_dir}/MANIFEST.sha256" >&2
  exit 2
fi
manifest_sha256="$(sha256sum "${data_dir}/MANIFEST.sha256" | awk '{print $1}')"
(
  cd "${data_dir}"
  sha256sum -c MANIFEST.sha256
)

submit_analysis() {
  local job_name="$1"
  local analysis="$2"
  local exclude_contrastive="$3"
  local analysis_output_dir="$4"
  local exports="ALL,DISSERTATION_REPO_DIR=${repo_dir},DISSERTATION_DATA_DIR=${data_dir},DISSERTATION_OUTPUT_DIR=${analysis_output_dir},ANALYSIS=${analysis},EXCLUDE_CONTRASTIVE=${exclude_contrastive},EXPECTED_GIT_COMMIT=${repo_commit},EXPECTED_MANIFEST_SHA256=${manifest_sha256}"
  local submitted
  submitted="$(sbatch --parsable --job-name="${job_name}" \
    --export="${exports}" "${analysis_jobscript}")"
  printf '%s\n' "${submitted%%;*}"
}

submit_diagnostic() {
  local job_name="$1"
  local fit_path="$2"
  if [[ ! -f "${fit_path}" ]]; then
    echo "Missing cached fit for diagnostics: ${fit_path}" >&2
    exit 2
  fi
  local exports="ALL,DISSERTATION_REPO_DIR=${repo_dir},BRMS_FIT_PATH=${fit_path},BRMS_DIAGNOSTIC_OUTPUT_DIR=${output_dir}/diagnostics,EXPECTED_GIT_COMMIT=${repo_commit}"
  local submitted
  submitted="$(sbatch --parsable --job-name="${job_name}" \
    --export="${exports}" "${diagnostic_jobscript}")"
  printf '%s\n' "${submitted%%;*}"
}

bridge_primary_id="$(submit_analysis \
  rq2_reading_cnmt rq2_reading_cnmt_bridge false "${output_dir}")"
bridge_exclude_id="$(submit_analysis \
  rq2_reading_cnmt_excl rq2_reading_cnmt_bridge true \
  "${sensitivity_output_dir}")"
rq3_long_primary_id="$(submit_analysis \
  rq3_coef_long rq3_coefficient_long_refit false "${output_dir}")"
rq3_long_exclude_id="$(submit_analysis \
  rq3_coef_long_excl rq3_coefficient_long_refit true \
  "${sensitivity_output_dir}")"

rq1_exclude_diag_id="$(submit_diagnostic \
  rq1_coef_excl_diag \
  "${data_dir}/brm_cache/rq1_coef_maximal_v2_exclude_contrastive.rds")"
rq2_exclude_diag_id="$(submit_diagnostic \
  rq2_joint_excl_diag \
  "${data_dir}/brm_cache/rq2_joint_maximal_v4_exclude_contrastive.rds")"
assumption_exports="ALL,DISSERTATION_REPO_DIR=${repo_dir},DISSERTATION_DATA_DIR=${data_dir},DISSERTATION_OUTPUT_DIR=${output_dir},EXPECTED_GIT_COMMIT=${repo_commit}"
assumption_id="$(sbatch --parsable --job-name=model_assumptions \
  --dependency="afterok:${rq3_long_primary_id}" \
  --export="${assumption_exports}" "${assumption_jobscript}")"
assumption_id="${assumption_id%%;*}"
printf 'Git commit: %s\n' "${repo_commit}"
printf 'Input manifest SHA-256: %s\n' "${manifest_sha256}"
printf '%-28s %s\n' \
  rq2_reading_cnmt "${bridge_primary_id}" \
  rq2_reading_cnmt_exclude "${bridge_exclude_id}" \
  rq3_coef_long "${rq3_long_primary_id}" \
  rq3_coef_long_exclude "${rq3_long_exclude_id}" \
  rq1_coef_exclude_diagnostics "${rq1_exclude_diag_id}" \
  rq2_joint_exclude_diagnostics "${rq2_exclude_diag_id}" \
  model_assumption_checks "${assumption_id}"
