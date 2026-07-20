#!/bin/bash

# Submit after csf3_check.sbatch has completed successfully.
# Usage: hpc/submit_core_jobs.sh /path/to/data [/path/to/output]

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 DATA_DIR [OUTPUT_DIR]" >&2
  exit 2
fi

data_dir="$(cd "$1" && pwd)"
output_dir="${2:-${data_dir}/results}"
mkdir -p "${output_dir}"
jobscript_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
jobscript="${jobscript_dir}/csf3_analysis.sbatch"
shared_export="ALL,DISSERTATION_DATA_DIR=${data_dir},DISSERTATION_OUTPUT_DIR=${output_dir}"

submit() {
  local analysis="$1"
  local dependency="${2:-}"
  local args=(--parsable --job-name="${analysis}" --export="${shared_export},ANALYSIS=${analysis}")
  if [[ -n "${dependency}" ]]; then
    args+=(--dependency="afterok:${dependency}")
  fi
  local submitted
  submitted="$(sbatch "${args[@]}" "${jobscript}")"
  printf '%s\n' "${submitted%%;*}"
}

rq1_coef_id="$(submit rq1_coef)"
rq1_cv_id="$(submit rq1_cv)"
rq2_joint_id="$(submit rq2_joint)"
rq2_cv_id="$(submit rq2_cv)"
rq3_cv_id="$(submit rq3_cv)"
rq1_robustness_id="$(submit rq1_robustness "${rq1_cv_id}")"
rq2_beyond_id="$(submit rq2_beyond "${rq1_cv_id}")"
rq2_stoplight_id="$(submit rq2_joint_stoplight "${rq2_joint_id}")"

printf '%-24s %s\n' \
  rq1_coef "${rq1_coef_id}" \
  rq1_cv "${rq1_cv_id}" \
  rq1_robustness "${rq1_robustness_id}" \
  rq2_joint "${rq2_joint_id}" \
  rq2_joint_stoplight "${rq2_stoplight_id}" \
  rq2_beyond "${rq2_beyond_id}" \
  rq2_cv "${rq2_cv_id}" \
  rq3_cv "${rq3_cv_id}"
