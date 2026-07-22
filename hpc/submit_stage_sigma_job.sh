#!/bin/bash

# Submit the targeted RQ2 stage-specific residual-scale sensitivity model.
# Usage: hpc/submit_stage_sigma_job.sh DATA_DIR [OUTPUT_DIR]

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 DATA_DIR [OUTPUT_DIR]" >&2
  exit 2
fi

data_dir="$(cd "$1" && pwd)"
output_dir="${2:-${data_dir}/results}"
jobscript_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${jobscript_dir}/.." && pwd)"
analysis_jobscript="${jobscript_dir}/csf3_analysis.sbatch"
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

mkdir -p "${output_dir}" "${data_dir}/brm_cache"
manifest_sha256="$(sha256sum "${data_dir}/MANIFEST.sha256" | awk '{print $1}')"
(
  cd "${data_dir}"
  sha256sum -c MANIFEST.sha256
)

exports="ALL,DISSERTATION_REPO_DIR=${repo_dir},DISSERTATION_DATA_DIR=${data_dir},DISSERTATION_OUTPUT_DIR=${output_dir},ANALYSIS=rq2_stage_sigma,EXCLUDE_CONTRASTIVE=false,EXPECTED_GIT_COMMIT=${repo_commit},EXPECTED_MANIFEST_SHA256=${manifest_sha256}"
submitted="$(sbatch --parsable --job-name=rq2_stage_sigma \
  --export="${exports}" "${analysis_jobscript}")"
job_id="${submitted%%;*}"

printf 'Git commit: %s\n' "${repo_commit}"
printf 'Input manifest SHA-256: %s\n' "${manifest_sha256}"
printf 'rq2_stage_sigma %s\n' "${job_id}"
