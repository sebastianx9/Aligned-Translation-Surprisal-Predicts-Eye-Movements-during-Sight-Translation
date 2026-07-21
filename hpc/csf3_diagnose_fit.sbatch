#!/bin/bash --login
#SBATCH -p multicore
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --job-name=brms-diagnostics
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail
module purge
module load apps/gcc/R/4.4.1

: "${BRMS_FIT_PATH:?Set BRMS_FIT_PATH before sbatch}"
export R_LIBS_USER="${R_LIBS_USER:-${HOME}/R/csf3-4.4}"

repo_dir="${DISSERTATION_REPO_DIR:-${SLURM_SUBMIT_DIR}}"
output_dir="${BRMS_DIAGNOSTIC_OUTPUT_DIR:-$(dirname "${BRMS_FIT_PATH}")}"

actual_git_commit="$(git -C "${repo_dir}" rev-parse HEAD)"
echo "Git commit: ${actual_git_commit}"
if [[ -n "$(git -C "${repo_dir}" status --porcelain --untracked-files=no)" ]]; then
  echo "Repository working tree is not clean; refusing non-reproducible diagnostics." >&2
  git -C "${repo_dir}" status --short >&2
  exit 2
fi
if [[ -n "${EXPECTED_GIT_COMMIT:-}" && \
      "${actual_git_commit}" != "${EXPECTED_GIT_COMMIT}" ]]; then
  echo "Repository changed after submission: expected ${EXPECTED_GIT_COMMIT}, got ${actual_git_commit}" >&2
  exit 2
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript "${repo_dir}/hpc/diagnose_brms_fit.R" \
  "--fit-path=${BRMS_FIT_PATH}" \
  "--output-dir=${output_dir}"
