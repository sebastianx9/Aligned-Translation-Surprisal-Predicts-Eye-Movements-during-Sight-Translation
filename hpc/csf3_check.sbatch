#!/bin/bash --login
#SBATCH -p multicore
#SBATCH -n 2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=dissertation-check
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail
module purge
module load apps/gcc/R/4.4.1

: "${DISSERTATION_DATA_DIR:?Set DISSERTATION_DATA_DIR before sbatch}"
export R_LIBS_USER="${R_LIBS_USER:-${HOME}/R/csf3-4.4}"

jobscript_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${jobscript_dir}/.." && pwd)"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript "${repo_dir}/hpc/check_r_environment.R" \
  "--data-dir=${DISSERTATION_DATA_DIR}"
