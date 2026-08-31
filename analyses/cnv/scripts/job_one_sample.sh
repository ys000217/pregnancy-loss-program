#!/usr/bin/env bash
# One sample, run inside a scheduled job (not a submitter).
# Usage: bash scripts/job_one_sample.sh SAMPLE_ID
set -euo pipefail
SAMPLE="${1:?usage: job_one_sample.sh SAMPLE_ID}"

# Login shells on Jhinno often lack conda until profile is sourced.
if [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
  # shellcheck disable=SC1091
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "${HOME}/.bashrc" ]]; then
  # shellcheck disable=SC1091
  source "${HOME}/.bashrc"
fi
conda activate cnv10x

CNV_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "${CNV_DIR}"
bash scripts/run_sample.sh "${SAMPLE}"
