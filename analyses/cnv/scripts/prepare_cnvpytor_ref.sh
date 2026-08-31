#!/usr/bin/env bash
# One-time on a login/compute node before the cohort (GC scan of the FASTA is slow).
#   conda activate cnv10x
#   bash scripts/prepare_cnvpytor_ref.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/scripts/lib_source_config.sh"
source "${ROOT}/scripts/lib_cnvpytor.sh"
cnvpytor_ensure_ref
echo "conf: ${cnvpytor_conf}"
echo "gc:   ${cnvpytor_gc}"
echo "Next: ONLY=0002C bash scripts/submit_per_sample.sh"
