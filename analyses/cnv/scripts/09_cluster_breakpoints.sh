#!/usr/bin/env bash
# Cohort breakpoint clustering of LARGE_HIGH beds. Run after per-sample merge.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "${ROOT}/scripts/lib_source_config.sh"

INDIR="${INDIR:-${WORKDIR}/merged}"
OUTDIR="${OUTDIR:-${WORKDIR}/cohort}"

python3 "${ROOT}/scripts/09_cluster_breakpoints.py" \
  --indir "${INDIR}" \
  --glob "*.cnv.high.bed" \
  --bp-pad "${CLUSTER_BP_PAD:-100000}" \
  --bp-pad-large "${CLUSTER_BP_PAD_LARGE:-500000}" \
  --ro "${MERGE_RO:-0.50}" \
  -o "${OUTDIR}"
