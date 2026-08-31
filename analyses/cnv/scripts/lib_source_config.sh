#!/usr/bin/env bash
# Source config.sh after stripping Windows CR so cluster bash does not see `$'\r'`.
# Export _CNV_ROOT so config.sh can resolve ref/ without BASH_SOURCE (stdin source).
export _CNV_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1090
source /dev/stdin <<< "$(tr -d '\r' < "${_CNV_ROOT}/config.sh")"
# Belt-and-suspenders if an older config.sh was synced without merge knobs
export MAX_CNV_EVENT="${MAX_CNV_EVENT:-10000000}"
export MERGE_MASK_FRAC="${MERGE_MASK_FRAC:-0.50}"
if [[ -z "${HARD_MASK_BED:-}" || ! -s "${HARD_MASK_BED}" ]]; then
  if [[ -s "${_CNV_ROOT}/ref/hard_mask.grch38.refseq.bed" ]]; then
    export HARD_MASK_BED="${_CNV_ROOT}/ref/hard_mask.grch38.refseq.bed"
  fi
fi
