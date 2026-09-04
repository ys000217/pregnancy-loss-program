#!/usr/bin/env bash
# Build REF_FASTA indices once on a login/compute node BEFORE submitting the cohort.
# Single-sample jobs (01_wgs_align.sh) only check that these exist; they do not create them
# (avoids many jobs racing to write the same public index).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "${ROOT}/scripts/lib_source_config.sh"

if [[ ! -s "${REF_FASTA}" ]]; then
  echo "ERROR REF_FASTA not found: ${REF_FASTA}" >&2
  exit 1
fi

need_bwa=0
if [[ ! -s "${REF_FASTA}.bwt.2bit.64" && ! -s "${REF_FASTA}.0123" ]]; then
  need_bwa=1
fi
if [[ "${need_bwa}" -eq 1 ]]; then
  echo "Building bwa-mem2 index: ${REF_FASTA}"
  bwa-mem2 index "${REF_FASTA}"
else
  echo "OK bwa-mem2 index present"
fi

if [[ ! -s "${REF_FASTA}.fai" ]]; then
  echo "Building samtools faidx: ${REF_FASTA}"
  samtools faidx "${REF_FASTA}"
else
  echo "OK ${REF_FASTA}.fai present"
fi

echo "Reference indices ready for ${REF_FASTA}"
