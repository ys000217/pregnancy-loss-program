#!/usr/bin/env bash
# Median depth + breadth. At 10x, report median not mean; fail if median < MIN_MEDIAN_DEPTH.
set -euo pipefail
source "$(dirname "$0")/lib_source_config.sh"

SAMPLE="${1:?usage: 02_qc_coverage.sh SAMPLE ont.bam wgs.bam}"
ONT_BAM="${2:?ONT BAM}"
WGS_BAM="${3:?WGS BAM}"

qcdir="${WORKDIR}/qc"
mkdir -p "${qcdir}"
if [[ ! -s "${ONT_BAM}" || ! -s "${WGS_BAM}" ]]; then
  echo "ERROR BAM missing: ${ONT_BAM} ${WGS_BAM}" >&2
  exit 1
fi

run_mosdepth() {
  local prefix="$1"
  local bam="$2"
  mosdepth --by "${BIN_PRIMARY}" -Q "${MOSDEPTH_MAPQ}" -t 4 "${qcdir}/${prefix}" "${bam}"
}

run_mosdepth "${SAMPLE}.ont" "${ONT_BAM}"
run_mosdepth "${SAMPLE}.wgs" "${WGS_BAM}"

python3 "$(dirname "$0")/02_qc_stats.py" "${SAMPLE}" "${qcdir}" "${MIN_MEDIAN_DEPTH}" "${MIN_BREADTH_5X}" "${QC_ABORT:-1}" "${QC_ABORT_MEDIAN:-3}" "${QC_ABORT_BREADTH:-0.40}"
