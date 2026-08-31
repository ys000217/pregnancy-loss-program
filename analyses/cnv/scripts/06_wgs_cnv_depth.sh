#!/usr/bin/env bash
# Illumina depth CNV with CNVpytor at 100 kb and 500 kb.
set -euo pipefail
source "$(dirname "$0")/lib_source_config.sh"
source "$(dirname "$0")/lib_cnvpytor.sh"

SAMPLE="${1:?usage: 06_wgs_cnv_depth.sh SAMPLE wgs.bam}"
WGS_BAM="${2:?WGS BAM}"

outdir="${WORKDIR}/wgs_cnv"
mkdir -p "${outdir}"
root="${outdir}/${SAMPLE}.pytor"
rm -f "${root}"

cnvpytor_ensure_ref
mapfile -t chroms < <(cnvpytor_primary_chroms)
cnvpytor_rd_his_call "${root}" "${WGS_BAM}" "${outdir}" "${SAMPLE}" "${chroms[@]}"

echo "WGS depth CNV: ${outdir}/${SAMPLE}.cnvpytor.${BIN_PRIMARY}.tsv"
