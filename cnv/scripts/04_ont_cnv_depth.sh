#!/usr/bin/env bash
# ONT depth CNV: Spectre optional, CNVpytor required.
set -euo pipefail
source "$(dirname "$0")/lib_source_config.sh"
source "$(dirname "$0")/lib_cnvpytor.sh"

SAMPLE="${1:?usage: 04_ont_cnv_depth.sh SAMPLE ont.bam}"
ONT_BAM="${2:?ONT BAM}"

outdir="${WORKDIR}/ont_cnv"
mkdir -p "${outdir}"

call_cnvpytor() {
  local bam="$1"
  local root="${outdir}/${SAMPLE}.pytor"
  cnvpytor_ensure_ref
  mapfile -t chroms < <(cnvpytor_primary_chroms)
  cnvpytor_rd_his_call "${root}" "${bam}" "${outdir}" "${SAMPLE}" "${chroms[@]}"
}

if command -v spectre >/dev/null 2>&1; then
  cov_prefix="${WORKDIR}/qc/${SAMPLE}.ont"
  if [[ ! -s "${cov_prefix}.regions.bed.gz" ]]; then
    mosdepth --by "${BIN_PRIMARY}" -Q "${MOSDEPTH_MAPQ}" -t 4 "${cov_prefix}" "${ONT_BAM}"
  fi
  spectre CNVCaller \
    --coverage "${cov_prefix}.regions.bed.gz" \
    --sample-id "${SAMPLE}" \
    --output-dir "${outdir}/spectre" \
    --reference "${REF_FASTA}" \
    --min-cnv-len "${MIN_DEPTH_CNV}" \
    || echo "WARN Spectre failed; continuing with CNVpytor" >&2
fi

call_cnvpytor "${ONT_BAM}"
echo "ONT depth CNV: ${outdir}/${SAMPLE}.cnvpytor.${BIN_PRIMARY}.tsv"
