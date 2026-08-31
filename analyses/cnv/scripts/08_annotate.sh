#!/usr/bin/env bash
# Annotate HIGH-confidence paired CNVs.
set -euo pipefail
source "$(dirname "$0")/lib_source_config.sh"

SAMPLE="${1:?usage: 08_annotate.sh SAMPLE}"
bed="${WORKDIR}/merged/${SAMPLE}.cnv.high.bed"
outdir="${WORKDIR}/annot"
mkdir -p "${outdir}"

if [[ ! -f "${bed}" ]]; then
  echo "ERROR missing ${bed}" >&2
  exit 1
fi

awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$4}' "${bed}" > "${outdir}/${SAMPLE}.high.bed"
n=$(awk 'END{print NR+0}' "${outdir}/${SAMPLE}.high.bed")
if [[ "${n}" -eq 0 ]]; then
  echo "WARN no LARGE_HIGH events for ${SAMPLE}; skip AnnotSV" >&2
  exit 0
fi

if command -v AnnotSV >/dev/null 2>&1; then
  AnnotSV \
    -SVinputFile "${outdir}/${SAMPLE}.high.bed" \
    -genomeBuild "${ANNOTSV_GENOME}" \
    -annotationMode both \
    -outputDir "${outdir}" \
    -outputFile "${SAMPLE}.annotsv"
  echo "AnnotSV: ${outdir}/${SAMPLE}.annotsv.tsv"
else
  echo "WARN AnnotSV not on PATH; wrote ${outdir}/${SAMPLE}.high.bed only" >&2
fi
