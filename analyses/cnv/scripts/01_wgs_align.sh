#!/usr/bin/env bash
# Align paired-end Illumina FASTQ to the same GRCh38 used for ONT BAMs.
# R1/R2 may be comma-separated lane files (Novogene _L1_ / _L2_). Concatenate then fastp.
set -euo pipefail
source "$(dirname "$0")/lib_source_config.sh"

SAMPLE="${1:?usage: 01_wgs_align.sh SAMPLE}"
R1_LIST="${2:?R1 FASTQ (comma-separated lanes OK)}"
R2_LIST="${3:?R2 FASTQ (comma-separated lanes OK)}"

outdir="${WORKDIR}/wgs_bam"
tmpdir="${WORKDIR}/tmp/${SAMPLE}"
mkdir -p "${outdir}" "${tmpdir}" "${WORKDIR}/qc"

split_csv() {
  local s="$1"
  s="${s//,/ }"
  # shellcheck disable=SC2206
  local arr=($s)
  printf '%s\n' "${arr[@]}"
}

mapfile -t R1S < <(split_csv "${R1_LIST}")
mapfile -t R2S < <(split_csv "${R2_LIST}")
if [[ ${#R1S[@]} -ne ${#R2S[@]} ]]; then
  echo "ERROR R1/R2 lane count mismatch" >&2
  exit 1
fi
for f in "${R1S[@]}" "${R2S[@]}"; do
  if [[ ! -s "${f}" ]]; then
    echo "ERROR missing FASTQ: ${f}" >&2
    exit 1
  fi
done
if [[ ! -s "${REF_FASTA}" ]]; then
  echo "ERROR REF_FASTA not found: ${REF_FASTA}" >&2
  exit 1
fi
# Do not build indexes here — many jsub jobs would race on the same files.
# Run once before the cohort: bash scripts/prepare_ref_index.sh
if [[ ! -s "${REF_FASTA}.bwt.2bit.64" && ! -s "${REF_FASTA}.0123" ]]; then
  echo "ERROR bwa-mem2 index missing for ${REF_FASTA}" >&2
  echo "  bash scripts/prepare_ref_index.sh" >&2
  exit 1
fi
if [[ ! -s "${REF_FASTA}.fai" ]]; then
  echo "ERROR samtools faidx missing: ${REF_FASTA}.fai" >&2
  echo "  bash scripts/prepare_ref_index.sh" >&2
  exit 1
fi

cat "${R1S[@]}" > "${tmpdir}/${SAMPLE}.r1.cat.fq.gz"
cat "${R2S[@]}" > "${tmpdir}/${SAMPLE}.r2.cat.fq.gz"

fastp \
  -i "${tmpdir}/${SAMPLE}.r1.cat.fq.gz" \
  -I "${tmpdir}/${SAMPLE}.r2.cat.fq.gz" \
  -o "${tmpdir}/${SAMPLE}.r1.fq.gz" \
  -O "${tmpdir}/${SAMPLE}.r2.fq.gz" \
  --thread 8 \
  --detect_adapter_for_pe \
  --html "${WORKDIR}/qc/${SAMPLE}.fastp.html" \
  --json "${WORKDIR}/qc/${SAMPLE}.fastp.json" \
  --disable_quality_filtering

bwa-mem2 mem -t "${THREADS}" \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
  "${REF_FASTA}" \
  "${tmpdir}/${SAMPLE}.r1.fq.gz" \
  "${tmpdir}/${SAMPLE}.r2.fq.gz" \
  | samtools sort -@ 4 -o "${tmpdir}/${SAMPLE}.sorted.bam"

if command -v sambamba >/dev/null 2>&1; then
  sambamba markdup -t "${THREADS}" \
    "${tmpdir}/${SAMPLE}.sorted.bam" \
    "${outdir}/${SAMPLE}.markdup.bam"
else
  # Proper samtools markdup needs name-order + fixmate -m before coordinate sort
  echo "WARN sambamba missing; using samtools fixmate + markdup" >&2
  samtools sort -n -@ 4 -o "${tmpdir}/${SAMPLE}.namesort.bam" "${tmpdir}/${SAMPLE}.sorted.bam"
  samtools fixmate -m -@ 4 \
    "${tmpdir}/${SAMPLE}.namesort.bam" \
    "${tmpdir}/${SAMPLE}.fixmate.bam"
  samtools sort -@ 4 -o "${tmpdir}/${SAMPLE}.coordsort.bam" "${tmpdir}/${SAMPLE}.fixmate.bam"
  samtools markdup -@ "${THREADS}" -s \
    "${tmpdir}/${SAMPLE}.coordsort.bam" \
    "${outdir}/${SAMPLE}.markdup.bam"
  samtools index -@ 4 "${outdir}/${SAMPLE}.markdup.bam"
  rm -f "${tmpdir}/${SAMPLE}.namesort.bam" \
        "${tmpdir}/${SAMPLE}.fixmate.bam" \
        "${tmpdir}/${SAMPLE}.coordsort.bam"
fi

if [[ ! -s "${outdir}/${SAMPLE}.markdup.bam.bai" && ! -s "${outdir}/${SAMPLE}.markdup.bam.csi" ]]; then
  samtools index -@ 4 "${outdir}/${SAMPLE}.markdup.bam"
fi

rm -f "${tmpdir}/${SAMPLE}.sorted.bam" \
      "${tmpdir}/${SAMPLE}.r1.fq.gz" \
      "${tmpdir}/${SAMPLE}.r2.fq.gz" \
      "${tmpdir}/${SAMPLE}.r1.cat.fq.gz" \
      "${tmpdir}/${SAMPLE}.r2.cat.fq.gz"

echo "WGS BAM: ${outdir}/${SAMPLE}.markdup.bam"
