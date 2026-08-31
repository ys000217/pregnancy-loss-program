#!/usr/bin/env bash
# Illumina breakpoint CNVs: DEL/DUP only. Manta if standalone install, else DELLY.
set -euo pipefail
source "$(dirname "$0")/lib_source_config.sh"

SAMPLE="${1:?usage: 05_wgs_sv.sh SAMPLE wgs.bam}"
WGS_BAM="${2:?WGS BAM}"

outdir="${WORKDIR}/wgs_sv"
mkdir -p "${outdir}"
raw_vcf="${outdir}/${SAMPLE}.wgs_sv.raw.vcf.gz"
cnv_vcf="${outdir}/${SAMPLE}.wgs_sv.cnv.vcf.gz"
legacy_vcf="${outdir}/${SAMPLE}.manta.cnv.vcf.gz"

if [[ -s "${cnv_vcf}" ]]; then
  echo "SKIP WGS SV (exists): ${cnv_vcf}"
  exit 0
fi

call_manta() {
  local run_dir="${WORKDIR}/wgs_sv/${SAMPLE}_manta"
  rm -rf "${run_dir}"
  configManta.py \
    --bam "${WGS_BAM}" \
    --referenceFasta "${REF_FASTA}" \
    --runDir "${run_dir}"
  "${run_dir}/runWorkflow.py" -j "${THREADS}"
  cp "${run_dir}/results/variants/diploidSV.vcf.gz" "${raw_vcf}"
  cp "${run_dir}/results/variants/diploidSV.vcf.gz.tbi" "${raw_vcf}.tbi"
}

call_delly() {
  local bcf="${outdir}/${SAMPLE}.delly.bcf"
  local usage
  usage="$(delly 2>&1 || true)"
  # Recent DELLY: `sr` (short-read SV). Older: `call`. `call` on new builds prints
  # "Unrecognized command call" and exit 1.
  if grep -qw sr <<<"${usage}"; then
    echo "DELLY subcommand: sr"
    delly sr -g "${REF_FASTA}" -o "${bcf}" "${WGS_BAM}"
  elif grep -qw call <<<"${usage}"; then
    echo "DELLY subcommand: call"
    delly call -g "${REF_FASTA}" -o "${bcf}" "${WGS_BAM}"
  else
    echo "ERROR delly has neither 'sr' nor 'call'. Help:" >&2
    echo "${usage}" >&2
    exit 1
  fi
  bcftools view -Oz -o "${raw_vcf}" "${bcf}"
  bcftools index -f "${raw_vcf}"
  rm -f "${bcf}"
}

if command -v configManta.py >/dev/null 2>&1; then
  echo "WGS SV caller: Manta"
  call_manta
elif command -v delly >/dev/null 2>&1; then
  echo "WGS SV caller: DELLY (Manta not on PATH)"
  call_delly
else
  echo "ERROR need configManta.py (standalone Manta) or delly in PATH" >&2
  exit 1
fi

bcftools view -i "(SVTYPE=\"DEL\" || SVTYPE=\"DUP\") && ABS(SVLEN)>=${MIN_SVLEN}" \
  "${raw_vcf}" \
  -Oz -o "${cnv_vcf}"
bcftools index -f "${cnv_vcf}"
ln -sf "$(basename "${cnv_vcf}")" "${legacy_vcf}"
ln -sf "$(basename "${cnv_vcf}").tbi" "${legacy_vcf}.tbi"

echo "WGS SV-CNV: ${cnv_vcf}"
