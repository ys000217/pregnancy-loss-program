#!/usr/bin/env bash
# ONT breakpoint CNVs: Sniffles2, then keep DEL/DUP only.
set -euo pipefail
source "$(dirname "$0")/lib_source_config.sh"

SAMPLE="${1:?usage: 03_ont_sv.sh SAMPLE ont.bam}"
ONT_BAM="${2:?ONT BAM}"

outdir="${WORKDIR}/ont_sv"
mkdir -p "${outdir}"

# Sniffles2 refuses to overwrite an existing VCF.
if [[ -s "${outdir}/${SAMPLE}.sniffles.cnv.vcf.gz" ]]; then
  echo "SKIP ONT SV (exists): ${outdir}/${SAMPLE}.sniffles.cnv.vcf.gz"
  exit 0
fi

sniffles_args=(
  --input "${ONT_BAM}"
  --vcf "${outdir}/${SAMPLE}.sniffles.raw.vcf.gz"
  --snf "${outdir}/${SAMPLE}.snf"
  --reference "${REF_FASTA}"
  --threads "${THREADS}"
  --minsvlen "${MIN_SVLEN}"
  --output-rnames
)
if [[ -n "${TR_BED}" && -s "${TR_BED}" ]]; then
  sniffles_args+=(--tandem-repeats "${TR_BED}")
else
  echo "WARN TR_BED unset or missing; Sniffles2 runs without --tandem-repeats (more repeat FPs)" >&2
fi

sniffles "${sniffles_args[@]}"

bcftools view -i "QUAL>=${MIN_QUAL_SV} && (SVTYPE=\"DEL\" || SVTYPE=\"DUP\") && ABS(SVLEN)>=${MIN_SVLEN}" \
  "${outdir}/${SAMPLE}.sniffles.raw.vcf.gz" \
  -Oz -o "${outdir}/${SAMPLE}.sniffles.cnv.vcf.gz"
bcftools index -f "${outdir}/${SAMPLE}.sniffles.cnv.vcf.gz"

echo "ONT SV-CNV: ${outdir}/${SAMPLE}.sniffles.cnv.vcf.gz"
