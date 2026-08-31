#!/usr/bin/env bash
# Discover ONT BAM location, reference build, and suggest config.sh edits.
# Usage: bash scripts/discover_paths.sh [sample_id]
# Example: bash scripts/discover_paths.sh 0002C
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/scripts/lib_source_config.sh"

SAMPLE="${1:-0002C}"

echo "=== search ONT merged BAM for sample ${SAMPLE} ==="
candidates=(
  "${ONT_BAM_ROOT}/${SAMPLE}/${SAMPLE}.merged.bam"
  "${SHARE_ROOT:-/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang}/ONT_processed_data/modbam_Dorodo/${SAMPLE}/${SAMPLE}.merged.bam"
)
found_bam=""
for p in "${candidates[@]}"; do
  if [[ -s "${p}" ]]; then
    found_bam="${p}"
    echo "FOUND ${p}"
    break
  else
    echo "miss  ${p}"
  fi
done

if [[ -z "${found_bam}" ]]; then
  echo
  echo "Try manual search:"
  echo "  find /data/jhinno/appform/data/share/5250028 -name '${SAMPLE}.merged.bam' 2>/dev/null | head"
  echo "  find ~/${PROJECT_ID} -name '${SAMPLE}.merged.bam' 2>/dev/null | head"
  exit 1
fi

ont_root="$(dirname "$(dirname "${found_bam}")")"
echo
echo "Suggested: export ONT_BAM_ROOT=${ont_root}"

echo
echo "=== BAM header (@SQ) ==="
samtools view -H "${found_bam}" | grep '^@SQ' | head -5

echo
echo "=== reference path in header (if any) ==="
samtools view -H "${found_bam}" | grep -E '^@SQ|^@PG' | head -10

echo
echo "=== config.sh paths ==="
for p in REF_FASTA TR_BED WORKDIR NGS_RAWDATA; do
  eval "v=\${$p-}"
  if [[ -e "${v}" ]]; then echo "OK   ${p}=${v}"; else echo "MISS ${p}=${v}"; fi
done

echo
echo "=== next steps ==="
echo "1. Put GRCh38.fa at REF_FASTA (must match @SQ names above)."
echo "2. Put tandem-repeat BED at TR_BED (Sniffles2; e.g. human_GRCh38_TR.bed)."
echo "3. mkdir -p ${WORKDIR}"
echo "4. Edit config.sh: ONT_BAM_ROOT=${ont_root}"
echo "5. python3 scripts/00_scan_fastq.py --ngs-root \"${NGS_RAWDATA}\" -o ${MANIFEST}"
