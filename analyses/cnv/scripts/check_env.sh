#!/usr/bin/env bash
# Check tools and (optionally) reference files on the cluster. Run after conda activate.
set -u

fail=0
need() {
  local cmd="$1"
  if command -v "${cmd}" >/dev/null 2>&1; then
    printf "OK   %-12s  %s\n" "${cmd}" "$(command -v "${cmd}")"
  else
    printf "MISS %-12s\n" "${cmd}"
    fail=1
  fi
}

echo "=== required ==="
need python3
need fastp
need bwa-mem2
need samtools
need mosdepth
need sniffles
need cnvpytor
need bcftools
need bedtools

echo
echo "=== WGS SV caller (need one of Manta or DELLY) ==="
wgs_sv_ok=0
if command -v configManta.py >/dev/null 2>&1; then
  printf "OK   %-12s  %s\n" configManta.py "$(command -v configManta.py)"
  wgs_sv_ok=1
fi
if command -v delly >/dev/null 2>&1; then
  printf "OK   %-12s  %s\n" delly "$(command -v delly)"
  wgs_sv_ok=1
fi
if [[ "${wgs_sv_ok}" -eq 0 ]]; then
  echo "MISS need configManta.py (standalone Manta) OR delly"
  fail=1
fi

echo
echo "=== markdup (one of these) ==="
if command -v sambamba >/dev/null 2>&1; then
  printf "OK   %-12s  %s\n" sambamba "$(command -v sambamba)"
elif samtools markdup -h >/dev/null 2>&1; then
  echo "OK   samtools markdup"
else
  echo "MISS sambamba and samtools markdup"
  fail=1
fi

echo
echo "=== optional ==="
if command -v spectre >/dev/null 2>&1; then
  printf "OK   %-12s  %s\n" spectre "$(command -v spectre)"
else
  echo "skip spectre (ONT depth will use CNVpytor only)"
fi
if command -v AnnotSV >/dev/null 2>&1; then
  printf "OK   %-12s  %s\n" AnnotSV "$(command -v AnnotSV)"
else
  echo "skip AnnotSV (HIGH BED will be unannotated)"
fi

echo
echo "=== versions ==="
python3 --version 2>/dev/null || true
fastp --version 2>/dev/null | head -1 || true
bwa-mem2 version 2>/dev/null | head -1 || true
samtools --version 2>/dev/null | head -1 || true
mosdepth --version 2>/dev/null || true
sniffles --version 2>/dev/null || true
bcftools --version 2>/dev/null | head -1 || true
bedtools --version 2>/dev/null || true

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
if [[ -f "${ROOT}/config.sh" ]]; then
  set +e
  # Strip CR in case this file was copied from Windows
  # shellcheck source=/dev/null
  source /dev/stdin <<< "$(tr -d '\r' < "${ROOT}/config.sh")"
  set -e
  echo
  echo "=== paths from config.sh ==="
  for p in REF_FASTA TR_BED ONT_BAM_ROOT NGS_RAWDATA WORKDIR; do
    eval "v=\${$p-}"
    if [[ -z "${v}" ]]; then
      echo "UNSET ${p}"
    elif [[ -e "${v}" ]]; then
      echo "OK    ${p}=${v}"
    else
      echo "MISS  ${p}=${v}"
    fi
  done
fi

echo
if [[ "${fail}" -ne 0 ]]; then
  echo "REQUIRED TOOLS MISSING. Create env with:"
  echo "  conda env create -f env/environment.yml && conda activate cnv10x"
  exit 1
fi
echo "Required tools are present. Next: lock REF_FASTA against an ONT BAM header, then build manifest.tsv."
