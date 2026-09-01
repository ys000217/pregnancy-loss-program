#!/usr/bin/env bash
# Re-run merge + annotate only (callsets already on disk). Useful after rule changes.
# Usage: bash scripts/remake_merge.sh SAMPLE
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "${ROOT}/scripts/lib_source_config.sh"

SAMPLE="${1:?usage: remake_merge.sh SAMPLE}"
SCRIPTS="${ROOT}/scripts"

if [[ ! -s "${MANIFEST}" ]]; then
  echo "ERROR manifest not found: ${MANIFEST}" >&2
  exit 1
fi

eval "$(python3 - "${SAMPLE}" "${MANIFEST}" <<'PY'
import csv, shlex, sys
sid, path = sys.argv[1], sys.argv[2]
with open(path, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row.get("ont_id") == sid:
            print("ont_id=" + shlex.quote(row["ont_id"]))
            print("sex=" + shlex.quote((row.get("sex") or "").strip()))
            sys.exit(0)
sys.exit(2)
PY
)" || { echo "ERROR ${SAMPLE} not in ${MANIFEST}"; exit 1; }

# Resolve hard mask (required for usable LARGE_HIGH)
if [[ -z "${HARD_MASK_BED:-}" || ! -s "${HARD_MASK_BED}" ]]; then
  for _mask_cand in \
    "${ROOT}/ref/hard_mask.grch38.refseq.bed" \
    "${_CNV_ROOT:-$ROOT}/ref/hard_mask.grch38.refseq.bed" \
    "${SHARE_ROOT}/cnv/ref/hard_mask.grch38.refseq.bed" \
    "${WORKDIR}/ref/hard_mask.grch38.refseq.bed"
  do
    if [[ -s "${_mask_cand}" ]]; then
      HARD_MASK_BED="${_mask_cand}"
      break
    fi
  done
  unset _mask_cand 2>/dev/null || true
fi
if [[ -z "${HARD_MASK_BED:-}" || ! -s "${HARD_MASK_BED}" ]]; then
  echo "ERROR hard mask missing — copy the BED before remake:" >&2
  echo "  mkdir -p ${ROOT}/ref" >&2
  echo "  # from Windows: D:\\ONT\\analyses\\cnv\\ref\\hard_mask.grch38.refseq.bed" >&2
  echo "  # to:           ${ROOT}/ref/hard_mask.grch38.refseq.bed" >&2
  exit 1
fi
export HARD_MASK_BED
echo "HARD_MASK_BED=${HARD_MASK_BED}"

merge_args=(
  --sample "${ont_id}"
  --ont-sv "${WORKDIR}/ont_sv/${ont_id}.sniffles.cnv.vcf.gz"
  --wgs-sv "${WORKDIR}/wgs_sv/${ont_id}.wgs_sv.cnv.vcf.gz"
  --ont-cnvpytor "${WORKDIR}/ont_cnv/${ont_id}.cnvpytor.${BIN_PRIMARY}.tsv"
  --wgs-cnvpytor "${WORKDIR}/wgs_cnv/${ont_id}.cnvpytor.${BIN_PRIMARY}.tsv"
  --ont-cnvpytor-large "${WORKDIR}/ont_cnv/${ont_id}.cnvpytor.${BIN_LARGE}.tsv"
  --wgs-cnvpytor-large "${WORKDIR}/wgs_cnv/${ont_id}.cnvpytor.${BIN_LARGE}.tsv"
  --ro "${MERGE_RO:-0.50}"
  --min-depth "${MIN_DEPTH_CNV:-100000}"
  --max-event "${MAX_CNV_EVENT:-10000000}"
  --mask-frac "${MERGE_MASK_FRAC:-0.50}"
  --cnvpytor-q0-max "${CNVPYTOR_Q0_MAX:-0.5}"
  --cnvpytor-pn-max "${CNVPYTOR_PN_MAX:-0.5}"
  --cnvpytor-eval-max "${CNVPYTOR_EVAL_MAX:-1e-4}"
  --sex "${sex:-}"
  --require-hard-mask
  --hard-mask "${HARD_MASK_BED}"
  -o "${WORKDIR}/merged"
)
if [[ "${CNVPYTOR_QC:-1}" != "1" ]]; then
  merge_args+=(--no-cnvpytor-qc)
fi
if [[ "${KEEP_SEX_CHROM:-0}" == "1" ]]; then
  merge_args+=(--keep-sex-chrom)
fi

python3 "${SCRIPTS}/07_merge_paired.py" "${merge_args[@]}"
bash "${SCRIPTS}/08_annotate.sh" "${ont_id}"
echo "remade ${ont_id}: LARGE_HIGH -> ${WORKDIR}/merged/${ont_id}.cnv.high.bed"
