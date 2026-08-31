#!/usr/bin/env bash
# Run the full 10x paired CNV pipeline for one sample in MANIFEST (manifest.from_ont.tsv).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "${ROOT}/scripts/lib_source_config.sh"

SAMPLE="${1:?usage: run_sample.sh SAMPLE}"
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
            bam = row.get("ont_bam") or ""
            print("ont_id=" + shlex.quote(row["ont_id"]))
            print("wgs_r1=" + shlex.quote(row["wgs_r1"]))
            print("wgs_r2=" + shlex.quote(row["wgs_r2"]))
            print("sex=" + shlex.quote((row.get("sex") or "").strip()))
            print("ONT_BAM=" + shlex.quote(bam))
            sys.exit(0)
sys.exit(2)
PY
)" || { echo "ERROR ${SAMPLE} not in ${MANIFEST}"; exit 1; }

if [[ ! -s "${ONT_BAM}" ]]; then
  echo "ERROR ONT BAM missing: ${ONT_BAM}" >&2
  exit 1
fi

WGS_BAM="${WORKDIR}/wgs_bam/${ont_id}.markdup.bam"

mkdir -p "${WORKDIR}"/{qc,wgs_bam,ont_sv,ont_cnv,wgs_sv,wgs_cnv,merged,annot,tmp}

echo "=== ${ont_id}  ONT=${ONT_BAM}"
echo "=== ${ont_id}  R1=${wgs_r1}"

if [[ ! -s "${WGS_BAM}" ]]; then
  bash "${SCRIPTS}/01_wgs_align.sh" "${ont_id}" "${wgs_r1}" "${wgs_r2}"
fi

bash "${SCRIPTS}/02_qc_coverage.sh" "${ont_id}" "${ONT_BAM}" "${WGS_BAM}"
bash "${SCRIPTS}/03_ont_sv.sh" "${ont_id}" "${ONT_BAM}"
bash "${SCRIPTS}/04_ont_cnv_depth.sh" "${ont_id}" "${ONT_BAM}"
bash "${SCRIPTS}/05_wgs_sv.sh" "${ont_id}" "${WGS_BAM}"
bash "${SCRIPTS}/06_wgs_cnv_depth.sh" "${ont_id}" "${WGS_BAM}"

merge_args=(
  --sample "${ont_id}"
  --ont-sv "${WORKDIR}/ont_sv/${ont_id}.sniffles.cnv.vcf.gz"
  --wgs-sv "${WORKDIR}/wgs_sv/${ont_id}.wgs_sv.cnv.vcf.gz"
  --ont-cnvpytor "${WORKDIR}/ont_cnv/${ont_id}.cnvpytor.${BIN_PRIMARY}.tsv"
  --wgs-cnvpytor "${WORKDIR}/wgs_cnv/${ont_id}.cnvpytor.${BIN_PRIMARY}.tsv"
  --ro "${MERGE_RO:-0.50}"
  --min-depth "${MIN_DEPTH_CNV:-100000}"
  --max-event "${MAX_CNV_EVENT:-10000000}"
  --mask-frac "${MERGE_MASK_FRAC:-0.50}"
  --sex "${sex:-}"
  -o "${WORKDIR}/merged"
)
if [[ -n "${HARD_MASK_BED:-}" && -s "${HARD_MASK_BED}" ]]; then
  merge_args+=(--hard-mask "${HARD_MASK_BED}")
else
  echo "WARN HARD_MASK_BED unset/missing; LARGE_HIGH will not apply genomic hard mask" >&2
fi
spectre_hit="$(find "${WORKDIR}/ont_cnv/spectre" -name "${ont_id}*.bed" -print -quit 2>/dev/null || true)"
if [[ -n "${spectre_hit}" ]]; then
  merge_args+=(--ont-spectre "${spectre_hit}")
fi

python3 "${SCRIPTS}/07_merge_paired.py" "${merge_args[@]}"
bash "${SCRIPTS}/08_annotate.sh" "${ont_id}"

echo "=== done ${ont_id}"
echo "LARGE_HIGH: ${WORKDIR}/merged/${ont_id}.cnv.high.bed"
echo "SHARED_SV:  ${WORKDIR}/merged/${ont_id}.cnv.shared_sv.bed"
