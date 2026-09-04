#!/usr/bin/env bash
# Submit one Jhinno jsub job per sample (same #JSUB style as jsub_abnormal_rpl_g8.sh).
#
#   ONLY=0002C bash scripts/submit_per_sample.sh   # one job first
#   bash scripts/submit_per_sample.sh                # all in MANIFEST
#   FORCE=1 ONLY=0002C bash scripts/submit_per_sample.sh  # rerun even if done/${sample}.done exists
#
# Skip uses done/${sample}.done (written only after merge+annotate finish), not cnv.high.bed
# alone — so a failed AnnotSV / mid-pipeline crash is not treated as complete.
#
# Before first cohort submit: bash scripts/prepare_ref_index.sh
#
# Cores: default 16 to match config.sh THREADS. DMR jobs used 56; WGS bwa+sniffles
# rarely needs that many, and 678 x 56 will stall the queue.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/scripts/lib_source_config.sh"

NCORES="${NCORES:-16}"
QUEUE="${QUEUE:-normal}"
HPC_LOGS="${HPC_LOGS:-${WORKDIR}/logs}"
JOBDIR="${WORKDIR}/job_scripts"
mkdir -p "${HPC_LOGS}" "${JOBDIR}"

MANIFEST="${MANIFEST:-${WORKDIR}/manifest.from_ont.tsv}"
if [[ ! -s "${MANIFEST}" ]]; then
  echo "ERROR missing ${MANIFEST}" >&2
  exit 1
fi
if ! command -v jsub >/dev/null 2>&1; then
  echo "ERROR jsub not on PATH" >&2
  exit 1
fi

# JSUB -J: letters, numbers, underscore only
job_name() {
  echo "cnv_$1" | tr -c 'A-Za-z0-9_' '_' | cut -c1-40
}

write_job() {
  local sample="$1"
  local jn
  jn="$(job_name "${sample}")"
  local job="${JOBDIR}/${sample}.jsub.sh"
  cat > "${job}" <<EOF
#!/bin/bash
#JSUB -q ${QUEUE}
#JSUB -n ${NCORES}
#JSUB -J ${jn}
#JSUB -o ${HPC_LOGS}/${jn}.out.%J
#JSUB -e ${HPC_LOGS}/${jn}.err.%J
#JSUB -cwd ${ROOT}
#JSUB -R span[hosts=1]

set -euo pipefail
echo "host=\$(hostname) sample=${sample} start=\$(date)"
if [[ -f "\${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
  source "\${HOME}/miniconda3/etc/profile.d/conda.sh"
fi
conda activate cnv10x
export THREADS=${NCORES}
bash ${ROOT}/scripts/run_sample.sh ${sample}
echo "host=\$(hostname) sample=${sample} done=\$(date)"
EOF
  chmod +x "${job}"
  echo "${job}"
}

submit_one() {
  local sample="$1"
  if [[ "${FORCE:-0}" != "1" && -f "${WORKDIR}/done/${sample}.done" ]]; then
    echo "SKIP ${sample} (done/${sample}.done exists; FORCE=1 to rerun)"
    return
  fi
  local job
  job="$(write_job "${sample}")"
  if [[ "${DRY_RUN:-0}" == "1" ]]; then
    echo "DRY jsub ${job}"
    return
  fi
  jsub "${job}"
  echo "SUBMITTED ${sample} -> ${job}"
}

n=0
max="${MAX_JOBS:-0}"
only="${ONLY:-}"

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r ont_id rest; do
  [[ -z "${ont_id}" ]] && continue
  if [[ -n "${only}" && "${ont_id}" != "${only}" ]]; then
    continue
  fi
  submit_one "${ont_id}"
  n=$((n + 1))
  if [[ "${max}" -gt 0 && "${n}" -ge "${max}" ]]; then
    echo "stopped after MAX_JOBS=${max}"
    break
  fi
done

echo "logs: ${HPC_LOGS}/cnv_*.out.<jobid>"
echo "job scripts: ${JOBDIR}"
