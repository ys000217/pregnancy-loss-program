#!/bin/bash
# ================================================================
# PMD Batch Detection — 56-Core Parallel (xargs, no extra deps)
# ================================================================
# Usage:
#   cd /data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/methlybed_modkit
#   bash setup.sh          # one-time: generate samples.txt, check env
#   jsub jsub_pmd_parallel.sh
#
# 56 cores × ~15 min/sample → ceil(655/56)=12 waves → ~3 h wall time
# ================================================================
#JSUB -q normal
#JSUB -n 56
#JSUB -J pmd_detect
#JSUB -o /data/home/4250160/hpc_logs/pmd_detect.out.%J
#JSUB -e /data/home/4250160/hpc_logs/pmd_detect.err.%J
#JSUB -cwd /data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/methlybed_modkit
#JSUB -R span[hosts=1]
# ================================================================

BASE_DIR="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/methlybed_modkit"
SCRIPT="${BASE_DIR}/call_pmd_dnmtools.py"
SAMPLE_LIST="${BASE_DIR}/samples.txt"
LOG_DIR="/data/home/4250160/hpc_logs/pmd"
N_PARALLEL=56

mkdir -p "${LOG_DIR}"

echo "=============================================="
echo "PMD Detection — 56-Core Parallel"
echo "Samples: $(wc -l < "${SAMPLE_LIST}")"
echo "Cores:   ${N_PARALLEL}"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=============================================="

# ── Per-sample worker (called via xargs) ──
run_one() {
    local SAMPLE="$1"
    local SAMPLE_DIR="${BASE_DIR}/${SAMPLE}"
    local BEDMETHYL="${SAMPLE_DIR}/${SAMPLE}.cpg.5mC.bed.gz"
    local OUTPUT="${SAMPLE_DIR}/${SAMPLE}_pmd.bed"
    local LOGFILE="${LOG_DIR}/${SAMPLE}_pmd.log"

    # ---- skip conditions ----
    if [ ! -f "${BEDMETHYL}" ]; then
        echo "[SKIP] ${SAMPLE}: bedMethyl not found"
        echo "SKIP:no_bedmethyl:${SAMPLE}" >> "${LOG_DIR}/summary.txt"
        return 0
    fi

    if [ -f "${OUTPUT}" ]; then
        echo "[SKIP] ${SAMPLE}: output already exists"
        echo "SKIP:exists:${SAMPLE}" >> "${LOG_DIR}/summary.txt"
        return 0
    fi

    # ---- run PMD detection ----
    echo "[START] ${SAMPLE} ($(date '+%H:%M:%S'))"

    if python "${SCRIPT}" "${BEDMETHYL}" -o "${OUTPUT}" -v > "${LOGFILE}" 2>&1; then
        N_PMD=$(grep -vc '^#' "${OUTPUT}" 2>/dev/null || echo "?")
        echo "[OK]    ${SAMPLE}: ${N_PMD} PMDs"
        echo "OK:${N_PMD}:${SAMPLE}" >> "${LOG_DIR}/summary.txt"
        return 0
    else
        echo "[FAIL]  ${SAMPLE}: see ${LOGFILE}"
        echo "FAIL:${SAMPLE}" >> "${LOG_DIR}/summary.txt"
        return 1   # xargs will continue despite this
    fi
}

export -f run_one
export BASE_DIR SCRIPT LOG_DIR

# ── Run: xargs -P launches up to 56 parallel workers ──
# Each worker picks the next sample from the list (FIFO scheduling).
cat "${SAMPLE_LIST}" | \
    xargs -P "${N_PARALLEL}" -I {} bash -c 'run_one "$@"' _ {}

# ── Summary ──
echo ""
echo "=============================================="
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=============================================="
if [ -f "${LOG_DIR}/summary.txt" ]; then
    echo "OK:        $(grep -c '^OK:'     "${LOG_DIR}/summary.txt" 2>/dev/null || echo 0)"
    echo "Failed:    $(grep -c '^FAIL:'   "${LOG_DIR}/summary.txt" 2>/dev/null || echo 0)"
    echo "Skipped:   $(grep -c '^SKIP:'   "${LOG_DIR}/summary.txt" 2>/dev/null || echo 0)"
fi
