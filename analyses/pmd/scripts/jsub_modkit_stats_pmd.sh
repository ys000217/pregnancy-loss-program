#!/bin/bash
#JSUB -q normal
#JSUB -n 28
#JSUB -J modkit_stats_pmd
#JSUB -o /data/home/4250160/hpc_logs/modkit_stats_pmd.%J.out
#JSUB -e /data/home/4250160/hpc_logs/modkit_stats_pmd.%J.err
#JSUB -cwd /data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/figure1/modkit_stats

set -euo pipefail

############################
# 基础路径
############################

BASE_DIR="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/figure1"

METHYL_DIR="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/methlybed_modkit"

ANNO_DIR="${BASE_DIR}/Anno_data"

OUT_BASE="${BASE_DIR}/modkit_stats"

SAMPLE_LIST="${OUT_BASE}/pmd_all_samples_g8910.tsv"

MODKIT="/data/home/4250160/miniconda3/envs/modkit/bin/modkit"

THREADS=20
IO_THREADS=8

############################
# 自检
############################

echo "=========================================="
echo "=== modkit stats : ALL samples PMD ==="
echo "=========================================="
echo "Sample list  : ${SAMPLE_LIST}"
echo "Methyl dir   : ${METHYL_DIR}"
echo "Anno dir     : ${ANNO_DIR}"
echo "Out base     : ${OUT_BASE}"
echo "Threads      : ${THREADS}"
echo "IO threads   : ${IO_THREADS}"
echo "Host         : $(hostname)"
echo "Start time   : $(date)"
echo "=========================================="
echo

[[ -s "${SAMPLE_LIST}" ]] || {
    echo "ERROR: sample list not found or empty:"
    echo "${SAMPLE_LIST}"
    echo "Run submit_modkitstats_pmd.sh first."
    exit 1
}

[[ -x "${MODKIT}" ]] || {
    echo "ERROR: modkit executable not found:"
    echo "${MODKIT}"
    exit 1
}

for GW in 8 9 10; do
    REGIONS="${ANNO_DIR}/control_g${GW}_consensus_pmd_frequency_ge70_refseq.bed"
    [[ -s "${REGIONS}" ]] || {
        echo "ERROR: RefSeq PMD regions missing for g${GW}:"
        echo "${REGIONS}"
        echo "Run submit_modkitstats_pmd.sh first."
        exit 1
    }
done

############################
# 逐样本计算（已有非空结果则跳过）
############################

N_TOTAL=0
N_DONE=0
N_SKIP=0
N_FAIL=0

# 跳过表头；不要用管道，以免子 shell
while IFS=$'\t' read -r SAMPLE GW GROUP1 GROUP4; do
    [[ "${SAMPLE}" == "sample_id" ]] && continue
    [[ -z "${SAMPLE}" ]] && continue

    N_TOTAL=$((N_TOTAL + 1))

    REGIONS="${ANNO_DIR}/control_g${GW}_consensus_pmd_frequency_ge70_refseq.bed"
    IN_BED="${METHYL_DIR}/${SAMPLE}/${SAMPLE}.cpg.5mC.bed.gz"
    OUT_DIR="${OUT_BASE}/${SAMPLE}"
    OUT_TSV="${OUT_DIR}/${SAMPLE}_g${GW}_pmd_stats.tsv"

    mkdir -p "${OUT_DIR}"

    # 已有结果则跳过（兼容旧名：*_RPL/SPL/case/CONTROL_g*_pmd_stats.tsv）
    EXISTING=""
    if [[ -s "${OUT_TSV}" ]]; then
        EXISTING="${OUT_TSV}"
    else
        for f in "${OUT_DIR}/${SAMPLE}_"*"_g${GW}_pmd_stats.tsv"; do
            [[ -e "${f}" ]] || continue
            case "${f}" in
                *nonpmd*) continue ;;
            esac
            if [[ -s "${f}" ]]; then
                EXISTING="${f}"
                break
            fi
        done
    fi

    if [[ -n "${EXISTING}" ]]; then
        echo "[SKIP] ${SAMPLE} g${GW} (exists: ${EXISTING})"
        N_SKIP=$((N_SKIP + 1))
        continue
    fi

    if [[ ! -f "${IN_BED}" ]]; then
        echo "[FAIL] ${SAMPLE} g${GW}: input BED missing: ${IN_BED}"
        N_FAIL=$((N_FAIL + 1))
        continue
    fi

    echo "------------------------------------------"
    echo "[RUN ] ${SAMPLE}  gw=${GW}  Group1=${GROUP1}  Group4=${GROUP4}"
    echo "  regions : ${REGIONS}"
    echo "  input   : ${IN_BED}"
    echo "  output  : ${OUT_TSV}"
    echo "  time    : $(date)"

    if "${MODKIT}" stats \
        --regions "${REGIONS}" \
        --threads "${THREADS}" \
        --io-threads "${IO_THREADS}" \
        --out-table "${OUT_TSV}" \
        "${IN_BED}"
    then
        echo "[DONE] ${SAMPLE} g${GW}"
        N_DONE=$((N_DONE + 1))
    else
        echo "[FAIL] ${SAMPLE} g${GW}: modkit stats failed"
        rm -f "${OUT_TSV}"
        N_FAIL=$((N_FAIL + 1))
    fi

done < "${SAMPLE_LIST}"

############################
# 汇总
############################

echo
echo "=========================================="
echo "ALL DONE"
echo "Total listed : ${N_TOTAL}"
echo "Computed     : ${N_DONE}"
echo "Skipped      : ${N_SKIP}"
echo "Failed       : ${N_FAIL}"
echo "Finish time  : $(date)"
echo "=========================================="

if [[ "${N_FAIL}" -gt 0 ]]; then
    exit 1
fi
