#!/usr/bin/env bash

set -euo pipefail

############################
# 全样本 PMD modkit stats：只提交 1 个 jsub
# 样本来源：clinical_649.tsv 中 Group3∈{8,9,10} 的全部样本
# （不区分 control / case / abnormal）
############################

BASE_DIR="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/figure1"

METHYL_DIR="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/methlybed_modkit"

ANNO_DIR="${BASE_DIR}/Anno_data"

MODKIT_STATS_DIR="${BASE_DIR}/modkit_stats"

CLINICAL="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/figure3/meQTL/prepare_methylation/clinical_649.tsv"

CHR_MAP="${ANNO_DIR}/chr_name_map_hg38.clean.tsv"

JOB_SCRIPT="${MODKIT_STATS_DIR}/jsub_modkit_stats_pmd.sh"

SAMPLE_LIST="${MODKIT_STATS_DIR}/pmd_all_samples_g8910.tsv"

############################
# 检查
############################

if [[ ! -f "${CLINICAL}" ]]; then
    echo "ERROR: clinical table not found:"
    echo "${CLINICAL}"
    exit 1
fi

if [[ ! -f "${JOB_SCRIPT}" ]]; then
    echo "ERROR: job script not found:"
    echo "${JOB_SCRIPT}"
    echo "Copy jsub_modkit_stats_pmd.sh to ${MODKIT_STATS_DIR}/ first."
    exit 1
fi

if [[ ! -f "${CHR_MAP}" ]]; then
    echo "ERROR: chromosome name map not found:"
    echo "${CHR_MAP}"
    exit 1
fi

############################
# UCSC → RefSeq：G8 / G9 / G10 三套 PMD 区域
############################

convert_pmd_bed() {
    local gw="$1"
    local ucsc="${ANNO_DIR}/control_g${gw}_consensus_pmd_frequency_ge70.bed"
    local refseq="${ANNO_DIR}/control_g${gw}_consensus_pmd_frequency_ge70_refseq.bed"

    if [[ ! -f "${ucsc}" ]]; then
        echo "ERROR: PMD regions BED not found: ${ucsc}"
        exit 1
    fi

    if [[ -s "${refseq}" ]]; then
        echo "RefSeq PMD g${gw} already exists: ${refseq} ($(wc -l < "${refseq}") intervals)"
        return 0
    fi

    echo "Converting PMD BED g${gw} → RefSeq..."
    local tmp="${refseq}.tmp.$$"
    awk -F'\t' '
        NR == FNR { map[$1] = $2; next }
        /^#/ || NF < 3 { next }
        ($1 in map) { print map[$1], $2, $3; next }
        { print "WARNING: unmapped chrom: " $1 > "/dev/stderr" }
    ' OFS='\t' "${CHR_MAP}" "${ucsc}" > "${tmp}"

    if [[ ! -s "${tmp}" ]]; then
        echo "ERROR: converted RefSeq PMD BED is empty for g${gw}"
        rm -f "${tmp}"
        exit 1
    fi

    mv "${tmp}" "${refseq}"
    echo "  wrote ${refseq} ($(wc -l < "${refseq}") intervals)"
}

for GW in 8 9 10; do
    convert_pmd_bed "${GW}"
done
echo

############################
# 从临床表构建全样本名单：sample_id  gw  Group1  Group4
############################

echo "Building sample list from clinical (Group3 in 8/9/10, all Group1/Group4)..."

awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        gsub(/\r/, "", $i)
        h[tolower($i)] = i
    }
    if (!("sample_id" in h) || !("group1" in h) || !("group3" in h) || !("group4" in h)) {
        print "ERROR: clinical header missing Sample_ID/Group1/Group3/Group4" > "/dev/stderr"
        exit 2
    }
    print "sample_id", "gw", "Group1", "Group4"
    next
}
{
    for (i = 1; i <= NF; i++) gsub(/\r/, "", $i)
    sid = $h["sample_id"]
    week = $h["group3"] + 0
    if (week != 8 && week != 9 && week != 10) next
    print sid, week, $h["group1"], $h["group4"]
}
' OFS='\t' "${CLINICAL}" > "${SAMPLE_LIST}"

N_ALL=$(($(wc -l < "${SAMPLE_LIST}") - 1))
N8=$(awk -F'\t' 'NR>1 && $2==8 {c++} END{print c+0}' "${SAMPLE_LIST}")
N9=$(awk -F'\t' 'NR>1 && $2==9 {c++} END{print c+0}' "${SAMPLE_LIST}")
N10=$(awk -F'\t' 'NR>1 && $2==10 {c++} END{print c+0}' "${SAMPLE_LIST}")

echo "Sample list : ${SAMPLE_LIST}"
echo "Total       : ${N_ALL}  (g8=${N8}, g9=${N9}, g10=${N10})"
echo

if [[ "${N_ALL}" -le 0 ]]; then
    echo "ERROR: no samples in list"
    exit 1
fi

############################
# 提交唯一一个 jsub
############################

echo "=========================================="
echo "=== all-sample PMD modkit submission ==="
echo "=========================================="
echo "Clinical     : ${CLINICAL}"
echo "Sample list  : ${SAMPLE_LIST}"
echo "Job script   : ${JOB_SCRIPT}"
echo "Submitting ONE job..."
echo "=========================================="
echo

jsub < "${JOB_SCRIPT}"

echo
echo "Submitted. Monitor logs under /data/home/4250160/hpc_logs/modkit_stats_pmd.*"
