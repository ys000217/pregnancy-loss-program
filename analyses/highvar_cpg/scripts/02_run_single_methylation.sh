#!/bin/bash

# =========================
# Serial methylation processing (LOGIN SAFE)
# ONE SAMPLE AT A TIME
# =========================

BASE="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/figure3/meQTL/prepare_methylation"

MAP="${BASE}/sample_file_map.tsv"
OUT="${BASE}/CpG_long_tmp_serial"
FINAL="${BASE}/CpG_long_table_serial.tsv"

mkdir -p "${OUT}"

echo "===================================="
echo "Start time: $(date)"
echo "Input: ${MAP}"
echo "Output dir: ${OUT}"
echo "===================================="

# =========================
# loop over all samples
# =========================

tail -n +2 "${MAP}" | while read -r sample file; do

    echo "------------------------------------"
    echo "Processing: ${sample}"
    echo "File: ${file}"
    echo "Time: $(date)"
    echo "------------------------------------"

    zcat "${file}" | awk -v sample="${sample}" 'BEGIN{OFS="\t"}
    {
        chr=$1;
        pos=$2;
        cov=$10;
        beta=$11/100;

        if(cov >= 5){
            print chr":"pos, sample, beta
        }
    }' > "${OUT}/${sample}.tsv"

    # 简单保护：避免IO瞬间打爆
    sleep 0.2

done

echo "Merging all samples..."

cat "${OUT}"/*.tsv > "${FINAL}"

echo "===================================="
echo "Done: $(date)"
echo "Final file: ${FINAL}"
echo "===================================="
