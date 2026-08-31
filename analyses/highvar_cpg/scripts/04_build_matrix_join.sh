#!/bin/bash
#JSUB -q normal
#JSUB -n 32
#JSUB -J build_cpg_matrix
#JSUB -o build_matrix.%J.out
#JSUB -e build_matrix.%J.err

set -euo pipefail

BASE="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/figure3/meQTL/prepare_methylation"

DIR="${BASE}/CpG_long_tmp_serial"
CPG_LIST="${BASE}/CpG_95pct.list"
OUT="${BASE}/CpG_matrix.tsv"

TMPDIR="${BASE}/tmp_matrix"
mkdir -p ${TMPDIR}

echo "===================================="
echo "Start: $(date)"
echo "CpG list: ${CPG_LIST}"
echo "Input dir: ${DIR}"
echo "Output: ${OUT}"
echo "===================================="

# =========================
# Step 1: 初始化 matrix
# =========================
echo "Initializing matrix..."

cp ${CPG_LIST} ${OUT}

# 确保排序（join 必须）
sort -k1,1 -o ${OUT} ${OUT}

# =========================
# Step 2: 逐样本加入
# =========================
for f in ${DIR}/*.tsv
do
    sample=$(basename $f .tsv)

    echo "------------------------------------"
    echo "Processing sample: ${sample}"
    echo "Time: $(date)"
    echo "------------------------------------"

    TMP_SAMPLE="${TMPDIR}/${sample}.tsv"

    # 过滤 + 提取 beta + 排序
    awk 'NR==FNR {keep[$1]; next}
         ($1 in keep) {print $1"\t"$3}' ${CPG_LIST} $f \
    | sort -k1,1 > ${TMP_SAMPLE}

    # join 到矩阵
    join -a1 -e NA -o auto ${OUT} ${TMP_SAMPLE} > ${OUT}.tmp

    mv ${OUT}.tmp ${OUT}
    rm -f ${TMP_SAMPLE}

done

echo "===================================="
echo "Done: $(date)"
echo "Matrix file: ${OUT}"
echo "===================================="
