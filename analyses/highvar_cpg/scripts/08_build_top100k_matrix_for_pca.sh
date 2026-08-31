#!/bin/bash
#JSUB -J cpg_top100k_matrix
#JSUB -q normal
#JSUB -n 32
#JSUB -o cpg_top100k_matrix.%J.out
#JSUB -e cpg_top100k_matrix.%J.err


BASE="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/figure3/meQTL/prepare_methylation"

DIR="${BASE}/CpG_long_tmp_serial"
VAR="${BASE}/CpG_variance.tsv"

TOP_LIST="${BASE}/CpG_top100k.list"
OUT="${BASE}/CpG_top100k_matrix.tsv"

TMPDIR="${BASE}/tmp_top100k"
mkdir -p ${TMPDIR}

echo "===================================="
echo "Start: $(date)"
echo "===================================="

# =========================
# Step 1: 选top 100k CpG（按variance）
# =========================
echo "Selecting top 100k CpGs..."

sort -k3,3nr ${VAR} | head -n 100000 | cut -f1 > ${TOP_LIST}

echo "Top CpG list saved: ${TOP_LIST}"

# =========================
# Step 2: 样本顺序（固定！）
# =========================
echo "Collecting sample list..."

ls ${DIR}/*.tsv | sed 's#.*/##; s/.tsv//' | sort > ${TMPDIR}/sample.list

# =========================
# Step 3: 初始化 matrix（CpG列）
# =========================
echo "Initializing matrix..."

cp ${TOP_LIST} ${TMPDIR}/matrix.body
sort -k1,1 -o ${TMPDIR}/matrix.body ${TMPDIR}/matrix.body

# =========================
# Step 4: 填充矩阵
# =========================
echo "Building matrix..."

while read sample
do
    f=${DIR}/${sample}.tsv
    echo "Processing ${sample} ... $(date)"

    TMP_SAMPLE=${TMPDIR}/${sample}.tmp

    awk 'NR==FNR {keep[$1]; next}
         ($1 in keep) {print $1"\t"$3}' ${TOP_LIST} ${f} \
    | sort -k1,1 > ${TMP_SAMPLE}

    join -a1 -e NA -o auto ${TMPDIR}/matrix.body ${TMP_SAMPLE} > ${TMPDIR}/matrix.tmp

    mv ${TMPDIR}/matrix.tmp ${TMPDIR}/matrix.body
    rm -f ${TMP_SAMPLE}

done < ${TMPDIR}/sample.list

# =========================
# Step 5: 添加 header
# =========================
echo "Adding header..."

echo -ne "CpG_ID" > ${TMPDIR}/header
while read s
do
    echo -ne "\t${s}" >> ${TMPDIR}/header
done < ${TMPDIR}/sample.list
echo "" >> ${TMPDIR}/header

cat ${TMPDIR}/header ${TMPDIR}/matrix.body > ${OUT}

echo "===================================="
echo "Done: $(date)"
echo "Matrix: ${OUT}"
echo "===================================="
