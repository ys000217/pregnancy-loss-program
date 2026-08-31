#!/bin/bash
#JSUB -q normal
#JSUB -n 8
#JSUB -J ext_648
#JSUB -o ext_648.%J.out
#JSUB -e ext_648.%J.err

set -euo pipefail

# --- 1. 路径定义 ---
BASE="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/figure3/meQTL/prepare_methylation"
# 使用包含 648 个样本的新列表
LIST="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/figure3/gwas/prepare_covariates/keep_samples_648.list"
TSV_DIR="${BASE}/CpG_long_tmp_serial"
TEMPLATE="${BASE}/CpG_final.list"

# 输出文件名
OUT_BED="./648samples_methymatrix.bed"
# 临时工作目录 (以时间戳命名，防止任务重叠)
TMP_WORK="./tmp_extract_648_$(date +%s)"

mkdir -p "${TMP_WORK}"

echo "===================================================="
echo "开始时间: $(date)"
echo "目标样本数: $(wc -l < "$LIST")"
echo "===================================================="

# --- 2. 准备位点基准 ---
echo "Step 1: 提取 3,059,870 个位点 ID 模板..."
cut -f1 "${TEMPLATE}" > "${TMP_WORK}/index.txt"

# --- 3. 提取样本数据 ---
echo "Step 2: 循环提取各样本 Beta 值 (查字典模式)..."
# 记录真正提取成功的样本 ID，确保后续表头对齐
TRUE_LIST="${TMP_WORK}/true_samples.list"
> "$TRUE_LIST"

while read id; do
    # 检查原始 .tsv 文件是否存在
    if [ ! -f "${TSV_DIR}/${id}.tsv" ]; then
        echo "警告: 找不到样本 ${id}.tsv，跳过该样本"
        continue
    fi

    # 查字典：以 index.txt 为准，查不到补 NA
    awk -F'\t' '
        NR==FNR { val[$1]=$3; next }
        {
            if ($1 in val) print val[$1];
            else print "NA"
        }' "${TSV_DIR}/${id}.tsv" "${TMP_WORK}/index.txt" > "${TMP_WORK}/${id}.col"
    
    echo "$id" >> "$TRUE_LIST"
done < "${LIST}"

# --- 4. 构建坐标列 ---
echo "Step 3: 生成坐标信息 (chr, start, end, phenotype_id)..."
awk -F'\t' 'BEGIN{OFS="\t"}{
    split($1, a, ":");
    chr=a[1]; pos=a[2];
    # 转换 NC 编号为 chr 格式 (针对人类 GRCh38)
    gsub("NC_0000","chr",chr); gsub("\\..*","",chr);
    if(chr=="chr23") chr="chrX"; if(chr=="chr24") chr="chrY";
    # 移除 chr01 -> chr1 等多余的零
    sub("chr0", "chr", chr);
    print chr, pos-1, pos, $1
}' "${TEMPLATE}" > "${TMP_WORK}/coords.txt"

# --- 5. 合并矩阵 ---
echo "Step 4: 合并表头与所有样本列..."

# 1. 生成表头 (只包含真正提取成功的样本)
echo -ne "#chr\tstart\tend\tphenotype_id" > "${TMP_WORK}/header.txt"
while read id; do
    echo -ne "\t${id}"
done < "$TRUE_LIST" >> "${TMP_WORK}/header.txt"
echo "" >> "${TMP_WORK}/header.txt"

# 2. 生成所有待合并文件的路径列表 (规避 Argument list too long)
# 这一步非常重要：确保坐标列在前，后续样本列按顺序跟随
COL_FILES="${TMP_WORK}/col_files.path"
echo "${TMP_WORK}/coords.txt" > "$COL_FILES"
awk '{print tmp "/" $1 ".col"}' tmp="${TMP_WORK}" "$TRUE_LIST" >> "$COL_FILES"

# 3. 使用 xargs 驱动 paste 进横向合并 (最安全且快的方法)
# xargs 会处理过长的参数列表
paste $(cat "$COL_FILES") > "${TMP_WORK}/body.txt"

# 4. 最终组装
cat "${TMP_WORK}/header.txt" "${TMP_WORK}/body.txt" > "${OUT_BED}"

echo "===================================================="
echo "任务完成: $(date)"
echo "生成文件: ${OUT_BED}"
echo "最终行数: $(wc -l < "${OUT_BED}")"
echo "最终列数: $(head -n 1 "${OUT_BED}" | awk '{print NF}')"
echo "===================================================="

# 清理临时文件 (如果需要调试可以注释掉下面这行)
rm -rf "${TMP_WORK}"
