#!/bin/bash
echo "=== 抬升断点中 第2/3列 为负或0 的数量 ==="
awk -F'\t' '($2<0 || $3<0 || $2==$3)' /tmp/grch38_bp.bed | wc -l
echo "=== 示例(负坐标或异常) ==="
awk -F'\t' '($2<0 || $3<0 || $2==$3)' /tmp/grch38_bp.bed | head -8
echo "=== 第3列(end)<=0 的记录 ==="
awk -F'\t' '$3<=0' /tmp/grch38_bp.bed | head -8
echo "=== 输出 VCF 中 POS 为负的记录 ==="
grep -v '^#' /mnt/d/ONT/clinical_649.GRCh38.correct.vcf | awk -F'\t' '$2<0 || $2 ~ /^-/' | head -5 | cut -c1-160
