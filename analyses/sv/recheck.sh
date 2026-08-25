#!/bin/bash
OUT=/mnt/d/ONT/clinical_649.GRCh38.correct.vcf
echo "=== 输出中 POS<0 或 END<0 的记录数 ==="
grep -v '^#' $OUT | awk -F'\t' '$2<0 || $8 ~ /END=-/' | wc -l
echo "=== 输出中 POS<0 的记录(前3条) ==="
grep -v '^#' $OUT | awk -F'\t' '$2<0' | head -3 | cut -c1-120
echo "=== bcftools 解析是否通过(统计 SVTYPE) ==="
/home/administrator/miniconda3/envs/liftover/bin/bcftools query -f '%INFO/SVTYPE\n' $OUT 2>&1 | sort | uniq -c | head
echo "=== query 成功输出的行数 ==="
/home/administrator/miniconda3/envs/liftover/bin/bcftools query -f '%INFO/SVTYPE\n' $OUT 2>/dev/null | wc -l
