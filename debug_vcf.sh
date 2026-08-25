#!/bin/bash
TOOLS=/home/administrator/miniconda3/envs/liftover/bin
OUT=/mnt/d/ONT/clinical_649.GRCh38.correct.vcf
echo "=== 数据记录数 ==="
grep -vc '^#' $OUT
echo "=== bcftools query SVTYPE 计数(带错误输出) ==="
$TOOLS/bcftools query -f '%INFO/SVTYPE\n' $OUT 2>&1 | head -5
echo "--- query 输出总行数 ---"
$TOOLS/bcftools query -f '%INFO/SVTYPE\n' $OUT 2>/dev/null | wc -l
echo "=== 前3条原始数据记录(截断) ==="
grep -v '^#' $OUT | head -3 | cut -c1-200
echo "=== 中间随机一条原始记录(第30000条) ==="
grep -v '^#' $OUT | sed -n '30000p' | cut -c1-300
