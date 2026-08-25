#!/bin/bash
TOOLS=/home/administrator/miniconda3/envs/liftover/bin
OUT=/mnt/d/ONT/clinical_649.GRCh38.correct.vcf
echo "=== 最终输出记录数 ==="
grep -vc '^#' $OUT
echo "=== SVTYPE 分布 ==="
$TOOLS/bcftools query -f '%INFO/SVTYPE\n' $OUT 2>/dev/null | sort | uniq -c
echo "=== 样本数 ==="
$TOOLS/bcftools query -l $OUT 2>/dev/null | wc -l
echo "=== 之前负坐标的 SV 是否还在(应为空) ==="
$TOOLS/bcftools query -f '%ID\n' $OUT 2>/dev/null | grep -E '370_Sniffles2.DEL.446S0|0027C_Sniffles2.INS.6S0|Sniffles2.DEL.ABM0' || echo "(已正确剔除)"
echo "=== 抽查 DEL 的 END>POS(有真实长度)前3条 ==="
$TOOLS/bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' $OUT 2>/dev/null | awk -F'\t' '$4=="DEL" && $3>$2' | head -3
