#!/bin/bash
TOOLS=/home/administrator/miniconda3/envs/liftover/bin
echo "=== 输出头部 contig(前2条) ==="
grep '^##contig' /mnt/d/ONT/clinical_649.GRCh38.correct.vcf | head -2
echo ""
echo "=== 输出中前5个 DEL(END>POS, 即有真实长度) ==="
$TOOLS/bcftools query -f '%ID\t%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' \
  /mnt/d/ONT/clinical_649.GRCh38.correct.vcf 2>/dev/null | \
  awk -F'\t' '$5=="DEL" && $4>$3' | head -5
echo ""
echo "=== 各 SVTYPE 数量 ==="
$TOOLS/bcftools query -f '%INFO/SVTYPE\n' /mnt/d/ONT/clinical_649.GRCh38.correct.vcf 2>/dev/null | sort | uniq -c
