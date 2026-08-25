#!/bin/bash
TOOLS=/home/administrator/miniconda3/envs/liftover/bin
OUT=/mnt/d/ONT/clinical_649.GRCh38.correct.vcf
echo "=== 之前负坐标的 SV 现在的坐标(应为正值) ==="
$TOOLS/bcftools query -f '%ID\t%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' $OUT 2>/dev/null | \
  grep -E '370_Sniffles2.DEL.446S0|0027C_Sniffles2.INS.6S0|Sniffles2.DEL.ABM0'
