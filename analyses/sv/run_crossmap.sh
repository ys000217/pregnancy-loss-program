#!/bin/bash
TOOLS=/home/administrator/miniconda3/envs/liftover/bin
echo "=== CrossMap vcf 开始 $(date) ==="
$TOOLS/CrossMap vcf \
  /mnt/e/genotype_data/liftover/CN1_to_GRCh38.chain \
  /mnt/d/ONT/clinical_649.vcf \
  /mnt/e/genotype_data/liftover/GRCh38_chr.fasta \
  /mnt/d/ONT/clinical_649.GRCh38.vcf
echo "=== CrossMap 结束 $(date), exit=$? ==="
echo "输出记录数:"
grep -vc '^#' /mnt/d/ONT/clinical_649.GRCh38.vcf 2>/dev/null || echo "无输出文件"
echo "unmap 行数:"
wc -l /mnt/d/ONT/clinical_649.GRCh38.vcf.unmap 2>/dev/null || echo "无 unmap 文件"
