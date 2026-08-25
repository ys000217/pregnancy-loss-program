#!/bin/bash
set -e
TOOLS=/home/administrator/miniconda3/envs/liftover/bin
echo "=== Step1: bcftools 筛选 649 临床样本 ==="
$TOOLS/bcftools view -S /mnt/d/ONT/clinical_649.samples.txt \
  -O v -o /mnt/d/ONT/clinical_649.vcf \
  /mnt/e/genotype_data/S957.1kGPont.merged.649clin_1027kgp.vcf
echo "样本数(应=649):"
$TOOLS/bcftools query -l /mnt/d/ONT/clinical_649.vcf | wc -l
echo "记录数:"
grep -vc '^#' /mnt/d/ONT/clinical_649.vcf

echo ""
echo "=== Step2: CrossMap vcf 做 CN1->GRCh38 liftover ==="
$TOOLS/CrossMap vcf \
  /mnt/e/genotype_data/liftover/CN1_to_GRCh38.chain \
  /mnt/d/ONT/clinical_649.vcf \
  /mnt/e/genotype_data/liftover/GRCh38_chr.fasta \
  /mnt/d/ONT/clinical_649.GRCh38.vcf 2>&1 | tail -5
echo "liftover 完成"
echo "输出记录数:"
grep -vc '^#' /mnt/d/ONT/clinical_649.GRCh38.vcf
echo "unmap 记录数(若有):"
ls -la /mnt/d/ONT/clinical_649.GRCh38.vcf.unmap 2>/dev/null && wc -l /mnt/d/ONT/clinical_649.GRCh38.vcf.unmap || echo "无 unmap 文件"
