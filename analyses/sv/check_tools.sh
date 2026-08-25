#!/bin/bash
echo "=== WSL 工具检查 ==="
for t in bcftools samtools CrossMap picard minimap2; do
  printf '== %s == ' "$t"
  command -v "$t" 2>/dev/null || echo "MISSING in PATH"
done
echo ""
echo "=== conda 环境 ==="
ls -d /home/*/miniconda3/envs/*/ 2>/dev/null
ls -d /root/miniconda3/envs/*/ 2>/dev/null
ls -d ~/miniconda3/envs/*/ 2>/dev/null
echo ""
echo "=== 各 env 里的关键工具 ==="
for d in /home/*/miniconda3/envs/*/bin /root/miniconda3/envs/*/bin ~/miniconda3/envs/*/bin; do
  [ -d "$d" ] && echo "[$d]" && ls "$d" 2>/dev/null | grep -iE '^(bcftools|samtools|CrossMap|picard|minimap2)'
done
echo ""
echo "=== 数据目录可访问性 ==="
ls -d /mnt/e/genotype_data 2>/dev/null && ls -d /mnt/d/ONT 2>/dev/null
