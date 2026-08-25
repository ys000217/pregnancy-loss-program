#!/bin/bash
export PATH="$HOME/miniconda3/envs/liftover/bin:$PATH"
echo "=== tool versions ==="
minimap2 --version
samtools --version | head -1
bcftools --version | head -1
CrossMap 2>&1 | head -1
echo "=== faidx CN1 ==="
samtools faidx /mnt/e/genotype_data/liftover/CN1_chr.fasta && echo CN1_FAIDX_OK
echo "=== faidx GRCh38 ==="
samtools faidx /mnt/e/genotype_data/liftover/GRCh38_chr.fasta && echo GRCH38_FAIDX_OK
echo "=== CN1 fai ==="
cat /mnt/e/genotype_data/liftover/CN1_chr.fasta.fai
echo "=== GRCh38 fai ==="
cat /mnt/e/genotype_data/liftover/GRCh38_chr.fasta.fai
