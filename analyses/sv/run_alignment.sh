#!/bin/bash
export PATH="$HOME/miniconda3/envs/liftover/bin:$PATH"
cd /mnt/e/genotype_data/liftover
echo "[align] start $(date)"
minimap2 -cx asm5 -c --secondary=no -t 16   CN1_chr.fasta GRCh38_chr.fasta   > CN1_to_GRCh38.paf   2> CN1_to_GRCh38.paf.log
echo "[align] done $(date)"
echo "=== PAF line count ==="
wc -l CN1_to_GRCh38.paf
echo "=== PAF size ==="
ls -la CN1_to_GRCh38.paf
echo "=== log tail ==="
tail -5 CN1_to_GRCh38.paf.log
