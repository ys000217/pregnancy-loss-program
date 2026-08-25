#!/bin/bash
cd /mnt/e/genotype_data/liftover
awk -F'	' '!seen[$4]++' GRCh38_breakpoints.bed > GRCh38_breakpoints.unique.bed
echo DONE
wc -l GRCh38_breakpoints.unique.bed
