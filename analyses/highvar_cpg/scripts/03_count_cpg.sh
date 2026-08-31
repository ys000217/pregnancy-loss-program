#!/bin/bash
#JSUB -q normal
#JSUB -n 56
#JSUB -J cpg_count
#JSUB -o cpg_count.%J.out
#JSUB -e cpg_count.%J.err

BASE="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/figure3/meQTL/prepare_methylation"
DIR="${BASE}/CpG_long_tmp_serial"
OUT="${BASE}/CpG_count.tsv"

echo "Start: $(date)"

cut -f1 ${DIR}/*.tsv \
  | sort --parallel=16 \
  | uniq -c \
  | awk '{printf "%s\t%s\n",$2,$1}' > ${OUT}

echo "Done: $(date)"
