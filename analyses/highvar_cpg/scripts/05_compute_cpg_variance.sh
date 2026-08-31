#!/bin/bash
#JSUB -q normal
#JSUB -n 32
#JSUB -J cpg_variance
#JSUB -o cpg_variance.%J.out
#JSUB -e cpg_variance.%J.err

set -euo pipefail

BASE="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/figure3/meQTL/prepare_methylation"

DIR="${BASE}/CpG_long_tmp_serial"
CPG_LIST="${BASE}/CpG_95pct.list"
OUT="${BASE}/CpG_variance.tsv"

echo "===================================="
echo "Start: $(date)"
echo "CpG list: ${CPG_LIST}"
echo "Input dir: ${DIR}"
echo "Output: ${OUT}"
echo "===================================="

awk -v list=${CPG_LIST} '
BEGIN {
    # 读入 CpG 白名单
    while ((getline < list) > 0) {
        keep[$1] = 1
    }
    close(list)
}

{
    cpg = $1
    if (cpg in keep) {
        val = $3 + 0   # 强制转数值
        sum[cpg] += val
        sumsq[cpg] += val * val
        count[cpg]++
    }
}

END {
    for (c in count) {
        mean = sum[c] / count[c]
        var = (sumsq[c] / count[c]) - (mean * mean)
        printf "%s\t%d\t%.6f\n", c, count[c], var
    }
}
' ${DIR}/*.tsv > ${OUT}

echo "===================================="
echo "Done: $(date)"
echo "Output: ${OUT}"
echo "===================================="
