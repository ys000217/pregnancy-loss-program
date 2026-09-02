#!/bin/bash
#JSUB -q normal
#JSUB -n 32
#JSUB -J w8to10_pool
#JSUB -o w8to10_pool.%J.out
#JSUB -e w8to10_pool.%J.err

# =============================================================================
# 8–10 周合并高变 CpG（job 12）——此前未做过
#
# 相对脚本 11：样本集合相同（all / noabn），但 8+9+10 周合成一层，
# 不在周内分别算方差。
#
#  1) all   : Group3∈{8,9,10} 全部样本（含 clinical abnormal）  expect ~514
#  2) noabn : 同上，剔除 clinical abnormal                     expect ~467
#
# 30 例 abnormal-like 不在 8/9/10 周，故 noabn 不再剔除 like。
# 合并层：覆盖 ≥90% → 方差 → top 10000 矩阵（不插补）
# =============================================================================

set -euo pipefail

BASE="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/figure3/meQTL/prepare_methylation"
DIR="${BASE}/CpG_long_tmp_serial"
CLINICAL="${BASE}/clinical_649.tsv"
POOL_ROOT="${BASE}/week8to10_pooled_hvar"
COVERAGE_FRAC="0.90"
TOP_N="10000"
RESUME="1"
BIN="W8_10"

mkdir -p "${POOL_ROOT}"
TMPDIR="${POOL_ROOT}/tmp_w8to10_$$"
mkdir -p "${TMPDIR}"

echo "===================================="
echo "Start: $(date)"
echo "BASE: ${BASE}"
echo "POOL_ROOT: ${POOL_ROOT}"
echo "Coverage: ${COVERAGE_FRAC}  TopN: ${TOP_N}  RESUME: ${RESUME}"
echo "===================================="

safe_top_n() {
    local infile="$1"
    local n="$2"
    local outfile="$3"
    sort -k3,3nr "${infile}" > "${TMPDIR}/var_sorted.tmp"
    awk -v n="${n}" 'NR <= n {print $1}' "${TMPDIR}/var_sorted.tmp" > "${outfile}"
    rm -f "${TMPDIR}/var_sorted.tmp"
}

if [ ! -f "${CLINICAL}" ]; then
    echo "ERROR: Clinical file not found: ${CLINICAL}" >&2
    exit 1
fi

INDEX="${TMPDIR}/tsv_index.tsv"
> "${INDEX}"
for f in "${DIR}"/*.tsv; do
    [ -f "${f}" ] || continue
    bn=$(basename "${f}" .tsv)
    bl=$(echo "${bn}" | tr 'A-Z' 'a-z')
    printf '%s\t%s\n' "${bl}" "${f}" >> "${INDEX}"
done

# dropabn=0: 8–10 周全部；dropabn=1: 剔除 clinical abnormal
prepare_cohort() {
    local tag="$1"
    local dropabn="$2"
    local nlo="$3"
    local nhi="$4"

    OUTDIR="${POOL_ROOT}/${tag}"
    mkdir -p "${OUTDIR}"
    OUT_SAMPLES="${OUTDIR}/samples.tsv"
    OUT_FILE_LIST="${OUTDIR}/sample_files.list"
    MISSING="${OUTDIR}/missing_samples.list"

    echo "======== cohort ${tag} (drop_abnormal=${dropabn}, pooled ${BIN}) ========"

    awk -F'\t' -v dropabn="${dropabn}" -v bin="${BIN}" '
    BEGIN { OFS="\t"; print "sample_id", "Class3", "Group1", "Group4", "gw_week", "ga_bin" }
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            gsub(/\r/, "", $i)
            h[tolower($i)] = i
        }
        if (!("sample_id" in h) || !("group1" in h) || !("group4" in h) || !("group3" in h)) {
            print "ERROR: clinical header missing Sample_ID/Group1/Group3/Group4" > "/dev/stderr"
            exit 2
        }
        next
    }
    {
        for (i = 1; i <= NF; i++) gsub(/\r/, "", $i)
        sid = $h["sample_id"]
        g1 = tolower($h["group1"])
        g4 = tolower($h["group4"])
        week = $h["group3"] + 0

        if (week < 8 || week > 10) next

        if (g1 == "control") class3 = "control"
        else if (g1 == "case" && g4 == "abnormal") {
            if (dropabn == 1) next
            class3 = "abnormal_case"
        }
        else if (g1 == "case") class3 = "normal_case"
        else next

        print sid, class3, $h["group1"], $h["group4"], week, bin
    }
    ' "${CLINICAL}" > "${OUT_SAMPLES}"

    local nkeep
    nkeep=$(tail -n +2 "${OUT_SAMPLES}" | wc -l | awk '{print $1}')
    echo "  n=${nkeep} (expect ${nlo}-${nhi})"
    if [ "${nkeep}" -lt "${nlo}" ] || [ "${nkeep}" -gt "${nhi}" ]; then
        echo "ERROR: unexpected n=${nkeep} for cohort ${tag}" >&2
        exit 1
    fi
    echo "  by Class3:"
    awk -F'\t' 'NR>1 {n[$2]++} END {for (c in n) print "    " c, n[c]}' "${OUT_SAMPLES}" | sort
    echo "  by week:"
    awk -F'\t' 'NR>1 {n[$5]++} END {for (w in n) print "    week" w, n[w]}' "${OUT_SAMPLES}" | sort

    > "${OUT_FILE_LIST}"
    > "${MISSING}"
    tail -n +2 "${OUT_SAMPLES}" | while IFS=$'\t' read -r sid class3 g1 g4 week bin; do
        f="${DIR}/${sid}.tsv"
        if [ -f "${f}" ]; then
            printf '%s\t%s\t%s\n' "${sid}" "${f}" "${bin}" >> "${OUT_FILE_LIST}"
            continue
        fi
        key=$(echo "${sid}" | tr 'A-Z' 'a-z')
        hit=$(awk -F'\t' -v k="${key}" '$1==k {print $2; exit}' "${INDEX}")
        if [ -n "${hit}" ] && [ -f "${hit}" ]; then
            printf '%s\t%s\t%s\n' "${sid}" "${hit}" "${bin}" >> "${OUT_FILE_LIST}"
        else
            echo "${sid}" >> "${MISSING}"
        fi
    done

    local nfound nmiss
    nfound=$(wc -l < "${OUT_FILE_LIST}" | awk '{print $1}')
    nmiss=$(wc -l < "${MISSING}" | awk '{print $1}')
    echo "  long tables: ${nfound}  missing: ${nmiss}"
    if [ "${nfound}" -lt 10 ]; then
        echo "ERROR: too few long tables for ${tag}" >&2
        exit 1
    fi
}

process_bin() {
    local bin="$1"
    local bdir="${OUTDIR}/${bin}"
    mkdir -p "${bdir}"

    local flist="${bdir}/sample_files.list"
    awk -F'\t' -v b="${bin}" '$3==b {print $1"\t"$2}' "${OUT_FILE_LIST}" > "${flist}"
    local n
    n=$(wc -l < "${flist}" | awk '{print $1}')
    echo "---- ${bin}: n=${n} ----"
    if [ "${n}" -lt 10 ]; then
        echo "WARNING: skip ${bin}, n<10"
        return 0
    fi

    local minc
    minc=$(awk -v n="${n}" -v f="${COVERAGE_FRAC}" 'BEGIN { c = int(n * f + 0.999); if (c < 2) c = 2; print c }')
    echo "  min coverage: ${minc} / ${n}"

    local countf="${bdir}/CpG_count.tsv"
    local white="${bdir}/CpG_90pct.list"
    local varf="${bdir}/CpG_variance.tsv"
    local toplist="${bdir}/CpG_top${TOP_N}.list"
    local matrix="${bdir}/CpG_matrix.tsv"
    local annot="${bdir}/sample_annotation.tsv"
    local stamp="${bdir}/run_stamp.txt"

    local need_count=1 need_white=1 need_var=1
    if [ "${RESUME}" = "1" ] && [ -f "${stamp}" ] && [ -f "${countf}" ] && [ -f "${white}" ] && [ -f "${varf}" ]; then
        local oldn
        oldn=$(awk -F'=' '/^N=/{print $2}' "${stamp}" | tail -n 1)
        if [ "${oldn}" = "${n}" ]; then
            echo "  RESUME: skip count/whitelist/variance"
            need_count=0
            need_white=0
            need_var=0
        fi
    fi

    if [ "${need_count}" = "1" ] || [ ! -f "${countf}" ]; then
        echo "  Count CpGs..."
        while IFS=$'\t' read -r sid fpath; do
            cut -f1 "${fpath}"
        done < "${flist}" \
            | sort --parallel=16 \
            | uniq -c \
            | awk '{printf "%s\t%s\n", $2, $1}' > "${countf}"
    fi

    if [ "${need_white}" = "1" ] || [ ! -f "${white}" ]; then
        echo "  Whitelist..."
        awk -v minc="${minc}" '$2 >= minc {print $1}' "${countf}" | sort -k1,1 > "${white}"
    fi
    local nw
    nw=$(wc -l < "${white}" | awk '{print $1}')
    echo "  whitelist CpGs: ${nw}"

    if [ "${need_var}" = "1" ] || [ ! -f "${varf}" ]; then
        echo "  Variance..."
        local fileargs="${TMPDIR}/${bin}.files.args"
        cut -f2 "${flist}" > "${fileargs}"
        awk -v list="${white}" '
        BEGIN {
            while ((getline < list) > 0) {
                gsub(/\r/, "", $1)
                keep[$1] = 1
            }
            close(list)
        }
        {
            cpg = $1
            if (!(cpg in keep)) next
            val = $3 + 0
            sum[cpg] += val
            sumsq[cpg] += val * val
            count[cpg]++
        }
        END {
            for (c in count) {
                mean = sum[c] / count[c]
                var = (sumsq[c] / count[c]) - (mean * mean)
                if (var < 0) var = 0
                printf "%s\t%d\t%.6f\n", c, count[c], var
            }
        }
        ' $(cat "${fileargs}") > "${varf}"
    fi

    {
        echo "N=${n}"
        echo "MIN_COUNT=${minc}"
        echo "DATE=$(date -Iseconds)"
    } > "${stamp}"

    echo "  Top ${TOP_N}..."
    safe_top_n "${varf}" "${TOP_N}" "${toplist}"
    local ntop
    ntop=$(wc -l < "${toplist}" | awk '{print $1}')

    echo "  Matrix (${ntop} CpGs x ${n} samples)..."
    cp "${toplist}" "${TMPDIR}/${bin}.body"
    sort -k1,1 -o "${TMPDIR}/${bin}.body" "${TMPDIR}/${bin}.body"

    while IFS=$'\t' read -r sid fpath; do
        echo "    + ${sid}"
        awk 'NR==FNR {gsub(/\r/,"",$1); keep[$1]; next}
             {gsub(/\r/,"",$1)}
             ($1 in keep) {print $1"\t"$3}' "${toplist}" "${fpath}" \
            | sort -k1,1 > "${TMPDIR}/${bin}.${sid}.tsv"
        join -a1 -e NA -o auto "${TMPDIR}/${bin}.body" "${TMPDIR}/${bin}.${sid}.tsv" > "${TMPDIR}/${bin}.tmp"
        mv "${TMPDIR}/${bin}.tmp" "${TMPDIR}/${bin}.body"
        rm -f "${TMPDIR}/${bin}.${sid}.tsv"
    done < "${flist}"

    echo -ne "CpG_ID" > "${matrix}"
    while IFS=$'\t' read -r sid fpath; do
        echo -ne "\t${sid}" >> "${matrix}"
    done < "${flist}"
    echo "" >> "${matrix}"
    cat "${TMPDIR}/${bin}.body" >> "${matrix}"
    rm -f "${TMPDIR}/${bin}.body"

    echo -e "sample_id\tClass3\tGroup1\tGroup4\tgw_week\tga_bin" > "${annot}"
    while IFS=$'\t' read -r sid fpath; do
        awk -F'\t' -v s="${sid}" '$1==s {print; exit}' "${OUT_SAMPLES}" >> "${annot}"
    done < "${flist}"

    echo "  Done ${bin}: ${matrix}"
}

# all: 339+101+74=514 ; noabn: 308+93+66=467
prepare_cohort "all" 0 500 530
process_bin "${BIN}"

prepare_cohort "noabn" 1 440 490
process_bin "${BIN}"

rm -rf "${TMPDIR}"

echo "===================================="
echo "Done: $(date)"
echo "all   : ${POOL_ROOT}/all/${BIN}   expect ~514"
echo "noabn : ${POOL_ROOT}/noabn/${BIN} expect ~467"
echo "===================================="
