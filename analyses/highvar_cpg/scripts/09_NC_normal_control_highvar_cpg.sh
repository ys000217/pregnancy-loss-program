#!/bin/bash
#JSUB -q normal
#JSUB -n 32
#JSUB -J nc_hvar_cpg
#JSUB -o nc_hvar_cpg.%J.out
#JSUB -e nc_hvar_cpg.%J.err

# =============================================================================
# 步骤 2（敏感性测试，已完成，不必为改叙述而重跑）
# 剔除临床 abnormal 后，在 control + normal_case 上重算高变 CpG。
# 目的：确认全队列 abnormal 枝上的 30 例 normal_case 在去掉 abnormal 后仍成簇
# （abnormal-like）。后续主流样本按孕周重算见 10_GA_stratified_*.sh。
# - 不使用全队列 CpG_95pct.list / CpG_variance.tsv / top100k matrix
# - 输入：02 产出的 CpG_long_tmp_serial/{sample}.tsv + 临床表
#
# 修复（相对 v1 / job 217523）:
# 1) sort|head 在 set -o pipefail 下触发 SIGPIPE → exit 141
# 2) clinical Windows 换行导致 Group4 匹配失败，abnormal 未剔除
# 3) 支持断点续跑：已完成且样本数正确的中间文件可跳过
# =============================================================================

set -euo pipefail

BASE="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/figure3/meQTL/prepare_methylation"
DIR="${BASE}/CpG_long_tmp_serial"
CLINICAL="${BASE}/clinical_649.tsv"

# NC 子集内：位点至少出现在多少比例的样本中（仅统计 normal+control）
NC_COVERAGE_FRAC="0.90"
# 输出方差最高的前 N 个位点（用于后续聚类/PCA）
TOP_N="10000"

# 断点续跑：1=若中间文件存在且样本分类正确则跳过重算；0=强制全重跑
RESUME="1"

OUT_NC_SAMPLES="${BASE}/NC_samples.tsv"
OUT_COUNT="${BASE}/CpG_count_NC.tsv"
OUT_WHITELIST="${BASE}/CpG_NC_90pct.list"
OUT_VAR="${BASE}/CpG_variance_NC.tsv"
OUT_TOP_LIST="${BASE}/CpG_topNC_${TOP_N}.list"
OUT_MATRIX="${BASE}/CpG_matrix_NC.tsv"
OUT_ANNOT="${BASE}/NC_matrix_sample_annotation.tsv"
OUT_FILE_LIST="${BASE}/NC_sample_files.list"

TMPDIR="${BASE}/tmp_nc_hvar_$$"
mkdir -p "${TMPDIR}"

echo "===================================="
echo "Start: $(date)"
echo "BASE: ${BASE}"
echo "Long table dir: ${DIR}"
echo "Clinical: ${CLINICAL}"
echo "NC coverage fraction: ${NC_COVERAGE_FRAC}"
echo "Top N: ${TOP_N}"
echo "RESUME: ${RESUME}"
echo "===================================="

# 安全取 top N，避免 sort|head 的 SIGPIPE (exit 141)
safe_top_n() {
    local infile="$1"
    local n="$2"
    local outfile="$3"
    # sort 写完再截断，不用 head 关管道
    sort -k3,3nr "${infile}" > "${TMPDIR}/var_sorted.tmp"
    awk -v n="${n}" 'NR <= n {print $1}' "${TMPDIR}/var_sorted.tmp" > "${outfile}"
    rm -f "${TMPDIR}/var_sorted.tmp"
}

# =============================================================================
# Step 1: 从临床表构建 normal_case + control 样本列表（剔除 abnormal）
# =============================================================================
echo "Step 1: Build NC sample list (exclude abnormal_case)..."

if [ ! -f "${CLINICAL}" ]; then
    echo "ERROR: Clinical file not found: ${CLINICAL}" >&2
    exit 1
fi

# 统一去掉 \r，避免 Group4==abnormal 匹配失败
awk -F'\t' '
BEGIN { OFS="\t"; print "sample_id", "Class3", "Group1", "Group4" }
NR == 1 {
    for (i = 1; i <= NF; i++) {
        gsub(/\r/, "", $i)
        h[tolower($i)] = i
    }
    if (!("sample_id" in h) || !("group1" in h) || !("group4" in h)) {
        print "ERROR: clinical header missing Sample_ID/Group1/Group4" > "/dev/stderr"
        exit 2
    }
    next
}
{
    for (i = 1; i <= NF; i++) gsub(/\r/, "", $i)
    g1 = tolower($h["group1"])
    g4 = tolower($h["group4"])
    sid = $h["sample_id"]

    if (g1 == "control") {
        class3 = "control"
    } else if (g1 == "case" && g4 == "abnormal") {
        next
    } else if (g1 == "case") {
        class3 = "normal_case"
    } else {
        next
    }
    print sid, class3, $h["group1"], $h["group4"]
}
' "${CLINICAL}" > "${OUT_NC_SAMPLES}"

NC_N=$(tail -n +2 "${OUT_NC_SAMPLES}" | wc -l | awk '{print $1}')
NC_CONTROL=$(awk -F'\t' '$2=="control"' "${OUT_NC_SAMPLES}" | wc -l | awk '{print $1}')
NC_NORMAL=$(awk -F'\t' '$2=="normal_case"' "${OUT_NC_SAMPLES}" | wc -l | awk '{print $1}')
NC_ABN_SKIPPED=$(awk -F'\t' 'BEGIN{c=0} NR>1{c++} END{print c}' "${CLINICAL}")
# 粗校验：剔除后应明显少于临床表总行数-1
echo "NC samples total: ${NC_N} (control=${NC_CONTROL}, normal_case=${NC_NORMAL})"

if [ "${NC_N}" -ge 640 ]; then
    echo "ERROR: NC sample count=${NC_N} looks like abnormal was NOT excluded (expect ~600)." >&2
    echo "Check clinical_649.tsv Group4 field / line endings." >&2
    exit 1
fi

# 匹配 long table 文件（大小写不敏感；不用 find|head，避免 SIGPIPE）
> "${OUT_FILE_LIST}"
MISSING="${TMPDIR}/nc_missing_samples.list"
> "${MISSING}"

# 预建小写文件名索引
INDEX="${TMPDIR}/tsv_index.tsv"
> "${INDEX}"
for f in "${DIR}"/*.tsv; do
    [ -f "${f}" ] || continue
    bn=$(basename "${f}" .tsv)
    # bash 小写（兼容大小写样本 ID）
    bl=$(echo "${bn}" | tr 'A-Z' 'a-z')
    printf '%s\t%s\n' "${bl}" "${f}" >> "${INDEX}"
done

tail -n +2 "${OUT_NC_SAMPLES}" | while IFS=$'\t' read -r sid class3 g1 g4; do
    f="${DIR}/${sid}.tsv"
    if [ -f "${f}" ]; then
        printf '%s\t%s\n' "${sid}" "${f}" >> "${OUT_FILE_LIST}"
        continue
    fi
    key=$(echo "${sid}" | tr 'A-Z' 'a-z')
    hit=$(awk -F'\t' -v k="${key}" '$1==k {print $2; exit}' "${INDEX}")
    if [ -n "${hit}" ] && [ -f "${hit}" ]; then
        printf '%s\t%s\n' "${sid}" "${hit}" >> "${OUT_FILE_LIST}"
    else
        echo "${sid}" >> "${MISSING}"
    fi
done

NC_FOUND=$(wc -l < "${OUT_FILE_LIST}" | awk '{print $1}')
NC_MISSING=$(wc -l < "${MISSING}" | awk '{print $1}')

echo "Long table matched: ${NC_FOUND} ; missing: ${NC_MISSING}"
if [ "${NC_MISSING}" -gt 0 ]; then
    echo "WARNING: missing sample tsv files (first 10):"
    head -n 10 "${MISSING}" || true
fi

if [ "${NC_FOUND}" -lt 10 ]; then
    echo "ERROR: Too few NC samples with long table files." >&2
    exit 1
fi

MIN_COUNT=$(awk -v n="${NC_FOUND}" -v f="${NC_COVERAGE_FRAC}" \
    'BEGIN { c = int(n * f + 0.999); if (c < 2) c = 2; print c }')
echo "Minimum NC sample coverage per CpG: ${MIN_COUNT} / ${NC_FOUND}"

# 若上次（job 217523）误把 abnormal 算进方差，强制重算 count/whitelist/variance
FORCE_RECOMPUTE="0"
if [ -f "${OUT_VAR}" ] && [ "${RESUME}" = "1" ]; then
    # 旧 run 在 649 样本上算的；若现 NC_FOUND 更小，旧结果不可用
    OLD_HINT=$(head -n 5 "${OUT_COUNT}" 2>/dev/null | wc -l || echo 0)
    if [ -f "${BASE}/NC_samples.tsv.bak_wrong649" ]; then
        FORCE_RECOMPUTE="1"
    fi
    # 用标记文件：若存在且样本数为旧错误值，删掉中间结果
    if [ -f "${OUT_COUNT}" ]; then
        # 无法从 count 反推样本数；用注释文件检查
        :
    fi
fi

# 明确：上次错误分类产生的中间文件应作废
# 用户若保留旧 CpG_count_NC / whitelist / variance，且当时为 649 样本，必须重算
# 用 OUT_FILE_LIST 行数写入 stamp，与 RESUME 对比
STAMP="${BASE}/NC_run_stamp.txt"
NEED_RECOUNT="1"
NEED_WHITELIST="1"
NEED_VAR="1"

if [ "${RESUME}" = "1" ] && [ -f "${STAMP}" ]; then
    OLD_N=$(awk -F'=' '/^NC_FOUND=/{print $2}' "${STAMP}" | tail -n 1)
    if [ -n "${OLD_N}" ] && [ "${OLD_N}" = "${NC_FOUND}" ] \
        && [ -f "${OUT_COUNT}" ] && [ -f "${OUT_WHITELIST}" ] && [ -f "${OUT_VAR}" ]; then
        echo "RESUME: stamp matches NC_FOUND=${NC_FOUND}; will skip Steps 2-4 if files exist."
        NEED_RECOUNT="0"
        NEED_WHITELIST="0"
        NEED_VAR="0"
    else
        echo "RESUME: stamp mismatch or missing intermediates → recompute Steps 2-4."
        echo "  stamp NC_FOUND=${OLD_N:-NA} ; current NC_FOUND=${NC_FOUND}"
    fi
else
    echo "Will compute Steps 2-4 (RESUME=${RESUME})."
fi

# =============================================================================
# Step 2: 仅在 NC 样本上统计每位点出现次数
# =============================================================================
if [ "${NEED_RECOUNT}" = "1" ] || [ ! -f "${OUT_COUNT}" ]; then
    echo "Step 2: Count CpG coverage in NC samples only..."
    while IFS=$'\t' read -r sid fpath; do
        cut -f1 "${fpath}"
    done < "${OUT_FILE_LIST}" \
        | sort --parallel=16 \
        | uniq -c \
        | awk '{printf "%s\t%s\n", $2, $1}' > "${OUT_COUNT}"
    echo "CpG count file: ${OUT_COUNT}"
else
    echo "Step 2: SKIP (reuse ${OUT_COUNT})"
fi
echo "Distinct CpGs observed in NC: $(wc -l < "${OUT_COUNT}")"

# =============================================================================
# Step 3: NC 覆盖度白名单
# =============================================================================
if [ "${NEED_WHITELIST}" = "1" ] || [ ! -f "${OUT_WHITELIST}" ]; then
    echo "Step 3: Build NC coverage whitelist (>= ${NC_COVERAGE_FRAC} of NC samples)..."
    awk -v minc="${MIN_COUNT}" '$2 >= minc {print $1}' "${OUT_COUNT}" \
        | sort -k1,1 > "${OUT_WHITELIST}"
else
    echo "Step 3: SKIP (reuse ${OUT_WHITELIST})"
fi
NC_WHITELIST_N=$(wc -l < "${OUT_WHITELIST}" | awk '{print $1}')
echo "NC whitelist CpGs: ${NC_WHITELIST_N}"
echo "Saved: ${OUT_WHITELIST}"

# =============================================================================
# Step 4: 在 NC 样本上计算方差
# =============================================================================
if [ "${NEED_VAR}" = "1" ] || [ ! -f "${OUT_VAR}" ]; then
    echo "Step 4: Compute variance within NC samples..."
    # 用文件列表喂 awk，避免命令行过长 / ARG_MAX
    FILEARGS="${TMPDIR}/nc_files.args"
    cut -f2 "${OUT_FILE_LIST}" > "${FILEARGS}"
    awk -v list="${OUT_WHITELIST}" '
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
    ' $(cat "${FILEARGS}") > "${OUT_VAR}"
else
    echo "Step 4: SKIP (reuse ${OUT_VAR})"
fi
echo "Variance file: ${OUT_VAR}"
echo "CpGs with variance: $(wc -l < "${OUT_VAR}")"

# 写入 stamp（供下次 RESUME）
{
    echo "NC_FOUND=${NC_FOUND}"
    echo "NC_CONTROL=${NC_CONTROL}"
    echo "NC_NORMAL=${NC_NORMAL}"
    echo "MIN_COUNT=${MIN_COUNT}"
    echo "DATE=$(date -Iseconds)"
} > "${STAMP}"

# =============================================================================
# Step 5: 取 NC 方差 top N（避免 SIGPIPE）
# =============================================================================
echo "Step 5: Select top ${TOP_N} variable CpGs (NC-specific)..."
safe_top_n "${OUT_VAR}" "${TOP_N}" "${OUT_TOP_LIST}"
TOP_ACTUAL=$(wc -l < "${OUT_TOP_LIST}" | awk '{print $1}')
echo "Top list saved: ${OUT_TOP_LIST} (n=${TOP_ACTUAL})"

# =============================================================================
# Step 6: 构建 NC 样本 × top CpG 矩阵（不插补，缺失为 NA）
# =============================================================================
echo "Step 6: Build NC methylation matrix..."

cp "${OUT_TOP_LIST}" "${TMPDIR}/matrix.body"
sort -k1,1 -o "${TMPDIR}/matrix.body" "${TMPDIR}/matrix.body"

while IFS=$'\t' read -r sid fpath; do
    echo "  Adding sample: ${sid}  $(date)"

    TMP_SAMPLE="${TMPDIR}/${sid}.tsv"

    awk 'NR==FNR {gsub(/\r/,"",$1); keep[$1]; next}
         {gsub(/\r/,"",$1)}
         ($1 in keep) {print $1"\t"$3}' "${OUT_TOP_LIST}" "${fpath}" \
        | sort -k1,1 > "${TMP_SAMPLE}"

    join -a1 -e NA -o auto "${TMPDIR}/matrix.body" "${TMP_SAMPLE}" > "${TMPDIR}/matrix.tmp"
    mv "${TMPDIR}/matrix.tmp" "${TMPDIR}/matrix.body"
    rm -f "${TMP_SAMPLE}"
done < "${OUT_FILE_LIST}"

# 表头（列顺序与 join 循环一致）
echo -ne "CpG_ID" > "${OUT_MATRIX}"
while IFS=$'\t' read -r sid fpath; do
    echo -ne "\t${sid}" >> "${OUT_MATRIX}"
done < "${OUT_FILE_LIST}"
echo "" >> "${OUT_MATRIX}"

cat "${TMPDIR}/matrix.body" >> "${OUT_MATRIX}"

# Class3 注释表
echo -e "sample_id\tClass3\tGroup1\tGroup4" > "${OUT_ANNOT}"
while IFS=$'\t' read -r sid fpath; do
    line=$(awk -F'\t' -v s="${sid}" '$1==s {print; exit}' "${OUT_NC_SAMPLES}")
    if [ -n "${line}" ]; then
        echo "${line}" | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4}'
    else
        echo -e "${sid}\tNA\tNA\tNA"
    fi
done < "${OUT_FILE_LIST}" >> "${OUT_ANNOT}"

# =============================================================================
# Summary
# =============================================================================
rm -rf "${TMPDIR}"

echo "===================================="
echo "Done: $(date)"
echo "Outputs:"
echo "  NC sample list       : ${OUT_NC_SAMPLES}"
echo "  NC sample files      : ${OUT_FILE_LIST}"
echo "  NC CpG count         : ${OUT_COUNT}"
echo "  NC whitelist         : ${OUT_WHITELIST}"
echo "  NC variance          : ${OUT_VAR}"
echo "  NC top CpG list      : ${OUT_TOP_LIST}"
echo "  NC matrix            : ${OUT_MATRIX}"
echo "  NC sample annotation : ${OUT_ANNOT}"
echo "Matrix rows (CpGs): ${TOP_ACTUAL}"
echo "Matrix cols (samples): ${NC_FOUND}"
echo "===================================="
