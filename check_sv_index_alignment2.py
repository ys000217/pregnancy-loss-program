#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
SV 编号对齐校验 (正确版)
========================
VCF 是 CN1 坐标, BED 是 CrossMap 升到 GRCh38 后的坐标, 故不能比位置。
正确做法: 比对坐标无关的 SVTYPE(以及染色体跨坐标系统大体一致)。
另查 CN1_breakpoints.bed(升前) 确认编号 = VCF 记录顺序。
"""
import os
import pandas as pd

VCF = r"E:\genotype_data\S957.1kGPont.merged.649clin_1027kgp.vcf"
BED = r"E:\genotype_data\liftover\GRCh38_breakpoints.unique.bed"
LIFT_DIR = r"E:\genotype_data\liftover"

# ---- 1. 读 VCF 全部数据记录, idx -> svtype ----
vcf_type = {}
idx = 0
with open(VCF, encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.split("\t", 8)         # 只切到 INFO(第8列), 快
        info = p[7]
        svtype = ""
        for field in info.split(";"):
            if field.startswith("SVTYPE="):
                svtype = field[7:]
                break
        vcf_type[idx] = svtype
        idx += 1
print("VCF 数据记录总数 =", idx)

# ---- 2. 读 BED, sv_idx -> 该SV的svtype集合 ----
bed = pd.read_csv(BED, sep="\t", header=None,
                  names=["chrom", "start", "end", "svid", "svtype", "side"])
bed["sv_idx"] = bed["svid"].str.extract(r"SV(\d+)")[0].astype(int)
bed_type = bed.groupby("sv_idx")["svtype"].agg(lambda s: set(s))

# ---- 3. 比对 SVTYPE ----
n_common = n_type_match = n_type_diff = 0
diff_examples = []
for k, btypes in bed_type.items():
    vt = vcf_type.get(k)
    if vt is None:
        continue
    n_common += 1
    if vt in btypes:
        n_type_match += 1
    else:
        n_type_diff += 1
        if len(diff_examples) < 20:
            diff_examples.append((k, vt, btypes))

print("BED 中的 SV 有 VCF 对应记录 =", n_common)
print("SVTYPE 一致 =", n_type_match, "| 不一致 =", n_type_diff)
print("不一致率 = %.4f%%" % (100.0 * n_type_diff / max(1, n_common)))
print("VCF 总记录 %d - BED 有对应 %d = 被 liftover 丢弃/缺失 %d 条"
      % (idx, n_common, idx - n_common))
if diff_examples:
    print("不一致样例(前20):", diff_examples)

# ---- 4. 查 CN1 升前文件 ----
print("\n=== liftover 目录 ===")
for fn in sorted(os.listdir(LIFT_DIR)):
    full = os.path.join(LIFT_DIR, fn)
    if os.path.isfile(full):
        print("  %-40s %d" % (fn, os.path.getsize(full)))
