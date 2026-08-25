#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
SV 编号对齐校验
===============
核对 BED 的 "SV<k>" 编号 是否 == VCF 数据记录的 0-based 顺序(即 chunk1 的 sv_idx)。
方法: 取 VCF 前 N 条数据记录(chrom/POS/SVTYPE), 与 BED 中 SV<0..N-1> 的
      (chrom/断点位置/svtype) 逐条比对。
"""
import pandas as pd

VCF = r"E:\genotype_data\S957.1kGPont.merged.649clin_1027kgp.vcf"
BED = r"E:\genotype_data\liftover\GRCh38_breakpoints.unique.bed"
N = 40   # 检查前 40 个 SV

# ================= 1. BED 概览 =================
bed = pd.read_csv(BED, sep="\t", header=None,
                  names=["chrom", "start", "end", "svid", "svtype", "side"])
bed["sv_idx"] = bed["svid"].str.extract(r"SV(\d+)")[0].astype(int)
print("=== BED 概览 ===")
print("rows =", len(bed), "| unique SV =", bed.sv_idx.nunique(),
      "| sv_idx range =", bed.sv_idx.min(), "-", bed.sv_idx.max())
print("svtype 分布:")
print(bed.svtype.value_counts().to_string())
print()

# ================= 2. VCF 前 N 条数据记录 =================
vcf_recs = {}   # idx -> (chrom, pos, svtype)
idx = 0
with open(VCF, encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("#"):
            continue
        if idx >= N:
            break
        p = line.split("\t")
        chrom, pos = p[0], int(p[1])
        svtype = ""
        for field in p[7].split(";"):
            if field.startswith("SVTYPE="):
                svtype = field[7:]
                break
        vcf_recs[idx] = (chrom, pos, svtype)
        idx += 1
print("=== VCF 前 %d 条数据记录 ===" % N)
for k in range(min(N, idx)):
    print("VCF rec %2d ->" % k, vcf_recs[k])
print()

# ================= 3. 逐条比对 =================
print("=== 比对 (BED SV<k> vs VCF rec k) ===")
print("%-4s %-22s %-22s %s" % ("k", "VCF(chrom,pos,type)", "BED(chrom,pos,type)", "判定"))
n_match = n_diff = 0
for k in range(min(N, idx)):
    vc, vp, vt = vcf_recs[k]
    sub = bed[bed.sv_idx == k]
    if sub.empty:
        print("%-4d %-22s %-22s %s" % (k, "(%s,%d,%s)" % (vc, vp, vt), "(无此编号)", "BED缺失!"))
        n_diff += 1
        continue
    # BED 该 SV 的所有断点
    bpos = [int(r.end) for _, r in sub.iterrows()]
    bchroms = set(sub.chrom)
    btypes = set(sub.svtype)
    # 判定: 染色体一致 且 svtype 一致 且 任一断点与 VCF POS 接近(<500bp)
    chrom_ok = vc in bchroms
    type_ok = vt in btypes
    pos_ok = any(abs(b - vp) < 500 for b in bpos)
    ok = chrom_ok and type_ok and pos_ok
    if ok:
        n_match += 1
    else:
        n_diff += 1
    bdesc = "; ".join("%s:%d(%s)" % (r.chrom, r.end, r.svtype) for _, r in sub.iterrows())
    print("%-4d %-22s %-22s %s" % (k, "(%s,%d,%s)" % (vc, vp, vt), bdesc,
                                   "MATCH" if ok else "!! 不一致"))
print()
print("匹配 %d / 不一致 %d (共 %d)" % (n_match, n_diff, n_match + n_diff))
