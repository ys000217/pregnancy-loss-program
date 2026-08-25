#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
================================================================================
胎盘 ONT 甲基化 × 种系结构变异 (SV) 分析 —— 合并主流程
================================================================================

复现论文
    Global DNA methylation differences involving germline structural variation
    impact gene expression in pediatric brain tumors (Nat Commun 2025, 16:4713)

设计要点（本版最终设定）
    1. 基因集      : 仅 protein_coding（gencode v46 basic, GRCh38，20,065 个）
    2. 窗口        : 上游100kb / 下游100kb / 基因体 / Roadmap启动子 / Roadmap增强子
                     （去掉易饱和的 ±1Mb 宽窗口；启动子/增强子用 E091 胎盘 ChromHMM 注释定义）
    3. 携带率过滤  : 窗口携带率限制在 5%–95%（随样本量动态），
                     保证「有/无 SV」两组各自 ≥5% 样本，避免饱和/过稀
    4. 自变量      : 窗口内"不同 SV 数量"（剂量-反应），替代二元"有无 SV"，
                     规避热点区 SV 过多导致 0/1 变量几乎人人 = 1 的问题
    5. 回归模型    : M值 ~ n_SV + Gestational_Week
                     (chunk4 仅保留孕周协变量；已去掉 Age 与 4 类细胞比例)
    6. 多重检验    : BH-FDR < 10%（在全量检验上），任一窗口显著即该位点入选

流程 8 个 chunk（每 chunk 一个阶段，可单独或按区间运行）
    CHUNK 0 : 解析 gencode GTF        -> genes_grch38.tsv          (一次性,可选)
    CHUNK 1 : VCF 提取 SV 载体矩阵    -> sv_carriers.tsv
    CHUNK 2 : 蛋白编码基因 4 窗口     -> genes_windows.tsv
    CHUNK 3 : 基因窗口 × SV 断点      -> gene_window_patients.tsv   (每患者 SV 计数)
    CHUNK 4 : 主回归(全基因组)        -> sv_methylation_results.tsv + sv_methylation_pvals.npy
    CHUNK 5 : BH-FDR 筛选             -> sv_methylation_sig_pairs.tsv + sv_methylation_selected_sites.tsv
    CHUNK 6 : SV × case/control 交互  -> interaction_results.tsv + interaction_significant.tsv
    CHUNK 7 : 汇报图 + 汇总表         -> figures/*.png + report_summary.tsv

用法
    python sv_methylation_pipeline.py            # 运行全部 chunk 0-7
    python sv_methylation_pipeline.py 2          # 只运行 chunk 2
    python sv_methylation_pipeline.py 2 6        # 运行 chunk 2-6

环境
    C:\ProgramData\anaconda3\python.exe
    (numpy / pandas / scipy / matplotlib)
================================================================================
"""

import os
import sys
import time
import bisect
import gzip
import re
import numpy as np
import pandas as pd
from collections import defaultdict
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================================
# 0. 全局配置（所有 chunk 共用；改这里即可切换输入/阈值）
# ============================================================================
# ---- 输入文件 ----
MATRIX   = r"E:\甲基化数据矩阵\甲基化CpG位点矩阵.txt"   # 11GB 甲基化矩阵(位点×样本)
VCF      = r"E:\genotype_data\S957.1kGPont.merged.649clin_1027kgp.vcf"  # 过滤后 SV VCF
BED      = r"E:\genotype_data\liftover\GRCh38_breakpoints.unique.bed"  # 已 lift 到 GRCh38 的断点
GTF_GZ   = r"D:\ONT\gencode.v46.basic.annotation.gtf.gz"               # gencode 注释(GRCh38)
COV      = r"D:\ONT\matrix_covariates.tsv"    # 样本协变量表(由矩阵71行头转置得到)
FID      = r"D:\ONT\matrix_fid.txt"           # 648 个样本 ID(矩阵列顺序)
CLINICAL = r"D:\ONT\clinical_649.tsv"         # 临床表(Group4 列标注 abnormal 异常样本)

# ---- 中间/输出文件 ----
GENES    = r"D:\ONT\genes_grch38.tsv"         # chunk0 输出
GWIN     = r"D:\ONT\genes_windows.tsv"        # chunk2 输出
CAR      = r"D:\ONT\sv_carriers.tsv"          # chunk1 输出
GPAT     = r"D:\ONT\gene_window_patients.tsv" # chunk3 输出
RESULTS  = r"D:\ONT\sv_methylation_results.tsv"      # chunk4 输出(全量检验结果)
PVALS    = r"D:\ONT\sv_methylation_pvals.npy"        # chunk4 输出(p 值数组)
SIGPAIRS = r"D:\ONT\sv_methylation_sig_pairs.tsv"    # chunk5 输出
SELSITES = r"D:\ONT\sv_methylation_selected_sites.tsv"
PSTAR    = r"D:\ONT\p_star.txt"               # chunk5 传给 chunk7 的 FDR 阈值
INTRES   = r"D:\ONT\interaction_results.tsv"  # chunk6 输出
INTSIG   = r"D:\ONT\interaction_significant.tsv"
FIGDIR   = r"D:\ONT\figures"

# ---- 分析参数 ----
N_HEADER   = 71        # 矩阵前 71 行为 "#" 开头的协变量头
WINDOWS    = [("up", "up_start", "up_end"), ("dn", "dn_start", "dn_end"),
              ("body", "body_start", "body_end")]   # 基因区间窗口
# ---- Roadmap 染色质状态窗口(启动子/增强子, 非区间, 见 chunk3) ----
ROADMAP          = r"D:\ONT\figure2\Roadmap_placenta\E091_Placenta_18state_hg38_chr1_22.bed"
PROMOTER_STATES  = {"TssA", "TssFlnk", "TssFlnkU", "TssFlnkD"}
ENHANCER_STATES  = {"EnhA1", "EnhA2", "EnhWk", "EnhG1", "EnhG2"}
STATE_WINDOWS    = ["prom", "enh"]   # prom=启动子状态, enh=增强子状态
STATE_LO         = 3                # 状态窗口(SV 稀少)的携带者下限
FDR        = 0.10      # 显著阈值
GENE_TYPE  = "protein_coding"   # 只保留蛋白编码基因

# ---- 异常样本处理 ----
# 47 个全局高甲基化样本(clinical Group4=='abnormal', 全部为 case)：
#   EXCLUDE_ABNORMAL=True  -> 主分析: 剔除这 47 人(剩 601 样本)
#   EXCLUDE_ABNORMAL=False -> 敏感性: 保留全部 648 样本, 并在回归中加 abnormal 二值协变量
# 携带率 5%–95% 边界随样本量动态计算(见 carrier_bounds)。
ABNORMAL_TAG     = "abnormal"
EXCLUDE_ABNORMAL = True    # 主分析: 剔除 47 异常样本(剩 601)

# ============================================================================
# 工具函数（各 chunk 共用）
# ============================================================================

def parse_site(s):
    """解析位点 ID 'NC_000006.12:105311423' -> (chrom, pos 1-based)。
    NC_ 前 3 位后的数字即染色体号(chr1-22)，'.'后是 NCBI 版本号，忽略。
    注意: 仅常染色体(chr1-22)返回 'chrN'; chrX/chrY/chrM(num>22)返回原始 accession,
    在 chunk4/6 中因与基因表的 'chrX/Y/M' 不匹配而被有意丢弃(本分析只做常染色体)。"""
    acc, pos = s.split(":")
    num = int(acc.split(".")[0][3:])
    return ("chr%d" % num if num <= 22 else acc), int(pos)


def mtransform(v):
    """甲基化比例(beta) -> M 值。
    严格按 Du et al. 2010 (BMC Bioinformatics 11:587): M = log2(beta / (1 - beta))。
    说明:
      * 论文公式在强度层面带 offset alpha=100: M = log2((Meth+alpha)/(Unmeth+alpha));
        本数据只有 beta(比例)、无原始甲基化/未甲基化强度, 故以 beta 裁剪 [1e-3, 1-1e-3]
        作为边界的等价处理(相当于给 logit 加一个小伪计数)。
      * 缺失值用位点均值填补 —— 这是本流水线固定设计矩阵的需要, 非论文定义的一部分。"""
    m = np.nanmean(v)
    v = np.where(np.isnan(v), m, v)
    v = np.clip(v, 1e-3, 1 - 1e-3)
    return np.log2(v / (1.0 - v))


def nearest_gene_tables(genes_windows_path):
    """为 CpG->基因映射预建"每条染色体按 TSS 排序"的索引。
    返回 (chrom_tss, chrom_gid, chrom_gname)：
        chrom_tss[c]  : 该染色体所有基因 TSS 升序数组
        chrom_gid[c]  : 与 TSS 对齐的 gene_id 列表
        chrom_gname[c]: 与 TSS 对齐的 gene_name 列表
    映射规则: 每个 CpG 取 TSS 距离最近的基因，若距离 ≤1Mb 则归属之，否则丢弃。"""
    genes = pd.read_csv(genes_windows_path, sep="\t")
    genes["tss"] = np.where(genes.strand == "-", genes.end, genes.start).astype(int)
    per = defaultdict(list)
    for _, g in genes.iterrows():
        per[g.chrom].append((int(g.tss), g.gene_id, g.gene_name))
    ct, cg, cn = {}, {}, {}
    for c, lst in per.items():
        lst.sort()
        ct[c] = np.array([x[0] for x in lst], dtype=np.int64)
        cg[c] = [x[1] for x in lst]
        cn[c] = [x[2] for x in lst]
    return ct, cg, cn


def abnormal_set():
    """从临床表读 Group4==ABNORMAL_TAG 的异常样本 ID 集合(字符串)。"""
    if not os.path.exists(CLINICAL):
        return set()
    c = pd.read_csv(CLINICAL, sep="\t")
    return set(c.loc[c["Group4"] == ABNORMAL_TAG, "Sample_ID"].astype(str))


def analysis_fids():
    """本次分析实际使用的样本 ID 列表(按是否剔除异常而定)。"""
    fids = [l.strip() for l in open(FID, encoding="utf-8") if l.strip()]
    if EXCLUDE_ABNORMAL:
        bad = abnormal_set()
        fids = [f for f in fids if f not in bad]
    return fids


def load_cov(fids):
    """读取协变量表并按 fids 顺序重排, 附 _abnormal 二值列(用于敏感性分析)。
    返回的 cov 行与 fids 严格对齐。"""
    cov = pd.read_csv(COV, sep="\t")
    cov["FID"] = cov["FID"].astype(str)
    cov = cov.set_index("FID").loc[fids].reset_index()
    cov["_abnormal"] = cov["FID"].isin(abnormal_set()).astype(float)
    return cov


def carrier_bounds(n):
    """携带率 5%–95% 对应携带者数(随样本量动态计算)。"""
    return int(np.ceil(0.05 * n)), int(np.floor(0.95 * n))


def load_roadmap_states():
    """加载 Roadmap E091 胎盘 ChromHMM 18-state 注释。
    返回 {chrom: (starts, ends, cat)}，cat: 0=other, 1=promoter, 2=enhancer。"""
    d = pd.read_csv(ROADMAP, sep="\t", header=None,
                    names=["chrom", "start", "end", "state"])
    def cat(s):
        if s in PROMOTER_STATES:
            return 1
        if s in ENHANCER_STATES:
            return 2
        return 0
    d["cat"] = d["state"].map(cat)
    idx = {}
    for chrom, g in d.groupby("chrom"):
        g = g.sort_values("start")
        idx[chrom] = (g["start"].to_numpy(), g["end"].to_numpy(), g["cat"].to_numpy())
    return idx


def roadmap_state(rm, chrom, pos):
    """返回断点 (chrom, pos) 的 Roadmap 状态: 'prom' / 'enh' / 'other'。"""
    arr = rm.get(chrom)
    if arr is None:
        return "other"
    starts, ends, cats = arr
    i = np.searchsorted(starts, pos, side="right") - 1
    if i >= 0 and ends[i] > pos:
        c = int(cats[i])
        if c == 1:
            return "prom"
        if c == 2:
            return "enh"
    return "other"


# ============================================================================
# CHUNK 0: 解析 gencode GTF -> 基因坐标表（一次性准备，可跳过）
# 输入  : GTF_GZ（gencode v46 basic, GRCh38）
# 输出  : GENES（genes_grch38.tsv）
# 说明  : 只抽 feature=="gene" 的行，保留 chrom/start/end/strand/gene_id/gene_name/gene_type
# ============================================================================
def chunk0_parse_gtf():
    if os.path.exists(GENES):
        print("[chunk0] %s 已存在，跳过" % GENES)
        return
    genes = []
    with gzip.open(GTF_GZ, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            attrs = p[8]
            mgid = re.search(r'gene_id "([^"]+)"', attrs)
            mgname = re.search(r'gene_name "([^"]+)"', attrs)
            mgtype = re.search(r'gene_type "([^"]+)"', attrs)
            if not (mgid and mgname and mgtype):   # 缺字段的 gene 行跳过
                continue
            genes.append((p[0], int(p[3]), int(p[4]), p[6],
                          mgid.group(1), mgname.group(1), mgtype.group(1)))
    with open(GENES, "w", encoding="utf-8") as w:
        w.write("chrom\tstart\tend\tstrand\tgene_id\tgene_name\tgene_type\n")
        for g in genes:
            w.write("%s\t%d\t%d\t%s\t%s\t%s\t%s\n" % g)
    print("[chunk0] wrote %s (%d genes)" % (GENES, len(genes)))


# ============================================================================
# CHUNK 1: 从 VCF 提取 患者×SV 载体矩阵
# 输入  : VCF（过滤后 1676 样本 = 649 临床 + 1027 1kGP）、FID（648 甲基化样本）
# 输出  : CAR（sv_carriers.tsv: sv_idx / svtype / n_carriers / carriers）
# 说明  : * 每个 VCF 数据记录按出现顺序编号为 SV<idx>（与断点 BED 的 id 一一对应）
#          * 剔除 TRA（易位），只留 DEL/INS/DUP/INV
#          * carrier 判定: 目标样本 GT(第一个冒号字段)含 "1" 即携带该 SV
# ============================================================================
def chunk1_extract_sv_carriers():
    fids = analysis_fids()   # 主分析时载体统计也只算保留样本(保证 n_carriers 分母一致)
    fid_set = set(fids)
    with open(VCF) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                break
    samples = header[9:]
    col_of_fid = {s: i for i, s in enumerate(samples) if s in fid_set}
    print("[chunk1] fids=%d found_in_vcf=%d" % (len(fids), len(col_of_fid)), flush=True)
    assert len(col_of_fid) == len(fids), \
        "部分 FID 不在 VCF 样本列中(%d/%d)，请检查 FID 与 VCF 样本 ID 一致性" % (len(col_of_fid), len(fids))
    target_cols = [col_of_fid[f] for f in fids]
    out = open(CAR, "w", encoding="utf-8")
    out.write("sv_idx\tsvtype\tn_carriers\tcarriers\n")
    idx = 0
    with open(VCF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            svtype = ""
            for field in p[7].split(";"):
                if field.startswith("SVTYPE="):
                    svtype = field[7:]
                    break
            if svtype == "TRA":          # 剔除易位
                idx += 1
                continue
            carriers = []
            for k, c in enumerate(target_cols):
                gt = p[9 + c]
                colon = gt.find(":")
                if colon > 0:
                    gt = gt[:colon]
                if "1" in gt:
                    carriers.append(fids[k])
            out.write("%d\t%s\t%d\t%s\n" % (idx, svtype, len(carriers), ",".join(carriers)))
            idx += 1
    out.close()
    print("[chunk1] done, non-TRA SVs =", idx)


# ============================================================================
# CHUNK 2: 蛋白编码基因 → 4 个窗口坐标
# 输入  : GENES（chunk0 产物）
# 输出  : GWIN（genes_windows.tsv）
# 说明  : * 只保留 protein_coding（排除 lncRNA、假基因、各类非编码 RNA）
#          * 窗口定义(strand-aware)：+链 TSS=start/TES=end；-链反之
#              up    = [TSS-100kb, TSS]          上游(启动子)
#              dn    = [TES, TES+100kb]           下游(3'调控)
#              body  = [start, end]               基因体
#              mb    = [start-1Mb, end+1Mb]       ±1Mb(本版不使用，仅保留备用)
# ============================================================================
def chunk2_gene_windows():
    genes = pd.read_csv(GENES, sep="\t")
    genes = genes[genes.gene_type == GENE_TYPE].copy()
    print("[chunk2] protein_coding genes =", len(genes))
    rows = []
    for _, g in genes.iterrows():
        chrom, start, end, strand = g.chrom, int(g.start), int(g.end), g.strand
        tss, tes = (end, start) if strand == "-" else (start, end)
        rows.append([chrom, start, end, strand, g.gene_id, g.gene_name, g.gene_type,
                     tss - 100000, tss, tes, tes + 100000, start, end,
                     start - 1000000, end + 1000000])
    df = pd.DataFrame(rows, columns=[
        "chrom", "start", "end", "strand", "gene_id", "gene_name", "gene_type",
        "up_start", "up_end", "dn_start", "dn_end", "body_start", "body_end",
        "mb_start", "mb_end"])
    df.to_csv(GWIN, sep="\t", index=False)
    print("[chunk2] wrote %s, n_genes = %d" % (GWIN, len(df)))


# ============================================================================
# CHUNK 3: 基因窗口 × SV 断点 → 每患者 SV 计数
# 输入  : GWIN + BED(断点) + CAR(载体)
# 输出  : GPAT（gene_window_patients.tsv）
#         列: gene_id / gene_name / gene_type / window / n_carriers / carriers
#         carriers = "fid:SV数,fid:SV数,..."（每个样本窗口内的不同 SV 个数）
# 说明  : * 区间窗口 WINDOWS(up/dn/body) + Roadmap 状态窗口(prom/enh)
#          * 断点落入窗口即算；同一 DEL 的 L/R 两断点用 (sv_idx,fid) 去重只计 1 次
#          * n_carriers = 携带≥1个 SV 的患者数（用于后续 5%-95% 过滤）
# ============================================================================
def chunk3_gene_window_patients():
    bed = pd.read_csv(BED, sep="\t", header=None,
                      names=["chrom", "start", "end", "svid", "svtype", "side"])
    bed["sv_idx"] = bed["svid"].str.extract(r"SV(\d+)")[0].astype(int)
    bed = bed[bed.svtype != "TRA"].copy()
    car = pd.read_csv(CAR, sep="\t")
    car["carriers"] = car["carriers"].fillna("")
    bed = bed.merge(car[["sv_idx", "carriers"]], on="sv_idx", how="inner")

    rm = load_roadmap_states()             # Roadmap 启动子/增强子注释
    per_chrom = defaultdict(list)          # chrom -> [(断点坐标, sv_idx, fid, state)]
    for _, r in bed.iterrows():
        if r.carriers == "":
            continue
        st = roadmap_state(rm, r.chrom, int(r.end))
        for f in r.carriers.split(","):
            per_chrom[r.chrom].append((int(r.end), int(r.sv_idx), f, st))
    cp, cs, cf, cst = {}, {}, {}, {}
    for c, lst in per_chrom.items():
        lst.sort()
        cp[c] = [x[0] for x in lst]
        cs[c] = [x[1] for x in lst]
        cf[c] = [x[2] for x in lst]
        cst[c] = [x[3] for x in lst]

    gw = pd.read_csv(GWIN, sep="\t")
    out = open(GPAT, "w", encoding="utf-8")
    out.write("gene_id\tgene_name\tgene_type\twindow\tn_carriers\tcarriers\n")
    for _, g in gw.iterrows():
        p = cp.get(g.chrom)
        if not p:
            continue
        sv, f, st = cs[g.chrom], cf[g.chrom], cst[g.chrom]
        # ---- 区间窗口(up/dn/body) ----
        for wname, ws, we in WINDOWS:
            a, b = int(g[ws]), int(g[we])
            lo = bisect.bisect_left(p, a)
            hi = bisect.bisect_right(p, b)
            if hi <= lo:
                continue
            pairs = set(zip(sv[lo:hi], f[lo:hi]))     # 去重同一 SV 的 L/R 断点
            counts = defaultdict(int)
            for sid, fid in pairs:
                counts[fid] += 1
            if not counts:
                continue
            enc = ",".join("%s:%d" % (fid, cnt) for fid, cnt in counts.items())
            out.write("%s\t%s\t%s\t%s\t%d\t%s\n" % (g.gene_id, g.gene_name, g.gene_type,
                                                     wname, len(counts), enc))
        # ---- Roadmap 状态窗口(prom/enh): 邻域 [up_start, dn_end] 内、对应状态的断点 ----
        a, b = int(g["up_start"]), int(g["dn_end"])
        lo = bisect.bisect_left(p, a)
        hi = bisect.bisect_right(p, b)
        prom_pairs, enh_pairs = set(), set()
        for k in range(lo, hi):
            sk = st[k]
            if sk == "prom":
                prom_pairs.add((sv[k], f[k]))
            elif sk == "enh":
                enh_pairs.add((sv[k], f[k]))
        for sname, pairs in (("prom", prom_pairs), ("enh", enh_pairs)):
            if not pairs:
                continue
            counts = defaultdict(int)
            for sid, fid in pairs:
                counts[fid] += 1
            enc = ",".join("%s:%d" % (fid, cnt) for fid, cnt in counts.items())
            out.write("%s\t%s\t%s\t%s\t%d\t%s\n" % (g.gene_id, g.gene_name, g.gene_type,
                                                     sname, len(counts), enc))
    out.close()
    print("[chunk3] wrote", GPAT)


# ============================================================================
# CHUNK 4: 主回归（全基因组，每 CpG × 每窗口）
# 输入  : MATRIX + COV + FID + GPAT + GWIN
# 输出  : RESULTS（sv_methylation_results.tsv）+ PVALS（p 值数组）
# 说明  : * 流式读取矩阵(跳过 71 行头)，每 CpG: 解析坐标 -> 映射最近蛋白编码基因(TSS≤1Mb)
#          * 对该基因的每个有效窗口(携带率5-95%)跑回归:
#              M值 ~ n_SV + Gestational_Week
#          * 关键优化: 同一 (基因,窗口) 的设计矩阵 X 对所有 CpG 相同，
#            先算 hat=(X'X)^-1 X' 缓存，每个 CpG 直接矩阵乘，避免 400 万次求逆
#          * 结果每行一个 (位点,窗口) 检验；p 值同时存入 numpy 数组供 FDR 用
# ============================================================================
def chunk4_regression():
    all_fids = [l.strip() for l in open(FID, encoding="utf-8") if l.strip()]
    fids = analysis_fids()
    n = len(fids)
    col_idx = [all_fids.index(f) for f in fids]   # 保留样本在矩阵列中的位置
    fid2idx = {f: i for i, f in enumerate(fids)}
    cov = load_cov(fids)
    Xcov = cov[["Gestational_Week"]].to_numpy(float)   # chunk4 只保留孕周协变量
    if not EXCLUDE_ABNORMAL:                           # 敏感性分析: 加 abnormal 协变量
        Xcov = np.column_stack([Xcov, cov[["_abnormal"]].to_numpy(float)])
    assert Xcov.shape[0] == n

    gwin = pd.read_csv(GPAT, sep="\t")
    lo, hi = carrier_bounds(n)
    is_state = gwin["window"].isin(STATE_WINDOWS)
    gwin = gwin[(is_state & (gwin.n_carriers >= STATE_LO)) |
                (~is_state & (gwin.n_carriers >= lo) & (gwin.n_carriers <= hi))].copy()
    gene_windows = defaultdict(list)
    for _, r in gwin.iterrows():
        mask = np.zeros(n)
        if isinstance(r.carriers, str) and r.carriers:
            for part in r.carriers.split(","):
                fid, cnt = part.rsplit(":", 1)
                if fid in fid2idx:
                    mask[fid2idx[fid]] = float(cnt)      # 每患者 SV 计数
        gene_windows[r.gene_id].append((r.window, mask, int(r.n_carriers)))
    print("[chunk4] genes_with_valid_windows=%d total_windows=%d"
          % (len(gene_windows), sum(len(v) for v in gene_windows.values())), flush=True)

    ct, cg, cn = nearest_gene_tables(GWIN)

    hat_cache = {}
    def get_hat(gid, win, mask):
        key = (gid, win)
        if key not in hat_cache:
            X = np.column_stack([np.ones(n), mask, Xcov])
            XtX_inv = np.linalg.inv(X.T @ X)
            hat_cache[key] = (X, XtX_inv @ X.T, np.diag(XtX_inv))
        return hat_cache[key]

    dof = n - (2 + Xcov.shape[1])   # 截距 + n_SV + 协变量数
    out = open(RESULTS, "w", encoding="utf-8")
    out.write("site\tchrom\tpos\tgene_id\tgene_name\twindow\teffect\tse\tt\tpval\tn_carriers\n")
    pvals = []
    t0 = time.time()
    n_site = n_tested = n_skip = 0
    reader = pd.read_csv(MATRIX, sep="\t", header=None, skiprows=N_HEADER, chunksize=20000,
                         na_values=["NA", ""], low_memory=False, dtype={0: str})
    for chunk in reader:
        ids = chunk.iloc[:, 0].values
        Y = chunk.iloc[:, 1:].to_numpy(float)[:, col_idx]
        for j in range(Y.shape[0]):
            site = ids[j]
            chrom, pos = parse_site(site)
            tss = ct.get(chrom)
            if tss is None or len(tss) == 0:
                n_skip += 1
                continue
            i = bisect.bisect_left(tss, pos)
            best, bestd = None, 1e18
            if i < len(tss) and int(tss[i]) - pos < bestd:
                best, bestd = i, int(tss[i]) - pos
            if i - 1 >= 0 and pos - int(tss[i - 1]) < bestd:
                best, bestd = i - 1, pos - int(tss[i - 1])
            if best is None or bestd > 1_000_000:
                n_skip += 1
                continue
            gid, gnm = cg[chrom][best], cn[chrom][best]
            wins = gene_windows.get(gid)
            if not wins:
                n_skip += 1
                continue
            M = mtransform(Y[j])
            for win, mask, ncar in wins:
                X, hat, diag = get_hat(gid, win, mask)
                beta = hat @ M
                resid = M - X @ beta
                sigma2 = float(resid @ resid) / dof
                se = float(np.sqrt(sigma2 * diag[1]))
                if se <= 0:            # 全相同位点(σ²=0)等退化情形, 跳过
                    continue
                tt = float(beta[1] / se)
                pv = float(2 * stats.t.sf(abs(tt), dof))
                out.write("%s\t%s\t%d\t%s\t%s\t%s\t%.6g\t%.6g\t%.6g\t%.6g\t%d\n"
                          % (site, chrom, pos, gid, gnm, win, float(beta[1]), se, tt, pv, ncar))
                pvals.append(pv)
                n_tested += 1
            n_site += 1
        if n_site % 200000 < 20000:
            print("[chunk4] sites=%d tested=%d skip=%d %.1fs"
                  % (n_site, n_tested, n_skip, time.time() - t0), flush=True)
    out.close()
    np.save(PVALS, np.array(pvals))
    print("[chunk4] DONE sites=%d tested=%d total_sec=%.1f"
          % (n_site, n_tested, time.time() - t0), flush=True)


# ============================================================================
# CHUNK 5: BH-FDR 筛选
# 输入  : PVALS + RESULTS
# 输出  : SIGPAIRS（显著"位点-窗口"对）、SELSITES（入选位点）、PSTAR（阈值供画图）
# 说明  : * Benjamini-Hochberg: q_i = min_{j>=i} p_j*m/j（排序后累计最小）
#          * q<0.10 对应 p<=p_star；一个位点任一窗口显著即入选
# ============================================================================
def chunk5_fdr():
    p = np.load(PVALS)
    m = len(p)
    order = np.argsort(p)
    sp = p[order]
    cummin = np.minimum.accumulate(sp[::-1] * m / np.arange(m, 0, -1))[::-1]
    q = np.empty(m)
    q[order] = cummin
    n_sig = int((q < FDR).sum())
    print("[chunk5] tests=%d  significant pairs=%d (FDR<%.0f%%)"
          % (m, n_sig, FDR * 100))
    if n_sig == 0:
        print("[chunk5] 无显著对，终止")
        return
    p_star = float(p[q < FDR].max())
    with open(PSTAR, "w") as fh:
        fh.write(repr(p_star))
    res = pd.read_csv(RESULTS, sep="\t")
    sig = res[res.pval <= p_star].copy()
    sig.to_csv(SIGPAIRS, sep="\t", index=False)
    sel = sig.groupby("site").agg(
        gene_id=("gene_id", "first"), gene_name=("gene_name", "first"),
        chrom=("chrom", "first"), pos=("pos", "first"),
        n_sig_windows=("window", "count"), min_p=("pval", "min"),
        windows=("window", lambda s: ",".join(sorted(s)))).reset_index().sort_values("min_p")
    sel.to_csv(SELSITES, sep="\t", index=False)
    print("[chunk5] p_star=%.3g  selected_sites=%d  genes=%d"
          % (p_star, len(sel), sel.gene_name.nunique()))
    print(sig.window.value_counts().to_string())


# ============================================================================
# CHUNK 6: SV × case/control 交互回归（仅在显著对上做）
# 输入  : SIGPAIRS + MATRIX + COV + GPAT
# 输出  : INTRES（interaction_results.tsv）+ INTSIG（显著交互）
# 说明  : * 模型: M ~ n_SV + Status + n_SV×Status + Gestational_Week（5 参数）
#          * 交互项系数=第4个(beta[3])，衡量"SV 剂量效应在 case 中是否不同"
#          * 可估性检查: case/control × 有/无SV 四格均 ≥2（防共线/奇异）
#          * 用 pinv 兜底奇异；重新流式读矩阵，只处理显著对对应的位点
# ============================================================================
def chunk6_interaction():
    all_fids = [l.strip() for l in open(FID, encoding="utf-8") if l.strip()]
    fids = analysis_fids()
    n = len(fids)
    col_idx = [all_fids.index(f) for f in fids]
    fid2idx = {f: i for i, f in enumerate(fids)}
    cov = load_cov(fids)
    status = cov["Status"].to_numpy(float)
    Xcov = cov[["Gestational_Week"]].to_numpy(float)   # chunk6 只保留孕周协变量
    if not EXCLUDE_ABNORMAL:                           # 敏感性分析: 加 abnormal 协变量
        Xcov = np.column_stack([Xcov, cov[["_abnormal"]].to_numpy(float)])

    gwin = pd.read_csv(GPAT, sep="\t")
    lo, hi = carrier_bounds(n)
    masks = {}
    for _, r in gwin.iterrows():
        if r.window in STATE_WINDOWS:
            if r.n_carriers < STATE_LO:
                continue
        elif r.n_carriers < lo or r.n_carriers > hi:
            continue
        mask = np.zeros(n)
        if isinstance(r.carriers, str) and r.carriers:
            for part in r.carriers.split(","):
                fid, cnt = part.rsplit(":", 1)
                if fid in fid2idx:
                    mask[fid2idx[fid]] = float(cnt)
        masks[(r.gene_id, r.window)] = mask

    sig = pd.read_csv(SIGPAIRS, sep="\t")
    site_info = defaultdict(list)
    for _, r in sig.iterrows():
        site_info[r.site].append((r.gene_id, r.gene_name, r.window, int(r.n_carriers)))
    print("[chunk6] sites to reprocess:", len(site_info), flush=True)

    dof = n - (4 + Xcov.shape[1])   # 截距 + n_SV + Status + 交互项 + 协变量数
    out = open(INTRES, "w", encoding="utf-8")
    out.write("site\tgene_id\tgene_name\twindow\tsv_effect\tsv_p\tint_effect\tint_p\t"
              "n_case_car\tn_case_non\tn_ctrl_car\tn_ctrl_non\tcase_car_frac\tctrl_car_frac\n")
    reader = pd.read_csv(MATRIX, sep="\t", header=None, skiprows=N_HEADER, chunksize=20000,
                         na_values=["NA", ""], low_memory=False, dtype={0: str})
    for chunk in reader:
        ids = chunk.iloc[:, 0].values
        Y = chunk.iloc[:, 1:].to_numpy(float)[:, col_idx]
        for j in range(Y.shape[0]):
            site = ids[j]
            if site not in site_info:
                continue
            M = mtransform(Y[j])
            for gid, gnm, win, ncar in site_info[site]:
                mask = masks.get((gid, win))
                if mask is None:
                    continue
                case_car = int(((mask >= 1) & (status == 1)).sum())
                case_non = int(((mask == 0) & (status == 1)).sum())
                ctrl_car = int(((mask >= 1) & (status == 0)).sum())
                ctrl_non = int(((mask == 0) & (status == 0)).sum())
                X = np.column_stack([np.ones(n), mask, status, mask * status, Xcov])
                XtX_inv = np.linalg.pinv(X.T @ X)
                beta = XtX_inv @ X.T @ M
                resid = M - X @ beta
                sigma2 = float(resid @ resid) / dof
                se = np.sqrt(sigma2 * np.diag(XtX_inv))
                p_sv = float(2 * stats.t.sf(abs(beta[1] / se[1]), dof)) if se[1] > 0 else float("nan")
                estimable = (case_car >= 2 and case_non >= 2 and ctrl_car >= 2 and ctrl_non >= 2)
                if estimable and se[3] > 0:
                    int_eff, p_int = float(beta[3]), float(2 * stats.t.sf(abs(beta[3] / se[3]), dof))
                else:
                    int_eff, p_int = float("nan"), float("nan")
                case_n = case_car + case_non
                ctrl_n = ctrl_car + ctrl_non
                case_frac = (case_car / case_n) if case_n > 0 else float("nan")
                ctrl_frac = (ctrl_car / ctrl_n) if ctrl_n > 0 else float("nan")
                out.write("%s\t%s\t%s\t%s\t%.6g\t%.6g\t%.6g\t%.6g\t%d\t%d\t%d\t%d\t%.6g\t%.6g\n"
                          % (site, gid, gnm, win, float(beta[1]), p_sv, int_eff, p_int,
                             case_car, case_non, ctrl_car, ctrl_non, case_frac, ctrl_frac))
    out.close()
    # ---- 保存显著交互(对 int_p 做 BH-FDR) ----
    d = pd.read_csv(INTRES, sep="\t")
    est = d[d.int_p.notna()]
    pp = est.int_p.values
    mm = len(pp)
    order = np.argsort(pp); ssp = pp[order]
    qq = np.empty(mm); qq[order] = np.minimum.accumulate(ssp[::-1] * mm / np.arange(mm, 0, -1))[::-1]
    if (qq < FDR).any():
        p_star = float(pp[qq < FDR].max())
        sigi = est[est.int_p <= p_star].sort_values("int_p")
        sigi.to_csv(INTSIG, sep="\t", index=False)
        print("[chunk6] significant interactions FDR<%.0f%%: %d"
              % (FDR * 100, len(sigi)))
    else:
        print("[chunk6] no significant interaction")


# ============================================================================
# CHUNK 7: 汇报图 + 汇总表
# 输入  : RESULTS / SIGPAIRS / SELSITES / PVALS / PSTAR / INTSIG
# 输出  : FIGDIR/*.png + report_summary.tsv
# 说明  : 6 张图: 曼哈顿 / 窗口柱状 / Top基因 / 携带者分布 / 火山 / QQ
# ============================================================================
def chunk7_figures():
    os.makedirs(FIGDIR, exist_ok=True)
    plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
    p_star = float(open(PSTAR).read()) if os.path.exists(PSTAR) else 1e-4

    res = pd.read_csv(RESULTS, sep="\t",
                      usecols=["site", "chrom", "pos", "gene_name", "window", "effect", "pval", "n_carriers"])
    sig = pd.read_csv(SIGPAIRS, sep="\t")
    sel = pd.read_csv(SELSITES, sep="\t")
    pvals = np.load(PVALS)

    # ---- Fig1 Manhattan (每 CpG 最小 p) ----
    per_site = res.groupby("site").agg(min_p=("pval", "min"),
                                       chrom=("chrom", "first"), pos=("pos", "first")).reset_index()
    chrom_order = ["chr%d" % i for i in range(1, 23)]
    chrom_size = {"chr1":248956422,"chr2":242193529,"chr3":198295559,"chr4":190214555,
     "chr5":181538259,"chr6":170805979,"chr7":159345973,"chr8":145138636,"chr9":138394717,
     "chr10":133797422,"chr11":135086622,"chr12":133275309,"chr13":114364328,"chr14":107043718,
     "chr15":101991189,"chr16":90338345,"chr17":83257441,"chr18":80373285,"chr19":58617616,
     "chr20":64444167,"chr21":46709983,"chr22":50818468}
    per_site = per_site[per_site.chrom.isin(chrom_order)]
    base, off = {}, 0
    for c in chrom_order:
        base[c] = off
        off += chrom_size[c]
    per_site["gpos"] = per_site.chrom.map(base) + per_site.pos
    per_site["neglog"] = -np.log10(per_site.min_p.clip(1e-300))
    bg = per_site[(per_site.min_p < 0.05) & (per_site.min_p > p_star)]
    hit = per_site[per_site.min_p <= p_star]
    fig, ax = plt.subplots(figsize=(14, 4.2))
    ax.scatter(bg.gpos, bg.neglog, s=1.2, c="#bbbbbb", rasterized=True)
    ax.scatter(hit.gpos, hit.neglog, s=2.5, c="#d62728", rasterized=True)
    ax.axhline(-np.log10(p_star), color="orange", ls="--", lw=1, label="FDR<10%% (p=%.1e)" % p_star)
    ax.set_xticks([base[c] + chrom_size[c] / 2 for c in chrom_order])
    ax.set_xticklabels(chrom_order, fontsize=8)
    ax.set_xlim(0, off); ax.set_ylabel("-log10(p)")
    ax.set_title("SV-associated differential methylation (per-CpG min-p)")
    ax.legend(loc="upper right")
    plt.tight_layout(); plt.savefig(os.path.join(FIGDIR, "fig1_manhattan.png"), dpi=150); plt.close()

    # ---- Fig2 窗口柱状 ----
    labels = {"up": "Upstream 100kb", "dn": "Downstream 100kb", "body": "Gene body",
              "prom": "Promoter", "enh": "Enhancer"}
    wc = sig.window.value_counts().reindex(["up", "dn", "body", "prom", "enh"])
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.bar([labels[w] for w in wc.index], wc.values,
           color=["#4c72b0", "#55a868", "#c44e52", "#dd8452", "#64b5cd"])
    for i, v in enumerate(wc.values):
        ax.text(i, v + 30, str(v), ha="center")
    ax.set_ylabel("# significant CpG-window pairs (FDR<10%)")
    ax.set_title("SV effect by window")
    plt.tight_layout(); plt.savefig(os.path.join(FIGDIR, "fig2_window_bar.png"), dpi=150); plt.close()

    # ---- Fig3 Top 基因 ----
    gc = sel.gene_name.value_counts().head(20)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.barh(range(len(gc))[::-1], gc.values[::-1], color="#4c72b0")
    ax.set_yticks(range(len(gc))[::-1]); ax.set_yticklabels(gc.index[::-1], fontsize=8)
    ax.set_xlabel("# SV-associated CpGs"); ax.set_title("Top 20 genes")
    plt.tight_layout(); plt.savefig(os.path.join(FIGDIR, "fig3_top_genes.png"), dpi=150); plt.close()

    # ---- Fig4 携带者数分布 ----
    bins = [3, 5, 10, 20, 50, 100, 200, 500, 1000]
    cnt, _ = np.histogram(sig.n_carriers, bins=bins)
    blab = ["3-5", "5-10", "10-20", "20-50", "50-100", "100-200", "200-500", "500-1000"]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(blab, cnt, color="#55a868")
    for i, v in enumerate(cnt):
        ax.text(i, v + 20, str(v), ha="center", fontsize=8)
    ax.set_xlabel("# SV carriers in window"); ax.set_ylabel("# significant pairs")
    ax.set_title("Most hits have many carriers (robust)")
    plt.tight_layout(); plt.savefig(os.path.join(FIGDIR, "fig4_ncarriers.png"), dpi=150); plt.close()

    # ---- Fig5 火山 ----
    colors = {"up": "#4c72b0", "dn": "#55a868", "body": "#c44e52", "prom": "#dd8452", "enh": "#64b5cd"}
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for w, grp in sig.groupby("window"):
        ax.scatter(grp.effect, -np.log10(grp.pval), s=6, alpha=0.6, label=labels[w], color=colors[w])
    ax.axhline(-np.log10(p_star), color="grey", ls="--", lw=1)
    ax.set_xlabel("Effect (M-value per SV)"); ax.set_ylabel("-log10(p)")
    ax.set_title("Volcano of significant pairs"); ax.legend(fontsize=8)
    plt.tight_layout(); plt.savefig(os.path.join(FIGDIR, "fig5_volcano.png"), dpi=150); plt.close()

    # ---- Fig6 QQ + lambda ----
    chi2 = stats.chi2.ppf(1 - pvals, 1)
    lam = np.median(chi2) / stats.chi2.ppf(0.5, 1)
    sub = np.random.default_rng(0).choice(len(pvals), min(150000, len(pvals)), replace=False)
    obs = -np.log10(np.sort(pvals[sub]))
    exp = -np.log10((np.arange(1, len(sub) + 1) - 0.5) / len(sub))
    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.scatter(exp, obs, s=1, c="#333333", rasterized=True)
    mx = max(exp.max(), obs.max()) * 1.05
    ax.plot([0, mx], [0, mx], color="red", lw=1)
    ax.set_xlabel("Expected -log10(p)"); ax.set_ylabel("Observed -log10(p)")
    ax.set_title("QQ plot (lambda_gc = %.2f)" % lam)
    plt.tight_layout(); plt.savefig(os.path.join(FIGDIR, "fig6_qq.png"), dpi=150); plt.close()

    # ---- 汇总表 ----
    summary = {"n_samples": len(analysis_fids()), "n_sites_tested": res.site.nunique(), "n_tests": len(pvals),
               "n_sig_pairs": len(sig), "n_selected_sites": sel.site.nunique(),
               "n_genes_hit": sel.gene_name.nunique(), "p_threshold": p_star,
               "lambda_gc": round(lam, 3), "frac_hyper": round((sig.effect > 0).mean(), 3),
               "frac_hypo": round((sig.effect < 0).mean(), 3)}
    sd = pd.DataFrame(list(summary.items()), columns=["metric", "value"])
    sd.to_csv(r"D:\ONT\report_summary.tsv", sep="\t", index=False)
    print("[chunk7] figures + report_summary.tsv saved")


# ============================================================================
# 入口：解析 chunk 区间并顺序执行
#   无参数       -> 全部 0-7
#   python xx.py 2      -> 只跑 chunk 2
#   python xx.py 2 6    -> 跑 chunk 2-6
# ============================================================================
CHUNKS = {0: chunk0_parse_gtf, 1: chunk1_extract_sv_carriers, 2: chunk2_gene_windows,
          3: chunk3_gene_window_patients, 4: chunk4_regression, 5: chunk5_fdr,
          6: chunk6_interaction, 7: chunk7_figures}

if __name__ == "__main__":
    args = [int(a) for a in sys.argv[1:]]
    lo = args[0] if len(args) >= 1 else 0
    hi = args[1] if len(args) >= 2 else lo if len(args) == 1 else 7
    for c in range(lo, hi + 1):
        if c not in CHUNKS:
            print("chunk %d 不存在" % c)
            continue
        print("\n" + "=" * 70)
        print("CHUNK %d 开始" % c)
        print("=" * 70)
        t = time.time()
        CHUNKS[c]()
        print("CHUNK %d 完成, %.1fs" % (c, time.time() - t))
