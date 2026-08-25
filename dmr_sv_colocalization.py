#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
DMR × SV 联合分析（阶段 0 + 1）
================================
目标: 找「abnormal 特异」且与「abnormal hyper-DMR(±100kb)」共定位的 SV。

口径:
  * DMR = abnormal vs 真 control 的差异甲基化区; effect_size = a-b = control-abnormal,
        故 effect_size < 0 表示 abnormal 更高(hyper)。
  * DMR 筛选: name=different 且 effect<0 且 num_sites>=3, 按 |effect| 取 top 10%。
  * SV 特异: 对每个 SV 做 abnormal(47) vs normal-case(448) 的 Fisher 精确检验。
  * 共定位: SV 断点落在 hyper-DMR 的 [start-100kb, end+100kb] 内。
  * 5 个 DMR 文件分开做(时间维度: 各自孕周/类型); spl_g9 作废不使用。
不读 11GB 甲基化矩阵, 只用 SV 携带者 + 断点 BED + 临床表 + DMR 文件。
"""
import re
import numpy as np
import pandas as pd
from bisect import bisect_right
from collections import defaultdict
from scipy import stats

# ================= 配置 =================
CLINICAL  = r"D:\ONT\clinical_649.tsv"
SV_CAR    = r"D:\ONT\sv_carriers.tsv"
SV_BED    = r"E:\genotype_data\liftover\GRCh38_breakpoints.unique.bed"
DMR_FILES = {
    "spl_g8":  r"D:\ONT\figure2\abnormal_spl_g8\segments_genome.bed",
    "spl_g10": r"D:\ONT\figure2\abnormal_spl_g10\segments_genome.bed",
    "rpl_g8":  r"D:\ONT\figure2\abnormal_rpl_g8\segments_genome.bed",
    "rpl_g9":  r"D:\ONT\figure2\abnormal_rpl_g9\segments_genome.bed",
    "rpl_g10": r"D:\ONT\figure2\abnormal_rpl_g10\segments_genome.bed",
}
TOP_FRAC  = 0.10        # DMR 取 top 10%
MIN_SITES = 3           # DMR 最少 CpG 数
WINDOW    = 100_000     # 共定位窗口 ±100kb
SV_P      = 0.05        # SV abnormal 富集 p 阈值(首筛, 未 FDR)
OUT       = r"D:\ONT\dmr_sv_colocalization.tsv"

# ================= 工具 =================
def chrom_map(acc):
    """NCBI accession -> chrN (仅常染色体)。"""
    m = re.match(r"NC_(\d+)\.", acc)
    if not m:
        return None
    n = int(m.group(1))
    return ("chr%d" % n) if 1 <= n <= 22 else None


def load_clinical():
    c = pd.read_csv(CLINICAL, sep="\t", dtype={"Sample_ID": str})
    abnormal    = set(c.loc[c.Group4 == "abnormal", "Sample_ID"])
    normal_case = set(c.loc[(c.Group1 == "case") & (c.Group4 != "abnormal"), "Sample_ID"])
    control     = set(c.loc[c.Group1 == "control", "Sample_ID"])
    return abnormal, normal_case, control


def load_dmrs():
    """返回 {name: {chrom: [(start, end, effect_size), ...]}}，仅 hyper top10%。"""
    out = {}
    for name, path in DMR_FILES.items():
        d = pd.read_csv(path, sep="\t", comment="#", header=None)
        d.columns = ["chrom", "start", "end", "name", "score", "num_sites",
                     "a_counts", "b_counts", "a_pct", "b_pct", "a_frac", "b_frac",
                     "effect_size", "cohen_h", "ch_low", "ch_high"]
        d = d[(d["name"] == "different") & (d["effect_size"] < 0)
              & (d["num_sites"] >= MIN_SITES)].copy()
        d["chrom"] = d["chrom"].map(chrom_map)
        d = d.dropna(subset=["chrom"])
        d = d.sort_values("effect_size")                 # 最负(最强 hyper)在前
        d = d.head(max(1, int(len(d) * TOP_FRAC)))
        per = defaultdict(list)
        for r in d.itertuples(index=False):
            per[r.chrom].append((int(r.start), int(r.end), float(r.effect_size)))
        for c in per:
            per[c].sort()
        out[name] = per
        print("[DMR %s] hyper top10%% = %d 条" % (name, len(d)))
    return out


def load_svs():
    """返回 (car_df, sv_breakpoints)。car_df 每行含 sv_idx/svtype/carrier_set。"""
    bed = pd.read_csv(SV_BED, sep="\t", header=None,
                      names=["chrom", "start", "end", "svid", "svtype", "side"])
    bed["sv_idx"] = bed["svid"].str.extract(r"SV(\d+)")[0].astype(int)
    bed = bed[bed["svtype"] != "TRA"]
    sv_bp = defaultdict(list)
    for r in bed.itertuples(index=False):
        sv_bp[r.sv_idx].append((r.chrom, int(r.start), int(r.end)))

    car = pd.read_csv(SV_CAR, sep="\t")
    car["carrier_set"] = car["carriers"].fillna("").map(
        lambda s: set(s.split(",")) if s else set())
    return car, sv_bp


def main():
    abnormal, normal_case, control = load_clinical()
    print("[clinical] abnormal=%d  normal_case=%d  control=%d"
          % (len(abnormal), len(normal_case), len(control)))

    dmrs = load_dmrs()

    car, sv_bp = load_svs()
    print("[SV] non-TRA 且有断点 = %d" % len(sv_bp))

    # ---- 诊断: 样本 ID 匹配 ----
    all_ids = set()
    for s in car["carrier_set"]:
        all_ids |= s
    print("[diag] abnormal 出现在 SV 载体中: %d/%d" % (len(abnormal & all_ids), len(abnormal)))
    print("[diag] normal_case 出现在 SV 载体中: %d/%d" % (len(normal_case & all_ids), len(normal_case)))
    print("[diag] control 出现在 SV 载体中: %d/%d" % (len(control & all_ids), len(control)))

    # ---- 每个 SV: 三组携带计数 + Fisher(abnormal vs normal-case) ----
    n_ab, n_nc, n_co = len(abnormal), len(normal_case), len(control)
    rows = []
    for _, r in car.iterrows():
        cs = r["carrier_set"]
        a  = len(cs & abnormal)
        nc = len(cs & normal_case)
        co = len(cs & control)
        # Fisher 2x2: [携带, 未携带] x [abnormal, normal-case]
        odds, p = stats.fisher_exact([[a, n_ab - a], [nc, n_nc - nc]])
        if p > SV_P:
            continue
        rows.append({
            "sv_idx": int(r["sv_idx"]), "svtype": r["svtype"],
            "n_abnormal": a, "n_normal": nc, "n_control": co,
            "freq_abnormal": a / n_ab, "freq_normal": nc / n_nc,
            "p_ab_vs_normal": p, "or_ab_vs_normal": odds,
        })

    sig = pd.DataFrame(rows)
    print("[SV] abnormal 富集(p<%.3g) = %d 个" % (SV_P, len(sig)))

    # ---- 共定位: 仅对显著 SV ----
    # 预建每个文件的 (chrom -> starts, ends, effects) numpy 便于查询
    dmr_idx = {}
    for name, per in dmrs.items():
        idx = {}
        for chrom, lst in per.items():
            idx[chrom] = (np.array([x[0] for x in lst]),
                          np.array([x[1] for x in lst]),
                          np.array([x[2] for x in lst]))
        dmr_idx[name] = idx

    hit_cols = []
    for _, s in sig.iterrows():
        sid = int(s["sv_idx"])
        hits = []
        for name, idx in dmr_idx.items():
            for chrom, bstart, bend in sv_bp[sid]:
                arr = idx.get(chrom)
                if arr is None:
                    continue
                starts, ends, effs = arr
                bp = bstart                       # 断点 0-based 起点
                i = bisect_right(starts, bp + WINDOW)   # starts[:i] <= bp+WINDOW
                if i == 0:
                    continue
                cand_ends = ends[:i]
                cand_effs = effs[:i]
                m = np.where(cand_ends >= bp - WINDOW)[0]
                for k in m:
                    hits.append((name, chrom, int(starts[k]), int(ends[k]), float(effs[k])))
        # 去重 + 摘要
        uniq = sorted(set(hits))
        hit_cols.append({
            "n_dmr_hits": len(uniq),
            "dmr_hits": "; ".join("%s:%s:%d-%d(ef=%.3f)" % (n, c, s, e, f)
                                  for n, c, s, e, f in uniq),
        })

    sig = sig.reset_index(drop=True)
    sig = pd.concat([sig, pd.DataFrame(hit_cols)], axis=1)
    sig["breakpoints"] = [ ";".join("%s:%d" % (c, st) for c, st, _ in sv_bp[int(i)])
                           for i in sig["sv_idx"] ]
    sig = sig.sort_values(["n_dmr_hits", "p_ab_vs_normal"], ascending=[False, True])

    sig.to_csv(OUT, sep="\t", index=False)
    n_hit = int((sig["n_dmr_hits"] > 0).sum())
    print("[结果] 显著 SV 共 %d, 其中与 hyper-DMR 共定位 %d 个" % (len(sig), n_hit))
    print("[保存] %s" % OUT)

    # ---- 打印 top 共定位候选 ----
    print("\n===== 与 DMR 共定位的候选 SV (top 30) =====")
    cols = ["sv_idx", "svtype", "n_abnormal", "n_normal", "n_control",
            "freq_abnormal", "freq_normal", "p_ab_vs_normal", "or_ab_vs_normal",
            "n_dmr_hits", "breakpoints", "dmr_hits"]
    top = sig[sig["n_dmr_hits"] > 0].head(30)
    pd.set_option("display.max_colwidth", 120)
    print(top[cols].to_string(index=False))


if __name__ == "__main__":
    main()
