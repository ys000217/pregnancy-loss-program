#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
受控富集检验: abnormal 富集的 SV 是否比背景 SV 更靠近/更强 hyper-DMR?
==================================================================
对每个 SV 计算「DMR 邻近强度分数」= 断点邻域(0/1k/10k/100k)内最强 hyper-DMR 的 |effect_size|。
然后比较: abnormal 富集(p<0.05) vs 背景(p>=0.05), 以及秩相关(-log10 p vs 分数)。
"""
import re
import numpy as np
import pandas as pd
from bisect import bisect_right
from collections import defaultdict
from scipy import stats

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
TOP_FRAC  = 0.10
MIN_SITES = 3
WINDOWS   = [0, 1_000, 10_000, 100_000]   # 0=断点落在DMR内

def chrom_map(acc):
    m = re.match(r"NC_(\d+)\.", acc)
    if not m:
        return None
    n = int(m.group(1))
    return ("chr%d" % n) if 1 <= n <= 22 else None


def load_clinical():
    c = pd.read_csv(CLINICAL, sep="\t", dtype={"Sample_ID": str})
    ab = set(c.loc[c.Group4 == "abnormal", "Sample_ID"])
    nc = set(c.loc[(c.Group1 == "case") & (c.Group4 != "abnormal"), "Sample_ID"])
    return ab, nc


def load_dmrs():
    """{name: {chrom: (starts, ends, |effect|)}}"""
    out = {}
    for name, path in DMR_FILES.items():
        d = pd.read_csv(path, sep="\t", comment="#", header=None)
        d.columns = ["chrom","start","end","name","score","num_sites",
                     "a_counts","b_counts","a_pct","b_pct","a_frac","b_frac",
                     "effect_size","cohen_h","ch_low","ch_high"]
        d = d[(d["name"] == "different") & (d["effect_size"] < 0)
              & (d["num_sites"] >= MIN_SITES)].copy()
        d["chrom"] = d["chrom"].map(chrom_map)
        d = d.dropna(subset=["chrom"])
        d = d.sort_values("effect_size").head(max(1, int(len(d) * TOP_FRAC)))
        per = defaultdict(list)
        for r in d.itertuples(index=False):
            per[r.chrom].append((int(r.start), int(r.end), abs(float(r.effect_size))))
        idx = {}
        for c, lst in per.items():
            lst.sort()
            idx[c] = (np.array([x[0] for x in lst]),
                      np.array([x[1] for x in lst]),
                      np.array([x[2] for x in lst]))
        out[name] = idx
    return out


def load_svs():
    bed = pd.read_csv(SV_BED, sep="\t", header=None,
                      names=["chrom","start","end","svid","svtype","side"])
    bed["sv_idx"] = bed["svid"].str.extract(r"SV(\d+)")[0].astype(int)
    bed = bed[bed["svtype"] != "TRA"]
    sv_bp = defaultdict(list)
    for r in bed.itertuples(index=False):
        sv_bp[r.sv_idx].append((r.chrom, int(r.start), int(r.end)))
    car = pd.read_csv(SV_CAR, sep="\t")
    car["cset"] = car["carriers"].fillna("").map(lambda s: set(s.split(",")) if s else set())
    return car, sv_bp


def max_eff_within(bp, starts, ends, effs, W):
    """断点 bp 邻域 [bp-W, bp+W] 内最强 hyper-DMR 的 |effect|, 无则 0。"""
    i = bisect_right(starts, bp + W)
    if i == 0:
        return 0.0
    m = ends[:i] >= bp - W
    if not m.any():
        return 0.0
    return float(effs[:i][m].max())


def main():
    ab, nc = load_clinical()
    dmrs = load_dmrs()
    car, sv_bp = load_svs()

    n_ab, n_nc = len(ab), len(nc)
    sv_rows = []
    for _, r in car.iterrows():
        cs = r["cset"]
        a = len(cs & ab)
        n = len(cs & nc)
        odds, p = stats.fisher_exact([[a, n_ab - a], [n, n_nc - n]])
        sv_rows.append((int(r["sv_idx"]), p))
    sv_p = pd.DataFrame(sv_rows, columns=["sv_idx", "p_ab"])

    # 每个 SV 的邻近强度分数 (对每个窗口)
    scores = {W: np.zeros(len(sv_p)) for W in WINDOWS}
    for k, (sid, _) in enumerate(sv_p.itertuples(index=False)):
        for W in WINDOWS:
            best = 0.0
            for name, idx in dmrs.items():
                for chrom, bstart, bend in sv_bp[sid]:
                    arr = idx.get(chrom)
                    if arr is None:
                        continue
                    e = max_eff_within(bstart, arr[0], arr[1], arr[2], W)
                    best = max(best, e)
            scores[W][k] = best

    sv_p["nlp"] = -np.log10(sv_p["p_ab"].clip(1e-300))
    is_enr = sv_p["p_ab"] < 0.05

    print("=== 受控富集检验 (abnormal 富集 SV vs 背景 SV) ===")
    print("SV 总数=%d, p<0.05 富集=%d, 背景=%d"
          % (len(sv_p), int(is_enr.sum()), int((~is_enr).sum())))
    print()
    print("%-8s %-10s %-10s %-10s %-14s %-14s %s" %
          ("窗口", "富集均值", "背景均值", "MWU_p", "Spearman_rho", "Spearman_p", "方向"))
    for W in WINDOWS:
        s = scores[W]
        enr = s[is_enr]
        bg = s[~is_enr]
        u, p_mwu = stats.mannwhitneyu(enr, bg, alternative="two-sided")
        rho, p_sp = stats.spearmanr(sv_p["nlp"].values, s)
        direction = "富集>背景" if enr.mean() > bg.mean() else "富集<背景"
        print("%-8s %-10.4f %-10.4f %-10.4g %-14.4g %-14.4g %s"
              % (("%d" % W if W else "0(内部)"), enr.mean(), bg.mean(),
                 p_mwu, rho, p_sp, direction))

    # 也看 top1% vs 其余
    print()
    print("=== top1% 富集 SV vs 其余 (按 -log10 p) ===")
    topk = int(0.01 * len(sv_p))
    order = np.argsort(-sv_p["nlp"].values)
    top_idx = order[:topk]
    rest_idx = order[topk:]
    for W in WINDOWS:
        s = scores[W]
        u, p = stats.mannwhitneyu(s[top_idx], s[rest_idx], alternative="two-sided")
        print("窗口 %-8s: top1%% 均值=%.4f  其余均值=%.4f  MWU_p=%.4g"
              % (("%d" % W if W else "0(内部)"), s[top_idx].mean(), s[rest_idx].mean(), p))


if __name__ == "__main__":
    main()
