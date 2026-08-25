import pandas as pd
import numpy as np
from scipy import stats

P_STAR = 7.524288202816965e-05   # BH-FDR<10% 阈值(来自 p_star.txt)
ROADMAP = r"D:\ONT\figure2\Roadmap_placenta\E091_Placenta_18state_hg38_chr1_22.bed"
RESULTS = r"D:\ONT\sv_methylation_results.tsv"

# ---- 1. Roadmap E091 启动子/增强子状态 ----
PROM = {"TssA", "TssFlnk", "TssFlnkU", "TssFlnkD"}
ENH = {"EnhA1", "EnhA2", "EnhWk", "EnhG1", "EnhG2"}
rm = pd.read_csv(ROADMAP, sep="\t", header=None, names=["chrom", "start", "end", "state"])
def cat(s):
    if s in PROM:
        return "prom"
    if s in ENH:
        return "enh"
    return "other"
rm["cat"] = rm["state"].map(cat)
idx = {}
for chrom, g in rm.groupby("chrom"):
    g = g.sort_values("start")
    idx[chrom] = (g["start"].to_numpy(), g["end"].to_numpy(), g["cat"].to_numpy())

# ---- 2. 每个被检验位点的最小 p + 所在染色质状态 ----
res = pd.read_csv(RESULTS, sep="\t", usecols=["site", "chrom", "pos", "pval"])
per = res.groupby("site").agg(min_p=("pval", "min"),
                              chrom=("chrom", "first"), pos=("pos", "first")).reset_index()
per["state"] = "other"
for chrom, arr in idx.items():
    m = per.chrom == chrom
    if not m.any():
        continue
    starts, ends, cats = arr
    pos = per.loc[m, "pos"].to_numpy()
    i = np.searchsorted(starts, pos, side="right") - 1
    ok = (i >= 0) & (ends[i] > pos)
    st = np.full(len(pos), "other", dtype=object)
    st[ok] = cats[i[ok]]
    per.loc[m, "state"] = st
per["hit"] = per.min_p <= P_STAR

n_hit = int(per.hit.sum())
n_bg = int((~per.hit).sum())
print("被检验位点: %d   显著位点(hit): %d   非显著背景: %d" % (len(per), n_hit, n_bg))

# ---- 3. Fisher 精确检验 + 富集倍数 ----
print("\n%-6s | %-24s | %-22s | %-10s | %-10s | %-8s" %
      ("状态", "显著位点在该状态/不在", "背景位点在该状态/不在", "富集倍数", "OR", "p值"))
for state in ("prom", "enh"):
    a = int((per.hit & (per.state == state)).sum())       # 显著 且 在状态内
    b = int((per.hit & (per.state != state)).sum())       # 显著 且 不在
    c = int((~per.hit & (per.state == state)).sum())      # 背景 且 在状态内
    d = int((~per.hit & (per.state != state)).sum())      # 背景 且 不在
    or_, pv = stats.fisher_exact([[a, b], [c, d]])
    frac_hit = a / (a + b) if (a + b) else 0.0
    frac_bg = c / (c + d) if (c + d) else 0.0
    fold = frac_hit / frac_bg if frac_bg > 0 else float("inf")
    print("%-6s | %-24s | %-22s | %-10.2f | %-10.3f | %.3g"
          % (state, "%d / %d" % (a, b), "%d / %d" % (c, d), fold, or_, pv))
    print("        显著位点中占比=%.4f   背景中占比=%.4f" % (frac_hit, frac_bg))
