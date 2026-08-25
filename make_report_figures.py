import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import os

FIGDIR = r"D:\ONT\figures"
os.makedirs(FIGDIR, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})

RES = r"D:\ONT\sv_ewas_results.tsv"
SIG = r"D:\ONT\sv_ewas_sig_pairs.tsv"
SEL = r"D:\ONT\sv_ewas_selected_sites.tsv"
PV = r"D:\ONT\sv_ewas_pvals.npy"
INT = r"D:\ONT\interaction_significant.tsv"

P_STAR = 1.05937e-4  # BH FDR<10% threshold (narrow + 5-95% + SV-count run)

# ---- load ----
res = pd.read_csv(RES, sep="\t", usecols=["site", "chrom", "pos", "gene_name", "window", "effect", "pval", "n_carriers"])
sig = pd.read_csv(SIG, sep="\t")
sel = pd.read_csv(SEL, sep="\t")
pvals = np.load(PV)

# ============ Fig 1: Manhattan (per-site min p) ============
per_site = res.groupby("site").agg(min_p=("pval", "min"), chrom=("chrom", "first"), pos=("pos", "first")).reset_index()
chrom_order = ["chr%d" % i for i in range(1, 23)]
chrom_size = {"chr1":248956422,"chr2":242193529,"chr3":198295559,"chr4":190214555,
 "chr5":181538259,"chr6":170805979,"chr7":159345973,"chr8":145138636,"chr9":138394717,
 "chr10":133797422,"chr11":135086622,"chr12":133275309,"chr13":114364328,"chr14":107043718,
 "chr15":101991189,"chr16":90338345,"chr17":83257441,"chr18":80373285,"chr19":58617616,
 "chr20":64444167,"chr21":46709983,"chr22":50818468}
per_site = per_site[per_site.chrom.isin(chrom_order)]
offset = 0
cum_pos = []
ticks = []
for c in chrom_order:
    cum_pos.append(offset)
    ticks.append((offset + chrom_size[c] / 2, c))
    offset += chrom_size[c]
# map chrom -> base
base = {}
off = 0
for c in chrom_order:
    base[c] = off
    off += chrom_size[c]
per_site["gpos"] = per_site["chrom"].map(base) + per_site["pos"]
per_site["neglog"] = -np.log10(per_site["min_p"].clip(1e-300))

bg = per_site[(per_site.min_p < 0.05) & (per_site.min_p > P_STAR)]
hit = per_site[per_site.min_p <= P_STAR]
fig, ax = plt.subplots(figsize=(14, 4.2))
ax.scatter(bg.gpos, bg.neglog, s=1.2, c="#bbbbbb", rasterized=True)
ax.scatter(hit.gpos, hit.neglog, s=2.5, c="#d62728", rasterized=True)
ax.axhline(-np.log10(P_STAR), color="orange", ls="--", lw=1, label="FDR<10%% (p=%.1e)" % P_STAR)
ax.set_xticks([t[0] for t in ticks]); ax.set_xticklabels([t[1] for t in ticks], fontsize=8)
ax.set_xlim(0, off); ax.set_ylabel("-log10(p)"); ax.set_title("SV-associated differential methylation (per-CpG min-p across 4 windows)")
ax.legend(loc="upper right")
plt.tight_layout(); plt.savefig(os.path.join(FIGDIR, "fig1_manhattan.png"), dpi=150); plt.close()

# ============ Fig 2: window bar ============
wc = sig.window.value_counts().reindex(["up", "dn", "body", "mb"])
labels = {"up": "Upstream 100kb", "dn": "Downstream 100kb", "body": "Gene body", "mb": "+/-1Mb"}
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.bar([labels[w] for w in wc.index], wc.values, color=["#4c72b0", "#55a868", "#c44e52", "#8172b2"])
for i, v in enumerate(wc.values):
    ax.text(i, v + 40, str(v), ha="center", fontsize=10)
ax.set_ylabel("# significant CpG-window pairs (FDR<10%)")
ax.set_title("SV effect is strongest in regulatory windows")
plt.tight_layout(); plt.savefig(os.path.join(FIGDIR, "fig2_window_bar.png"), dpi=150); plt.close()

# ============ Fig 3: top genes ============
gc = sel.gene_name.value_counts().head(20)
fig, ax = plt.subplots(figsize=(7, 6))
ax.barh(range(len(gc))[::-1], gc.values[::-1], color="#4c72b0")
ax.set_yticks(range(len(gc))[::-1]); ax.set_yticklabels(gc.index[::-1], fontsize=8)
ax.set_xlabel("# SV-associated CpGs")
ax.set_title("Top 20 genes by # SV-associated CpGs")
plt.tight_layout(); plt.savefig(os.path.join(FIGDIR, "fig3_top_genes.png"), dpi=150); plt.close()

# ============ Fig 4: n_carriers histogram ============
bins = [3, 5, 10, 20, 50, 100, 200, 500, 1000]
cnt, _ = np.histogram(sig.n_carriers, bins=bins)
blab = ["3-5", "5-10", "10-20", "20-50", "50-100", "100-200", "200-500", "500-1000"]
fig, ax = plt.subplots(figsize=(7, 4))
ax.bar(blab, cnt, color="#55a868")
for i, v in enumerate(cnt):
    ax.text(i, v + 20, str(v), ha="center", fontsize=8)
ax.set_xlabel("# SV carriers in window"); ax.set_ylabel("# significant pairs")
ax.set_title("Most hits come from windows with many SV carriers (robust)")
plt.tight_layout(); plt.savefig(os.path.join(FIGDIR, "fig4_ncarriers.png"), dpi=150); plt.close()

# ============ Fig 5: volcano (significant pairs) ============
fig, ax = plt.subplots(figsize=(6.5, 5))
colors = {"up": "#4c72b0", "dn": "#55a868", "body": "#c44e52", "mb": "#8172b2"}
for w, grp in sig.groupby("window"):
    ax.scatter(grp.effect, -np.log10(grp.pval), s=6, alpha=0.6, label=labels[w], color=colors[w])
ax.axhline(-np.log10(P_STAR), color="grey", ls="--", lw=1)
ax.set_xlabel("Effect size (M-value, SV carriers vs non-carriers)")
ax.set_ylabel("-log10(p)")
ax.set_title("Volcano plot of significant pairs (FDR<10%)")
ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(FIGDIR, "fig5_volcano.png"), dpi=150); plt.close()

# ============ Fig 6: QQ plot + lambda ============
chi2 = stats.chi2.ppf(1 - pvals, 1)
lam = np.median(chi2) / stats.chi2.ppf(0.5, 1)
sub = np.random.default_rng(0).choice(len(pvals), 150000, replace=False)
obs = -np.log10(np.sort(pvals[sub]))
exp = -np.log10((np.arange(1, len(sub) + 1) - 0.5) / len(sub))
fig, ax = plt.subplots(figsize=(5.5, 5))
ax.scatter(exp, obs, s=1, c="#333333", rasterized=True)
mx = max(exp.max(), obs.max()) * 1.05
ax.plot([0, mx], [0, mx], color="red", lw=1)
ax.set_xlabel("Expected -log10(p)"); ax.set_ylabel("Observed -log10(p)")
ax.set_title("QQ plot (lambda_gc = %.2f)" % lam)
plt.tight_layout(); plt.savefig(os.path.join(FIGDIR, "fig6_qq.png"), dpi=150); plt.close()

# ============ summary table ============
summary = {
    "n_samples": 648,
    "n_cpg_total": 3059870,
    "n_tests": len(pvals),
    "n_sig_pairs": len(sig),
    "n_selected_sites": sel.site.nunique(),
    "n_genes_hit": sel.gene_name.nunique(),
    "p_threshold": P_STAR,
    "lambda_gc": round(lam, 3),
    "frac_hyper": round((sig.effect > 0).mean(), 3),
    "frac_hypo": round((sig.effect < 0).mean(), 3),
}
sd = pd.DataFrame(list(summary.items()), columns=["metric", "value"])
sd.to_csv(r"D:\ONT\report_summary.tsv", sep="\t", index=False)
print(sd.to_string(index=False))
print("\nfigures saved to", FIGDIR)
