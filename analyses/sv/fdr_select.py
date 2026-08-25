import numpy as np, pandas as pd

p = np.load(r"D:\ONT\sv_ewas_pvals.npy")
m = len(p)
print("total tests:", m)

# BH-FDR q-values
order = np.argsort(p)
sp = p[order]
cummin = np.minimum.accumulate(sp[::-1] * m / np.arange(m, 0, -1))[::-1]
q = np.empty(m)
q[order] = cummin
n_sig = int((q < 0.10).sum())
print("significant (site,window) pairs at FDR<10%%: %d / %d (%.2f%%)" % (n_sig, m, 100.0 * n_sig / m))

if n_sig == 0:
    print("no significant pairs; min q =", cummin.min())
    raise SystemExit

p_star = float(p[q < 0.10].max())
print("pval threshold (BH FDR<10%%): %.6g" % p_star)

res = pd.read_csv(r"D:\ONT\sv_ewas_results.tsv", sep="\t")
sig = res[res.pval <= p_star].copy()
print("sig pair rows:", len(sig))
print("--- sig pairs by window ---")
print(sig.window.value_counts().to_string())

sel = sig.groupby("site").agg(
    gene_id=("gene_id", "first"), gene_name=("gene_name", "first"),
    chrom=("chrom", "first"), pos=("pos", "first"),
    n_sig_windows=("window", "count"), min_p=("pval", "min"),
    windows=("window", lambda s: ",".join(sorted(s)))
).reset_index().sort_values("min_p")
print("selected sites:", len(sel))
print("--- selected sites by window(s) ---")
print(sel.windows.value_counts().head(10).to_string())
print("--- top 20 genes by # selected sites ---")
print(sel.gene_name.value_counts().head(20).to_string())
print("--- top 10 selected sites ---")
print(sel.head(10).to_string(index=False))

sel.to_csv(r"D:\ONT\sv_ewas_selected_sites.tsv", sep="\t", index=False)
sig.to_csv(r"D:\ONT\sv_ewas_sig_pairs.tsv", sep="\t", index=False)
print("saved sv_ewas_selected_sites.tsv and sv_ewas_sig_pairs.tsv")
