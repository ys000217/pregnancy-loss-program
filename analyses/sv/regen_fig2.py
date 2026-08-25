import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
sig = pd.read_csv(r"D:\ONT\sv_methylation_sig_pairs.tsv", sep="\t")

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
plt.tight_layout()
plt.savefig(os.path.join(r"D:\ONT\figures", "fig2_window_bar.png"), dpi=150)
plt.close()
print("fig2_window_bar.png regenerated")
print("labels:", [labels[w] for w in wc.index])
