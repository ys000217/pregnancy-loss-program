import pandas as pd

# 1) KSR1P1 region SVs (chr1:25.6M-27.6M) + carriers
bed = pd.read_csv(r"E:\genotype_data\liftover\GRCh38_breakpoints.unique.bed", sep="\t",
                  header=None, names=["chrom","start","end","svid","svtype","side"])
bed["sv_idx"] = bed["svid"].str.extract(r"SV(\d+)")[0].astype(int)
car = pd.read_csv(r"D:\ONT\sv_carriers.tsv", sep="\t")
b = bed[(bed.chrom=="chr1") & (bed.end>=25_600_000) & (bed.end<=27_600_000)]
b = b.merge(car[["sv_idx","svtype","n_carriers"]], on="sv_idx", how="left")
print("KSR1P1 region (±1Mb chr1:25.6-27.6M) breakpoints:")
print(b[["end","svid","svtype_x","n_carriers"]].sort_values("end").to_string(index=False))

# 2) KSR1P1 significant pairs (effect direction)
sig = pd.read_csv(r"D:\ONT\sv_ewas_sig_pairs.tsv", sep="\t")
k = sig[sig.gene_name=="KSR1P1"].sort_values("pval")
print("\nKSR1P1 sig pairs:", len(k))
print(k[["site","window","effect","pval","n_carriers"]].head(8).to_string(index=False))
print("effect sign: + = %d, - = %d" % ((k.effect>0).sum(), (k.effect<0).sum()))
