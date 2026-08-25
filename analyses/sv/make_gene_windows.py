import pandas as pd

genes = pd.read_csv(r"D:\ONT\genes_grch38.tsv", sep="\t")
# keep protein_coding only (excludes lncRNA, all pseudogenes incl. snRNA/rRNA pseudogenes, and other ncRNA)
genes = genes[genes.gene_type == "protein_coding"].copy()
print("genes after filter (protein_coding only):", len(genes))
print(genes.gene_type.value_counts().to_string())
rows = []
for _, g in genes.iterrows():
    chrom, start, end, strand = g.chrom, int(g.start), int(g.end), g.strand
    # strand-aware TSS/TES
    if strand == "-":
        tss, tes = end, start
    else:
        tss, tes = start, end
    # 4 windows (0-based half-open coordinates for overlap)
    up_s, up_e = tss - 100000, tss          # upstream 100kb (5' side)
    dn_s, dn_e = tes, tes + 100000          # downstream 100kb (3' side)
    bd_s, bd_e = start, end                 # gene body
    mb_s, mb_e = start - 1000000, end + 1000000  # ±1Mb
    rows.append([chrom, start, end, strand, g.gene_id, g.gene_name, g.gene_type,
                 up_s, up_e, dn_s, dn_e, bd_s, bd_e, mb_s, mb_e])

df = pd.DataFrame(rows, columns=[
    "chrom", "start", "end", "strand", "gene_id", "gene_name", "gene_type",
    "up_start", "up_end", "dn_start", "dn_end", "body_start", "body_end", "mb_start", "mb_end"])
df.to_csv(r"D:\ONT\genes_windows.tsv", sep="\t", index=False)
print("wrote genes_windows.tsv, n_genes =", len(df))
print(df.head(3).to_string())
# sanity: check a few window spans
print("protein_coding genes:", (df.gene_type == "protein_coding").sum())
