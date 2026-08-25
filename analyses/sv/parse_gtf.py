import gzip, re
from collections import Counter

path = r"D:\ONT\gencode.v46.basic.annotation.gtf.gz"
out = r"D:\ONT\genes_grch38.tsv"

# header + first lines
with gzip.open(path, "rt") as f:
    head = [next(f) for _ in range(4)]
print("--- first lines ---")
for l in head:
    print(l.rstrip()[:130])

genes = []
with gzip.open(path, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "gene":
            continue
        chrom, start, end, strand, attrs = p[0], int(p[3]), int(p[4]), p[6], p[8]
        gid = re.search(r'gene_id "([^"]+)"', attrs).group(1)
        gname = re.search(r'gene_name "([^"]+)"', attrs).group(1)
        gtype = re.search(r'gene_type "([^"]+)"', attrs).group(1)
        genes.append((chrom, start, end, strand, gid, gname, gtype))

print("total gene entries:", len(genes))
print("chrom naming sample:", sorted(set(g[0] for g in genes))[:6])
types = Counter(g[6] for g in genes)
print("gene_type counts:", dict(types.most_common(8)))

with open(out, "w", encoding="utf-8") as w:
    w.write("chrom\tstart\tend\tstrand\tgene_id\tgene_name\tgene_type\n")
    for g in genes:
        w.write("%s\t%d\t%d\t%s\t%s\t%s\t%s\n" % g)
print("wrote", out, len(genes), "rows")
