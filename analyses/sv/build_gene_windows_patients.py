import pandas as pd
import bisect
from collections import defaultdict

BED = r"E:\genotype_data\liftover\GRCh38_breakpoints.unique.bed"
CAR = r"D:\ONT\sv_carriers.tsv"
GW = r"D:\ONT\genes_windows.tsv"
OUT = r"D:\ONT\gene_window_patients.tsv"

bed = pd.read_csv(BED, sep="\t", header=None,
                  names=["chrom", "start", "end", "svid", "svtype", "side"])
bed["sv_idx"] = bed["svid"].str.extract(r"SV(\d+)")[0].astype(int)
bed = bed[bed.svtype != "TRA"].copy()

car = pd.read_csv(CAR, sep="\t")
car["carriers"] = car["carriers"].fillna("")
bed = bed.merge(car[["sv_idx", "carriers"]], on="sv_idx", how="inner")

# (chrom, pos, sv_idx, fid) tuples
per_chrom = defaultdict(list)
for _, r in bed.iterrows():
    if r.carriers == "":
        continue
    for f in r.carriers.split(","):
        per_chrom[r.chrom].append((int(r.end), int(r.sv_idx), f))

chrom_p, chrom_sv, chrom_f = {}, {}, {}
for c, lst in per_chrom.items():
    lst.sort()
    chrom_p[c] = [x[0] for x in lst]
    chrom_sv[c] = [x[1] for x in lst]
    chrom_f[c] = [x[2] for x in lst]

gw = pd.read_csv(GW, sep="\t")
# narrow windows only: up / dn / body (drop +-1Mb)
windows = [("up", "up_start", "up_end"), ("dn", "dn_start", "dn_end"),
           ("body", "body_start", "body_end")]

out = open(OUT, "w", encoding="utf-8")
out.write("gene_id\tgene_name\tgene_type\twindow\tn_carriers\tcarriers\n")
ncar_dist = defaultdict(int)
for _, g in gw.iterrows():
    c = g.chrom
    p = chrom_p.get(c)
    if not p:
        continue
    sv, f = chrom_sv[c], chrom_f[c]
    for wname, ws, we in windows:
        a = int(g[ws]); b = int(g[we])
        lo = bisect.bisect_left(p, a)
        hi = bisect.bisect_right(p, b)
        if hi <= lo:
            continue
        pairs = set(zip(sv[lo:hi], f[lo:hi]))   # dedup L/R of same SV
        counts = defaultdict(int)
        for sid, fid in pairs:
            counts[fid] += 1
        if not counts:
            continue
        enc = ",".join("%s:%d" % (fid, cnt) for fid, cnt in counts.items())
        out.write("%s\t%s\t%s\t%s\t%d\t%s\n" % (g.gene_id, g.gene_name, g.gene_type, wname, len(counts), enc))
        ncar_dist[len(counts)] += 1
out.close()
print("wrote", OUT, flush=True)
print("carrier-rate distribution (n carriers: n windows):")
for k in sorted(ncar_dist):
    print("  %d: %d" % (k, ncar_dist[k]))
