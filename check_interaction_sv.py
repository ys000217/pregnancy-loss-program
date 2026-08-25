import pandas as pd

# load
cov = pd.read_csv(r"D:\ONT\matrix_covariates.tsv", sep="\t")
status = dict(zip(cov.FID, cov.Status))   # FID -> 1=case 0=control
gw = pd.read_csv(r"D:\ONT\genes_windows.tsv", sep="\t")
bed = pd.read_csv(r"E:\genotype_data\liftover\GRCh38_breakpoints.unique.bed", sep="\t",
                  header=None, names=["chrom","start","end","svid","svtype","side"])
bed["sv_idx"] = bed["svid"].str.extract(r"SV(\d+)")[0].astype(int)
bed = bed[bed.svtype != "TRA"]
car = pd.read_csv(r"D:\ONT\sv_carriers.tsv", sep="\t")

def case_ctrl_counts(fid_str):
    cc = 0; ct = 0
    if isinstance(fid_str, str) and fid_str:
        for f in fid_str.split(","):
            if status.get(f) == 1:
                cc += 1
            elif status.get(f) == 0:
                ct += 1
    return cc, ct

# top interaction genes (gene, window)
targets = [("ZNF718","mb"),("DIRAS2","mb"),("KCNAB2","mb"),("ZFP37","mb"),("LARGE1","mb"),("CWH43","dn")]

for gname, win in targets:
    g = gw[gw.gene_name == gname]
    if g.empty:
        print("gene not found:", gname); continue
    g = g.iloc[0]
    if win == "mb":
        a, b = g.mb_start, g.mb_end
    elif win == "up":
        a, b = g.up_start, g.up_end
    elif win == "dn":
        a, b = g.dn_start, g.dn_end
    else:
        a, b = g.body_start, g.body_end
    # SVs with breakpoint in window
    bd = bed[(bed.chrom == g.chrom) & (bed.end >= a) & (bed.end <= b)]
    bd = bd.merge(car[["sv_idx","carriers"]], on="sv_idx", how="left")
    bd["case"], bd["ctrl"] = zip(*bd.carriers.apply(case_ctrl_counts))
    bd = bd[bd.case + bd.ctrl > 0].sort_values("end")
    print("\n===== %s  %s window (%s:%d-%d) =====" % (gname, win, g.chrom, a, b))
    print("  SVs in window: %d   (pos | type | case_car | ctrl_car)" % len(bd))
    for _, r in bd.head(40).iterrows():
        print("    %d  %s  %d  %d" % (r.end, r.svtype, r.case, r.ctrl))
