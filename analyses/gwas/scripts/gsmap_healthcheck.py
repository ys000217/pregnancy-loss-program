#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gsMap feasibility health-check for the user's own RPL GWAS (PLINK .glm.logistic.hybrid).
- Computes mean chi^2 / lambda_GC (the decisive S-LDSC "polygenic signal" test).
- Maps CHR:POS -> rsID via the gsMap 1000G EUR LD reference (.bim), quantifying SNP overlap.
- Writes a gsMap-format sumstats file (SNP A1 A2 Z N) for the mapped SNPs.
"""
import os
import pandas as pd
import numpy as np

GWAS_DIR = "/mnt/d/ONT/figure3/gwas"
BIM_DIR = "/mnt/d/gsMap/gsMap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink"

COLS = ["CHROM", "POS", "ID", "REF", "ALT", "A1", "OMITTED", "TEST",
        "OBS_CT", "OR", "LOG_OR_SE", "Z_STAT", "P"]

FILES = [
    "combine_gwas_result_v1.Status.glm.logistic.hybrid",
    "gwas_result_RPL_v1.Status.glm.logistic.hybrid",
    "gwas_result_RPL_v3_simplified.Status.glm.logistic.hybrid",
]

# ---- build CHROM_POS -> rsID map from 1000G EUR .bim (gsMap LD reference) ----
print("Building rsID map from 1000G EUR .bim files ...")
id2rs = {}
for i in range(1, 23):
    bf = os.path.join(BIM_DIR, f"1000G.EUR.QC.{i}.bim")
    if not os.path.exists(bf):
        continue
    with open(bf) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            chrom, rs, pos = p[0], p[1], p[3]
            key = chrom + "_" + pos
            if key not in id2rs:  # keep first rsID at a multi-allelic position
                id2rs[key] = rs
print(f"  reference entries: {len(id2rs):,}\n")

summary_rows = []

for fn in FILES:
    path = os.path.join(GWAS_DIR, fn)
    if not os.path.exists(path):
        print(f"[skip] {fn} not found")
        continue
    df = pd.read_csv(path, sep="\t", comment="#", header=None,
                     names=COLS, low_memory=False, na_values=[".", "NA"])

    z = pd.to_numeric(df["Z_STAT"], errors="coerce")
    p = pd.to_numeric(df["P"], errors="coerce")
    valid = z.notna() & p.notna()

    chisq = (z[valid] ** 2).astype(float)
    mean_chisq = chisq.mean()
    lambda_gc = chisq.median() / 0.4549364
    max_chisq = chisq.max()
    n_sig = int((p[valid] < 5e-8).sum())
    n = df["OBS_CT"].median()

    print(f"=== {fn} ===")
    print(f"  total variants        : {len(df):,}")
    print(f"  valid Z (non-NaN)     : {int(valid.sum()):,}")
    print(f"  median OBS_CT (N)     : {n:.0f}")
    print(f"  mean chi^2            : {mean_chisq:.3f}   (S-LDSC needs >> 1; gsMap warns < 1.02)")
    print(f"  lambda_GC             : {lambda_gc:.3f}")
    print(f"  max chi^2             : {max_chisq:.3f}")
    print(f"  # p<5e-8              : {n_sig}")

    # map to rsID (overlap with gsMap LD reference)
    df["_key"] = df["CHROM"].astype(str) + "_" + df["POS"].astype(str)
    df["_rsid"] = df["_key"].map(id2rs)
    mapped = df["_rsid"].notna() & valid
    print(f"  mapped to 1000G EUR rsID : {int(mapped.sum()):,}  ({100*mapped.sum()/len(df):.1f}%)")

    summary_rows.append(dict(file=fn, mean_chisq=round(mean_chisq, 3),
                             lambda_gc=round(lambda_gc, 3), max_chisq=round(max_chisq, 3),
                             n_sig=n_sig, n_valid=int(valid.sum()), n_mapped=int(mapped.sum()),
                             n_total=len(df)))

    # ---- write gsMap-format sumstats for the primary (combine) file ----
    if fn == FILES[0]:
        out = df.loc[mapped, ["_rsid", "A1", "OMITTED", "Z_STAT", "OBS_CT"]].copy()
        out.columns = ["SNP", "A1", "A2", "Z", "N"]
        out["A1"] = out["A1"].str.upper()
        out["A2"] = out["A2"].str.upper()
        # drop strand-ambiguous A/T and C/G SNPs (same rule as gsMap format_sumstats)
        amb = out.apply(lambda r: {r["A1"], r["A2"]} in ({"A", "T"}, {"C", "G"}), axis=1)
        out = out[~amb]
        out["Z"] = out["Z"].astype(float)
        out["N"] = out["N"].astype(int)
        outfile = os.path.join(GWAS_DIR, "combine_gwas_result_v1.gsmap.sumstats.gz")
        out.to_csv(outfile, sep="\t", index=False, compression="gzip")
        print(f"\n  wrote gsMap-format sumstats -> {outfile}  (n={len(out):,} after strand filter)\n")

print("\n================ SUMMARY ================")
for r in summary_rows:
    verdict = "NO polygenic signal (mean chi^2 ~ 1)" if r["mean_chisq"] < 1.05 else "borderline"
    print(f"{r['file']}\n   mean_chi2={r['mean_chisq']}  lambda={r['lambda_gc']}  "
          f"n_valid={r['n_valid']:,}  n_mapped={r['n_mapped']:,}  -> {verdict}")
