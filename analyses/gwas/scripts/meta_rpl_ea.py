#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Meta-analysis: Chinese RPL (user data, hg38, N~608) + Japanese RPL (Sonehara 2024, hg19).
Both are East Asian recurrent pregnancy loss case/control GWAS.

Steps:
 1. liftover user variants hg38 -> hg19 (using the pre-computed user_hg19.bed).
 2. Build HapMap3 (chr, BP) -> rsID map from gsMap LDSC weights.
 3. Stream Sonehara; for each SNP meta-combine (sample-size weighted, sqrt(N_eff))
    with the user study when the same chr:pos + alleles match; harmonize effect alleles.
 4. Write gsMap-format sumstats (SNP A1 A2 Z N) restricted to HapMap3 SNPs.
 5. Report mean chi^2 / lambda_GC / max chi^2 and the OLA1 chr2 locus combined Z.
"""
import os
import gzip
import math
import pandas as pd
import numpy as np

USER = "/mnt/d/ONT/figure3/gwas/combine_gwas_result_v1.Status.glm.logistic.hybrid"
USER_BED = "/mnt/d/gsMap/RPL_GWAS/user_hg19.bed"           # cols: chr, start0, end, name=CHR:POS_hg38
SONEHARA = "/mnt/d/gsMap/RPL_GWAS/japanese_rpl/hum0197.v20.gwas.v1/GWASsummary_RPL_Japanese_SoneharaNatCommun2024.txt"
WEIGHTS_DIR = "/mnt/d/gsMap/gsMap_resource/LDSC_resource/weights_hm3_no_hla"
OUT = "/mnt/d/gsMap/RPL_GWAS/RPL_meta_EA.sumstats.gz"
OUT_REPORT = "/mnt/d/gsMap/RPL_GWAS/RPL_meta_EA.report.txt"

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
VALID = {"AC", "AG", "CA", "CT", "GA", "GT", "TC", "TG"}   # strand-unambiguous pairs (A1+A2)

def comp(a):
    return str(a).translate(COMP).upper()

# ---------------- 1. liftover map: name "CHR:POS_hg38" -> hg19 pos ----------------
print("[1] loading liftover map ...")
lift = {}
with open(USER_BED) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4:
            name = p[3]                 # e.g. "1:55164"  (hg38 CHR:POS)
            hg19_pos = int(p[1]) + 1    # bed start is 0-based
            lift[name] = hg19_pos
print(f"    lifted variants available: {len(lift):,}")

# ---------------- 2. load user study (hg38) ----------------
print("[2] loading user study (Chinese RPL, hg38) ...")
COLS = ["CHROM","POS","ID","REF","ALT","A1","OMITTED","TEST","OBS_CT","OR","LOG_OR_SE","Z_STAT","P"]
udf = pd.read_csv(USER, sep="\t", comment="#", header=None, names=COLS,
                  low_memory=False, na_values=[".", "NA"])
udf["CHROM"] = udf["CHROM"].astype(str).str.strip().str.replace("chr", "", regex=False)
for c in ["POS", "OBS_CT", "Z_STAT"]:
    udf[c] = pd.to_numeric(udf[c], errors="coerce")
udf = udf.dropna(subset=["Z_STAT", "OBS_CT", "POS"])
udf = udf[udf["CHROM"].isin([str(i) for i in range(1, 23)])].copy()   # autosomes only
udf["_name"] = udf["CHROM"] + ":" + udf["POS"].astype(int).astype(str)
udf["_pos19"] = udf["_name"].map(lift)
um = udf.dropna(subset=["_pos19"]).copy()
um["_pos19"] = um["_pos19"].astype(int)
um["_chr"] = um["CHROM"]
um["EFFECT"] = um["A1"].astype(str).str.upper()
um["OTHER"] = um["OMITTED"].astype(str).str.upper()
um["z"] = um["Z_STAT"].astype(float)
um["N"] = um["OBS_CT"].astype(float)
print(f"    user variants with Z: {len(udf):,} ; lifted to hg19: {len(um):,}")

# user dict keyed by "chr:pos19" -> (effect, other, z, N)
user_dict = {}
_keys = um["_chr"].astype(str) + ":" + um["_pos19"].astype(str)
for k, e, o, z, n in zip(_keys, um["EFFECT"], um["OTHER"], um["z"], um["N"]):
    if k not in user_dict:            # keep first occurrence
        user_dict[k] = (e, o, z, n)
print(f"    user dict entries: {len(user_dict):,}")

# ---------------- 3. HapMap3 (chr,BP) -> rsID ----------------
print("[3] building HapMap3 rsID map from LDSC weights ...")
rsid_map = {}
for i in range(1, 23):
    wf = os.path.join(WEIGHTS_DIR, f"weights.{i}.l2.ldscore.gz")
    if not os.path.exists(wf):
        continue
    w = pd.read_csv(wf, sep="\t", usecols=["CHR", "SNP", "BP"])
    for _, r in w.iterrows():
        k = f"{int(r.CHR)}:{int(r.BP)}"
        if k not in rsid_map:
            rsid_map[k] = str(r.SNP)
print(f"    HapMap3 rsIDs: {len(rsid_map):,}")

# ---------------- 4. stream Sonehara & meta-combine ----------------
print("[4] streaming Sonehara (Japanese RPL, hg19) and meta-combining ...")
N_CASE, N_CTRL = 1728, 24353
N_EFF_SONE = 4.0 * N_CASE * N_CTRL / (N_CASE + N_CTRL)
w_sone = math.sqrt(N_EFF_SONE)
print(f"    Sonehara N_eff = {N_EFF_SONE:.1f}  (w = {w_sone:.1f})")

use_cols = ["CHR", "POS", "Allele1", "Allele2", "BETA", "SE"]
n_out = 0
n_matched = 0
n_sone_only = 0
n_amb_drop = 0
n_no_rsid = 0
z2_list = []
matched_keys = set()
ola1_rows = []     # chr2:174,914,728-175,114,728 (hg19)

chunk_iter = pd.read_csv(SONEHARA, sep="\t", usecols=use_cols, chunksize=2_000_000)
first = True
with gzip.open(OUT, "wt") as out:
    out.write("SNP\tA1\tA2\tZ\tN\n")
    for chunk in chunk_iter:
        chunk["CHR"] = pd.to_numeric(chunk["CHR"], errors="coerce")
        chunk["POS"] = pd.to_numeric(chunk["POS"], errors="coerce")
        chunk = chunk.dropna(subset=["CHR", "POS", "BETA", "SE"])
        chunk = chunk[(chunk["CHR"] >= 1) & (chunk["CHR"] <= 22)]   # autosomes only
        chunk["z"] = (chunk["BETA"] / chunk["SE"]).astype(float)
        chunk["EFFECT"] = chunk["Allele2"].astype(str).str.upper()
        chunk["OTHER"] = chunk["Allele1"].astype(str).str.upper()
        chunk["key"] = chunk["CHR"].astype(int).astype(str) + ":" + chunk["POS"].astype(int).astype(str)

        rows = []
        for key, eff, oth, z in zip(chunk["key"], chunk["EFFECT"], chunk["OTHER"], chunk["z"]):
            eff_z = z
            n_eff = N_EFF_SONE
            matched = False
            u = user_dict.get(key)
            if u is not None:
                u_eff, u_oth, u_z, u_n = u
                w_user = math.sqrt(float(u_n))
                # harmonize Sonehara z to the user effect allele
                if {eff, oth} == {u_eff, u_oth}:
                    if eff == u_eff:
                        s_z = z
                    else:
                        s_z = -z
                elif {comp(eff), comp(oth)} == {u_eff, u_oth}:
                    if comp(eff) == u_eff:
                        s_z = z
                    else:
                        s_z = -z
                else:
                    s_z = None          # multi-allelic / mismatch -> treat as separate
                if s_z is not None:
                    eff_z = (u_z * w_user + s_z * w_sone) / math.sqrt(w_user ** 2 + w_sone ** 2)
                    n_eff = float(u_n) + N_EFF_SONE
                    matched = True
                    matched_keys.add(key)
                    # keep the USER effect/other alleles as the output convention
                    eff, oth = u_eff, u_oth

            # OLA1 hg19 window diagnostic
            if key.startswith("2:") and 174_914_728 <= int(key.split(":")[1]) <= 175_114_728:
                ola1_rows.append((key, eff, oth, eff_z, n_eff, "matched" if matched else "sone_only"))

            # strand-ambiguous drop
            if (eff + oth) not in VALID:
                n_amb_drop += 1
                continue
            rsid = rsid_map.get(key)
            if rsid is None:
                n_no_rsid += 1
                continue

            rows.append(f"{rsid}\t{eff}\t{oth}\t{eff_z:.6f}\t{n_eff:.1f}")
            z2_list.append(eff_z * eff_z)
            if matched:
                n_matched += 1
            else:
                n_sone_only += 1

        out.write("\n".join(rows) + ("\n" if rows else ""))
        n_out += len(rows)
        if first:
            first = False
print(f"    output SNPs (HapMap3): {n_out:,}")
print(f"      meta-matched         : {n_matched:,}")
print(f"      Sonehara-only        : {n_sone_only:,}")
print(f"      strand-ambiguous drop: {n_amb_drop:,}")
print(f"      no HapMap3 rsID      : {n_no_rsid:,}")

# user-only SNPs (in HapMap3 but not matched to any Sonehara row)
n_user_only = 0
with gzip.open(OUT, "at") as out:
    for key, (eff, oth, z, n) in user_dict.items():
        if key in matched_keys:
            continue
        if (eff + oth) not in VALID:
            continue
        rsid = rsid_map.get(key)
        if rsid is None:
            continue
        out.write(f"{rsid}\t{eff}\t{oth}\t{z:.6f}\t{n:.1f}\n")
        z2_list.append(z * z)
        n_user_only += 1
print(f"      user-only           : {n_user_only:,}")

# ---------------- 5. summary stats ----------------
z2 = np.array(z2_list)
n_total = len(z2)
mean_chisq = float(z2.mean())
lambda_gc = float(np.median(z2) / 0.4549)
max_chisq = float(z2.max())

print("\n" + "=" * 64)
print("META-ANALYSIS SUMMARY (East Asian RPL: China + Japan)")
print("=" * 64)
print(f"total HapMap3 SNPs written : {n_total:,}")
print(f"mean chi^2                 : {mean_chisq:.6f}")
print(f"lambda_GC                  : {lambda_gc:.4f}")
print(f"max chi^2                  : {max_chisq:.2f}")
print(f"fuel (mean chi2 - 1)       : {mean_chisq - 1:.6f}")

# OLA1 locus (their top hit) — does it survive meta-analysis?
print("\nOLA1 chr2 hg19 174,914,728-175,114,728 (matched/meta z, top by |z|):")
ola_sorted = sorted(ola1_rows, key=lambda t: abs(t[3]), reverse=True)[:15]
for key, eff, oth, z, n, src in ola_sorted:
    print(f"  {key:>14s}  {eff}/{oth}  Z={z:+.3f}  N={n:.0f}  [{src}]")

with open(OUT_REPORT, "w") as fh:
    fh.write("RPL East Asian meta (China N~608 + Japan Sonehara N_eff~6454)\n")
    fh.write(f"mean_chi2={mean_chisq:.6f}\nlambda_gc={lambda_gc:.4f}\nmax_chi2={max_chisq:.2f}\n")
    fh.write(f"n_total={n_total} n_matched={n_matched} n_sone_only={n_sone_only} n_user_only={n_user_only}\n")
    fh.write("OLA1 top (key eff/oth Z N src):\n")
    for key, eff, oth, z, n, src in ola_sorted:
        fh.write(f"{key} {eff}/{oth} {z:.3f} {n:.0f} {src}\n")

print(f"\nwrote {OUT}")
print(f"wrote {OUT_REPORT}")
