#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify the meta-analysis harmonization on the user's OLA1 top hit and a few SNPs."""
import math
import pandas as pd

USER = "/mnt/d/ONT/figure3/gwas/combine_gwas_result_v1.Status.glm.logistic.hybrid"
USER_BED = "/mnt/d/gsMap/RPL_GWAS/user_hg19.bed"
SONEHARA = "/mnt/d/gsMap/RPL_GWAS/japanese_rpl/hum0197.v20.gwas.v1/GWASsummary_RPL_Japanese_SoneharaNatCommun2024.txt"

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
def comp(a): return str(a).translate(COMP).upper()

# liftover map
lift = {}
with open(USER_BED) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4:
            lift[p[3]] = int(p[1]) + 1

# user data
COLS = ["CHROM","POS","ID","REF","ALT","A1","OMITTED","TEST","OBS_CT","OR","LOG_OR_SE","Z_STAT","P"]
u = pd.read_csv(USER, sep="\t", comment="#", header=None, names=COLS, low_memory=False, na_values=[".","NA"])
u["CHROM"] = u["CHROM"].astype(str).str.replace("chr","",regex=False)
u["POS"] = pd.to_numeric(u["POS"], errors="coerce")
u["Z_STAT"] = pd.to_numeric(u["Z_STAT"], errors="coerce")
u = u.dropna(subset=["Z_STAT"])

# their top hit (from handoff: chr2:174095142, p=7.3e-7)
top = u[(u["CHROM"]=="2") & (u["POS"]==174095142)]
print("=== user top hit chr2:174095142 ===")
if len(top):
    r = top.iloc[0]
    name = f"2:174095142"
    pos19 = lift.get(name)
    print(f"  hg38 chr2:174095142 -> hg19 {pos19}")
    print(f"  REF={r['REF']} ALT={r['ALT']} A1(effect)={r['A1']} OMITTED={r['OMITTED']}")
    print(f"  Z_STAT={r['Z_STAT']:.4f}  OBS_CT={r['OBS_CT']}  P={r['P']}")
    # Sonehara at that hg19 pos
    if pos19:
        s = pd.read_csv(SONEHARA, sep="\t", usecols=["CHR","POS","Allele1","Allele2","BETA","SE","N_case","N_ctrl"])
        s = s[(s["CHR"]==2) & (s["POS"]==pos19)]
        if len(s):
            sr = s.iloc[0]
            s_z = sr["BETA"]/sr["SE"]
            print(f"  Sonehara @ hg19 {pos19}: Allele1(ref)={sr['Allele1']} Allele2(eff)={sr['Allele2']} z={s_z:.4f} p={sr['N_case']}/{sr['N_ctrl']}")
            u_eff = str(r["A1"]).upper(); u_oth = str(r["OMITTED"]).upper()
            s_eff = str(sr["Allele2"]).upper(); s_oth = str(sr["Allele1"]).upper()
            # harmonize
            if {s_eff, s_oth} == {u_eff, u_oth}:
                sz = s_z if s_eff == u_eff else -s_z
            elif {comp(s_eff), comp(s_oth)} == {u_eff, u_oth}:
                sz = s_z if comp(s_eff) == u_eff else -s_z
            else:
                sz = None
            print(f"  harmonized Sonehara z (aligned to {u_eff}) = {sz}")
            if sz is not None:
                N_EFF_S = 4.0*sr["N_case"]*sr["N_ctrl"]/(sr["N_case"]+sr["N_ctrl"])
                w_u = math.sqrt(float(r["OBS_CT"])); w_s = math.sqrt(N_EFF_S)
                z_meta = (r["Z_STAT"]*w_u + sz*w_s)/math.sqrt(w_u**2 + w_s**2)
                print(f"  N_eff_sone={N_EFF_S:.1f}  w_user={w_u:.1f} w_sone={w_s:.1f}")
                print(f"  >>> combined Z = {z_meta:.4f}  (user z={r['Z_STAT']:.3f} diluted)")
        else:
            print(f"  (Sonehara has no SNP at hg19 {pos19})")
else:
    print("  top hit not found at exact pos; checking nearest...")
    sub = u[(u["CHROM"]=="2") & (u["POS"].between(174094000, 174097000))]
    print(sub[["POS","A1","Z_STAT","P"]].head())

# their top hits from Table2 (top 5)
print("\n=== user top-5 hits (from Table2) liftover + Sonehara z ===")
top5 = ["2:174095142","2:174193638","2:174134460","2:174094826","2:174190584"]
s = pd.read_csv(SONEHARA, sep="\t", usecols=["CHR","POS","Allele1","Allele2","BETA","SE"])
for name in top5:
    chrom, pos = name.split(":")
    pos19 = lift.get(name)
    r = u[(u["CHROM"]==chrom) & (u["POS"]==int(pos))]
    if not len(r) or pos19 is None:
        print(f"  {name}: no liftover/user row"); continue
    r = r.iloc[0]
    ss = s[(s["CHR"]==int(chrom)) & (s["POS"]==pos19)]
    s_z = ss.iloc[0]["BETA"]/ss.iloc[0]["SE"] if len(ss) else float("nan")
    print(f"  {name} -> hg19 {pos19}: user z={r['Z_STAT']:+.3f} (p={r['P']:.1e}) | Sonehara z={s_z:+.3f}")
