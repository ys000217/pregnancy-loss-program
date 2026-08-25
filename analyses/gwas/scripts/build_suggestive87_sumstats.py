#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Incorporate the 87 suggestive GWAS hits (P < 1e-4) into a gsMap-format sumstats.

Strategy (same constraints as prior OLA1 boost experiment):
  1. Lift each hit CHR:POS from hg38 -> hg19.
  2. Map hg19 CHR:POS -> rsID via 1000G EUR .bim (gsMap LD panel).
  3. Start from RPL_combine_gsmap.sumstats.gz background.
  4. For matched rsIDs, set |Z|=8 (chisq=64 < gsMap chisq_max=80),
     preserving the original GWAS Z sign.

Outputs:
  - RPL_suggestive87_boost.sumstats.gz
  - RPL_suggestive87_mapping_report.txt
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pandas as pd

WORK = Path("/mnt/d/ONT/figure3/gwas")
HITS = WORK / "GWAS_suggestive_hits_1e-4.csv"
BIM_DIR = Path("/mnt/d/gsMap/gsMap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink")
SUMSTATS = Path("/mnt/d/gsMap/RPL_GWAS/RPL_combine_gsmap.sumstats.gz")
OUT_SUM = Path("/mnt/d/gsMap/RPL_GWAS/RPL_suggestive87_boost.sumstats.gz")
OUT_REPORT = WORK / "RPL_suggestive87_mapping_report.txt"
LIFTOVER = "/mnt/d/gsMap/tools/liftOver"
CHAIN = "/mnt/d/gsMap/tools/hg38ToHg19.over.chain.gz"
BOOST_Z = 8.0


def main() -> None:
    hits = pd.read_csv(HITS)
    print(f"[1] suggestive hits: {len(hits)}")

    bed_in = WORK / "_sug87_hg38.bed"
    bed_out = WORK / "_sug87_hg19.bed"
    bed_unmap = WORK / "_sug87_unmapped.bed"
    with open(bed_in, "w", encoding="utf-8") as fh:
        for _, r in hits.iterrows():
            chrom = str(r["CHROM"]).replace("chr", "")
            pos = int(r["POS"])
            name = f"{chrom}:{pos}"
            # bed: 0-based half-open single-base
            fh.write(f"chr{chrom}\t{pos - 1}\t{pos}\t{name}\n")

    r = subprocess.run(
        [LIFTOVER, str(bed_in), CHAIN, str(bed_out), str(bed_unmap)],
        capture_output=True,
        text=True,
    )
    print("[liftOver]", r.stderr.strip() or "ok")

    lift = {}
    with open(bed_out, encoding="utf-8") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                chrom19 = p[0].replace("chr", "")
                pos19 = int(p[1]) + 1
                name = p[3]
                lift[name] = (chrom19, pos19)
    print(f"[2] lifted to hg19: {len(lift)} / {len(hits)}")

    print("[3] building bim CHR:POS -> rsID map ...")
    pos2rs: dict[str, str] = {}
    for i in range(1, 23):
        bf = BIM_DIR / f"1000G.EUR.QC.{i}.bim"
        if not bf.exists():
            continue
        with open(bf, encoding="utf-8") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 6:
                    continue
                key = f"{p[0]}:{p[3]}"
                if key not in pos2rs:
                    pos2rs[key] = p[1]
    print(f"    bim entries: {len(pos2rs):,}")

    hits = hits.copy()
    hits["name"] = hits["CHROM"].astype(str).str.replace("^chr", "", regex=True) + ":" + hits["POS"].astype(int).astype(str)
    hits["hg19_chrom"] = hits["name"].map(lambda n: lift.get(n, (None, None))[0])
    hits["hg19_pos"] = hits["name"].map(
        lambda n: (None if lift.get(n, (None, None))[1] is None else int(lift.get(n)[1]))
    )
    hits["rsid"] = [
        pos2rs.get(f"{c}:{int(p)}") if c is not None and pd.notna(p) else None
        for c, p in zip(hits["hg19_chrom"], hits["hg19_pos"])
    ]
    n_rs = hits["rsid"].notna().sum()
    print(f"[4] mapped to rsID: {n_rs} / {len(hits)}")

    df = pd.read_csv(SUMSTATS, sep="\t", compression="gzip")
    print(f"[5] background sumstats: {len(df):,}")

    # preserve GWAS Z sign when boosting
    sign = hits.set_index("rsid")["Z_STAT"] if "Z_STAT" in hits.columns else None
    if sign is None:
        # fall back: negative if OR<1 else positive
        hits["_sign"] = hits["OR"].apply(lambda x: -1.0 if pd.notna(x) and x < 1 else 1.0)
        sign = hits.dropna(subset=["rsid"]).set_index("rsid")["_sign"]
    else:
        sign = sign.dropna()
        sign = sign.apply(lambda z: 1.0 if float(z) >= 0 else -1.0)

    target_rs = set(hits["rsid"].dropna())
    in_sum = df["SNP"].isin(target_rs)
    n_in = int(in_sum.sum())
    print(f"[6] suggestive rsIDs present in sumstats: {n_in}")

    sign_map = sign.to_dict()
    df["Z_new"] = df["Z"]
    boosted_idx = df.index[in_sum]
    signs = df.loc[boosted_idx, "SNP"].map(lambda s: float(sign_map.get(s, 1.0)))
    df.loc[boosted_idx, "Z_new"] = signs.values * BOOST_Z

    n_changed = int((df["Z_new"] != df["Z"]).sum())
    out = df[["SNP", "A1", "A2", "Z_new", "N"]].rename(columns={"Z_new": "Z"})
    out.to_csv(OUT_SUM, sep="\t", index=False, compression="gzip")
    print(f"[7] boosted {n_changed} SNPs to |Z|={BOOST_Z} -> {OUT_SUM}")

    mean_chi2 = float((out["Z"] ** 2).mean())
    in_sum_rs = set(df.loc[in_sum, "SNP"])
    report = [
        f"n_hits={len(hits)}",
        f"n_lifted={len(lift)}",
        f"n_rsid={n_rs}",
        f"n_in_sumstats={n_in}",
        f"n_boosted={n_changed}",
        f"boost_abs_Z={BOOST_Z}",
        f"mean_chi2_after={mean_chi2:.6f}",
        f"fuel_mean_chi2_minus_1={mean_chi2 - 1:.6f}",
        "",
        "mapped_hits:",
    ]
    for _, row in hits.iterrows():
        report.append(
            f"{row['SNP']}\thg19={row['hg19_chrom']}:{row['hg19_pos']}\t"
            f"rsid={row['rsid']}\tP={row['P']}\tin_sumstats={row['rsid'] in in_sum_rs}"
        )
    OUT_REPORT.write_text("\n".join(report) + "\n", encoding="utf-8")
    # also copy sumstats next to report for convenience
    import shutil

    shutil.copy2(OUT_SUM, WORK / "gsmap_RPL_results" / "RPL_suggestive87_boost.sumstats.gz")
    print(f"[8] report -> {OUT_REPORT}")


if __name__ == "__main__":
    main()
