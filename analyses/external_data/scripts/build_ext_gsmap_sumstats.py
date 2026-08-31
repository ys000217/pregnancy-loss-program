#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Split analyses/external_data GWAS hits into two classes and inject into gsMap sumstats:

  A) miscarriage / RPL  (Laisk, Sonehara, Fan, Liu2024, ...)
  B) pregnancy phenotypes (Liu2026 gestational traits; not miscarriage case-control)

Background = RPL_combine_gsmap.sumstats.gz (polygenic scaffold).
Hits present in background are boosted to |Z|=8.
Hits absent from background but present in 1000G EUR HapMap3 BIM are appended
with alleles from BIM and |Z|=8 (so LDSC can see them if they fall in weights).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/mnt/d/ONT/analyses/external_data")
HITS = ROOT / "metadata/gwas_hits.tsv"
OUT_DIR = ROOT / "results/gsmap_external"
SUMSTATS = Path("/mnt/d/gsMap/RPL_GWAS/RPL_combine_gsmap.sumstats.gz")
BIM_DIR = Path("/mnt/d/gsMap/gsMap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink")
WEIGHTS = Path("/mnt/d/gsMap/gsMap_resource/LDSC_resource/weights_hm3_no_hla")
BOOST_Z = 8.0
N_DEFAULT = 10000


def class_of(row: pd.Series) -> str:
    if str(row.get("study_short", "")).startswith("Liu2026") or str(
        row.get("hit_id", "")
    ).startswith("LIU2026"):
        return "pregnancy"
    return "miscarriage"


def effect_sign(row: pd.Series) -> float:
    ev = row.get("effect_value")
    try:
        v = float(ev)
        if np.isfinite(v):
            metric = str(row.get("effect_metric", "")).lower()
            if "or" in metric:
                return 1.0 if v >= 1 else -1.0
            return 1.0 if v >= 0 else -1.0
    except (TypeError, ValueError):
        pass
    return 1.0


def load_bim_alleles(rsids: set[str]) -> dict[str, tuple[str, str]]:
    """SNP -> (A1, A2) from 1000G EUR Phase3 BIM (A1=allele1 col5, A2=allele2 col6)."""
    out: dict[str, tuple[str, str]] = {}
    for i in range(1, 23):
        path = BIM_DIR / f"1000G.EUR.QC.{i}.bim"
        with open(path, encoding="utf-8") as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                snp = p[1]
                if snp in rsids and snp not in out:
                    out[snp] = (p[4], p[5])  # A1, A2
    return out


def load_hm3_weight_snps(rsids: set[str]) -> set[str]:
    found: set[str] = set()
    for i in range(1, 23):
        p = WEIGHTS / f"weights.{i}.l2.ldscore.gz"
        df = pd.read_csv(p, sep="\t", compression="gzip", usecols=["SNP"])
        found |= rsids & set(df["SNP"])
    return found


def boost_class(
    name: str,
    hits: pd.DataFrame,
    bg: pd.DataFrame,
    bim: dict[str, tuple[str, str]],
    in_weights: set[str],
) -> dict:
    hits = hits.copy()
    hits["rsid"] = hits["rsid"].astype(str).str.strip()
    chrom = hits["chrom"].astype(str).str.replace("^chr", "", regex=True)
    auto = chrom.isin([str(i) for i in range(1, 23)])
    hits = hits.loc[auto].drop_duplicates(subset=["rsid"])
    hits["_sign"] = hits.apply(effect_sign, axis=1)

    target = set(hits["rsid"])
    sign_map = hits.set_index("rsid")["_sign"].to_dict()
    bg_snps = set(bg["SNP"])

    out = bg.copy()
    # boost existing
    in_bg = out["SNP"].isin(target)
    n_in = int(in_bg.sum())
    idx = out.index[in_bg]
    signs = out.loc[idx, "SNP"].map(lambda s: float(sign_map.get(s, 1.0)))
    out.loc[idx, "Z"] = signs.values * BOOST_Z
    n_boosted_existing = n_in

    # append missing that are in BIM
    missing = [s for s in target if s not in bg_snps and s in bim]
    rows = []
    for snp in missing:
        a1, a2 = bim[snp]
        rows.append(
            {
                "SNP": snp,
                "A1": a1,
                "A2": a2,
                "Z": float(sign_map.get(snp, 1.0)) * BOOST_Z,
                "N": N_DEFAULT,
            }
        )
    n_appended = len(rows)
    if rows:
        out = pd.concat([out, pd.DataFrame(rows)], ignore_index=True)

    if out["N"].isna().any():
        out["N"] = out["N"].fillna(N_DEFAULT).astype(int)

    final = out[["SNP", "A1", "A2", "Z", "N"]]
    sum_path = OUT_DIR / f"EXT_{name}_boost.sumstats.gz"
    final.to_csv(sum_path, sep="\t", index=False, compression="gzip")
    gsmap_copy = Path(f"/mnt/d/gsMap/RPL_GWAS/EXT_{name}_boost.sumstats.gz")
    final.to_csv(gsmap_copy, sep="\t", index=False, compression="gzip")

    hit_list = OUT_DIR / f"EXT_{name}_hits.tsv"
    hits.to_csv(hit_list, sep="\t", index=False)

    in_bim = sorted(s for s in target if s in bim)
    in_w = sorted(s for s in target if s in in_weights)
    not_mappable = sorted(s for s in target if s not in bim)

    mean_chi2 = float((final["Z"] ** 2).mean())
    report = {
        "class": name,
        "n_hits_table": int(len(hits)),
        "n_unique_rsid": len(target),
        "n_in_sumstats_before": n_in,
        "n_boosted_existing": n_boosted_existing,
        "n_appended_from_bim": n_appended,
        "n_injectable": n_boosted_existing + n_appended,
        "n_in_hm3_weights": len(in_w),
        "rsids_in_bim": in_bim,
        "rsids_in_weights": in_w,
        "rsids_not_in_bim": not_mappable,
        "mean_chi2": mean_chi2,
        "fuel": mean_chi2 - 1,
        "sumstats": str(sum_path),
        "phenotypes": sorted(hits["phenotype"].astype(str).unique().tolist()),
        "studies": sorted(hits["study_short"].astype(str).unique().tolist()),
    }
    return report


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(HITS, sep="\t")
    df["_class"] = df.apply(class_of, axis=1)
    print("class counts:")
    print(df.groupby("_class").size())
    print("by study:")
    print(df.groupby(["_class", "study_short"]).size())

    all_rs = set(df["rsid"].astype(str).str.strip())
    print(f"\nloading BIM alleles for {len(all_rs)} rsIDs...")
    bim = load_bim_alleles(all_rs)
    print(f"  in 1000G EUR BIM (chr1-22): {len(bim)}")
    in_w = load_hm3_weight_snps(all_rs)
    print(f"  in HapMap3 LDSC weights: {len(in_w)}")

    bg = pd.read_csv(SUMSTATS, sep="\t", compression="gzip")
    print(f"background sumstats: {len(bg):,}")

    reports = []
    for name in ["miscarriage", "pregnancy"]:
        sub = df[df["_class"] == name]
        r = boost_class(name, sub, bg, bim, in_w)
        reports.append(r)
        print(
            f"\n[{name}] hits={r['n_hits_table']} "
            f"existing={r['n_boosted_existing']} appended={r['n_appended_from_bim']} "
            f"injectable={r['n_injectable']} in_weights={r['n_in_hm3_weights']} "
            f"mean_chi2={r['mean_chi2']:.4f} fuel={r['fuel']:.4f}"
        )
        print(f"  in_weights rsIDs: {r['rsids_in_weights']}")
        print(f"  not_in_bim: {r['rsids_not_in_bim']}")
        print(f"  studies={r['studies']}")

    rep_path = OUT_DIR / "EXT_two_class_build_report.txt"
    lines = []
    for r in reports:
        lines.append(
            f"{r['class']}\tn_hits={r['n_hits_table']}\t"
            f"n_boosted_existing={r['n_boosted_existing']}\t"
            f"n_appended={r['n_appended_from_bim']}\t"
            f"n_injectable={r['n_injectable']}\t"
            f"n_in_hm3_weights={r['n_in_hm3_weights']}\t"
            f"mean_chi2={r['mean_chi2']:.6f}\tfuel={r['fuel']:.6f}"
        )
        lines.append(f"  weights: {','.join(r['rsids_in_weights']) or '(none)'}")
        lines.append(f"  not_in_bim: {','.join(r['rsids_not_in_bim']) or '(none)'}")
    rep_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nwrote {rep_path}")


if __name__ == "__main__":
    main()
