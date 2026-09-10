#!/usr/bin/env python3
"""Aggregate gene-level Arnadottir matches and Fisher tests."""
from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from scan_arnadottir_genes_annovar import CANON, T9  # noqa: E402

PHENO = ROOT / "tables" / "suspect_abn" / "phenotype.tsv"
SITES = ROOT / "tables" / "external_arnadottir2025" / "arnadottir_gene_sites_in_648vcf.tsv"
LONG = ROOT / "tables" / "external_arnadottir2025" / "arnadottir_gene_carriers_long.tsv"
OUT = ROOT / "tables" / "external_arnadottir2025"


def fisher_two_sided(a: int, b: int, c: int, d: int) -> float:
    n = a + b + c + d
    if n == 0:
        return 1.0
    row1 = a + b
    col1 = a + c
    lo = max(0, row1 + col1 - n)
    hi = min(row1, col1)

    def dens(x: int) -> float:
        return math.comb(col1, x) * math.comb(n - col1, row1 - x) / math.comb(n, row1)

    p_obs = dens(a)
    s = 0.0
    for x in range(lo, hi + 1):
        px = dens(x)
        if px <= p_obs + 1e-15:
            s += px
    return min(1.0, s)


def fisher_or(a: int, b: int, c: int, d: int) -> float:
    if 0 in (a, b, c, d):
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    return (a / b) / (c / d)


def bh_fdr(pvals: list[float]) -> list[float]:
    m = len(pvals)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for k in range(m, 0, -1):
        i = order[k - 1]
        prev = min(prev, pvals[i] * m / k)
        q[i] = min(1.0, prev)
    return q


def load_pheno() -> dict[str, str]:
    out = {}
    with PHENO.open(encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for rec in r:
            out[rec["VCF_Sample_ID"]] = rec["Condition"].strip().lower()
    return out


def is_plp(clnsig: str) -> bool:
    s = (clnsig or "").replace(" ", "_")
    return "Pathogenic" in s or "Likely_pathogenic" in s


def site_class(rec: dict) -> dict[str, bool]:
    cons = rec.get("consequence", "")
    rare = rec.get("pop_af_lt1pct") == "yes" and rec.get("clinvar_ok") == "yes"
    return {
        "coding_rare": rare and rec.get("in_vcf") in {True, "TRUE", "True"} and int(rec.get("n_carriers") or 0) > 0,
        "lof_splice": cons in {"lof", "splicing"} and rec.get("pop_af_lt1pct") == "yes"
        and rec.get("clinvar_ok") == "yes"
        and rec.get("in_vcf") in {True, "TRUE", "True"}
        and int(rec.get("n_carriers") or 0) > 0,
        "plp": is_plp(rec.get("clnsig", ""))
        and rec.get("in_vcf") in {True, "TRUE", "True"}
        and int(rec.get("n_carriers") or 0) > 0,
    }


def main() -> None:
    pheno = load_pheno()
    n_ab = sum(1 for v in pheno.values() if v == "abnormal")
    n_nm = sum(1 for v in pheno.values() if v == "normal")
    n_ct = sum(1 for v in pheno.values() if v == "control")

    sites = list(csv.DictReader(SITES.open(encoding="utf-8"), delimiter="\t"))
    for rec in sites:
        rec["in_vcf"] = rec.get("in_vcf", "")
        rec["n_carriers"] = int(float(rec["n_carriers"] or 0))

    long_rows = list(csv.DictReader(LONG.open(encoding="utf-8"), delimiter="\t"))

    # gene -> class -> set of samples
    sets: dict[str, dict[str, dict[str, set[str]]]] = defaultdict(
        lambda: {k: {"abnormal": set(), "normal": set(), "control": set()} for k in ("coding_rare", "lof_splice", "plp")}
    )
    site_keep = []
    for rec in sites:
        flags = site_class(rec)
        g = rec["canon_gene"]
        ids = [x for x in (rec.get("carrier_ids") or "").split(",") if x]
        for sid in ids:
            cond = pheno.get(sid) or pheno.get(sid.lower())
            if cond not in {"abnormal", "normal", "control"}:
                continue
            for cl, ok in flags.items():
                if ok:
                    sets[g][cl][cond].add(sid)
        if any(flags.values()) or rec["n_carriers"] > 0:
            rec2 = dict(rec)
            rec2["class_coding_rare"] = "yes" if flags["coding_rare"] else "no"
            rec2["class_lof_splice"] = "yes" if flags["lof_splice"] else "no"
            rec2["class_plp"] = "yes" if flags["plp"] else "no"
            site_keep.append(rec2)

    rows_out = []
    for cl in ("coding_rare", "lof_splice", "plp"):
        pvals = []
        block = []
        for g in sorted(CANON):
            ab = len(sets[g][cl]["abnormal"])
            nm = len(sets[g][cl]["normal"])
            ct = len(sets[g][cl]["control"])
            a, b = ab, n_ab - ab
            c, d = nm, n_nm - nm
            p = fisher_two_sided(a, b, c, d)
            pvals.append(p)
            block.append(
                {
                    "class": cl,
                    "gene": g,
                    "euploid_likely_causal": "yes" if g in T9 else "no",
                    "n_ab_carrier": ab,
                    "n_ab_total": n_ab,
                    "n_normal_carrier": nm,
                    "n_normal_total": n_nm,
                    "n_control_carrier": ct,
                    "n_control_total": n_ct,
                    "or_ab_vs_normal": round(fisher_or(a, b, c, d), 4),
                    "p_fisher_ab_vs_normal": p,
                }
            )
        q = bh_fdr(pvals)
        for rec, qi in zip(block, q):
            rec["fdr_bh"] = qi
            rows_out.append(rec)

    gpath = OUT / "arnadottir_gene_burden.tsv"
    with gpath.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)

    spath = OUT / "arnadottir_gene_sites_classified.tsv"
    if site_keep:
        with spath.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(site_keep[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(site_keep)

    print("N_PHENO", n_ab, n_nm, n_ct)
    print("SITES", len(sites), "with_vcf", sum(1 for x in sites if str(x.get("in_vcf")).upper() == "TRUE"))
    print("WROTE", gpath)
    # print FDR<0.05
    hits = [r for r in rows_out if r["fdr_bh"] < 0.05]
    print("FDR05", len(hits))
    for r in hits:
        print(r["class"], r["gene"], r["n_ab_carrier"], r["n_normal_carrier"], r["p_fisher_ab_vs_normal"], r["fdr_bh"])


if __name__ == "__main__":
    main()
