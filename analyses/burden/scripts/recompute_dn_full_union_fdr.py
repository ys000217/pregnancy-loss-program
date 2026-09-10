#!/usr/bin/env python3
"""Recompute full-union FDR from existing DN site/gene count tables (no VCF)."""
from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_dn_ab_vs_normal import (  # noqa: E402
    OUT,
    apply_full_union_fdr,
    fisher_one_sided_row1_greater,
    fisher_or,
    fisher_two_sided,
)

N_AB = 77
N_NORM = 418


def recompute_row(n_ab: int, n_norm: int) -> dict:
    a, b = n_ab, N_AB - n_ab
    c, d = n_norm, N_NORM - n_norm
    return {
        "rate_abnormal": n_ab / N_AB,
        "rate_normal": n_norm / N_NORM,
        "or_ab_vs_normal": fisher_or(a, b, c, d),
        "p_fisher_two_sided": fisher_two_sided(a, b, c, d),
        "p_ab_greater": fisher_one_sided_row1_greater(a, b, c, d),
        "p_norm_greater": fisher_one_sided_row1_greater(c, d, a, b),
    }


def load_tsv(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    site_path = OUT / "dn_site_patterns_ab_vs_normal.tsv"
    gene_path = OUT / "dn_gene_burden_ab_vs_normal.tsv"
    sites = load_tsv(site_path)
    genes = load_tsv(gene_path)

    # drop obsolete directional columns if present
    drop_keys = {
        "fdr_ab_enriched",
        "fdr_norm_enriched",
        "in_ab_fdr_family",
        "in_norm_fdr_family",
        "fdr",
        "p_fisher",
    }

    new_sites = []
    for r in sites:
        n_ab = int(r["n_abnormal"])
        n_norm = int(r["n_normal"])
        stats = recompute_row(n_ab, n_norm)
        out = {k: v for k, v in r.items() if k not in drop_keys and k not in stats}
        # also drop old fdr_* that we will replace
        for k in list(out):
            if k.startswith("fdr_"):
                del out[k]
        out.update(stats)
        out["n_abnormal"] = n_ab
        out["n_normal"] = n_norm
        new_sites.append(out)

    apply_full_union_fdr(new_sites, n_ab_key="n_abnormal", n_norm_key="n_normal")
    new_sites.sort(
        key=lambda r: (
            min(r["fdr_two_sided"], r["fdr_ab_greater"], r["fdr_norm_greater"]),
            r["p_fisher_two_sided"],
        )
    )
    write_tsv(site_path, new_sites)

    new_genes = []
    for r in genes:
        n_ab = int(r["n_abnormal_carriers"])
        n_norm = int(r["n_normal_carriers"])
        stats = recompute_row(n_ab, n_norm)
        out = {k: v for k, v in r.items() if k not in drop_keys and k not in stats}
        for k in list(out):
            if k.startswith("fdr_"):
                del out[k]
        out.update(stats)
        out["n_abnormal_carriers"] = n_ab
        out["n_normal_carriers"] = n_norm
        new_genes.append(out)

    apply_full_union_fdr(
        new_genes, n_ab_key="n_abnormal_carriers", n_norm_key="n_normal_carriers"
    )
    new_genes.sort(
        key=lambda r: (
            min(r["fdr_two_sided"], r["fdr_ab_greater"], r["fdr_norm_greater"]),
            r["p_fisher_two_sided"],
        )
    )
    write_tsv(gene_path, new_genes)

    site_pat = Counter(r["pattern"] for r in new_sites)
    gene_pat = Counter(r["pattern"] for r in new_genes)

    def n05(rows, key):
        return sum(1 for r in rows if float(r[key]) < 0.05)

    sum_path = OUT / "dn_ab_vs_normal_summary.tsv"
    # keep previous non-FDR metrics if present
    prev = {}
    if sum_path.exists():
        with open(sum_path, encoding="utf-8") as f:
            next(f)
            for line in f:
                k, v = line.rstrip("\n").split("\t", 1)
                prev[k] = v

    with open(sum_path, "w", encoding="utf-8") as f:
        f.write("metric\tvalue\n")
        f.write(f"phenotype_file\t{prev.get('phenotype_file', 'suspect_abn/phenotype.tsv')}\n")
        f.write(f"n_abnormal\t{N_AB}\n")
        f.write(f"n_normal\t{N_NORM}\n")
        f.write("fdr_method\tfull_union_BH_two_sided_and_two_one_sided\n")
        f.write(
            "fdr_note\t"
            "BH on full union separately for two-sided, one-sided ab higher, "
            "one-sided normal higher; direction does not shrink m\n"
        )
        f.write(f"site_fdr_m\t{len(new_sites)}\n")
        f.write(f"gene_fdr_m\t{len(new_genes)}\n")
        f.write(f"site_n_fdr05_two\t{n05(new_sites, 'fdr_two_sided')}\n")
        f.write(f"site_n_fdr05_ab\t{n05(new_sites, 'fdr_ab_greater')}\n")
        f.write(f"site_n_fdr05_norm\t{n05(new_sites, 'fdr_norm_greater')}\n")
        f.write(f"gene_n_fdr05_two\t{n05(new_genes, 'fdr_two_sided')}\n")
        f.write(f"gene_n_fdr05_ab\t{n05(new_genes, 'fdr_ab_greater')}\n")
        f.write(f"gene_n_fdr05_norm\t{n05(new_genes, 'fdr_norm_greater')}\n")
        f.write(f"n_sites_with_ab_or_norm_carrier\t{len(new_sites)}\n")
        f.write(f"n_genes\t{len(new_genes)}\n")
        for k, v in site_pat.most_common():
            f.write(f"site_pattern_{k}\t{v}\n")
        for k, v in gene_pat.most_common():
            f.write(f"gene_pattern_{k}\t{v}\n")
        for k in [
            "mw_p_n_dn_sites",
            "mw_p_n_dn_genes",
            "abnormal_sites_median",
            "normal_sites_median",
            "normal_frac_sites_gt_abnormal_median",
            "n_carrier_rows",
        ]:
            if k in prev:
                f.write(f"{k}\t{prev[k]}\n")

    tp = OUT / "dn_ab_vs_normal_talking_points.txt"
    with open(tp, "w", encoding="utf-8") as f:
        f.write("DN-background SNV: abnormal vs normal (raw columns)\n")
        f.write("=" * 60 + "\n")
        f.write(
            "Valid inference: full-union BH for (1) two-sided difference, "
            "(2) one-sided ab higher, (3) one-sided normal higher.\n"
        )
        f.write("Directional subset FDR was removed (selective inference / invalid).\n")
        f.write(f"Phenotype: abnormal={N_AB}, normal={N_NORM}\n")
        f.write(f"FDR m: sites={len(new_sites)}, genes={len(new_genes)}\n")
        f.write(
            f"Site FDR<0.05 two/ab/norm="
            f"{n05(new_sites,'fdr_two_sided')}/{n05(new_sites,'fdr_ab_greater')}/{n05(new_sites,'fdr_norm_greater')}\n"
        )
        f.write(
            f"Gene FDR<0.05 two/ab/norm="
            f"{n05(new_genes,'fdr_two_sided')}/{n05(new_genes,'fdr_ab_greater')}/{n05(new_genes,'fdr_norm_greater')}\n"
        )
        f.write(f"Site patterns: {dict(site_pat)}\n")
        f.write(f"Gene patterns: {dict(gene_pat)}\n")

    # short significant lists
    for label, rows, keys in [
        (
            "sites",
            new_sites,
            ("fdr_two_sided", "fdr_ab_greater", "fdr_norm_greater"),
        ),
        (
            "genes",
            new_genes,
            ("fdr_two_sided", "fdr_ab_greater", "fdr_norm_greater"),
        ),
    ]:
        hits = []
        for r in rows:
            flags = []
            if float(r["fdr_two_sided"]) < 0.05:
                flags.append("two")
            if float(r["fdr_ab_greater"]) < 0.05:
                flags.append("ab")
            if float(r["fdr_norm_greater"]) < 0.05:
                flags.append("norm")
            if flags:
                hits.append((r, flags))
        outp = OUT / f"dn_{label}_fdr05_hits.tsv"
        if hits:
            base = dict(hits[0][0])
            base["sig_tests"] = ""
            with open(outp, "w", encoding="utf-8", newline="") as f:
                w = csv.DictWriter(
                    f, fieldnames=list(base.keys()), delimiter="\t"
                )
                w.writeheader()
                for r, flags in hits:
                    rr = dict(r)
                    rr["sig_tests"] = ",".join(flags)
                    w.writerow(rr)
        else:
            outp.write_text("no_hits\n", encoding="utf-8")
        print(f"{label} FDR<0.05 hits: {len(hits)} -> {outp}")

    print("Site patterns", dict(site_pat))
    print("Gene patterns", dict(gene_pat))
    print("DONE")


if __name__ == "__main__":
    main()
