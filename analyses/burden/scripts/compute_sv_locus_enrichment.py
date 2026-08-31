#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Genome-wide per-locus Fisher enrichment for all PASS SVs."""

from __future__ import annotations

import csv
import os
import sys
from collections import Counter
from pathlib import Path

# Reuse helpers from compute_burden.py
sys.path.insert(0, str(Path(__file__).resolve().parent))
from compute_burden import (  # noqa: E402
    OUT_DIR,
    PHENO_PATH,
    SV_VCF,
    bh_fdr,
    fisher_exact_2x2,
    has_alt_allele,
    load_phenotype,
    open_text,
    parse_info_field,
)

OUT_DIR.mkdir(parents=True, exist_ok=True)

FDR_THRESHOLD = float(os.environ.get("SV_LOCUS_FDR", "0.05"))


def classify_enrichment_pattern(row: dict) -> str:
    """Assign locus to mutually exclusive enrichment pattern."""
    ab, norm, ctrl = row["abnormal_rate"], row["normal_rate"], row["control_rate"]
    f_ab_ctrl = row["fdr_abnormal_vs_control"]
    f_ab_norm = row["fdr_abnormal_vs_normal"]
    f_norm_ctrl = row["fdr_normal_vs_control"]

    if is_abnormal_specific(row):
        return "abnormal_specific"
    if f_ab_ctrl < FDR_THRESHOLD and ab > ctrl and norm > ctrl and f_norm_ctrl < FDR_THRESHOLD:
        return "case_vs_control"
    if f_norm_ctrl < FDR_THRESHOLD and norm > ctrl:
        return "normal_vs_control"
    if f_ab_ctrl < FDR_THRESHOLD and ab > ctrl:
        return "abnormal_vs_control_only"
    return "none"


def is_abnormal_specific(row: dict) -> bool:
    """
    Abnormal-specific: enriched in abnormal vs control AND vs normal.

    Rejects loci where normal cases carry the SV at a similar rate (case-wide signal).
    """
    if row["fdr_abnormal_vs_control"] >= FDR_THRESHOLD:
        return False
    if row["abnormal_rate"] <= row["control_rate"]:
        return False
    if row["abnormal_rate"] <= row["normal_rate"]:
        return False
    if row["fdr_abnormal_vs_normal"] >= FDR_THRESHOLD:
        return False
    return True


def annotate_patterns(rows: list[dict]) -> list[dict]:
    for row in rows:
        row["delta_abnormal_vs_normal"] = row["abnormal_rate"] - row["normal_rate"]
        row["delta_abnormal_vs_control"] = row["abnormal_rate"] - row["control_rate"]
        row["abnormal_specific"] = int(is_abnormal_specific(row))
        row["enrichment_pattern"] = classify_enrichment_pattern(row)
    return rows


def group_samples(samples: list[str], pheno_map: dict[str, str]) -> dict[str, list[str]]:
    return {
        "abnormal": [s for s in samples if pheno_map.get(s) == "abnormal"],
        "normal": [s for s in samples if pheno_map.get(s) == "normal"],
        "control": [s for s in samples if pheno_map.get(s) == "control"],
    }


def fisher_pair(carriers: dict[str, int], group_a: list[str], group_b: list[str]) -> dict:
    a = sum(carriers.get(s, 0) for s in group_a)
    b = len(group_a) - a
    c = sum(carriers.get(s, 0) for s in group_b)
    d = len(group_b) - c
    na, nb = len(group_a), len(group_b)
    return {
        "carrier_a": a,
        "total_a": na,
        "carrier_b": c,
        "total_b": nb,
        "rate_a": a / na if na else 0.0,
        "rate_b": c / nb if nb else 0.0,
        "p_value": fisher_exact_2x2(a, b, c, d),
    }


def collect_pass_sv_loci(samples: list[str], vcf_path: Path) -> list[dict]:
    loci: list[dict] = []

    with open_text(vcf_path) as f:
        vcf_samples = None
        for line in f:
            if line.startswith("#CHROM"):
                vcf_samples = line.strip().split("\t")[9:]
                break
        if vcf_samples is None:
            raise RuntimeError("SV VCF header missing")

        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[6] != "PASS":
                continue
            svtype = parse_info_field(parts[7], "SVTYPE") or "UNK"
            sv_end = parse_info_field(parts[7], "END") or parts[1]
            coord_key = f"{parts[0]}:{parts[1]}:{sv_end}:{svtype}"
            gt_map = dict(zip(vcf_samples, parts[9:]))
            carriers = {}
            for sample in samples:
                gt = gt_map.get(sample, "./.")
                if has_alt_allele(gt):
                    carriers[sample] = 1
            if not carriers:
                continue
            loci.append(
                {
                    "coord_key": coord_key,
                    "chrom": parts[0],
                    "pos": parts[1],
                    "end": sv_end,
                    "svtype": svtype,
                    "sv_id": parts[2],
                    "cohort_carrier": len(carriers),
                    "carriers": carriers,
                }
            )
    return loci


def build_enrichment_rows(loci: list[dict], groups: dict[str, list[str]]) -> list[dict]:
    rows: list[dict] = []
    n_cohort = sum(len(v) for v in groups.values())

    for loc in loci:
        carriers = loc["carriers"]
        ab_vs_ctrl = fisher_pair(carriers, groups["abnormal"], groups["control"])
        ab_vs_norm = fisher_pair(carriers, groups["abnormal"], groups["normal"])
        norm_vs_ctrl = fisher_pair(carriers, groups["normal"], groups["control"])
        rows.append(
            {
                "coord_key": loc["coord_key"],
                "chrom": loc["chrom"],
                "pos": loc["pos"],
                "end": loc["end"],
                "svtype": loc["svtype"],
                "sv_id": loc["sv_id"],
                "cohort_carrier": loc["cohort_carrier"],
                "cohort_rate": loc["cohort_carrier"] / n_cohort if n_cohort else 0.0,
                "abnormal_carrier": ab_vs_ctrl["carrier_a"],
                "abnormal_total": ab_vs_ctrl["total_a"],
                "abnormal_rate": ab_vs_ctrl["rate_a"],
                "normal_carrier": ab_vs_norm["carrier_b"],
                "normal_total": ab_vs_norm["total_b"],
                "normal_rate": ab_vs_norm["rate_b"],
                "control_carrier": ab_vs_ctrl["carrier_b"],
                "control_total": ab_vs_ctrl["total_b"],
                "control_rate": ab_vs_ctrl["rate_b"],
                "p_abnormal_vs_control": ab_vs_ctrl["p_value"],
                "p_abnormal_vs_normal": ab_vs_norm["p_value"],
                "p_normal_vs_control": norm_vs_ctrl["p_value"],
            }
        )

    for p_col, fdr_col in [
        ("p_abnormal_vs_control", "fdr_abnormal_vs_control"),
        ("p_abnormal_vs_normal", "fdr_abnormal_vs_normal"),
        ("p_normal_vs_control", "fdr_normal_vs_control"),
    ]:
        qvals = bh_fdr([r[p_col] for r in rows])
        for row, q in zip(rows, qvals):
            row[fdr_col] = q

    rows.sort(key=lambda r: r["p_abnormal_vs_control"])
    return annotate_patterns(rows)


def write_tsv(rows: list[dict], path: Path) -> None:
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def write_top_hits(rows: list[dict], p_col: str, fdr_col: str, path: Path, n: int = 50) -> None:
    top = sorted(rows, key=lambda r: r[p_col])[:n]
    write_tsv(top, path)


def main() -> None:
    print("[1/3] Loading phenotype...")
    pheno_rows = load_phenotype(PHENO_PATH)
    samples = [r["VCF_Sample_ID"] for r in pheno_rows]
    pheno_map = {r["VCF_Sample_ID"]: r["Condition"] for r in pheno_rows}
    groups = group_samples(samples, pheno_map)
    print(
        f"  samples: {len(samples)} "
        f"(abnormal={len(groups['abnormal'])}, normal={len(groups['normal'])}, control={len(groups['control'])})"
    )

    print("[2/3] Scanning PASS SV loci...")
    loci = collect_pass_sv_loci(samples, SV_VCF)
    print(f"  unique PASS SV loci with carriers: {len(loci):,}")

    print("[3/3] Fisher enrichment (abnormal vs normal vs control)...")
    rows = build_enrichment_rows(loci, groups)

    ab_specific = [r for r in rows if r["abnormal_specific"] == 1]
    ab_specific.sort(key=lambda r: r["p_abnormal_vs_control"])

    all_path = OUT_DIR / "sv_locus_enrichment_all_pass.tsv"
    write_tsv(rows, all_path)
    print(f"  wrote {all_path}")

    ab_spec_path = OUT_DIR / "sv_locus_enrichment_abnormal_specific.tsv"
    write_tsv(ab_specific, ab_spec_path)
    print(f"  abnormal-specific loci: {len(ab_specific)}")

    sig_rows = [
        r
        for r in rows
        if min(
            r["fdr_abnormal_vs_control"],
            r["fdr_abnormal_vs_normal"],
            r["fdr_normal_vs_control"],
        )
        < FDR_THRESHOLD
    ]
    sig_path = OUT_DIR / "sv_locus_enrichment_fdr05_any_comparison.tsv"
    write_tsv(sig_rows, sig_path)
    print(f"  FDR<{FDR_THRESHOLD} in any pairwise comparison: {len(sig_rows)}")

    write_top_hits(
        ab_specific if ab_specific else rows,
        "p_abnormal_vs_control",
        "fdr_abnormal_vs_control",
        OUT_DIR / "sv_locus_enrichment_top50_abnormal_vs_control.tsv",
    )
    write_top_hits(
        rows,
        "p_abnormal_vs_normal",
        "fdr_abnormal_vs_normal",
        OUT_DIR / "sv_locus_enrichment_top50_abnormal_vs_normal.tsv",
    )
    write_top_hits(
        rows,
        "p_normal_vs_control",
        "fdr_normal_vs_control",
        OUT_DIR / "sv_locus_enrichment_top50_normal_vs_control.tsv",
    )

    summary_path = OUT_DIR / "sv_locus_enrichment_summary.txt"
    pattern_counts = Counter(r["enrichment_pattern"] for r in rows)
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"pass_sv_loci\t{len(loci)}\n")
        f.write(f"fdr_threshold\t{FDR_THRESHOLD}\n")
        f.write(f"abnormal_specific_loci\t{len(ab_specific)}\n")
        f.write(
            "abnormal_specific_rule\t"
            "fdr(ab vs control)<thr AND ab_rate>control AND ab_rate>normal "
            "AND fdr(ab vs normal)<thr\n"
        )
        f.write(f"fdr05_any_pairwise\t{len(sig_rows)}\n")
        for pat in sorted(pattern_counts):
            f.write(f"pattern_{pat}\t{pattern_counts[pat]}\n")
        for label, p_col, fdr_col in [
            ("abnormal_vs_control", "p_abnormal_vs_control", "fdr_abnormal_vs_control"),
            ("abnormal_vs_normal", "p_abnormal_vs_normal", "fdr_abnormal_vs_normal"),
            ("normal_vs_control", "p_normal_vs_control", "fdr_normal_vs_control"),
        ]:
            n_fdr = sum(1 for r in rows if r[fdr_col] < FDR_THRESHOLD)
            f.write(f"fdr05_{label}\t{n_fdr}\n")
        if ab_specific:
            top = ab_specific[0]
            f.write(
                f"top_abnormal_specific\t{top['coord_key']}\t{top['sv_id']}\t"
                f"ab={top['abnormal_rate']:.3f}\tnormal={top['normal_rate']:.3f}\t"
                f"ctrl={top['control_rate']:.3f}\t"
                f"p_ab_ctrl={top['p_abnormal_vs_control']:.4g}\t"
                f"p_ab_norm={top['p_abnormal_vs_normal']:.4g}\n"
            )
    print(f"  wrote {summary_path}")
    print("Done.")


if __name__ == "__main__":
    main()
