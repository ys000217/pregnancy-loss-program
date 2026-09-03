#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Storyline figures: genome-wide burden, per-locus SV, EpiFactors, gene maps.

Uses Nature-style defaults (Arial 5–7 pt, 180 mm, Okabe–Ito).
Run:
  analyses/burden/.venv/Scripts/python.exe scripts/plot_narrative.py
"""

from __future__ import annotations

import csv
import json
import sys
import math
import re
import shutil
import time
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))
import epifactors_catalog as epi  # noqa: E402

from compute_burden import (  # noqa: E402
    PHENO_PATH as DEFAULT_PHENO,
    SV_ANNO,
    keep_gw_8_10,
    load_phenotype,
    load_plp_sv_records,
)
import compute_sv_locus_enrichment as locsv  # noqa: E402

MODULE = Path(__file__).resolve().parents[1]
TAB = Path(epi.os.environ.get("BURDEN_OUT", str(MODULE / "tables")))
# gw8_10: original 8–10 week narrative. all_suspect: full 648, clustering suspects → abnormal.
NARRATIVE = epi.os.environ.get("BURDEN_NARRATIVE", "gw8_10")
if NARRATIVE == "all_suspect":
    PLOT = MODULE / "plots" / "narrative_all_suspect"
    LOCUS_TAB = TAB / "suspect_abn"
    EPI_TAB = LOCUS_TAB / "epifactors"
    PHENO_PATH = LOCUS_TAB / "phenotype.tsv"
    FILTER_GW810 = False
else:
    PLOT = Path(epi.os.environ.get("BURDEN_PLOT", str(MODULE / "plots"))) / "narrative"
    LOCUS_TAB = TAB / "gw8_10"
    EPI_TAB = TAB / "gw8_10" / "epifactors"
    PHENO_PATH = DEFAULT_PHENO
    FILTER_GW810 = True
EPI_PLOT = PLOT / "epifactors"
CODING_MB = 85.77348

GROUPS = epi.GROUPS
GROUP_COLOR = epi.GROUP_COLOR
GROUP_LABEL = epi.GROUP_LABEL
MM = epi.MM

N_FIG4_GENES = 6
ENSEMBL = "https://rest.ensembl.org"


def rank_rare_snv_concentration(snv_long: list[dict], min_ab: int = 3) -> tuple[list[str], list[dict]]:
    """HHI of rare-SNV sites among abnormal carriers (higher = more concentrated on few sites)."""
    sites = defaultdict(lambda: defaultdict(set))
    people = defaultdict(set)
    for r in snv_long:
        if r.get("Condition") != "abnormal":
            continue
        if str(r.get("pop_rare")) not in {"1", "1.0", 1}:
            continue
        key = (r["chrom"], r["pos"], r["ref"], r["alt"])
        sites[r["gene"]][key].add(r["Sample_ID"])
        people[r["gene"]].add(r["Sample_ID"])
    ranked = []
    for gene, by_site in sites.items():
        n_ab = len(people[gene])
        if n_ab < min_ab:
            continue
        counts = [len(s) for s in by_site.values()]
        total = sum(counts)
        hhi = sum((c / total) ** 2 for c in counts) if total else 0.0
        ranked.append(
            {
                "gene": gene,
                "hhi": hhi,
                "n_abnormal": n_ab,
                "n_sites": len(by_site),
                "max_site_abnormal": max(counts),
                "max_site_share": max(counts) / total if total else 0.0,
            }
        )
    ranked.sort(key=lambda r: (-r["hhi"], -r["n_abnormal"], -r["max_site_abnormal"], r["gene"]))
    genes = [r["gene"] for r in ranked[:N_FIG4_GENES]]
    return genes, ranked


def rank_functional_sv_genes(sv_long: list[dict]) -> tuple[list[str], list[dict]]:
    """Genes with exon/splice or ACMG>=4 SV in abnormal; ranked by abnormal functional carriers."""
    people = defaultdict(lambda: {g: set() for g in GROUPS})
    people_exon = defaultdict(lambda: {g: set() for g in GROUPS})
    people_acmg = defaultdict(lambda: {g: set() for g in GROUPS})
    for r in sv_long:
        try:
            acmg = int(float(r.get("acmg_class") or 0))
        except ValueError:
            acmg = 0
        exon = inum(r.get("exon_or_splice")) == 1
        if not exon and acmg < 4:
            continue
        cond = r.get("Condition", "")
        if cond not in GROUPS:
            continue
        gene = r["gene"]
        people[gene][cond].add(r["Sample_ID"])
        if exon:
            people_exon[gene][cond].add(r["Sample_ID"])
        if acmg >= 4:
            people_acmg[gene][cond].add(r["Sample_ID"])
    ranked = []
    n = group_sizes()
    for gene, by in people.items():
        n_ab = len(by["abnormal"])
        if n_ab < 1:
            continue
        ranked.append(
            {
                "gene": gene,
                "func_abnormal_n": n_ab,
                "func_normal_n": len(by["normal"]),
                "func_control_n": len(by["control"]),
                "exon_abnormal_n": len(people_exon[gene]["abnormal"]),
                "acmg4_abnormal_n": len(people_acmg[gene]["abnormal"]),
                "func_abnormal_rate": n_ab / n["abnormal"],
                "func_normal_rate": len(by["normal"]) / n["normal"],
                "func_control_rate": len(by["control"]) / n["control"],
            }
        )
    ranked.sort(
        key=lambda r: (
            -r["exon_abnormal_n"],
            -r["acmg4_abnormal_n"],
            -r["func_abnormal_n"],
            r["gene"],
        )
    )
    genes = [r["gene"] for r in ranked[:8]]
    return genes, ranked


def group_sizes() -> dict[str, int]:
    rows = load_phenotype(PHENO_PATH)
    if FILTER_GW810:
        rows = [r for r in rows if keep_gw_8_10(r.get("Gestational_Week", ""))]
    c = Counter(r["Condition"] for r in rows)
    return {g: int(c.get(g, 0)) for g in GROUPS}


def apply_pheno_condition(rows: list[dict], id_field: str = "Sample_ID") -> list[dict]:
    pheno_map = {r["VCF_Sample_ID"]: r["Condition"] for r in load_phenotype(PHENO_PATH)}
    out = []
    for r in rows:
        sid = r.get(id_field, "")
        if sid in pheno_map:
            r = dict(r)
            r["Condition"] = pheno_map[sid]
        out.append(r)
    return out


def fnum(x, default=0.0) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def inum(x, default=0) -> int:
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return default


def sv_span(rec: dict) -> tuple[int, int]:
    a = inum(rec.get("start"))
    b = inum(rec.get("end") or rec.get("start"))
    return (min(a, b), max(a, b))


def sv_is_insertion(rec: dict) -> bool:
    svt = str(rec.get("svtype") or "").upper()
    s, e = sv_span(rec)
    return svt == "INS" or s == e or inum(rec.get("svlen")) == 0


def overlaps_interval(a0: int, a1: int, b0: int, b1: int, pad: int = 2) -> bool:
    return a0 <= b1 + pad and a1 >= b0 - pad


def sv_overlaps_drawn_exons(rec: dict, model: dict, pad: int = 2) -> bool:
    """True only if the SV interval hits an exon box drawn on this track (MANE/canonical)."""
    s, e = sv_span(rec)
    for ex in model.get("exons") or []:
        if overlaps_interval(s, e, int(ex["start"]), int(ex["end"]), pad=pad):
            return True
    return False


def save_fig(fig, stem: str, pad_inches: float = 0.0) -> None:
    PLOT.mkdir(parents=True, exist_ok=True)
    png = PLOT / f"{stem}.png"
    kw = {"facecolor": "white"}
    if pad_inches > 0:
        kw["bbox_inches"] = "tight"
        kw["pad_inches"] = pad_inches
    fig.savefig(png, **kw)
    pdf = PLOT / f"{stem}.pdf"
    try:
        fig.savefig(pdf, **kw)
    except PermissionError:
        alt = PLOT / f"{stem}_new.pdf"
        fig.savefig(alt, **kw)
        print(f"  PDF locked, wrote {alt.name} instead")
    print(f"  wrote {stem}")


def panel_label(ax, letter: str) -> None:
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes, fontsize=8, fontweight="bold",
            va="bottom", ha="left", color="black", clip_on=False)


# ---------------------------------------------------------------------------
# Fig. 1 genome-wide burden
# ---------------------------------------------------------------------------
def plot_burden() -> None:
    rows = apply_pheno_condition(epi.read_tsv(TAB / "sample_burden.tsv"))
    if FILTER_GW810:
        rows = [r for r in rows if keep_gw_8_10(r.get("Gestational_Week", ""))]
    for r in rows:
        r["SNV_mut_per_Mb"] = fnum(r["SNV_nonsyn_count"]) / CODING_MB
    by = {g: [r for r in rows if r["Condition"] == g] for g in GROUPS}
    fig, axes = plt.subplots(1, 3, figsize=(180 * MM, 64 * MM))
    specs = [
        ("SNV_mut_per_Mb", "Nonsynonymous SNV\n(mut/Mb)"),
        ("SV_total", "Alt-carrying PASS SV count"),
        ("SV_plp_count", "P/LP SV count\n(ACMG 4 or 5)"),
    ]
    ptab = []
    n_sizes = {g: len(by[g]) for g in GROUPS}
    for ax, (col, ylab), letter in zip(axes, specs, "abc"):
        data = {g: [fnum(r[col]) for r in by[g]] for g in GROUPS}
        epi.box_with_points(ax, data, ylab, group_ns=n_sizes)
        ps = epi.annotate_mw(ax, data)
        ptab.append({"panel": letter, "metric": col, **ps})
        panel_label(ax, letter)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.82, bottom=0.16, wspace=0.45)
    save_fig(fig, "01_genomewide_burden")
    plt.close(fig)
    epi.write_tsv(
        PLOT / "source_data" / "01_burden.tsv",
        [{"Sample_ID": r["Sample_ID"], "Condition": r["Condition"],
          "SNV_mut_per_Mb": fnum(r["SNV_nonsyn_count"]) / CODING_MB,
          "SV_total": r["SV_total"], "SV_plp_count": r["SV_plp_count"]}
         for r in rows],
        ["Sample_ID", "Condition", "SNV_mut_per_Mb", "SV_total", "SV_plp_count"],
    )
    epi.write_tsv(
        PLOT / "source_data" / "01_mannwhitney.tsv",
        ptab,
        ["panel", "metric", "abnormal_vs_normal", "abnormal_vs_control", "normal_vs_control"],
    )


# ---------------------------------------------------------------------------
# Fig. 2 per-locus SV (FDR ab vs control; 8–10 week denominators)
# ---------------------------------------------------------------------------
def load_plp_locus_rows() -> list[dict]:
    cache = LOCUS_TAB / "sv_plp_locus_rates.tsv"
    if cache.exists():
        return epi.read_tsv(cache)
    print("  scanning VCF for AnnotSV P/LP loci...")
    pheno_rows = load_phenotype(PHENO_PATH)
    if FILTER_GW810:
        pheno_rows = [r for r in pheno_rows if keep_gw_8_10(r.get("Gestational_Week", ""))]
    samples = [r["VCF_Sample_ID"] for r in pheno_rows]
    pheno_map = {r["VCF_Sample_ID"]: r["Condition"] for r in pheno_rows}
    groups = locsv.group_samples(samples, pheno_map)
    recs = load_plp_sv_records(SV_ANNO)
    coord_keys = {r["coord_key"] for r in recs}
    sv_ids = {r.get("sv_id") or "" for r in recs}
    loci = locsv.collect_pass_sv_loci_matching(samples, locsv.SV_VCF, coord_keys, sv_ids)
    rows = locsv.build_enrichment_rows(loci, groups)
    drop = {"carriers"}
    slim = [{k: v for k, v in r.items() if k not in drop} for r in rows]
    if slim:
        locsv.write_tsv(slim, cache)
    print(f"  P/LP PASS loci with a carrier: {len(slim)}")
    return slim


def norm_coord_key(raw: str) -> str:
    return str(raw or "").lower().replace("chr", "").strip()


def load_epifactors_genes() -> set[str]:
    panel = EPI_TAB / "gene_panel.tsv"
    if panel.exists():
        return {r["gene"].upper() for r in epi.read_tsv(panel)}
    return {g.upper() for g in epi.GENES}


def lookup_annotsv_loci(coord_keys: set[str]) -> dict[str, dict]:
    """Best AnnotSV row per coord key (full mode preferred)."""
    want = {norm_coord_key(k) for k in coord_keys}
    best: dict[str, dict] = {}
    from compute_burden import parse_acmg_class, parse_max_pop_af

    with open(SV_ANNO, newline="", encoding="utf-8", errors="replace") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = norm_coord_key(
                f"{row.get('SV_chrom')}:{row.get('SV_start')}:{row.get('SV_end')}:{row.get('SV_type')}"
            )
            if key not in want:
                continue
            acmg = parse_acmg_class(row.get("ACMG_class"))
            mode = row.get("Annotation_mode", "")
            score = fnum(row.get("AnnotSV_ranking_score"))
            prev = best.get(key)
            if prev is None or (mode == "full" and prev.get("mode") != "full") or (
                mode == prev.get("mode") and score > fnum(prev.get("ranking_score"))
            ):
                genes = [g.strip().upper() for g in str(row.get("Gene_name") or "").split(";") if g.strip()]
                best[key] = {
                    "coord_key": key,
                    "acmg_class": acmg if acmg is not None else "",
                    "acmg_label": row.get("ACMG_class", ""),
                    "ranking_score": row.get("AnnotSV_ranking_score", ""),
                    "gene_names": ";".join(genes),
                    "mode": mode,
                    "max_pop_af": parse_max_pop_af(row) if parse_max_pop_af(row) is not None else "",
                }
    return best


def summarize_fig2_loci(plot_rows: list[dict], plp_all: list[dict], plp_hi: list[dict]) -> None:
    epi_genes = load_epifactors_genes()
    keys = {r["coord_key"] for r in plot_rows}
    annot = lookup_annotsv_loci(keys)

    fig2a_rows = []
    for r in plot_rows:
        meta = annot.get(norm_coord_key(r["coord_key"]), {})
        genes = [g for g in str(meta.get("gene_names", "")).split(";") if g]
        epi_hits = sorted({g.upper() for g in genes if g.upper() in epi_genes})
        fig2a_rows.append(
            {
                "coord_key": r["coord_key"],
                "svtype": r.get("svtype", ""),
                "abnormal_rate": r.get("abnormal_rate", ""),
                "control_rate": r.get("control_rate", ""),
                "fdr_abnormal_vs_control": r.get("fdr_abnormal_vs_control", ""),
                "enrichment_pattern": r.get("enrichment_pattern", ""),
                "acmg_class": meta.get("acmg_class", ""),
                "acmg_label": meta.get("acmg_label", ""),
                "ranking_score": meta.get("ranking_score", ""),
                "annotsv_genes": meta.get("gene_names", ""),
                "epifactors_genes": ";".join(epi_hits),
                "epifactors_hit": int(bool(epi_hits)),
            }
        )

    def acmg_bucket(v) -> str:
        try:
            x = int(float(v))
        except (TypeError, ValueError):
            return "missing"
        if x >= 5:
            return "5_LP"
        if x == 4:
            return "4_P"
        if x == 3:
            return "3_VUS"
        if x == 2:
            return "2_LB"
        if x == 1:
            return "1_B"
        return "0"

    acmg_ct = Counter(acmg_bucket(r["acmg_class"]) for r in fig2a_rows)
    summary_a = [
        {"panel": "fig2a", "metric": "n_loci", "value": len(fig2a_rows)},
        {"panel": "fig2a", "metric": "n_with_acmg_class", "value": sum(1 for r in fig2a_rows if r["acmg_class"] != "")},
        {"panel": "fig2a", "metric": "n_with_ranking_score", "value": sum(1 for r in fig2a_rows if r["ranking_score"] != "")},
        {"panel": "fig2a", "metric": "n_epifactors_hit", "value": sum(int(r["epifactors_hit"]) for r in fig2a_rows)},
        {"panel": "fig2a", "metric": "n_acmg_ge4", "value": sum(1 for r in fig2a_rows if acmg_bucket(r["acmg_class"]) in {"4_P", "5_LP"})},
    ]
    for k, v in sorted(acmg_ct.items()):
        summary_a.append({"panel": "fig2a", "metric": f"acmg_{k}", "value": v})

    plp_detail = []
    plp_keys = {norm_coord_key(r["coord_key"]) for r in plp_all}
    hi_keys = {norm_coord_key(r["coord_key"]) for r in plp_hi}
    plp_annot = lookup_annotsv_loci(plp_keys) if plp_keys else {}
    for r in plp_all:
        meta = plp_annot.get(norm_coord_key(r["coord_key"]), {})
        genes = [g for g in str(meta.get("gene_names") or r.get("gene") or "").split(";") if g]
        epi_hits = sorted({g.upper() for g in genes if g.upper() in epi_genes})
        plp_detail.append(
            {
                "coord_key": r["coord_key"],
                "svtype": r.get("svtype", ""),
                "acmg_class": r.get("acmg_class", meta.get("acmg_class", "")),
                "acmg_label": r.get("acmg_label", meta.get("acmg_label", "")),
                "ranking_score": meta.get("ranking_score", r.get("ranking_score", "")),
                "annotsv_genes": meta.get("gene_names", r.get("gene", "")),
                "epifactors_genes": ";".join(epi_hits),
                "epifactors_hit": int(bool(epi_hits)),
                "abnormal_rate": r.get("abnormal_rate", ""),
                "normal_rate": r.get("normal_rate", ""),
                "control_rate": r.get("control_rate", ""),
                "abnormal_highest": int(norm_coord_key(r["coord_key"]) in hi_keys),
            }
        )
    summary_c = [
        {"panel": "fig2c", "metric": "n_plp_loci", "value": len(plp_all)},
        {"panel": "fig2c", "metric": "n_plp_abnormal_highest", "value": len(plp_hi)},
        {"panel": "fig2c", "metric": "n_plp_epifactors_hit", "value": sum(int(r["epifactors_hit"]) for r in plp_detail)},
        {
            "panel": "fig2c",
            "metric": "n_plp_abnormal_highest_epifactors_hit",
            "value": sum(int(r["epifactors_hit"]) for r in plp_detail if int(r["abnormal_highest"]) == 1),
        },
    ]

    epi.write_tsv(PLOT / "source_data" / "02_fig2a_locus_annotation.tsv", fig2a_rows, list(fig2a_rows[0].keys()) if fig2a_rows else ["coord_key"])
    epi.write_tsv(PLOT / "source_data" / "02_fig2c_plp_annotation.tsv", plp_detail, list(plp_detail[0].keys()) if plp_detail else ["coord_key"])
    epi.write_tsv(
        PLOT / "source_data" / "02_fig2_annotation_summary.tsv",
        summary_a + summary_c,
        ["panel", "metric", "value"],
    )
    print(
        f"  fig2a annot: {len(fig2a_rows)} loci, "
        f"ACMG≥4={summary_a[4]['value']}, EpiFactors={summary_a[3]['value']}"
    )
    print(
        f"  fig2c annot: {len(plp_all)} P/LP, EpiFactors={summary_c[2]['value']}, "
        f"ab-highest EpiFactors={summary_c[3]['value']}"
    )


def plot_sv_loci() -> None:
    path = LOCUS_TAB / "sv_locus_enrichment_fdr05_any_comparison.tsv"
    rows = epi.read_tsv(path)
    ab_ctrl = [r for r in rows if fnum(r["fdr_abnormal_vs_control"]) < 0.05]
    plot_rows = ab_ctrl if ab_ctrl else rows
    n_plot = len(plot_rows)
    pal = {
        "abnormal_specific": "#D55E00",
        "case_vs_control": "#0072B2",
        "abnormal_vs_control_only": "#E69F00",
        "normal_vs_control": "#009E73",
        "none": "#999999",
    }

    plp_all = load_plp_locus_rows()
    plp_rows = [
        r for r in plp_all
        if fnum(r["abnormal_rate"]) > fnum(r["normal_rate"])
        and fnum(r["abnormal_rate"]) > fnum(r["control_rate"])
    ]

    fig, axes = plt.subplots(1, 3, figsize=(180 * MM, 62 * MM))

    ax = axes[0]
    xs = [fnum(r["control_rate"]) * 100 for r in plot_rows]
    ys = [fnum(r["abnormal_rate"]) * 100 for r in plot_rows]
    cs = [pal.get(r["enrichment_pattern"], "#999999") for r in plot_rows]
    ax.plot([0, 100], [0, 100], color="#888888", lw=0.5, ls="--")
    ax.scatter(xs, ys, c=cs, s=9, linewidths=0, alpha=0.85)
    ax.set_xlabel("Control carrier rate (%)")
    ax.set_ylabel("Abnormal carrier rate (%)")
    ax.set_xlim(-2, 102)
    ax.set_ylim(-2, 102)
    panel_label(ax, "a")

    ax = axes[1]
    rates = [fnum(r["abnormal_rate"]) * 100 for r in plot_rows]
    ax.hist(rates, bins=np.arange(0, 105, 10), color="#0072B2", edgecolor="white", linewidth=0.4)
    ax.axvline(5, color="#D55E00", lw=0.7, ls="--")
    ax.set_xlabel(f"Abnormal carrier rate among\nthe {n_plot} FDR loci (%)")
    ax.set_ylabel("Number of SV loci")
    ax.annotate(
        "5%", xy=(5, 0.2), xytext=(14, 8),
        fontsize=5.5, color="#D55E00",
        arrowprops={"arrowstyle": "-", "color": "#D55E00", "lw": 0.4},
    )
    panel_label(ax, "b")

    ax = axes[2]
    ax.plot([0, 100], [0, 100], color="#888888", lw=0.5, ls="--")
    xs = [fnum(r["control_rate"]) * 100 for r in plp_rows]
    ys = [fnum(r["abnormal_rate"]) * 100 for r in plp_rows]
    ax.scatter(xs, ys, c="#D55E00", s=12, linewidths=0, alpha=0.85, zorder=3)
    ax.set_xlabel("Control carrier rate (%)")
    ax.set_ylabel("Abnormal carrier rate (%)")
    ax.set_xlim(-2, 102)
    ax.set_ylim(-2, 102)
    ax.set_title("P/LP, abnormal highest", fontsize=6, pad=2)
    ax.text(
        0.04, 0.96,
        f"n={len(plp_rows)} / {len(plp_all)} P/LP sites",
        transform=ax.transAxes, va="top", ha="left", fontsize=5.5,
    )
    panel_label(ax, "c")

    fig.subplots_adjust(left=0.07, right=0.98, top=0.86, bottom=0.22, wspace=0.42)
    save_fig(fig, "02_sv_locus_enrichment")
    plt.close(fig)

    src = [{
        "coord_key": r["coord_key"], "svtype": r["svtype"],
        "abnormal_rate": r["abnormal_rate"], "normal_rate": r["normal_rate"],
        "control_rate": r["control_rate"], "fdr_abnormal_vs_control": r["fdr_abnormal_vs_control"],
        "enrichment_pattern": r["enrichment_pattern"],
        "cohort_carrier": r.get("cohort_carrier", ""),
    } for r in plot_rows]
    fields = [
        "coord_key", "svtype", "abnormal_rate", "normal_rate", "control_rate",
        "fdr_abnormal_vs_control", "enrichment_pattern", "cohort_carrier",
    ]
    epi.write_tsv(PLOT / "source_data" / "02_fdr05_plotted.tsv", src, fields)
    epi.write_tsv(
        PLOT / "source_data" / "02_fdr05_ab_vs_control.tsv",
        [r for r in src if fnum(r["fdr_abnormal_vs_control"]) < 0.05],
        fields,
    )
    if plp_all:
        epi.write_tsv(
            PLOT / "source_data" / "02_plp_locus_rates.tsv",
            plp_all,
            [k for k in plp_all[0].keys() if k != "carriers"],
        )
    if plp_rows:
        epi.write_tsv(
            PLOT / "source_data" / "02_plp_abnormal_highest.tsv",
            plp_rows,
            [k for k in plp_rows[0].keys() if k != "carriers"],
        )
    summarize_fig2_loci(plot_rows, plp_all, plp_rows)


# ---------------------------------------------------------------------------
# Gene models (Ensembl MANE / canonical) + rare SNV lollipops
# ---------------------------------------------------------------------------
def ensembl_get(url: str, retries: int = 4) -> dict | list:
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Content-Type": "application/json", "User-Agent": "burden-narrative/1.0"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except Exception as exc:
            last_err = exc
            time.sleep(1.5 * (attempt + 1))
    raise last_err


def pick_transcript(gene: dict) -> dict:
    txs = gene.get("Transcript") or []
    for t in txs:
        mane = t.get("MANE") or []
        if any((m.get("type") or "").lower().find("select") >= 0 for m in mane if isinstance(m, dict)):
            return t
        if str(t.get("source", "")).lower().find("mane") >= 0:
            return t
    for t in txs:
        if t.get("is_canonical") in {1, True, "1"}:
            return t
    coding = [t for t in txs if t.get("biotype") == "protein_coding"]
    pool = coding or txs
    return max(pool, key=lambda t: abs(int(t["end"]) - int(t["start"])))


def fetch_gene_model(symbol: str, cache: Path) -> dict:
    cache.mkdir(parents=True, exist_ok=True)
    fp = cache / f"{symbol}.json"
    if fp.exists():
        return json.loads(fp.read_text(encoding="utf-8"))
    url = f"{ENSEMBL}/lookup/symbol/homo_sapiens/{symbol}?expand=1;content-type=application/json"
    gene = ensembl_get(url)
    tx = pick_transcript(gene)
    exons = sorted(tx.get("Exon") or [], key=lambda e: int(e["start"]))
    cds = None
    trans = tx.get("Translation") or {}
    if trans:
        cds = (int(trans.get("start", tx["start"])), int(trans.get("end", tx["end"])))
    model = {
        "symbol": symbol,
        "chrom": str(gene.get("seq_region_name")),
        "gene_start": int(gene["start"]),
        "gene_end": int(gene["end"]),
        "strand": int(gene.get("strand") or 1),
        "tx_id": tx.get("id", ""),
        "tx_start": int(tx["start"]),
        "tx_end": int(tx["end"]),
        "exons": [{"start": int(e["start"]), "end": int(e["end"])} for e in exons],
        "cds": cds,
    }
    fp.write_text(json.dumps(model), encoding="utf-8")
    time.sleep(0.15)
    return model


def parse_aa_for_gene(aa: str, gene: str) -> tuple[str, str]:
    """Return (exon_label, protein_change) from ANNOVAR AAChange."""
    gene = gene.upper()
    for part in str(aa or "").split(","):
        bits = part.split(":")
        if not bits or bits[0].upper() != gene:
            continue
        exon = next((b for b in bits if b.lower().startswith("exon")), "")
        prot = next((b for b in bits if b.startswith("p.")), "")
        return exon, prot
    return "", ""


def collect_rare_sites(gene: str, snv_long: list[dict]) -> list[dict]:
    sites = defaultdict(lambda: {"abnormal": 0, "normal": 0, "control": 0, "meta": None})
    for r in snv_long:
        if r["gene"] != gene:
            continue
        if str(r.get("pop_rare")) not in {"1", "1.0", 1}:
            continue
        key = (r["chrom"], r["pos"], r["ref"], r["alt"])
        sites[key][r.get("Condition", "")] += 1
        if sites[key]["meta"] is None:
            exon, prot = parse_aa_for_gene(r.get("aa_change", ""), gene)
            sites[key]["meta"] = {
                "chrom": r["chrom"],
                "pos": inum(r["pos"]),
                "consequence": r.get("consequence", ""),
                "exon": exon,
                "protein": prot,
                "variant": r.get("variant_pos", ""),
            }
    out = []
    for rec in sites.values():
        m = rec["meta"] or {}
        m.update({g: rec[g] for g in GROUPS})
        out.append(m)
    out.sort(key=lambda x: x["pos"])
    return out


def draw_gene_track(ax, model: dict, sites: list[dict], ylabel: str) -> None:
    g0, g1 = model["tx_start"], model["tx_end"]
    span = max(g1 - g0, 1)
    ax.plot([g0, g1], [0.35, 0.35], color="#222222", lw=0.7, solid_capstyle="butt", zorder=1)
    cds = model.get("cds")
    for i, ex in enumerate(model["exons"]):
        x0, x1 = ex["start"], ex["end"]
        coding = False
        if cds:
            a, b = max(x0, cds[0]), min(x1, cds[1])
            coding = a < b
        h = 0.28 if coding else 0.14
        ax.add_patch(Rectangle((x0, 0.35 - h / 2), max(x1 - x0, span * 0.0008), h,
                               facecolor="#000000", edgecolor="none", zorder=2))
    # lollipops: height encodes abnormal carrier count; colour = group mix
    ymax = max((s["abnormal"] for s in sites), default=1)
    for s in sites:
        x = s["pos"]
        if not (g0 - span * 0.02 <= x <= g1 + span * 0.02):
            continue
        n_ab = s["abnormal"]
        if n_ab > 0:
            h = 0.55 + 1.15 * (n_ab / ymax)
            ax.plot([x, x], [0.50, h], color=GROUP_COLOR["abnormal"], lw=0.6, zorder=3)
            ax.scatter([x], [h], s=10, color=GROUP_COLOR["abnormal"], zorder=4, linewidths=0)
        else:
            ax.scatter([x], [0.08], s=7, color="#0072B2", zorder=4, linewidths=0, marker="v", alpha=0.7)
    ax.set_xlim(g0 - span * 0.03, g1 + span * 0.03)
    ax.set_ylim(-0.15, 2.05)
    ax.set_yticks([])
    ax.set_ylabel(ylabel, fontstyle="italic", rotation=0, ha="right", va="center", labelpad=18)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    # Mb tick
    ax.set_xticks([g0, g1])
    ax.set_xticklabels([f"{g0/1e6:.2f} Mb", f"{g1/1e6:.2f} Mb"])
    strand = "→" if model["strand"] > 0 else "←"
    n_ab_sites = sum(1 for s in sites if s["abnormal"] > 0)
    ax.text(0.01, 0.92, f"{strand}  chr{model['chrom']}  {n_ab_sites} rare sites in abnormal",
            transform=ax.transAxes, fontsize=5.5, va="top")


def collect_sv_events(gene: str, sv_long: list[dict]) -> list[dict]:
    """Unique overlapping SVs per gene; carrier counts by group."""
    sites = defaultdict(lambda: {"abnormal": 0, "normal": 0, "control": 0, "meta": None})
    for r in sv_long:
        if r["gene"] != gene:
            continue
        key = (r.get("chrom"), r.get("start"), r.get("end"), r.get("svtype"), r.get("sv_id"))
        cond = r.get("Condition", "")
        if cond in GROUPS:
            sites[key][cond] += 1
        if sites[key]["meta"] is None:
            sites[key]["meta"] = {
                "chrom": r.get("chrom"),
                "start": inum(r.get("start")),
                "end": inum(r.get("end") or r.get("start")),
                "svtype": r.get("svtype", ""),
                "svlen": r.get("svlen", ""),
                "location": r.get("location", ""),
                "exon_or_splice": inum(r.get("exon_or_splice")),
                "acmg_class": r.get("acmg_class", ""),
                "variant": r.get("variant_pos", ""),
            }
    out = []
    for rec in sites.values():
        m = rec["meta"] or {}
        m.update({g: rec[g] for g in GROUPS})
        out.append(m)
    out.sort(key=lambda x: (x["start"], x["end"]))
    return out


def draw_sv_track(ax, model: dict, events: list[dict], ylabel: str) -> None:
    g0, g1 = model["tx_start"], model["tx_end"]
    span = max(g1 - g0, 1)
    ax.plot([g0, g1], [0.35, 0.35], color="#222222", lw=0.7, solid_capstyle="butt", zorder=1)
    cds = model.get("cds")
    for ex in model["exons"]:
        x0, x1 = ex["start"], ex["end"]
        coding = False
        if cds:
            a, b = max(x0, cds[0]), min(x1, cds[1])
            coding = a < b
        h = 0.28 if coding else 0.14
        ax.add_patch(Rectangle((x0, 0.35 - h / 2), max(x1 - x0, 1), h,
                               facecolor="#000000", edgecolor="none", zorder=2))
    vis = []
    for e in events:
        s, t = sv_span(e)
        if t < g0 or s > g1:
            continue
        vis.append(e)
    ymax = max((e["abnormal"] for e in vis), default=1) or 1
    for e in vis:
        s, t = sv_span(e)
        exonish = int(e.get("overlap_exon") or 0) == 1
        if e["abnormal"] > 0:
            y = 0.55 + 1.15 * (e["abnormal"] / ymax)
            color = "#D55E00" if exonish else "#E69F00"
            if sv_is_insertion(e):
                ax.vlines(s, 0.48, y, colors=color, lw=1.0, zorder=4)
                ax.plot(s, y, marker="o", color=color, markersize=3.2, zorder=5)
            else:
                ax.add_patch(Rectangle((s, 0.48), max(t - s, 1), y - 0.48,
                                       facecolor=color, edgecolor="none", alpha=0.85, zorder=3))
        else:
            color = "#0072B2"
            if sv_is_insertion(e):
                ax.plot(s, 0.11, marker="v", color=color, markersize=3.2, zorder=4)
            else:
                ax.add_patch(Rectangle((s, 0.02), max(t - s, 1), 0.18,
                                       facecolor=color, edgecolor="none", alpha=0.45, zorder=3))
    ax.set_xlim(g0 - span * 0.03, g1 + span * 0.03)
    ax.set_ylim(-0.15, 2.05)
    ax.set_yticks([])
    ax.set_ylabel(ylabel, fontstyle="italic", rotation=0, ha="right", va="center", labelpad=18)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.set_xticks([g0, g1])
    ax.set_xticklabels([f"{g0/1e6:.2f} Mb", f"{g1/1e6:.2f} Mb"])
    strand = "→" if model["strand"] > 0 else "←"
    n_ab = sum(1 for e in vis if e["abnormal"] > 0)
    n_ex = sum(1 for e in vis if e["abnormal"] > 0 and int(e.get("overlap_exon") or 0) == 1)
    ax.text(
        0.01, 0.92,
        f"{strand}  chr{model['chrom']}  {n_ab} SV in abnormal ({n_ex} overlap drawn exons)",
        transform=ax.transAxes, fontsize=5.5, va="top",
    )


def top_genes_from_rank(path: Path, n: int = N_FIG4_GENES) -> tuple[list[str], list[dict]]:
    if not path.exists():
        return [], []
    ranked = epi.read_tsv(path)
    genes = [r["gene"] for r in ranked[:n]]
    return genes, ranked


def plot_gene_sv_maps() -> None:
    rank_path = EPI_PLOT / "source_data" / "fig3c_sv_gene_rank.tsv"
    if not rank_path.exists():
        rank_path = PLOT / "source_data" / "03_sv_gene_rank.tsv"
    genes, ranked = top_genes_from_rank(rank_path, N_FIG4_GENES)
    epi.write_tsv(
        PLOT / "source_data" / "05_gene_rank_from_fig3c.tsv",
        ranked,
        epi.GENE_RANK_FIELDS if ranked else ["gene"],
    )
    if not genes:
        print("  no SV genes from fig3c rank")
        return
    sv_long = epi.read_tsv(EPI_TAB / "sv_variants.tsv")
    cache = TAB / "epifactors" / "ensembl_cache"
    models = {}
    for g in genes:
        print(f"  Ensembl SV {g}")
        models[g] = fetch_gene_model(g, cache)

    n_genes = len(genes)
    fig_h = (32 * n_genes + 16) * MM
    fig, axes = plt.subplots(n_genes, 1, figsize=(180 * MM, fig_h), sharex=False)
    if n_genes == 1:
        axes = [axes]
    rows = []
    for ax, gene in zip(axes, genes):
        events = collect_sv_events(gene, sv_long)
        for e in events:
            e["annotsv_exon"] = int(e.get("exon_or_splice") or 0)
            e["overlap_exon"] = int(sv_overlaps_drawn_exons(e, models[gene]))
        draw_sv_track(ax, models[gene], events, gene)
        for e in events:
            rows.append({"gene": gene, **e})
    axes[-1].set_xlabel("Genomic position (GRCh38)")
    has_ab_exon_ins = any(
        e["abnormal"] > 0 and int(e.get("overlap_exon") or 0) == 1 and sv_is_insertion(e) for e in rows
    )
    has_ab_exon_span = any(
        e["abnormal"] > 0 and int(e.get("overlap_exon") or 0) == 1 and not sv_is_insertion(e) for e in rows
    )
    has_ab_intron = any(
        e["abnormal"] > 0 and int(e.get("overlap_exon") or 0) != 1 for e in rows
    )
    has_other = any(e["abnormal"] <= 0 for e in rows)
    handles = []
    if has_ab_exon_ins:
        handles.append(Line2D(
            [0], [0], color="#D55E00", marker="o", lw=1.1, markersize=4,
            label="INS overlapping a drawn exon (height = n abnormal)",
        ))
    if has_ab_exon_span:
        handles.append(Patch(
            facecolor="#D55E00", edgecolor="none", alpha=0.85,
            label="DEL/DUP overlapping a drawn exon",
        ))
    if has_ab_intron:
        handles.append(Patch(
            facecolor="#E69F00", edgecolor="none", alpha=0.85,
            label="SV in abnormal, not overlapping drawn exons",
        ))
    if has_other:
        handles.append(Patch(
            facecolor="#0072B2", edgecolor="none", alpha=0.45,
            label="SV only in normal/control",
        ))
    handles.append(Line2D(
        [0], [0], color="black", marker="s", lw=0, markersize=5,
        label="Drawn transcript exon (taller = CDS)",
    ))
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=2,
        frameon=False,
        fontsize=5.5,
        bbox_to_anchor=(0.55, 0.995),
        borderaxespad=0,
    )
    fig.subplots_adjust(left=0.14, right=0.98, top=0.90, bottom=0.07, hspace=0.50)
    save_fig(fig, "05_candidate_gene_sv_maps", pad_inches=0.10)
    plt.close(fig)
    fields = [
        "gene", "chrom", "start", "end", "svtype", "svlen", "location",
        "exon_or_splice", "annotsv_exon", "overlap_exon", "acmg_class", "variant",
    ] + list(GROUPS)
    epi.write_tsv(PLOT / "source_data" / "05_sv_events.tsv", rows, fields)


def plot_gene_maps() -> None:
    rank_path = EPI_PLOT / "source_data" / "fig3d_snv_gene_rank.tsv"
    if not rank_path.exists():
        rank_path = PLOT / "source_data" / "03_snv_gene_rank.tsv"
    genes, ranked = top_genes_from_rank(rank_path, N_FIG4_GENES)
    epi.write_tsv(
        PLOT / "source_data" / "04_gene_rank_from_fig3d.tsv",
        ranked,
        epi.GENE_RANK_FIELDS if ranked else ["gene"],
    )
    if not genes:
        print("  no SNV genes from fig3d rank")
        return
    snv_long = epi.read_tsv(EPI_TAB / "snv_variants.tsv")
    cache = TAB / "epifactors" / "ensembl_cache"
    models = {}
    for g in genes:
        print(f"  Ensembl SNV {g}")
        models[g] = fetch_gene_model(g, cache)

    fig, axes = plt.subplots(len(genes), 1, figsize=(180 * MM, 28 * MM * max(len(genes), 2)), sharex=False)
    if len(genes) == 1:
        axes = [axes]
    site_rows = []
    for ax, gene in zip(axes, genes):
        sites = collect_rare_sites(gene, snv_long)
        draw_gene_track(ax, models[gene], sites, gene)
        for s in sites:
            site_rows.append({"gene": gene, **s})
    axes[-1].set_xlabel("Genomic position (GRCh38)")
    handles = [
        Line2D([0], [0], color=GROUP_COLOR["abnormal"], marker="o", lw=0.6, label="Rare SNV in abnormal (height = n carriers)", markersize=4),
        Line2D([0], [0], color="#0072B2", marker="v", lw=0, label="Rare SNV only in normal/control", markersize=4),
        Line2D([0], [0], color="black", marker="s", lw=0, label="Exon (taller = CDS)", markersize=4),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.55, 0.995))
    fig.subplots_adjust(left=0.14, right=0.98, top=0.94, bottom=0.05, hspace=0.55)
    save_fig(fig, "04_candidate_gene_maps")
    plt.close(fig)

    fields = ["gene", "chrom", "pos", "consequence", "exon", "protein", "variant"] + list(GROUPS)
    epi.write_tsv(PLOT / "source_data" / "04_rare_sites.tsv", site_rows, fields)

    # compact table of abnormal protein changes for the talk
    talk_rows = []
    for s in site_rows:
        if s["abnormal"] <= 0:
            continue
        talk_rows.append(s)
    epi.write_tsv(PLOT / "source_data" / "04_abnormal_rare_sites.tsv", talk_rows, fields)


def copy_epifactors_overview() -> None:
    src_pdf = EPI_PLOT / "fig1_epifactors_overview.pdf"
    if src_pdf.exists():
        shutil.copy2(src_pdf, PLOT / "03_epifactors_overview.pdf")
        png = EPI_PLOT / "fig1_epifactors_overview.png"
        if png.exists():
            shutil.copy2(png, PLOT / "03_epifactors_overview.png")
    mw = EPI_PLOT / "source_data" / "fig1ab_mannwhitney.tsv"
    if mw.exists():
        shutil.copy2(mw, PLOT / "source_data" / "03_mannwhitney.tsv")
    for src_name, dst_name in (
        ("fig3c_sv_gene_rank.tsv", "03_sv_gene_rank.tsv"),
        ("fig3d_snv_gene_rank.tsv", "03_snv_gene_rank.tsv"),
        ("fig3c_sv_site_rates.tsv", "03_sv_site_rates.tsv"),
    ):
        src = EPI_PLOT / "source_data" / src_name
        if src.exists():
            shutil.copy2(src, PLOT / "source_data" / dst_name)


def write_legends() -> None:
    p1 = {r["metric"]: r for r in epi.read_tsv(PLOT / "source_data" / "01_mannwhitney.tsv")}
    p3 = {}
    p3path = EPI_PLOT / "source_data" / "fig1ab_mannwhitney.tsv"
    alt = PLOT / "source_data" / "03_mannwhitney.tsv"
    if alt.exists():
        p3rows = epi.read_tsv(alt)
    elif p3path.exists():
        p3rows = epi.read_tsv(p3path)
    else:
        p3rows = []
    for r in p3rows:
        p3[r["panel"]] = r

    def pv(row, key):
        return epi.fmt_p(fnum(row[key])) if row else "NA"

    s_snv = p1.get("SNV_mut_per_Mb", {})
    s_sv = p1.get("SV_total", {})
    s_plp = p1.get("SV_plp_count", {})
    a3 = p3.get("a_rare_snv", {})
    b3 = p3.get("b_sv_any", {})

    n_by = Counter()
    burden_src = PLOT / "source_data" / "01_burden.tsv"
    if burden_src.exists():
        for r in epi.read_tsv(burden_src):
            n_by[r["Condition"]] += 1
    n_ab = n_by.get("abnormal", 0)
    n_nm = n_by.get("normal", 0)
    n_ct = n_by.get("control", 0)
    fdr_ab_path = PLOT / "source_data" / "02_fdr05_ab_vs_control.tsv"
    fdr_plot_path = PLOT / "source_data" / "02_fdr05_plotted.tsv"
    n_fdr = len(epi.read_tsv(fdr_ab_path)) if fdr_ab_path.exists() else 0
    n_fdr_plot = len(epi.read_tsv(fdr_plot_path)) if fdr_plot_path.exists() else 0
    plp_all_path = PLOT / "source_data" / "02_plp_locus_rates.tsv"
    plp_hi_path = PLOT / "source_data" / "02_plp_abnormal_highest.tsv"
    n_plp_all = len(epi.read_tsv(plp_all_path)) if plp_all_path.exists() else 0
    n_plp_hi = len(epi.read_tsv(plp_hi_path)) if plp_hi_path.exists() else 0
    fig2_sum = {r["metric"]: r["value"] for r in epi.read_tsv(PLOT / "source_data" / "02_fig2_annotation_summary.tsv")} if (PLOT / "source_data" / "02_fig2_annotation_summary.tsv").exists() else {}
    gw_note = "8–10 周分母" if FILTER_GW810 else "全队列分母"

    cohort = (
        "本套图为全部 648 例（不限孕周）。"
        "abnormal = 原临床 abnormal + 高变 CpG 聚类中 30 例疑似 abnormal（原 Class3=normal_case，已从 normal 划出）。"
        if not FILTER_GW810
        else "本套图仅保留孕周 8–10 周样本（g8/g9/g10，含 8+/9+/10+）。"
    )
    text = f"""图注与解读（按汇报顺序）

{cohort}
检验均为双侧 Mann–Whitney U（未做多重校正）。n = {n_ab} abnormal / {n_nm} normal / {n_ct} control。

========== Fig. 1 全基因组负荷 ==========
a, 规范编码外显子（85.77 Mb）上的非同义 SNV 速率（mut/Mb）。
   检验：Ab–N p={pv(s_snv,'abnormal_vs_normal')}；Ab–C p={pv(s_snv,'abnormal_vs_control')}；N–C p={pv(s_snv,'normal_vs_control')}。
   解读：abnormal 与 normal 没有差别，高甲基化这一刀并不对应更高的编码突变负担。
   abnormal/normal 相对 control 若显著，方向是略低而不是更高，排除“突变越多→高甲基化”。

b, 每人携带的 PASS SV 条数（只计 0/1、1/1 等带 alt 的基因型，不计 0/0）。
   检验：Ab–N p={pv(s_sv,'abnormal_vs_normal')}；Ab–C p={pv(s_sv,'abnormal_vs_control')}；N–C p={pv(s_sv,'normal_vs_control')}。
   解读：中位数大约每人 3–4 千条。Fig.2 检验的是队列里出现过的不重复位点（约 5.2 万），不是每人 5 万条。
   同一条常见 SV 会在几百人中重复出现，所以“每人几千”和“位点库五万”可以同时成立。

c, 每人携带的 AnnotSV P/LP（ACMG 4 或 5）条数（同样只计 alt）。
   检验：Ab–N p={pv(s_plp,'abnormal_vs_normal')}；Ab–C p={pv(s_plp,'abnormal_vs_control')}；N–C p={pv(s_plp,'normal_vs_control')}。
   解读：三组无显著差别；abnormal 并不携带更多注释致病 SV。

========== Fig. 2 逐位点 SV（{gw_note}） ==========
a, FDR<0.05 的 PASS SV 位点（共 {n_fdr_plot} 个；FDR ab vs ctrl = {n_fdr}）。
   附表 02_fig2a_locus_annotation.tsv：AnnotSV ACMG 等级、ranking score、是否命中 EpiFactors 基因。
   统计：ACMG≥4 = {fig2_sum.get('n_acmg_ge4', 'NA')} / {fig2_sum.get('n_loci', 'NA')}；EpiFactors 命中 = {fig2_sum.get('n_epifactors_hit', 'NA')}。
b, 上述位点在 abnormal 中的携带率分布。虚线=5%。
c, AnnotSV P/LP（ACMG 4–5）中 abnormal 携带率同时高于 normal 与 control 的位点（共 {n_plp_hi}/{n_plp_all}）。
   附表 02_fig2c_plp_annotation.tsv：P/LP 位点 EpiFactors 命中 = {fig2_sum.get('n_plp_epifactors_hit', 'NA')}；其中 abnormal 最高 = {fig2_sum.get('n_plp_abnormal_highest_epifactors_hit', 'NA')}。

========== Fig. 3 EpiFactors 基因库 ==========
a–b, 每人 rare SNV / panel SV 条数（箱线图 x 轴标注 n）。
   检验：a Ab–N p={pv(a3,'abnormal_vs_normal')}；Ab–C p={pv(a3,'abnormal_vs_control')}；N–C p={pv(a3,'normal_vs_control')}。
   b Ab–N p={pv(b3,'abnormal_vs_normal')}；Ab–C p={pv(b3,'abnormal_vs_control')}；N–C p={pv(b3,'normal_vs_control')}。

c, EpiFactors 基因：按 abnormal vs control 携带者 Fisher 原始 p 升序（第一行=携带率差最大），展示 top 12 基因的 SV 携带比例。
d, 同上排序规则，展示 rare SNV 携带比例（layout 与 c 相同）。

========== Fig. 4 rare SNV 基因图（来自 Fig.3d top {N_FIG4_GENES}） ==========
不再按 HHI 集中度选基因；直接取 Fig.3d 原始 p 最靠前的 {N_FIG4_GENES} 个基因画 lollipop。
黑块=外显子；橙棒=abnormal rare SNV；蓝=仅 normal/control。

========== Fig. 5 panel SV 基因图（来自 Fig.3c top {N_FIG4_GENES}） ==========
取 Fig.3c 原始 p 最靠前的 {N_FIG4_GENES} 个基因。深橙=与所画外显子相交的 SV；浅橙=内含子；圆点=INS。
附表 03_sv_gene_rank.tsv / 03_snv_gene_rank.tsv 为 Fig.3c/d 完整排序。
"""
    (PLOT / "figure_legends.txt").write_text(text, encoding="utf-8")
    (PLOT / "talking_points.txt").write_text(text, encoding="utf-8")


def main() -> None:
    only = sys.argv[2] if len(sys.argv) >= 3 and sys.argv[1] == "--only" else None
    epi.apply_nature_style()
    epi.PLOT_DIR = EPI_PLOT
    (PLOT / "source_data").mkdir(parents=True, exist_ok=True)
    if only in {None, "1"}:
        print("[1] burden")
        plot_burden()
    if only in {None, "2"}:
        print("[2] SV loci")
        plot_sv_loci()
    if only in {None, "3"}:
        print("[3] EpiFactors overview")
        epi.make_plots(EPI_TAB)
        copy_epifactors_overview()
    if only in {None, "4"}:
        print("[4] gene maps")
        plot_gene_maps()
    if only in {None, "5"}:
        print("[5] gene SV maps")
        plot_gene_sv_maps()
    if only is None:
        write_legends()
        if EPI_PLOT.exists() and EPI_PLOT.resolve() != PLOT.resolve():
            shutil.rmtree(EPI_PLOT)
            print(f"removed old plot dir {EPI_PLOT}")
    print(f"done -> {PLOT}")


if __name__ == "__main__":
    main()
