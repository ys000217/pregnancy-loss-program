#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Nature-style Manhattan + QQ for RPL GWAS (suggestive line P < 1e-4).

Follows Nature Portfolio artwork defaults (Arial/Helvetica, 5–7 pt,
Okabe–Ito colours, print-size canvas, PDF vector + PNG companion).
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# ---- paths ----
MODULE = Path(__file__).resolve().parents[1]
RESULTS = MODULE / "results"
# local hybrid (not in git); fall back to figure3 path
CANDIDATES = [
    Path(r"D:\ONT\figure3\gwas\combine_gwas_result_v1.Status.glm.logistic.hybrid"),
    MODULE.parent.parent / "figure3" / "gwas" / "combine_gwas_result_v1.Status.glm.logistic.hybrid",
]
HITS_CSV = RESULTS / "GWAS_suggestive_hits_1e-4.csv"

OUT_PDF = RESULTS / "Fig_GWAS_Manhattan_QQ_suggestive_1e-4.pdf"
OUT_PNG = RESULTS / "Fig_GWAS_Manhattan_QQ_suggestive_1e-4.png"
OUT_SRC_A = RESULTS / "Fig_GWAS_Manhattan_QQ_suggestive_1e-4_source_a.csv"
OUT_SRC_B = RESULTS / "Fig_GWAS_Manhattan_QQ_suggestive_1e-4_source_b.csv"
OUT_SUMMARY = RESULTS / "Fig_GWAS_Manhattan_QQ_suggestive_1e-4_summary.txt"

# ---- thresholds ----
SIG_P = 5e-8
SUG_P = 1e-4
YLIM = (0, 8.5)

# Okabe–Ito
OI = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermilion": "#D55E00",
    "purple": "#CC79A7",
}

MM = 1 / 25.4


def apply_nature_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.labelsize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.unicode_minus": False,
            "figure.dpi": 300,
            "savefig.dpi": 600,
            "savefig.bbox": None,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def chrom_to_int(x) -> float:
    x = str(x).strip().replace("chr", "")
    if x == "X":
        return 23.0
    if x == "Y":
        return 24.0
    if x in ("MT", "M"):
        return 25.0
    try:
        return float(x)
    except ValueError:
        return np.nan


def lambda_gc(p: np.ndarray) -> float:
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    chisq = stats.chi2.ppf(1.0 - p, df=1)
    return float(np.median(chisq) / stats.chi2.ppf(0.5, df=1))


def load_gwas(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "A1",
            "OMITTED",
            "TEST",
            "OBS_CT",
            "OR",
            "LOG_OR_SE",
            "Z_STAT",
            "P",
        ],
        low_memory=False,
        na_values=[".", "NA"],
    )
    df["CHR_NUM"] = df["CHROM"].map(chrom_to_int)
    for c in ["POS", "OBS_CT", "OR", "LOG_OR_SE", "Z_STAT", "P"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    gwas = df.dropna(subset=["P", "CHR_NUM", "POS"]).copy()
    gwas = gwas[(gwas["P"] > 0) & (gwas["P"] <= 1)]
    gwas = gwas[gwas["CHR_NUM"].between(1, 22)]
    gwas["NEGLOG10P"] = -np.log10(gwas["P"].clip(lower=1e-300))
    gwas["SNP"] = (
        gwas["CHROM"].astype(str).str.replace("^chr", "", regex=True)
        + ":"
        + gwas["POS"].astype(int).astype(str)
    )
    return gwas


def build_manhattan_coords(gwas: pd.DataFrame):
    chroms = sorted(gwas["CHR_NUM"].unique())
    offset = 0.0
    ticks, tick_labels = [], []
    xs = np.empty(len(gwas), dtype=float)
    chrom_colors = np.empty(len(gwas), dtype=object)
    palette = [OI["blue"], OI["sky"]]
    for i, c in enumerate(chroms):
        mask = gwas["CHR_NUM"].to_numpy() == c
        pos = gwas.loc[mask, "POS"].to_numpy()
        xs[mask] = pos + offset
        chrom_colors[mask] = palette[i % 2]
        ticks.append(offset + (pos.max() + pos.min()) / 2.0)
        # show every chrom label; Nature 6 pt is readable at 180 mm
        tick_labels.append(str(int(c)))
        offset += float(pos.max()) + 1e7
    return xs, chrom_colors, ticks, tick_labels, offset


def thin_qq(p: np.ndarray, max_points: int = 80_000):
    p = np.sort(p[np.isfinite(p) & (p > 0) & (p <= 1)])
    n = len(p)
    exp = -np.log10(np.arange(1, n + 1) / (n + 1))
    obs = -np.log10(p)
    if n > max_points:
        # keep denser sampling in the extreme tail
        idx = np.unique(
            np.concatenate(
                [
                    np.linspace(0, n - 1, max_points - 5_000, dtype=int),
                    np.arange(n - 5_000, n),
                ]
            )
        )
        exp, obs = exp[idx], obs[idx]
    return exp, obs, n


def main() -> None:
    apply_nature_style()
    RESULTS.mkdir(parents=True, exist_ok=True)

    inp = next((p for p in CANDIDATES if p.is_file()), None)
    if inp is None:
        raise FileNotFoundError("GWAS hybrid not found in expected paths")

    print(f"[1] reading {inp}")
    gwas = load_gwas(inp)
    lam = lambda_gc(gwas["P"].to_numpy())
    n_sig = int((gwas["P"] < SIG_P).sum())
    n_sug = int((gwas["P"] < SUG_P).sum())
    print(f"    n={len(gwas):,}  lambda_GC={lam:.4f}  GWS={n_sig}  suggestive={n_sug}")

    sug = gwas["P"] < SUG_P
    xs, chrom_colors, ticks, tick_labels, x_max = build_manhattan_coords(gwas)

    # ---- figure: double-column, panels a (Manhattan) + b (QQ) ----
    # 180 mm wide × ~70 mm tall
    fig = plt.figure(figsize=(180 * MM, 72 * MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.35, 1.0], wspace=0.28)
    ax_m = fig.add_subplot(gs[0, 0])
    ax_q = fig.add_subplot(gs[0, 1])

    # a) Manhattan — background SNPs
    bg = ~sug
    ax_m.scatter(
        xs[bg],
        gwas.loc[bg, "NEGLOG10P"].to_numpy(),
        c=list(chrom_colors[bg]),
        s=2.5,
        linewidths=0,
        rasterized=True,
        zorder=1,
    )
    # suggestive SNPs on top
    ax_m.scatter(
        xs[sug],
        gwas.loc[sug, "NEGLOG10P"].to_numpy(),
        c=OI["vermilion"],
        s=6,
        linewidths=0,
        zorder=3,
        label=f"Suggestive (P < {SUG_P:g}; n={n_sug})",
    )
    y_sug = -np.log10(SUG_P)
    y_sig = -np.log10(SIG_P)
    ax_m.axhline(y_sug, color=OI["blue"], linestyle="--", linewidth=0.7, zorder=2)
    ax_m.axhline(y_sig, color=OI["black"], linestyle=":", linewidth=0.7, zorder=2)
    # line labels (right side, small)
    ax_m.text(
        x_max * 0.995,
        y_sug + 0.18,
        f"Suggestive {SUG_P:g}",
        ha="right",
        va="bottom",
        fontsize=5.5,
        color=OI["blue"],
    )
    ax_m.text(
        x_max * 0.995,
        min(y_sig, YLIM[1] - 0.35) - 0.05,
        f"Genome-wide {SIG_P:g}",
        ha="right",
        va="top",
        fontsize=5.5,
        color=OI["black"],
    )

    ax_m.set_xlim(0, x_max)
    ax_m.set_ylim(*YLIM)
    ax_m.set_xticks(ticks)
    # label odd chromosomes + last to reduce clutter at print size
    shown = []
    for lab, t in zip(tick_labels, ticks):
        ci = int(lab)
        if ci % 2 == 1 or ci in (2, 22):
            shown.append((t, lab))
    ax_m.set_xticks([t for t, _ in shown])
    ax_m.set_xticklabels([lab for _, lab in shown])
    ax_m.set_xlabel("Chromosome")
    ax_m.set_ylabel(r"$-\log_{10}(P)$")
    ax_m.legend(loc="upper left", frameon=False, handletextpad=0.3, borderpad=0.2)

    # OLA1 peak callout (chr2 cluster)
    ola = gwas[(gwas["CHR_NUM"] == 2) & (gwas["POS"].between(174_070_000, 174_220_000))]
    if len(ola):
        peak = ola.loc[ola["P"].idxmin()]
        # x of peak among full gwas order
        peak_i = gwas.index.get_loc(peak.name)
        if isinstance(peak_i, slice):
            peak_i = peak_i.start
        ax_m.annotate(
            "OLA1",
            xy=(xs[peak_i], peak["NEGLOG10P"]),
            xytext=(xs[peak_i] - x_max * 0.06, peak["NEGLOG10P"] + 0.9),
            fontsize=6,
            fontstyle="italic",
            arrowprops=dict(arrowstyle="-", color=OI["black"], lw=0.5),
            ha="right",
            va="bottom",
        )

    # b) QQ
    exp, obs, n_qq = thin_qq(gwas["P"].to_numpy())
    lim = float(max(exp.max(), obs.max(), 1.0) * 1.02)
    ax_q.scatter(exp, obs, s=3, c=OI["blue"], linewidths=0, rasterized=True, zorder=2)
    ax_q.plot([0, lim], [0, lim], color=OI["black"], linestyle="--", linewidth=0.7, zorder=1)
    ax_q.set_xlim(0, lim)
    ax_q.set_ylim(0, lim)
    ax_q.set_xlabel(r"Expected $-\log_{10}(P)$")
    ax_q.set_ylabel(r"Observed $-\log_{10}(P)$")
    ax_q.text(
        0.04,
        0.96,
        rf"$\lambda_{{\mathrm{{GC}}}}={lam:.3f}$",
        transform=ax_q.transAxes,
        ha="left",
        va="top",
        fontsize=6,
    )

    # panel labels
    for ax, lab in ((ax_m, "a"), (ax_q, "b")):
        ax.text(
            -0.08 if ax is ax_m else -0.18,
            1.05,
            lab,
            transform=ax.transAxes,
            fontsize=8,
            fontweight="bold",
            va="bottom",
            ha="right",
        )

    fig.subplots_adjust(left=0.06, right=0.985, bottom=0.16, top=0.90)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG)
    plt.close(fig)
    print(f"[2] wrote {OUT_PDF.name} / {OUT_PNG.name}")

    # source data — panel a: suggestive hits + chrom midpoints; panel b: thinned QQ
    src_a = gwas.loc[sug, ["SNP", "CHROM", "POS", "REF", "ALT", "A1", "OR", "Z_STAT", "P", "NEGLOG10P"]].copy()
    src_a.insert(0, "x_plot", xs[sug])
    src_a.to_csv(OUT_SRC_A, index=False)
    pd.DataFrame({"expected_neglog10p": exp, "observed_neglog10p": obs}).to_csv(OUT_SRC_B, index=False)

    # refresh hits table in results/
    hits = gwas.loc[sug].sort_values("P").copy()
    hits["Level"] = np.where(hits["P"] < SIG_P, "Genome-wide Significant", "Suggestive")
    se = hits["LOG_OR_SE"]
    log_or = np.log(hits["OR"].clip(lower=1e-300))
    hits["OR_L95"] = np.exp(log_or - 1.96 * se)
    hits["OR_U95"] = np.exp(log_or + 1.96 * se)
    hits_out = RESULTS / "GWAS_suggestive_hits_1e-4.csv"
    hits[
        [
            "SNP",
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "A1",
            "OMITTED",
            "OBS_CT",
            "OR",
            "Z_STAT",
            "P",
            "Level",
        ]
    ].to_csv(hits_out, index=False)

    summary = "\n".join(
        [
            "Nature-style GWAS figure (suggestive line)",
            f"input: {inp}",
            f"n_variants: {len(gwas)}",
            f"lambda_GC: {lam:.6f}",
            f"genomewide_P: {SIG_P}",
            f"suggestive_P: {SUG_P}",
            f"n_genomewide: {n_sig}",
            f"n_suggestive: {n_sug}",
            f"figure_pdf: {OUT_PDF.name}",
            f"figure_png: {OUT_PNG.name}",
            f"source_a: {OUT_SRC_A.name}",
            f"source_b: {OUT_SRC_B.name}",
            "style: Nature Portfolio (180 mm; Arial 6–7 pt; Okabe–Ito)",
        ]
    )
    OUT_SUMMARY.write_text(summary + "\n", encoding="utf-8")
    print(f"[3] source data + summary written under {RESULTS}")


if __name__ == "__main__":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    main()
