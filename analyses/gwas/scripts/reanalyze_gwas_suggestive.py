#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Re-analyze PLINK GWAS logistic hybrid output with a small-sample suggestive line.

Thresholds (discussion with PI / small N ≈ 600):
  - Genome-wide significance : P < 5e-8   (field convention; report only)
  - Suggestive line (primary) : P < 1e-4   (looser than prior 1e-5 for small N)

Also computes genomic inflation factor λ_GC (median χ² method), mirroring
`pca/曼哈顿质控图.Rmd`. Sensitivity across alternate GWAS runs is deferred.
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# ---- paths ----
GWAS_DIR = Path(__file__).resolve().parent
INPUT = GWAS_DIR / "combine_gwas_result_v1.Status.glm.logistic.hybrid"
OUT_HITS = GWAS_DIR / "GWAS_suggestive_hits_1e-4.csv"
OUT_SUMMARY = GWAS_DIR / "GWAS_reanalysis_summary_1e-4.txt"
OUT_MANHATTAN = GWAS_DIR / "combine_GWAS_Manhattan_suggestive_1e-4.png"
OUT_QQ = GWAS_DIR / "combine_GWAS_QQ_suggestive_1e-4.png"

# ---- thresholds ----
SIG_P = 5e-8          # genome-wide
SUG_P = 1e-4          # suggestive (new, looser than previous 1e-5)
PLOT_YLIM = (0, 10)


def chrom_to_int(x: str) -> float:
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
    """Genomic inflation factor from median of χ²(1) transformed P."""
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    chisq = stats.chi2.ppf(1.0 - p, df=1)
    return float(np.median(chisq) / stats.chi2.ppf(0.5, df=1))


def manhattan(df: pd.DataFrame, out: Path) -> None:
    plot_df = df.dropna(subset=["CHR_NUM", "POS", "P"]).copy()
    plot_df = plot_df[plot_df["CHR_NUM"].between(1, 22)]
    plot_df["NEGLOG10P"] = -np.log10(plot_df["P"].clip(lower=1e-300))

    chroms = sorted(plot_df["CHR_NUM"].unique())
    # cumulative offset for x-axis
    offset = 0.0
    ticks = []
    tick_labels = []
    xs = np.empty(len(plot_df))
    colors = []
    palette = ["#27408B", "#87CEEB"]
    for i, c in enumerate(chroms):
        mask = plot_df["CHR_NUM"].values == c
        pos = plot_df.loc[mask, "POS"].values
        xs[mask] = pos + offset
        colors.extend([palette[i % 2]] * int(mask.sum()))
        mid = offset + (pos.max() + pos.min()) / 2.0
        ticks.append(mid)
        tick_labels.append(str(int(c)))
        offset += pos.max() + 1e7  # gap between chroms

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.scatter(xs, plot_df["NEGLOG10P"].values, c=colors, s=3, linewidths=0)
    ax.axhline(-np.log10(SUG_P), color="blue", linestyle="--", linewidth=1, label=f"suggestive {SUG_P:g}")
    ax.axhline(-np.log10(SIG_P), color="red", linestyle="--", linewidth=1, label=f"genome-wide {SIG_P:g}")
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels)
    ax.set_ylim(*PLOT_YLIM)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$-\log_{10}(P)$")
    ax.set_title("Manhattan Plot (combine_gwas_result_v1)")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def qq_plot(p: np.ndarray, lam: float, out: Path) -> None:
    p = np.sort(p[np.isfinite(p) & (p > 0) & (p <= 1)])
    n = len(p)
    exp = -np.log10(np.arange(1, n + 1) / (n + 1))
    obs = -np.log10(p)
    # thin for speed/plot size
    if n > 200_000:
        step = n // 200_000
        exp, obs = exp[::step], obs[::step]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(exp, obs, s=4, c="#27408B", linewidths=0)
    lim = max(exp.max(), obs.max(), 1)
    ax.plot([0, lim], [0, lim], "r--", linewidth=1)
    ax.set_xlabel(r"Expected $-\log_{10}(P)$")
    ax.set_ylabel(r"Observed $-\log_{10}(P)$")
    ax.set_title(f"QQ Plot (lambda_GC = {lam:.4f})")
    ax.set_xlim(0, lim * 1.02)
    ax.set_ylim(0, lim * 1.02)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main() -> None:
    print(f"[1] reading {INPUT}")
    df = pd.read_csv(
        INPUT,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "CHROM", "POS", "ID", "REF", "ALT", "A1", "OMITTED", "TEST",
            "OBS_CT", "OR", "LOG_OR_SE", "Z_STAT", "P",
        ],
        low_memory=False,
        na_values=[".", "NA"],
    )
    # header row may have been skipped via comment='#'; re-read with header if needed
    # PLINK files start with #CHROM — pandas comment='#' drops it; columns are correct.

    df["CHR_NUM"] = df["CHROM"].map(chrom_to_int)
    for c in ["POS", "OBS_CT", "OR", "LOG_OR_SE", "Z_STAT", "P"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    gwas = df.dropna(subset=["P", "CHR_NUM", "POS"]).copy()
    gwas = gwas[(gwas["P"] > 0) & (gwas["P"] <= 1)]
    gwas["SNP"] = gwas["CHROM"].astype(str).str.replace("^chr", "", regex=True) + ":" + gwas["POS"].astype(int).astype(str)

    n_total = len(gwas)
    lam = lambda_gc(gwas["P"].to_numpy())
    n_sig = int((gwas["P"] < SIG_P).sum())
    n_sug = int((gwas["P"] < SUG_P).sum())

    print(f"[2] variants with valid P: {n_total:,}")
    print(f"    lambda_GC = {lam:.4f}")
    print(f"    P < {SIG_P:g} (genome-wide): {n_sig}")
    print(f"    P < {SUG_P:g} (suggestive)  : {n_sug}")

    hits = gwas.loc[gwas["P"] < SUG_P].copy()
    hits = hits.sort_values("P")
    hits["Level"] = np.where(hits["P"] < SIG_P, "Genome-wide Significant", "Suggestive")
    # CI on OR scale
    se = hits["LOG_OR_SE"]
    log_or = np.log(hits["OR"].clip(lower=1e-300))
    hits["OR_L95"] = np.exp(log_or - 1.96 * se)
    hits["OR_U95"] = np.exp(log_or + 1.96 * se)
    hits["OR_95_CI"] = (
        hits["OR"].map(lambda x: f"{x:.4g}" if pd.notna(x) else "NA")
        + " ("
        + hits["OR_L95"].map(lambda x: f"{x:.4g}" if pd.notna(x) else "NA")
        + "-"
        + hits["OR_U95"].map(lambda x: f"{x:.4g}" if pd.notna(x) else "NA")
        + ")"
    )

    out_cols = [
        "SNP", "CHROM", "POS", "REF", "ALT", "A1", "OMITTED",
        "OBS_CT", "OR", "OR_95_CI", "Z_STAT", "P", "Level",
    ]
    hits[out_cols].to_csv(OUT_HITS, index=False)
    print(f"[3] wrote hits -> {OUT_HITS}")

    # locus summary
    auto = hits[hits["CHR_NUM"].between(1, 22)]
    chr_counts = auto.groupby(auto["CHR_NUM"].astype(int)).size().to_dict()
    if len(auto):
        locus_min = int(auto["POS"].min())
        locus_max = int(auto["POS"].max())
        main_chr = int(auto["CHR_NUM"].mode().iloc[0])
        on_main = auto[auto["CHR_NUM"] == main_chr]
        main_span = (int(on_main["POS"].min()), int(on_main["POS"].max()))
    else:
        main_chr, main_span, locus_min, locus_max = None, (None, None), None, None

    summary_lines = [
        "GWAS reanalysis summary",
        f"input: {INPUT.name}",
        f"n_variants: {n_total}",
        f"lambda_GC: {lam:.6f}",
        f"genomewide_P: {SIG_P}",
        f"suggestive_P: {SUG_P}",
        f"n_genomewide: {n_sig}",
        f"n_suggestive: {n_sug}",
        f"chr_distribution: {chr_counts}",
        f"main_chr: {main_chr}",
        f"main_chr_span_bp: {main_span[0]}-{main_span[1]}",
        f"all_hits_pos_range: {locus_min}-{locus_max}",
        f"hits_csv: {OUT_HITS.name}",
        f"manhattan: {OUT_MANHATTAN.name}",
        f"qq: {OUT_QQ.name}",
    ]
    OUT_SUMMARY.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print(f"[4] wrote summary -> {OUT_SUMMARY}")

    print("[5] plotting Manhattan / QQ ...")
    manhattan(gwas, OUT_MANHATTAN)
    qq_plot(gwas["P"].to_numpy(), lam, OUT_QQ)
    print(f"    -> {OUT_MANHATTAN.name}")
    print(f"    -> {OUT_QQ.name}")

    print("\n--- top 20 suggestive hits ---")
    print(
        hits[out_cols]
        .head(20)
        .to_string(index=False)
    )


if __name__ == "__main__":
    # avoid unicode issues on some Windows consoles
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    main()
