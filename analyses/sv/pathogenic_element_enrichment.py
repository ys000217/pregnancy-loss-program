#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
显著 CpG–window pair 在基因窗口 / Roadmap 7 类元件上的分布与富集

设计
----
1. Roadmap 18-state → 7 类：Prom / Enh / Tx / Biv / Repr / Het-ZNF / Quies
2. 富集「背景」两层含义同时使用：
   a) 计数宇宙：只在效应显著的 pair 内计数（致病栏则只在致病断点相关显著 pair）
   b) 期望对照：p0 来自同口径的非显著 pair 分布
      - 全部显著栏：随机抽样非显著 pair（默认 5000）
      - 致病栏：致病断点相关的非显著 pair（不足则全用；过多则抽样）
3. fig2a / fig2b：每栏画 观测(显著集合) vs 期望(非显著对照) 的比例，* = 富集 FDR<0.05

输出（D:\ONT\figures\ 与 D:\ONT\）
  - fig2a_window_gene_bar.png
  - fig2b_window_roadmap_bar.png
  - pair_7class_enrichment.tsv
  - pair_gene_window_enrichment.tsv
"""

from __future__ import print_function

import os
import shutil
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- paths ----
ANNOTSV = r"D:\ONT\clinical_649.GRCh38.annotsv.tsv"
ROADMAP = r"D:\ONT\figure2\Roadmap_placenta\E091_Placenta_18state_hg38_chr1_22.bed"
GWIN = r"D:\ONT\genes_windows.tsv"  # local large/intermediate; not in git
SIGPAIRS = r"D:\ONT\analyses\sv\results\sv_methylation_sig_pairs.tsv"
RESULTS = r"D:\ONT\sv_methylation_results.tsv"  # full matrix of tests; local only
PSTAR = r"D:\ONT\analyses\sv\results\p_star.txt"
FIGDIR = r"D:\ONT\analyses\sv\results\plots"
OUTDIR = r"D:\ONT\analyses\sv\results"

STATE_TO_CLASS = {}
for s in ("TssA", "TssFlnk", "TssFlnkU", "TssFlnkD"):
    STATE_TO_CLASS[s] = "Prom"
for s in ("EnhA1", "EnhA2", "EnhWk", "EnhG1", "EnhG2"):
    STATE_TO_CLASS[s] = "Enh"
for s in ("Tx", "TxWk"):
    STATE_TO_CLASS[s] = "Tx"
for s in ("TssBiv", "BivFlnk", "EnhBiv"):
    STATE_TO_CLASS[s] = "Biv"
for s in ("ReprPC", "ReprPCWk"):
    STATE_TO_CLASS[s] = "Repr"
for s in ("Het", "ZNF/Rpts"):
    STATE_TO_CLASS[s] = "Het-ZNF"
STATE_TO_CLASS["Quies"] = "Quies"

CLASS_ORDER = ["Prom", "Enh", "Tx", "Biv", "Repr", "Het-ZNF", "Quies"]
CLASS_COLORS = {
    "Prom": "#dd8452", "Enh": "#64b5cd", "Tx": "#55a868",
    "Biv": "#9a60b4", "Repr": "#c44e52", "Het-ZNF": "#8c6d31",
    "Quies": "#999999",
}
GENE_WINS = ["up", "dn", "body"]
GENE_LABELS = {"up": "Upstream 100kb", "dn": "Downstream 100kb", "body": "Gene body"}
GENE_COLORS = {"up": "#4c72b0", "dn": "#55a868", "body": "#c44e52"}
PROM_STATES = {"TssA", "TssFlnk", "TssFlnkU", "TssFlnkD"}
ENH_STATES = {"EnhA1", "EnhA2", "EnhWk", "EnhG1", "EnhG2"}
N_BG_SAMPLE = 5000
RNG_SEED = 42
Q_STAR = 0.05


def load_roadmap():
    rm = pd.read_csv(ROADMAP, sep="\t", header=None,
                     names=["chrom", "start", "end", "state"])
    rm["length"] = (rm["end"] - rm["start"]).clip(lower=0)
    idx = {}
    for chrom, g in rm.groupby("chrom"):
        g = g.sort_values("start")
        idx[chrom] = (g["start"].to_numpy(), g["end"].to_numpy(),
                      g["state"].to_numpy())
    return idx


def annotate_positions(chroms, positions, idx):
    out = np.full(len(chroms), "NA", dtype=object)
    df = pd.DataFrame({"chrom": chroms, "pos": positions})
    for chrom, arr in idx.items():
        m = df["chrom"].values == chrom
        if not m.any():
            continue
        starts, ends, states = arr
        pos = df.loc[m, "pos"].to_numpy()
        i = np.searchsorted(starts, pos, side="right") - 1
        ok = (i >= 0) & (ends[i] > pos)
        st = np.full(len(pos), "NA", dtype=object)
        st[ok] = states[i[ok]]
        out[np.where(m)[0]] = st
    return out


def state_to_class(states):
    return [STATE_TO_CLASS.get(s, "NA") for s in states]


def extract_pathogenic_svs(path=ANNOTSV):
    usecols = ["AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_length",
               "SV_type", "Annotation_mode", "AnnotSV_ranking_score", "ACMG_class"]
    rows = []
    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=200000,
                             dtype={"ACMG_class": str, "SV_chrom": str}):
        m = (chunk["Annotation_mode"] == "full") & (chunk["ACMG_class"].isin(["4", "5"]))
        if m.any():
            rows.append(chunk.loc[m])
    if not rows:
        raise RuntimeError("未在 AnnotSV 中找到 ACMG_class 4/5 的 full 行")
    df = pd.concat(rows, ignore_index=True)
    df["chrom"] = df["SV_chrom"].astype(str).map(
        lambda c: c if c.startswith("chr") else ("chr" + c))
    df["start"] = df["SV_start"].astype(int)
    df["end"] = df["SV_end"].astype(int)
    df["svtype"] = df["SV_type"].astype(str)
    return df


def pathogenic_breakpoints(patho):
    bps = []
    for _, r in patho.iterrows():
        chrom, a, b, t = r.chrom, int(r.start), int(r.end), r.svtype.upper()
        if t in ("DEL", "DUP", "INV", "CNV"):
            bps.append((chrom, a, t, r.AnnotSV_ID, r.ACMG_class, "L"))
            if b != a:
                bps.append((chrom, b, t, r.AnnotSV_ID, r.ACMG_class, "R"))
        elif t in ("INS", "INS:ME", "MEI"):
            bps.append((chrom, a, t, r.AnnotSV_ID, r.ACMG_class, "S"))
        else:
            bps.append((chrom, a, t, r.AnnotSV_ID, r.ACMG_class, "L"))
            if b != a:
                bps.append((chrom, b, t, r.AnnotSV_ID, r.ACMG_class, "R"))
    return pd.DataFrame(bps, columns=["chrom", "pos", "svtype", "AnnotSV_ID",
                                      "ACMG_class", "side"])


def build_patho_keys(gw, bp_by_chrom, rm_idx):
    """预计算所有与致病断点相交的 (gene_id, window)。"""
    keys = set()
    for _, g in gw.iterrows():
        chrom = str(g.chrom)
        bps = bp_by_chrom.get(chrom)
        if bps is None or len(bps) == 0:
            continue
        pos = bps["pos"].to_numpy()
        gid = str(g.gene_id)
        for wname, cs, ce in (("up", "up_start", "up_end"),
                              ("dn", "dn_start", "dn_end"),
                              ("body", "body_start", "body_end")):
            a, b = int(g[cs]), int(g[ce])
            if ((pos >= a) & (pos <= b)).any():
                keys.add((gid, wname))
        a, b = int(g["up_start"]), int(g["dn_end"])
        m = (pos >= a) & (pos <= b)
        if not m.any():
            continue
        sub = pos[m]
        states = annotate_positions(np.array([chrom] * len(sub)), sub, rm_idx)
        if any(s in PROM_STATES for s in states):
            keys.add((gid, "prom"))
        if any(s in ENH_STATES for s in states):
            keys.add((gid, "enh"))
    return keys


def frac_dict(series_or_counts, categories):
    """类别比例；只在 categories 内归一化（排除 NA）。"""
    if isinstance(series_or_counts, pd.Series):
        vc = series_or_counts.value_counts()
    else:
        vc = pd.Series(series_or_counts)
    tot = float(sum(int(vc.get(c, 0)) for c in categories))
    if tot <= 0:
        return {c: 0.0 for c in categories}, 0
    return {c: int(vc.get(c, 0)) / tot for c in categories}, int(tot)


def enrichment_within_set(counts, expected_frac, categories):
    """计数宇宙内相对 expected_frac（期望对照）做单侧二项富集 + BH-FDR。"""
    n = int(sum(int(counts.get(c, 0)) for c in categories))
    rows = []
    for cat in categories:
        obs = int(counts.get(cat, 0))
        p0 = float(expected_frac.get(cat, 0.0))
        if n == 0 or p0 <= 0:
            rows.append(dict(category=cat, n_obs=obs, n_total=n, frac_obs=np.nan,
                             frac_expected=p0, fold=np.nan, pval=np.nan))
            continue
        try:
            pval = stats.binomtest(obs, n, p0, alternative="greater").pvalue
        except AttributeError:
            pval = stats.binom_test(obs, n, p0, alternative="greater")
        frac_obs = obs / n
        fold = frac_obs / p0
        rows.append(dict(category=cat, n_obs=obs, n_total=n, frac_obs=frac_obs,
                         frac_expected=p0, fold=fold, pval=pval))
    out = pd.DataFrame(rows)
    p = out["pval"].to_numpy(dtype=float)
    m = np.isfinite(p)
    q = np.full(len(p), np.nan)
    if m.any():
        pv = p[m]
        order = np.argsort(pv)
        ranked = pv[order]
        mm = len(ranked)
        qv = np.minimum.accumulate((ranked * mm / np.arange(1, mm + 1))[::-1])[::-1]
        q_m = np.empty(mm)
        q_m[order] = qv
        q[np.where(m)[0]] = q_m
    out["qval"] = q
    out["enriched"] = out["qval"] < Q_STAR
    return out


def plot_obs_vs_expected(enr_a, enr_b, categories, labels, colors,
                         title_a, title_b, n_bg_a, n_bg_b, outpath, note):
    """两栏：观测比例(显著集合) vs 期望比例(非显著对照)；* = 富集。"""
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), sharey=True)
    panels = [(axes[0], enr_a, title_a, n_bg_a),
              (axes[1], enr_b, title_b, n_bg_b)]
    w = 0.38
    for ax, enr, title, n_bg in panels:
        x = np.arange(len(categories))
        fro = [float(enr.set_index("category").loc[c, "frac_obs"])
               if c in set(enr["category"]) and np.isfinite(
                   enr.set_index("category").loc[c, "frac_obs"]) else 0.0
               for c in categories]
        fre = [float(enr.set_index("category").loc[c, "frac_expected"])
               if c in set(enr["category"]) else 0.0 for c in categories]
        n_obs = int(enr["n_total"].iloc[0]) if len(enr) else 0
        cols = [colors.get(c, "#888888") for c in categories]
        ax.bar(x - w / 2, fro, width=w, color=cols, label="Observed (sig set)")
        ax.bar(x + w / 2, fre, width=w, color="#cccccc", edgecolor="#666666",
               label="Expected (non-sig bg)")
        enr_map = {r.category: r for _, r in enr.iterrows()}
        for i, c in enumerate(categories):
            star = "*" if (c in enr_map and bool(enr_map[c].enriched)) else ""
            ax.text(i - w / 2, fro[i] + 0.01, star, ha="center", fontsize=12,
                    color="#d62728", fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels([labels.get(c, c) for c in categories],
                           rotation=35, ha="right", fontsize=8)
        ax.set_ylabel("Fraction of pairs")
        ax.set_title("%s\nn_sig=%d | n_bg=%d" % (title, n_obs, n_bg), fontsize=10)
        ax.set_ylim(0, max(max(fro + fre + [0.05]) * 1.25, 0.1))
        ax.legend(fontsize=7, loc="upper right")
    fig.suptitle(note, fontsize=9, y=1.03)
    plt.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print("[fig] wrote", outpath)


def sample_nonsig_backgrounds(sig, patho_keys, rm_idx):
    """
    返回:
      bg_all: 抽样非显著 pair（全基因组口径）
      bg_patho: 致病断点相关非显著 pair
    """
    p_star = (float(open(PSTAR).read().strip()) if os.path.exists(PSTAR)
              else float(sig["pval"].max()))
    sig_keys = set(zip(sig["site"].astype(str), sig["window"].astype(str)))
    usecols = ["site", "chrom", "pos", "gene_id", "window", "pval"]

    all_parts = []
    patho_parts = []
    print("scanning RESULTS for non-sig backgrounds ...", flush=True)
    for ci, chunk in enumerate(pd.read_csv(RESULTS, sep="\t", usecols=usecols,
                                           chunksize=300000)):
        chunk = chunk[chunk["pval"] > p_star]
        if chunk.empty:
            continue
        sk = list(zip(chunk["site"].astype(str), chunk["window"].astype(str)))
        chunk = chunk.copy()
        chunk = chunk[[k not in sig_keys for k in sk]]
        if chunk.empty:
            continue
        pk = list(zip(chunk["gene_id"].astype(str), chunk["window"].astype(str)))
        is_patho = [k in patho_keys for k in pk]
        patho_parts.append(chunk.loc[is_patho])
        # 每块最多抽 1500，最后再压到 N_BG_SAMPLE（近似均匀）
        n_take = min(len(chunk), 1500)
        all_parts.append(chunk.sample(n=n_take, random_state=RNG_SEED + ci))
        if (ci + 1) % 10 == 0:
            print("  chunks=%d ..." % (ci + 1), flush=True)

    if all_parts:
        bg_all = pd.concat(all_parts, ignore_index=True)
        if len(bg_all) > N_BG_SAMPLE:
            bg_all = bg_all.sample(n=N_BG_SAMPLE, random_state=RNG_SEED)
    else:
        bg_all = pd.DataFrame(columns=usecols)

    if patho_parts:
        bg_patho = pd.concat(patho_parts, ignore_index=True)
        if len(bg_patho) > N_BG_SAMPLE:
            bg_patho = bg_patho.sample(n=N_BG_SAMPLE, random_state=RNG_SEED)
    else:
        bg_patho = bg_all.iloc[0:0].copy()

    for name, df in (("bg_all", bg_all), ("bg_patho", bg_patho)):
        if len(df) == 0:
            print("%s n=0" % name, flush=True)
            continue
        states = annotate_positions(df["chrom"].astype(str).values,
                                    df["pos"].astype(int).values, rm_idx)
        df["state"] = states
        df["class7"] = state_to_class(states)
        print("%s n=%d" % (name, len(df)), flush=True)

    return bg_all, bg_patho


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    idx = load_roadmap()

    print("=== Pathogenic SVs (AnnotSV ACMG 4+5) ===", flush=True)
    patho = extract_pathogenic_svs()
    print("n_pathogenic_SV =", len(patho), flush=True)
    patho.to_csv(os.path.join(OUTDIR, "pathogenic_sv_acmg45.tsv"),
                 sep="\t", index=False)

    bps = pathogenic_breakpoints(patho)
    bps["state"] = annotate_positions(bps["chrom"].values, bps["pos"].values, idx)
    bps["class7"] = state_to_class(bps["state"])
    bps.to_csv(os.path.join(OUTDIR, "pathogenic_sv_breakpoints_chromatin.tsv"),
               sep="\t", index=False)
    bp_by_chrom = {c: g for c, g in bps.groupby("chrom")}

    gw = pd.read_csv(GWIN, sep="\t")
    print("building patho gene×window keys ...", flush=True)
    patho_keys = build_patho_keys(gw, bp_by_chrom, idx)
    print("patho-related gene×window keys =", len(patho_keys), flush=True)

    print("\n=== Significant pairs ===", flush=True)
    sig = pd.read_csv(SIGPAIRS, sep="\t")
    if not {"chrom", "pos"}.issubset(sig.columns):
        res_mini = pd.read_csv(RESULTS, sep="\t",
                               usecols=["site", "chrom", "pos"])
        sig = sig.merge(res_mini.drop_duplicates("site"), on="site", how="left")
    sig = sig.dropna(subset=["chrom", "pos"]).copy()
    sig["pos"] = sig["pos"].astype(int)
    sig["chrom"] = sig["chrom"].astype(str)
    sig["state"] = annotate_positions(sig["chrom"].values, sig["pos"].values, idx)
    sig["class7"] = state_to_class(sig["state"])
    sig["patho_related"] = [
        (str(g), str(w)) in patho_keys
        for g, w in zip(sig["gene_id"], sig["window"])
    ]
    print("significant =", len(sig),
          " patho_related =", int(sum(sig["patho_related"])), flush=True)
    sig.to_csv(os.path.join(OUTDIR, "sig_pairs_with_chromatin.tsv"),
               sep="\t", index=False)

    sig_all = sig
    sig_patho = sig[sig["patho_related"]].copy()

    # ---- 非显著期望对照 ----
    print("\n=== Non-sig expected backgrounds ===", flush=True)
    bg_all, bg_patho = sample_nonsig_backgrounds(sig, patho_keys, idx)
    p0_rm_all, n_bg_rm_all = frac_dict(pd.Series(bg_all["class7"]) if len(bg_all) else pd.Series(dtype=object),
                                       CLASS_ORDER)
    p0_rm_patho, n_bg_rm_patho = frac_dict(
        pd.Series(bg_patho["class7"]) if len(bg_patho) else pd.Series(dtype=object),
        CLASS_ORDER)

    bg_gene_all = bg_all[bg_all["window"].isin(GENE_WINS)] if len(bg_all) else bg_all
    bg_gene_patho = bg_patho[bg_patho["window"].isin(GENE_WINS)] if len(bg_patho) else bg_patho
    p0_g_all, n_bg_g_all = frac_dict(
        bg_gene_all["window"] if len(bg_gene_all) else pd.Series(dtype=object), GENE_WINS)
    p0_g_patho, n_bg_g_patho = frac_dict(
        bg_gene_patho["window"] if len(bg_gene_patho) else pd.Series(dtype=object), GENE_WINS)

    # 保存背景分布
    pd.DataFrame([
        dict(set="nonsig_all_sample", category=c, frac=p0_rm_all[c], n_bg=n_bg_rm_all)
        for c in CLASS_ORDER
    ] + [
        dict(set="nonsig_patho_related", category=c, frac=p0_rm_patho[c], n_bg=n_bg_rm_patho)
        for c in CLASS_ORDER
    ]).to_csv(os.path.join(OUTDIR, "pair_7class_expected_background.tsv"),
              sep="\t", index=False)

    # ---- Roadmap 7 类 ----
    print("\n=== Roadmap 7-class enrichment (expected = non-sig bg) ===", flush=True)
    vc_all = pd.Series(sig_all["class7"]).value_counts()
    vc_patho = pd.Series(sig_patho["class7"]).value_counts() if len(sig_patho) else pd.Series(dtype=int)
    enr_rm_all = enrichment_within_set(vc_all, p0_rm_all, CLASS_ORDER)
    enr_rm_all.insert(0, "set", "all_significant_pairs")
    enr_rm_all.insert(1, "expected_from", "nonsig_random_sample")
    enr_rm_patho = enrichment_within_set(vc_patho, p0_rm_patho, CLASS_ORDER)
    enr_rm_patho.insert(0, "set", "patho_breakpoint_related_pairs")
    enr_rm_patho.insert(1, "expected_from", "nonsig_patho_related")
    enr_rm = pd.concat([enr_rm_all, enr_rm_patho], ignore_index=True)
    enr_rm.to_csv(os.path.join(OUTDIR, "pair_7class_enrichment.tsv"),
                  sep="\t", index=False)
    print(enr_rm_all[["category", "n_obs", "frac_obs", "frac_expected", "fold", "qval", "enriched"]]
          .to_string(index=False), flush=True)
    print("--- patho ---")
    print(enr_rm_patho[["category", "n_obs", "frac_obs", "frac_expected", "fold", "qval", "enriched"]]
          .to_string(index=False), flush=True)

    plot_obs_vs_expected(
        enr_rm_all, enr_rm_patho, CLASS_ORDER, {c: c for c in CLASS_ORDER}, CLASS_COLORS,
        "All significant pairs",
        "Pathogenic-SV breakpoint–related pairs",
        n_bg_rm_all, n_bg_rm_patho,
        os.path.join(FIGDIR, "fig2b_window_roadmap_bar.png"),
        "Roadmap 7 classes: observed vs non-sig expected (* FDR<0.05 enrichment)",
    )

    # ---- 基因窗口 ----
    print("\n=== Gene-window enrichment (expected = non-sig bg) ===", flush=True)
    gene_all = sig_all[sig_all["window"].isin(GENE_WINS)]
    gene_patho = sig_patho[sig_patho["window"].isin(GENE_WINS)]
    vc_g_all = gene_all["window"].value_counts()
    vc_g_patho = gene_patho["window"].value_counts() if len(gene_patho) else pd.Series(dtype=int)
    enr_g_all = enrichment_within_set(vc_g_all, p0_g_all, GENE_WINS)
    enr_g_all.insert(0, "set", "all_significant_pairs")
    enr_g_all.insert(1, "expected_from", "nonsig_random_sample")
    enr_g_patho = enrichment_within_set(vc_g_patho, p0_g_patho, GENE_WINS)
    enr_g_patho.insert(0, "set", "patho_breakpoint_related_pairs")
    enr_g_patho.insert(1, "expected_from", "nonsig_patho_related")
    enr_g = pd.concat([enr_g_all, enr_g_patho], ignore_index=True)
    enr_g.to_csv(os.path.join(OUTDIR, "pair_gene_window_enrichment.tsv"),
                 sep="\t", index=False)
    print(enr_g_all[["category", "n_obs", "frac_obs", "frac_expected", "fold", "qval", "enriched"]]
          .to_string(index=False), flush=True)
    print("--- patho ---")
    print(enr_g_patho[["category", "n_obs", "frac_obs", "frac_expected", "fold", "qval", "enriched"]]
          .to_string(index=False), flush=True)

    plot_obs_vs_expected(
        enr_g_all, enr_g_patho, GENE_WINS, GENE_LABELS, GENE_COLORS,
        "All significant pairs",
        "Pathogenic-SV breakpoint–related pairs",
        n_bg_g_all, n_bg_g_patho,
        os.path.join(FIGDIR, "fig2a_window_gene_bar.png"),
        "Gene windows: observed vs non-sig expected (* FDR<0.05 enrichment)",
    )
    shutil.copyfile(os.path.join(FIGDIR, "fig2a_window_gene_bar.png"),
                    os.path.join(FIGDIR, "fig2_window_bar.png"))

    # 分布对照补充图
    if len(bg_all):
        sc = pd.Series(sig_all["class7"]).value_counts(normalize=True)
        bc = pd.Series(bg_all["class7"]).value_counts(normalize=True)
        x = np.arange(len(CLASS_ORDER))
        w = 0.4
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.bar(x - w / 2, [sc.get(c, 0) for c in CLASS_ORDER], width=w,
               label="Significant", color="#d62728")
        ax.bar(x + w / 2, [bc.get(c, 0) for c in CLASS_ORDER], width=w,
               label="Non-sig expected (n=%d)" % len(bg_all), color="#999999")
        ax.set_xticks(x)
        ax.set_xticklabels(CLASS_ORDER, rotation=30, ha="right")
        ax.set_ylabel("Fraction")
        ax.set_title("Roadmap 7-class: sig vs non-sig expected background")
        ax.legend(fontsize=8)
        plt.tight_layout()
        p = os.path.join(FIGDIR, "fig_pair_7class_dist_vs_nonsig.png")
        plt.savefig(p, dpi=150)
        plt.close()
        print("[fig] wrote", p, flush=True)

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
