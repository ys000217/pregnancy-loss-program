#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
β 值 vs M 值 分布诊断
=====================
采样甲基化矩阵, 对比:
  1. β 值分布 (histogram)
  2. M 值分布 (logit 变换后)
  3. 均值-方差关系 (评估同方差性: β 的方差随均值呈倒U, M 的方差较平坦)
输出 3 张 PNG + 屏幕摘要, 用于判断哪种因变量适合线性回归。
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MATRIX   = r"E:\甲基化数据矩阵\EWAS_INPUT_NO_HEADER.txt"
N_HEADER = 71
N_SAMPLE = 30000   # 采样前 N 个 CpG(矩阵前段)

def mtransform_2d(B):
    """与 sv_methylation_pipeline.mtransform 相同(自然对数 logit), 向量化版本。
    按位点(行)均值填补缺失, 0/1 裁剪到 [1e-3, 1-1e-3]。"""
    rowmean = np.nanmean(B, axis=1, keepdims=True)
    B = np.where(np.isnan(B), rowmean, B)
    B = np.clip(B, 1e-3, 1.0 - 1e-3)
    return np.log(B / (1.0 - B))

def main():
    print("[read] sampling %d CpGs from matrix..." % N_SAMPLE, flush=True)
    df = pd.read_csv(MATRIX, sep="\t", header=None, skiprows=N_HEADER,
                     nrows=N_SAMPLE, na_values=["NA", ""], low_memory=False,
                     dtype={0: str})
    beta = df.iloc[:, 1:].to_numpy(float)          # (N_SAMPLE, n_samples)
    n_sites, n_samp = beta.shape
    print("[read] done: %d sites x %d samples" % (n_sites, n_samp), flush=True)

    # ---- 每 CpG 均值/方差 ----
    bmean = np.nanmean(beta, axis=1)
    bvar  = np.nanvar(beta, axis=1)

    # ---- M 变换 ----
    M = mtransform_2d(beta)
    mmean = M.mean(axis=1)
    mvar  = M.var(axis=1)

    bflat = beta[~np.isnan(beta)]
    mflat = M.ravel()

    # ================= 摘要统计 =================
    def summ(x, name):
        q = np.percentile(x, [1, 5, 50, 95, 99])
        from scipy import stats as _st
        return ("%-6s | N=%-9d | skew=%.2f | kurt=%.2f | p1=%.3g p5=%.3g med=%.3g p95=%.3g p99=%.3g"
                % (name, len(x), _st.skew(x), _st.kurtosis(x),
                   q[0], q[1], q[2], q[3], q[4]))

    print("\n=== 分布摘要 ===")
    print(summ(bflat, "beta"))
    print(summ(mflat, "M"))
    print("beta 边界聚集: <0.05 占 %.1f%%, >0.95 占 %.1f%%"
          % (100 * (bflat < 0.05).mean(), 100 * (bflat > 0.95).mean()))
    print("M 极端: |M|>3 占 %.2f%%, |M|>5 占 %.2f%%"
          % (100 * (np.abs(mflat) > 3).mean(), 100 * (np.abs(mflat) > 5).mean()))

    # ---- 均值-方差分箱(同方差性诊断) ----
    def binned_var(mean_vals, var_vals, lo, hi, nbins=20):
        edges = np.linspace(lo, hi, nbins + 1)
        centers, bvars = [], []
        for k in range(nbins):
            m = (mean_vals >= edges[k]) & (mean_vals < edges[k + 1])
            if m.sum() < 20:
                continue
            centers.append((edges[k] + edges[k + 1]) / 2)
            bvars.append(np.median(var_vals[m]))
        return np.array(centers), np.array(bvars)

    bc, bv = binned_var(bmean, bvar, 0.0, 1.0)
    mc, mv = binned_var(mmean, mvar, -8.0, 8.0)
    print("\n=== 方差异质性(分箱方差 max/min 比, 越小越同方差) ===")
    if len(bv) > 1 and bv.min() > 0:
        print("beta: 分箱方差范围 %.4g ~ %.4g, 比=%.1f" % (bv.min(), bv.max(), bv.max() / bv.min()))
    if len(mv) > 1 and mv.min() > 0:
        print("M   : 分箱方差范围 %.4g ~ %.4g, 比=%.1f" % (mv.min(), mv.max(), mv.max() / mv.min()))

    # ================= 图 =================
    plt.rcParams.update({"font.size": 11})

    # Fig1: beta 分布
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(bflat, bins=100, color="#4c72b0", alpha=0.85, density=True)
    ax.set_xlabel("β value (methylation proportion)")
    ax.set_ylabel("Density")
    ax.set_title("Distribution of β values (n=%d sites × %d samples)" % (n_sites, n_samp))
    plt.tight_layout(); plt.savefig(r"D:\ONT\diagnosis_beta_dist.png", dpi=130); plt.close()

    # Fig2: M 分布
    fig, ax = plt.subplots(figsize=(7, 4.5))
    xlim = np.percentile(mflat, [0.1, 99.9])
    ax.hist(mflat, bins=120, color="#c44e52", alpha=0.85, density=True)
    ax.set_xlim(xlim)
    ax.set_xlabel("M value = ln(β/(1−β))")
    ax.set_ylabel("Density")
    ax.set_title("Distribution of M values (logit-transformed)")
    plt.tight_layout(); plt.savefig(r"D:\ONT\diagnosis_M_dist.png", dpi=130); plt.close()

    # Fig3: 均值-方差 对比 (核心诊断)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].scatter(bmean, bvar, s=1, c="#4c72b0", alpha=0.35, rasterized=True)
    axes[0].plot(bc, bv, color="black", lw=2, label="binned median var")
    axes[0].set_xlabel("mean β (per CpG)")
    axes[0].set_ylabel("variance β (per CpG)")
    axes[0].set_title("β: variance peaks at 0.5 (heteroscedastic)")
    axes[0].legend()
    axes[1].scatter(mmean, mvar, s=1, c="#c44e52", alpha=0.35, rasterized=True)
    axes[1].plot(mc, mv, color="black", lw=2, label="binned median var")
    axes[1].set_xlabel("mean M (per CpG)")
    axes[1].set_ylabel("variance M (per CpG)")
    axes[1].set_title("M: variance ~ flat (homoscedastic)")
    axes[1].legend()
    plt.tight_layout(); plt.savefig(r"D:\ONT\diagnosis_meanvar.png", dpi=130); plt.close()

    print("\n[save] diagnosis_beta_dist.png / diagnosis_M_dist.png / diagnosis_meanvar.png")

if __name__ == "__main__":
    main()
