#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
同方差诊断(正确口径): 方差 vs 自变量 n_SV
=========================================
横轴 = 自变量 n_SV(窗口内 SV 剂量), 纵轴 = 该 n_SV 水平下样本的残差方差。
对 β 与 M 分别拟合  y ~ n_SV + Gestational_Week, 取残差, 按 n_SV 分箱算方差。
"""
import sys
sys.path.insert(0, r'D:\ONT')
import numpy as np
import pandas as pd
from collections import defaultdict
from bisect import bisect_left
import sv_methylation_pipeline as p
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fids = [l.strip() for l in open(p.FID, encoding="utf-8") if l.strip()]
n = len(fids)
fid2idx = {f: i for i, f in enumerate(fids)}

cov = pd.read_csv(p.COV, sep="\t")
gw = cov["Gestational_Week"].to_numpy(float)

# ---- 每个 (gene, window) 的 n_SV mask ----
gpat = pd.read_csv(p.GPAT, sep="\t")
masks = {}
for _, r in gpat.iterrows():
    m = np.zeros(n)
    if isinstance(r.carriers, str) and r.carriers:
        for part in r.carriers.split(","):
            fid, cnt = part.rsplit(":", 1)
            if fid in fid2idx:
                m[fid2idx[fid]] = float(cnt)
    masks[(r.gene_id, r.window)] = m

# ---- 只保留 n_SV 分布较好的 (gene, window): 有 0/1/2/3+ 各档 ----
def spread_ok(m):
    vals = np.floor(m).astype(int)
    return (vals == 0).sum() > 50 and (vals == 1).sum() > 30 and (vals >= 2).sum() > 20

good_windows = {}
for k, m in masks.items():
    if spread_ok(m):
        good_windows[k] = m

# 基因 -> 其"好窗口"列表
valid_genes = defaultdict(list)
for (gid, win), m in good_windows.items():
    valid_genes[gid].append(win)

print("good gene-windows =", len(good_windows), "genes =", len(valid_genes), flush=True)

# ---- 读矩阵采样, 把 CpG 映射到基因 ----
ct, cg, cn = p.nearest_gene_tables(p.GWIN)
collected = defaultdict(list)   # gene_id -> [beta row, ...]
reader = pd.read_csv(p.MATRIX, sep="\t", header=None, skiprows=p.N_HEADER,
                     chunksize=20000, na_values=["NA", ""], low_memory=False, dtype={0: str})
n_read = 0
LIMIT = 80000
for chunk in reader:
    ids = chunk.iloc[:, 0].values
    B = chunk.iloc[:, 1:].to_numpy(float)
    for j in range(B.shape[0]):
        chrom, pos = p.parse_site(ids[j])
        tss = ct.get(chrom)
        if tss is None or len(tss) == 0:
            continue
        i = bisect_left(tss, pos)
        best, bestd = None, 1e18
        if i < len(tss) and int(tss[i]) - pos < bestd:
            best, bestd = i, int(tss[i]) - pos
        if i - 1 >= 0 and pos - int(tss[i - 1]) < bestd:
            best, bestd = i - 1, pos - int(tss[i - 1])
        if best is not None and bestd <= 1_000_000:
            gid = cg[chrom][best]
            if gid in valid_genes:
                collected[gid].append(B[j])
    n_read += len(ids)
    if n_read >= LIMIT:
        break

print("read rows =", n_read, "genes with CpGs =", len(collected), flush=True)

# ---- 对每个 (gene, window, CpG) 拟合, 按 n_SV 分箱累计残差方差 ----
def mtransform(v):
    mm = np.nanmean(v)
    v = np.where(np.isnan(v), mm, v)
    v = np.clip(v, 1e-3, 1.0 - 1e-3)
    return np.log(v / (1.0 - v))


def impute_mean(v):
    """与 mtransform 相同的缺失填补(只补均值, 不变换不裁剪), 供 β 公平对比。"""
    mm = np.nanmean(v)
    return np.where(np.isnan(v), mm, v)

def fit_resid(y, X):
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    return y - X @ beta

# 分箱: n_SV = 0, 1, 2, >=3
def bins_of(m):
    v = np.floor(m).astype(int)
    b = np.zeros_like(v)
    b[v >= 3] = 3
    b[(v >= 1) & (v < 3)] = v[(v >= 1) & (v < 3)]
    return b

bin_names = ["0", "1", "2", "3+"]
beta_var = {k: [] for k in range(4)}
M_var    = {k: [] for k in range(4)}

used_windows = 0
for gid, rows in collected.items():
    if len(rows) < 20:
        continue
    for win in valid_genes[gid][:1]:     # 每基因取 1 个窗口即可
        m = masks[(gid, win)]
        X = np.column_stack([np.ones(n), m, gw])
        bins = bins_of(m)
        for y in rows:
            if np.isnan(y).all():
                continue
            yb = impute_mean(y)     # β: 填补缺失后直接回归
            Mb = mtransform(y)      # M: 填补+裁剪+logit
            for name, yy in (("beta", yb), ("M", Mb)):
                r = fit_resid(yy, X)
                for k in range(4):
                    sel = bins == k
                    if sel.sum() < 5:
                        continue
                    var = float(np.var(r[sel]))
                    (beta_var if name == "beta" else M_var)[k].append(var)
        used_windows += 1
        break  # 每个基因只处理一次(一个窗口)

print("used gene-windows =", used_windows, flush=True)

# ---- 汇总并画图 ----
def agg(d):
    return [np.median(v) if v else np.nan for v in (d[0], d[1], d[2], d[3])]

b_med = agg(beta_var)
m_med = agg(M_var)

print("=== 各 n_SV 箱残差方差(中位数) ===")
print("bin      :", bin_names)
print("beta var :", ["%.4g" % x for x in b_med])
print("M    var :", ["%.4g" % x for x in m_med])
print("估计数/箱 :", ["%d" % len(beta_var[k]) for k in range(4)])

fig, ax = plt.subplots(figsize=(7.5, 4.8))
x = np.arange(4)
ax.plot(x, b_med, "o-", color="#4c72b0", lw=2, label="β (heteroscedastic?)")
ax.plot(x, m_med, "s-", color="#c44e52", lw=2, label="M (homoscedastic?)")
ax.set_xticks(x); ax.set_xticklabels(["n_SV=0", "n_SV=1", "n_SV=2", "n_SV≥3"])
ax.set_xlabel("independent variable  n_SV  (SV dose in window)")
ax.set_ylabel("residual variance (median over CpG-window fits)")
ax.set_title("Heteroscedasticity check: variance of residuals vs n_SV")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(r"D:\ONT\diagnosis_var_vs_nsv.png", dpi=130); plt.close()
print("[save] diagnosis_var_vs_nsv.png")
