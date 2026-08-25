#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""量化: 边界堆积到底导致多少 CpG 退化(不变量/近不变量), 用于评估对回归的实际影响。"""
import numpy as np
import pandas as pd

MATRIX = r"E:\甲基化数据矩阵\EWAS_INPUT_NO_HEADER.txt"
N = 30000

df = pd.read_csv(MATRIX, sep="\t", header=None, skiprows=71, nrows=N,
                 na_values=["NA", ""], low_memory=False, dtype={0: str})
B = df.iloc[:, 1:].to_numpy(float)          # (N, 648)
n_sites, n_samp = B.shape

bvar = np.nanvar(B, axis=1)
bmean = np.nanmean(B, axis=1)

# 每 CpG 落在精确 0/1 边界的样本占比
boundary = (B == 0.0) | (B == 1.0)
per_site_boundary = boundary.sum(axis=1) / n_samp

fully_inv = (bvar <= 1e-12)                 # 完全无变异
near_inv  = (bvar < 1e-4)                   # 近无变异(几乎测不出)
heavy_bnd = (per_site_boundary >= 0.95)     # >=95% 样本都在 0/1 边界
usable     = (bvar >= 0.01)                 # 有足够变异可用于回归

print("sampled sites  = %d" % n_sites)
print("fully invariant (var=0)           : %.2f%%" % (100 * fully_inv.mean()))
print("near-invariant   (var<1e-4)        : %.2f%%" % (100 * near_inv.mean()))
print(">=95%% samples at exact 0/1        : %.2f%%" % (100 * heavy_bnd.mean()))
print("usable (var>=0.01)                 : %.2f%%" % (100 * usable.mean()))

# 双峰程度: 有变异的 CpG 中, 两侧都接近边界的比例
vary = bvar >= 0.01
both_bnd = (boundary.sum(axis=1) / n_samp) >= 0.5
print("among usable, >=50%% samples at 0/1 : %.2f%%" % (100 * (both_bnd[vary].mean() if vary.any() else 0)))
