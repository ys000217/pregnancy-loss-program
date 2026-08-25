#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Randomly subset the East Asian RPL meta sumstats to ~250k HapMap3 SNPs
(for gsMap runtime), verify mean chi^2 is preserved, write new sumstats."""
import gzip
import numpy as np
import pandas as pd

SRC = "/mnt/d/gsMap/RPL_GWAS/RPL_meta_EA.sumstats.gz"
OUT = "/mnt/d/gsMap/RPL_GWAS/RPL_meta_EA_250k.sumstats.gz"
SEED = 42
TARGET = 250_000

df = pd.read_csv(SRC, sep="\t")
print(f"full sumstats SNPs: {len(df):,}")
full_mean = float((df["Z"] ** 2).mean())
print(f"full mean chi^2   : {full_mean:.6f}")

n = min(TARGET, len(df))
sub = df.sample(n=n, random_state=SEED, replace=False)
sub_mean = float((sub["Z"] ** 2).mean())
sub_median = float(np.median(sub["Z"] ** 2))
print(f"subset SNPs       : {len(sub):,}")
print(f"subset mean chi^2 : {sub_mean:.6f}")
print(f"subset lambda_GC  : {sub_median / 0.4549:.4f}")

sub.to_csv(OUT, sep="\t", index=False, compression="gzip")
print(f"wrote {OUT}")
