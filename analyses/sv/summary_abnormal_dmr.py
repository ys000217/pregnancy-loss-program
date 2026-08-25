#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""摘要 abnormal DMR 文件: 行数 / different(DMR)数 / effect_size 方向分布。"""
import glob, os
import pandas as pd

files = sorted(glob.glob(r"D:\ONT\figure2\abnormal_*\segments_genome.bed"))
print("%-24s %10s %10s %10s %10s" % ("dir", "total", "different", "effect>0", "effect<0"))
for f in files:
    d = pd.read_csv(f, sep="\t", comment="#", header=None)
    # 列: 0=chrom 1=start 2=end 3=name 4=score 5=num_sites ... 12=effect_size
    diff = d[d[3] == "different"]
    pos = int((diff[12] > 0).sum())
    neg = int((diff[12] < 0).sum())
    name = os.path.basename(os.path.dirname(f))
    print("%-24s %10d %10d %10d %10d" % (name, len(d), len(diff), pos, neg))
