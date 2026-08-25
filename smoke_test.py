import numpy as np
import pandas as pd

MATRIX = r"E:\甲基化数据矩阵\EWAS_INPUT_NO_HEADER.txt"
COV = r"D:\ONT\matrix_covariates.tsv"
N_HEADER = 71

# 1) manually read first 2000 data rows (fast, no full-file scan)
rows = []
n = 0
with open(MATRIX, encoding="utf-8") as f:
    for i, line in enumerate(f):
        if i < N_HEADER:
            continue
        rows.append(line.rstrip("\n").split("\t"))
        n += 1
        if n >= 2000:
            break
print("read rows:", len(rows), "ncols:", len(rows[0]))

sids = [r[0] for r in rows]
vals = np.array([[x if x not in ("NA", "") else "nan" for x in r[1:]] for r in rows], dtype=float)
print("vals shape:", vals.shape)
print("site ids first/last:", sids[0], sids[-1])
print("value min/max/mean:", np.nanmin(vals), np.nanmax(vals), np.nanmean(vals))
na_per_site = np.isnan(vals).sum(axis=1)
print("NA per site: min=%d max=%d mean=%.1f" % (na_per_site.min(), na_per_site.max(), na_per_site.mean()))
print("sites 0 NA:", (na_per_site == 0).sum(), "/ 2000 ; <=5 NA:", (na_per_site <= 5).sum(), "/ 2000")

# chromosome coverage from site ids (NC_ accession)
from collections import Counter
accs = [s.split(":")[0] for s in sids]
print("distinct accessions:", sorted(set(accs))[:5], "... total", len(set(accs)))

# 2) load covariates, check NA
cov = pd.read_csv(COV, sep="\t")
print("cov shape:", cov.shape)
print("cov cols:", list(cov.columns))
print("cov NA counts:")
print(cov.isna().sum()[cov.isna().sum() > 0].to_string())
print("Status: 1=%d 0=%d" % ((cov.Status == 1).sum(), (cov.Status == 0).sum()))
