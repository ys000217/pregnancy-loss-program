import numpy as np
import pandas as pd
from collections import defaultdict
from scipy import stats

MATRIX = r"E:\甲基化数据矩阵\EWAS_INPUT_NO_HEADER.txt"
COV = r"D:\ONT\matrix_covariates.tsv"
FID = r"D:\ONT\matrix_fid.txt"
GWIN = r"D:\ONT\gene_window_patients.tsv"
SIG = r"D:\ONT\sv_ewas_sig_pairs.tsv"
OUT = r"D:\ONT\interaction_results.tsv"
N_HEADER = 71

fids = [l.strip() for l in open(FID, encoding="utf-8") if l.strip()]
n = len(fids)
fid2idx = {f: i for i, f in enumerate(fids)}
cov = pd.read_csv(COV, sep="\t")
status = cov["Status"].to_numpy(float)
cov_names = ["Age", "Gestational_Week", "Endothelial", "Hofbauer", "nRBC", "Stromal"]
Xcov = cov[cov_names].to_numpy(float)

gwin = pd.read_csv(GWIN, sep="\t")
masks = {}
for _, r in gwin.iterrows():
    if r.n_carriers < 33 or r.n_carriers > 615:
        continue
    mask = np.zeros(n)
    if isinstance(r.carriers, str) and r.carriers:
        for part in r.carriers.split(","):
            fid, cnt = part.rsplit(":", 1)
            if fid in fid2idx:
                mask[fid2idx[fid]] = float(cnt)
    masks[(r.gene_id, r.window)] = mask

sig = pd.read_csv(SIG, sep="\t")
site_info = defaultdict(list)
for _, r in sig.iterrows():
    site_info[r.site].append((r.gene_id, r.gene_name, r.window, int(r.n_carriers)))
print("sites to reprocess:", len(site_info), flush=True)

def mtransform(v):
    m = np.nanmean(v)
    v = np.where(np.isnan(v), m, v)
    v = np.clip(v, 1e-3, 1 - 1e-3)
    return np.log(v / (1.0 - v))

dof = n - 10
out = open(OUT, "w", encoding="utf-8")
out.write("site\tgene_id\tgene_name\twindow\tsv_effect\tsv_p\tint_effect\tint_p\tn_case_car\tn_case_non\tn_ctrl_car\tn_ctrl_non\n")
n_hit = 0
n_skip = 0
reader = pd.read_csv(MATRIX, sep="\t", header=None, skiprows=N_HEADER, chunksize=20000,
                     na_values=["NA", ""], low_memory=False, dtype={0: str})
for chunk in reader:
    ids = chunk.iloc[:, 0].values
    Y = chunk.iloc[:, 1:].to_numpy(float)
    for j in range(Y.shape[0]):
        site = ids[j]
        if site not in site_info:
            continue
        M = mtransform(Y[j])
        for gid, gnm, win, ncar in site_info[site]:
            mask = masks.get((gid, win))
            if mask is None:
                continue
            case_car = int(((mask >= 1) & (status == 1)).sum())
            case_non = int(((mask == 0) & (status == 1)).sum())
            ctrl_car = int(((mask >= 1) & (status == 0)).sum())
            ctrl_non = int(((mask == 0) & (status == 0)).sum())
            X = np.column_stack([np.ones(n), mask, status, mask * status, Xcov])
            XtX_inv = np.linalg.pinv(X.T @ X)
            beta = XtX_inv @ X.T @ M
            resid = M - X @ beta
            sigma2 = float(resid @ resid) / dof
            se = np.sqrt(sigma2 * np.diag(XtX_inv))
            t_sv = beta[1] / se[1] if se[1] > 0 else 0.0
            p_sv = float(2 * stats.t.sf(abs(t_sv), dof))
            estimable = (case_car >= 2 and case_non >= 2 and ctrl_car >= 2 and ctrl_non >= 2)
            if estimable:
                t_int = beta[3] / se[3] if se[3] > 0 else 0.0
                p_int = float(2 * stats.t.sf(abs(t_int), dof))
                int_eff = float(beta[3])
            else:
                p_int = float("nan")
                int_eff = float("nan")
                n_skip += 1
            out.write("%s\t%s\t%s\t%s\t%.6g\t%.6g\t%.6g\t%.6g\t%d\t%d\t%d\t%d\n"
                      % (site, gid, gnm, win, float(beta[1]), p_sv, int_eff, p_int,
                         case_car, case_non, ctrl_car, ctrl_non))
            n_hit += 1
out.close()
print("DONE interaction pairs=%d (skipped non-estimable=%d)" % (n_hit, n_skip), flush=True)
