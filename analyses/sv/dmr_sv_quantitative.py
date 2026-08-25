#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
阶段2 定量验证: SV 携带是否预测 hyper-DMR 的甲基化(读 11GB 矩阵)
================================================================
对最强 K 个 hyper-DMR:
  1) 抽 DMR 内 CpG, 算每样本的平均 M 值;
  2) 算每样本在 DMR±100kb 内的 n_SV(携带的不同 SV 数);
  3) 回归 DMR_M ~ n_SV + Gestational_Week (+ abnormal):
       (a) 全样本, 控制 abnormal;
       (b) 仅 non-abnormal case(448), 看 SV 单独是否预测高甲基化(最干净的因果检验)。
"""
import re
import numpy as np
import pandas as pd
from bisect import bisect_left, bisect_right
from collections import defaultdict
from scipy import stats

MATRIX   = r"E:\甲基化数据矩阵\甲基化CpG位点矩阵.txt"
CLINICAL = r"D:\ONT\clinical_649.tsv"
COV      = r"D:\ONT\matrix_covariates.tsv"
FID_FILE = r"D:\ONT\matrix_fid.txt"
SV_CAR   = r"D:\ONT\sv_carriers.tsv"
SV_BED   = r"E:\genotype_data\liftover\GRCh38_breakpoints.unique.bed"
DMR_FILES = {
    "spl_g8":  r"D:\ONT\figure2\abnormal_spl_g8\segments_genome.bed",
    "spl_g10": r"D:\ONT\figure2\abnormal_spl_g10\segments_genome.bed",
    "rpl_g8":  r"D:\ONT\figure2\abnormal_rpl_g8\segments_genome.bed",
    "rpl_g9":  r"D:\ONT\figure2\abnormal_rpl_g9\segments_genome.bed",
    "rpl_g10": r"D:\ONT\figure2\abnormal_rpl_g10\segments_genome.bed",
}
N_HEADER  = 71
TOP_FRAC  = 0.10
MIN_SITES = 3
K_DMR     = 500        # 最强 K 个合并 DMR
WINDOW    = 100_000
OUT       = r"D:\ONT\dmr_sv_quantitative.tsv"


def chrom_map(acc):
    m = re.match(r"NC_(\d+)\.", acc)
    if not m:
        return None
    n = int(m.group(1))
    return ("chr%d" % n) if 1 <= n <= 22 else None


def mtransform(v):
    m = np.nanmean(v)
    v = np.where(np.isnan(v), m, v)
    v = np.clip(v, 1e-3, 1.0 - 1e-3)
    return np.log2(v / (1.0 - v))


def ols(y, X):
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    sigma2 = float(resid @ resid) / dof if dof > 0 else np.nan
    se = np.sqrt(sigma2 * np.diag(XtX_inv))
    return beta, se, dof


def load_dmr_topk():
    """合并 5 文件的 hyper-DMR(去重叠), 取 |effect| 最强 K 个。"""
    intervals = []
    for path in DMR_FILES.values():
        d = pd.read_csv(path, sep="\t", comment="#", header=None)
        d.columns = ["chrom","start","end","name","score","num_sites",
                     "a_counts","b_counts","a_pct","b_pct","a_frac","b_frac",
                     "effect_size","cohen_h","ch_low","ch_high"]
        d = d[(d["name"] == "different") & (d["effect_size"] < 0)
              & (d["num_sites"] >= MIN_SITES)].copy()
        d["chrom"] = d["chrom"].map(chrom_map)
        d = d.dropna(subset=["chrom"])
        d = d.sort_values("effect_size").head(max(1, int(len(d) * TOP_FRAC)))
        for r in d.itertuples(index=False):
            intervals.append((r.chrom, int(r.start), int(r.end), abs(float(r.effect_size))))
    # 按 chrom 合并重叠
    by = defaultdict(list)
    for c, s, e, f in intervals:
        by[c].append((s, e, f))
    merged = []
    for c, lst in by.items():
        lst.sort()
        cs, ce, cf = lst[0]
        for s, e, f in lst[1:]:
            if s <= ce:
                ce = max(ce, e); cf = max(cf, f)
            else:
                merged.append((c, cs, ce, cf)); cs, ce, cf = s, e, f
        merged.append((c, cs, ce, cf))
    merged.sort(key=lambda x: -x[3])
    return merged[:K_DMR]


def main():
    # ---- 样本: 顺序 + 分组 + GW ----
    fids = [l.strip() for l in open(FID_FILE, encoding="utf-8") if l.strip()]
    n = len(fids)
    clin = pd.read_csv(CLINICAL, sep="\t", dtype={"Sample_ID": str})
    clin_map = clin.set_index("Sample_ID")
    abnormal = set(clin.loc[clin.Group4 == "abnormal", "Sample_ID"])
    normal_case = set(clin.loc[(clin.Group1 == "case") & (clin.Group4 != "abnormal"), "Sample_ID"])
    cov = pd.read_csv(COV, sep="\t")
    cov["FID"] = cov["FID"].astype(str)
    gw_map = dict(zip(cov["FID"], cov["Gestational_Week"].astype(float)))

    is_abnormal = np.array([1 if f in abnormal else 0 for f in fids], dtype=float)
    is_normal_case = np.array([1 if f in normal_case else 0 for f in fids], dtype=float)
    gw = np.array([gw_map.get(f, np.nan) for f in fids], dtype=float)
    gw = np.where(np.isnan(gw), np.nanmean(gw), gw)
    print("[samples] n=%d abnormal=%d normal_case=%d control=%d"
          % (n, int(is_abnormal.sum()), int(is_normal_case.sum()), int((1-is_abnormal-is_normal_case).sum())), flush=True)

    # ---- top-K DMR ----
    top = load_dmr_topk()
    print("[DMR] merged top %d" % len(top), flush=True)
    # 索引: chrom -> (starts, ends, dmr_id)
    idx = defaultdict(lambda: ([], [], []))
    for di, (c, s, e, f) in enumerate(top):
        idx[c][0].append(s); idx[c][1].append(e); idx[c][2].append(di)
    for c in idx:
        order = np.argsort(idx[c][0])
        idx[c] = (np.array(idx[c][0])[order], np.array(idx[c][1])[order], np.array(idx[c][2])[order])

    # ---- per-sample SV 断点 (chrom -> [(pos, sv_idx)]) ----
    bed = pd.read_csv(SV_BED, sep="\t", header=None,
                      names=["chrom","start","end","svid","svtype","side"])
    bed["sv_idx"] = bed["svid"].str.extract(r"SV(\d+)")[0].astype(int)
    bed = bed[bed["svtype"] != "TRA"]
    sv_bp = defaultdict(list)
    for r in bed.itertuples(index=False):
        sv_bp[r.sv_idx].append((r.chrom, int(r.start), int(r.end)))
    car = pd.read_csv(SV_CAR, sep="\t")
    car["cset"] = car["carriers"].fillna("").map(lambda s: set(s.split(",")) if s else set())
    fid2idx = {f: i for i, f in enumerate(fids)}
    sample_bp = [defaultdict(list) for _ in range(n)]   # per sample: chrom -> [(pos, sv_idx)]
    for _, r in car.iterrows():
        sid = int(r["sv_idx"])
        bps = sv_bp.get(sid)
        if not bps:
            continue
        for f in r["cset"]:
            i = fid2idx.get(f)
            if i is None:
                continue
            for chrom, s, e in bps:
                sample_bp[i][chrom].append((s, sid))
    for i in range(n):
        for c in sample_bp[i]:
            sample_bp[i][c].sort()
    print("[SV] per-sample breakpoints built", flush=True)

    # ---- 流式读矩阵, 抽 DMR 内 CpG, 累加每样本 M ----
    sumM = np.zeros((len(top), n))
    cnt = np.zeros(len(top), dtype=int)
    n_read = 0
    reader = pd.read_csv(MATRIX, sep="\t", header=None, skiprows=N_HEADER,
                         chunksize=20000, na_values=["NA", ""], low_memory=False,
                         dtype={0: str})
    for chunk in reader:
        ids = chunk.iloc[:, 0].values
        B = chunk.iloc[:, 1:].to_numpy(float)
        # 解析 site -> (chrom, pos)
        acc_pos = pd.Series(ids).str.split(":", expand=True)
        acc = acc_pos[0]; pos = acc_pos[1].astype(int).values
        num = acc.str.split(".", expand=True)[0].str[3:].astype(int).values
        chrom = np.array([("chr%d" % x) if x <= 22 else "" for x in num], dtype=object)
        for c, (starts, ends, dmr_ids) in idx.items():
            mask = chrom == c
            if not mask.any():
                continue
            p = pos[mask]
            j = np.searchsorted(starts, p, side="right") - 1
            ok = (j >= 0) & (ends[j] >= p)
            if not ok.any():
                continue
            p_ok = p[ok]; j_ok = j[ok]
            for pp, jj, rowi in zip(p_ok, j_ok, np.where(mask)[0][ok]):
                di = int(dmr_ids[jj])
                Mv = mtransform(B[rowi])
                sumM[di] += Mv
                cnt[di] += 1
        n_read += len(ids)
        if n_read % 200000 < 20000:
            print("[matrix] read %d rows" % n_read, flush=True)

    meanM = np.zeros((len(top), n))
    for di in range(len(top)):
        if cnt[di] > 0:
            meanM[di] = sumM[di] / cnt[di]
    print("[matrix] done, DMR with CpGs = %d/%d" % (int((cnt > 0).sum()), len(top)), flush=True)

    # ---- per-DMR per-sample n_SV ----
    n_sv = np.zeros((len(top), n), dtype=int)
    for di, (c, s, e, f) in enumerate(top):
        for i in range(n):
            lst = sample_bp[i].get(c)
            if not lst:
                continue
            pos_lst = [x[0] for x in lst]
            lo = bisect_left(pos_lst, s - WINDOW)
            hi = bisect_right(pos_lst, e + WINDOW)
            if hi <= lo:
                continue
            n_sv[di, i] = len({x[1] for x in lst[lo:hi]})

    # ---- 回归 ----
    rows = []
    for di, (c, s, e, f) in enumerate(top):
        if cnt[di] == 0:
            continue
        y = meanM[di]
        x = n_sv[di].astype(float)
        # (a) 全样本, 控制 abnormal
        Xa = np.column_stack([np.ones(n), x, gw, is_abnormal])
        ba, sea, dofa = ols(y, Xa)
        pa = 2 * stats.t.sf(abs(ba[1] / sea[1]), dofa) if sea[1] > 0 else np.nan
        # (b) 仅 normal-case
        sel = is_normal_case == 1
        if sel.sum() > 10:
            yb = y[sel]; xb = x[sel]; gwb = gw[sel]
            Xb = np.column_stack([np.ones(len(yb)), xb, gwb])
            bb, seb, dofb = ols(yb, Xb)
            pb = 2 * stats.t.sf(abs(bb[1] / seb[1]), dofb) if seb[1] > 0 else np.nan
            eff_b = bb[1]
        else:
            pb = eff_b = np.nan
        rows.append(dict(chrom=c, start=s, end=e, effect=f, n_cpg=cnt[di],
                         n_sv_total=int(x.sum()), n_samples_with_sv=int((x > 0).sum()),
                         beta_all=ba[1], p_all=pa,
                         beta_normal=eff_b, p_normal=pb))

    res = pd.DataFrame(rows)
    res.to_csv(OUT, sep="\t", index=False)
    print("[回归] DMR 数 = %d" % len(res), flush=True)

    # ---- 汇总 ----
    def summarize(col, name):
        v = res[col].dropna()
        if len(v) == 0:
            print("  %s: 无有效值" % name); return
        pos = (v > 0).mean()
        sig = (res[col] < 0.05).sum() if col in res else 0
        print("  %s: n=%d, 正效应占比=%.1f%%, 中位 p=%.3g, p<0.05 个数=%d"
              % (name, len(v), 100*pos, np.median(res[col].dropna()), int((res[col] < 0.05).sum())))

    print("\n=== n_SV 效应汇总 ===")
    summarize("p_all", "全样本(控制abnormal)")
    summarize("p_normal", "仅normal-case")

    # BH-FDR on p_normal
    pn = res["p_normal"].dropna().values
    if len(pn) > 0:
        order = np.argsort(pn); sp = pn[order]
        q = np.empty(len(pn)); q[order] = np.minimum.accumulate(sp[::-1] * len(pn) / np.arange(len(pn), 0, -1))[::-1]
        n_sig = int((q < 0.10).sum())
        print("  [仅normal-case] BH-FDR<10%% 显著 DMR 数 = %d" % n_sig)


if __name__ == "__main__":
    main()
