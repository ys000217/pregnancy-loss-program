#!/usr/bin/env python3
"""DN-background SNV: abnormal vs normal analysis using raw columns only.

Universe: loose dual-filter candidates with control==0 and (abnormal or normal
carrier). Derived columns (clinvar_tiers/stratum/score_*) are not used in tests.

Phenotype default: 47 label-abnormal + 30 CpG/NC-cluster suspects -> 77 abnormal
(tables/suspect_abn/phenotype.tsv).

Enrichment FDR (valid full-union families):
  - two-sided: H1 rates differ; BH on all union sites / all hit genes
  - one-sided ab: H1 ab rate higher; BH on same full union
  - one-sided normal: H1 normal rate higher; BH on same full union
Direction is used only for post-hoc pattern labels, not to shrink m.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from snv_loose_filter_rules import _parse_float  # noqa: E402

VCF = Path(r"D:/ONT/figure2/WGS_ONT_Intersection_648samples.vcf.gz")
PHENO = ROOT / "tables" / "suspect_abn" / "phenotype.tsv"
PHENO_FALLBACK = Path(r"D:/ONT/figure2/sample_phenotype_648.tsv")
CAND = ROOT / "tables" / "loose_dual_filter" / "loose_candidates.tsv"
OUT = ROOT / "tables" / "loose_dual_filter"
PLOT_DIR = ROOT / "plots" / "dn_ab_vs_normal"
GENE_LEN = ROOT / "tables" / "loose_dual_filter" / "gene_cds_length_bp.tsv"

OUT.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)


def norm_chr(c: str) -> str:
    c = str(c).strip()
    return c[3:] if c.lower().startswith("chr") else c


def parse_af(x) -> float | None:
    if x is None or str(x).strip() in {"", ".", "NA", "na", "nan", "-"}:
        return None
    try:
        return float(x)
    except ValueError:
        return _parse_float(x)


def load_pheno(path: Path) -> dict[str, str]:
    pheno: dict[str, str] = {}
    with open(path, encoding="utf-8") as f:
        h = f.readline().rstrip("\n").split("\t")
        si, ci = h.index("VCF_Sample_ID"), h.index("Condition")
        for line in f:
            p = line.rstrip("\n").split("\t")
            c = p[ci].strip().lower()
            if c in {"abnormal", "case"}:
                c = "abnormal"
            elif c in {"normal", "normal_case"}:
                c = "normal"
            pheno[p[si]] = c
    return pheno


def fisher_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact via enumerating tables with same margins."""
    # table [[a,b],[c,d]]; row1=a+b, col1=a+c
    n = a + b + c + d
    if n == 0:
        return 1.0
    row1 = a + b
    col1 = a + c
    lo = max(0, row1 + col1 - n)
    hi = min(row1, col1)

    def dens(x: int) -> float:
        # comb(col1,x)*comb(n-col1, row1-x) / comb(n, row1)
        return (
            math.comb(col1, x)
            * math.comb(n - col1, row1 - x)
            / math.comb(n, row1)
        )

    p_obs = dens(a)
    s = 0.0
    for x in range(lo, hi + 1):
        p = dens(x)
        if p <= p_obs + 1e-15:
            s += p
    return min(1.0, s)


def fisher_one_sided_row1_greater(a: int, b: int, c: int, d: int) -> float:
    """One-sided Fisher: P(X>=a) for higher carrier rate in row1 [[a,b],[c,d]]."""
    n = a + b + c + d
    if n == 0:
        return 1.0
    k = a + b
    K = a + c
    hi = min(k, K)
    total = math.comb(n, k)
    if total == 0:
        return 1.0
    s = 0
    for x in range(a, hi + 1):
        s += math.comb(K, x) * math.comb(n - K, k - x)
    return min(1.0, s / total)


def fisher_or(a: int, b: int, c: int, d: int) -> float:
    """OR with Haldane-Anscombe if zero cell."""
    if 0 in (a, b, c, d):
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    return (a / b) / (c / d)


def bh_fdr(pvals: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR q-values aligned to input order."""
    m = len(pvals)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for k in range(m, 0, -1):
        i = order[k - 1]
        val = pvals[i] * m / k
        prev = min(prev, val)
        q[i] = min(1.0, prev)
    return q


def mann_whitney_u(x: list[float], y: list[float]) -> tuple[float, float]:
    """Two-sided MW U via normal approx with tie correction; returns (U, p)."""
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return float("nan"), 1.0
    combined = [(v, 0) for v in x] + [(v, 1) for v in y]
    combined.sort(key=lambda t: t[0])
    ranks = [0.0] * len(combined)
    i = 0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]:
            j += 1
        avg = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[k] = avg
        i = j
    r1 = sum(ranks[k] for k, t in enumerate(combined) if t[1] == 0)
    u1 = r1 - n1 * (n1 + 1) / 2.0
    u2 = n1 * n2 - u1
    u = min(u1, u2)
    mu = n1 * n2 / 2.0
    # tie correction
    tie_sum = 0
    i = 0
    vals = [t[0] for t in combined]
    while i < len(vals):
        j = i
        while j < len(vals) and vals[j] == vals[i]:
            j += 1
        t = j - i
        if t > 1:
            tie_sum += t * t * t - t
        i = j
    n = n1 + n2
    sigma2 = n1 * n2 / 12.0 * ((n + 1) - tie_sum / (n * (n - 1))) if n > 1 else 0.0
    if sigma2 <= 0:
        return u, 1.0
    sigma = math.sqrt(sigma2)
    z = (u - mu) / sigma
    # two-sided normal survival
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return u, float(p)


def gene_tokens(field: str) -> list[str]:
    out = []
    for g in str(field or "").replace(";", ",").split(","):
        g = g.strip()
        if g and g not in {".", "NA", "nan"}:
            out.append(g)
    return out


def load_union_sites(path: Path, universe: str = "control0") -> dict[tuple, dict]:
    sites: dict[tuple, dict] = {}
    with open(path, encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            a = int(row["n_abnormal"])
            m = int(row["n_normal"])
            c = int(row["n_control"])
            if universe == "control0" and c != 0:
                continue
            if a == 0 and m == 0:
                continue
            key = (norm_chr(row["chrom"]), row["pos"], row["ref"], row["alt"])
            sites[key] = {
                "chrom": row["chrom"],
                "pos": row["pos"],
                "ref": row["ref"],
                "alt": row["alt"],
                "gene": row.get("gene", ""),
                "func_class": row.get("consequence", ""),
                "clnsig": row.get("clnsig", ""),
                "gnomad_af": row.get("gnomad_af", ""),
                "sift": row.get("sift", ""),
                "polyphen_hdiv": row.get("polyphen_hdiv", ""),
                "cadd": row.get("cadd", ""),
                "revel": row.get("revel", ""),
                "n_abnormal_pre": a,
                "n_normal_pre": m,
                "n_control_pre": c,
            }
    return sites


def scan_vcf_carriers(
    vcf: Path,
    sites: dict[tuple, dict],
    pheno: dict[str, str],
) -> dict[tuple, dict[str, str]]:
    """Return key -> {sample_id: gt} for carriers in phenotype set."""
    carriers: dict[tuple, dict[str, str]] = defaultdict(dict)
    with gzip.open(vcf, "rt", encoding="utf-8", errors="replace") as f:
        samp_idx = None
        n_hit = 0
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                vs = line.strip().split("\t")[9:]
                samp_idx = [(j, s) for j, s in enumerate(vs) if s in pheno]
                print(f"VCF phenotype samples indexed: {len(samp_idx)}", flush=True)
                continue
            p = line.rstrip("\n").split("\t")
            if p[6] != "PASS":
                continue
            key = (norm_chr(p[0]), p[1], p[3], p[4])
            if key not in sites:
                continue
            n_hit += 1
            assert samp_idx is not None
            gts = p[9:]
            for j, s in samp_idx:
                gt = gts[j].split(":")[0]
                if gt in {".", "./.", ".|."}:
                    continue
                alleles = gt.replace("|", "/").split("/")
                if any(a not in {"0", "."} for a in alleles):
                    carriers[key][s] = gt
            if n_hit % 20000 == 0:
                print(f"  ... matched union sites: {n_hit}", flush=True)
    print(f"VCF matched union sites: {n_hit}", flush=True)
    return carriers


def classify_pattern(
    n_ab: int,
    n_norm: int,
    rate_ab: float,
    rate_norm: float,
    fdr_two: float,
    fdr_ab: float,
    fdr_norm: float,
    alpha: float = 0.05,
) -> str:
    """Post-hoc label only; FDR families are always the full union."""
    if fdr_ab < alpha and rate_ab > rate_norm:
        return "abnormal_enriched"
    if fdr_norm < alpha and rate_norm > rate_ab:
        return "normal_enriched"
    if fdr_two < alpha:
        return "differential_two_sided"
    if n_ab > 0 and n_norm > 0:
        return "shared"
    if n_ab > 0 and n_norm == 0:
        return "abnormal_only"
    if n_norm > 0 and n_ab == 0:
        return "normal_only"
    return "other"


def apply_full_union_fdr(
    rows: list[dict],
    *,
    n_ab_key: str = "n_abnormal",
    n_norm_key: str = "n_normal",
) -> None:
    """BH on the full row set for two-sided and each one-sided test (m = n rows)."""
    fdr_two = bh_fdr([r["p_fisher_two_sided"] for r in rows])
    fdr_ab = bh_fdr([r["p_ab_greater"] for r in rows])
    fdr_norm = bh_fdr([r["p_norm_greater"] for r in rows])
    for r, a, b, c in zip(rows, fdr_two, fdr_ab, fdr_norm):
        r["fdr_two_sided"] = a
        r["fdr_ab_greater"] = b
        r["fdr_norm_greater"] = c
        r["pattern"] = classify_pattern(
            int(r[n_ab_key]),
            int(r[n_norm_key]),
            float(r["rate_abnormal"]),
            float(r["rate_normal"]),
            a,
            b,
            c,
        )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="DN-background abnormal vs normal")
    p.add_argument(
        "--cand",
        default=str(CAND),
        help="loose candidate TSV",
    )
    p.add_argument("--out", default=str(OUT), help="output table directory")
    p.add_argument(
        "--universe",
        choices=["control0", "freq_gate"],
        default="control0",
        help="control0: require n_control==0; freq_gate: keep pre-gated candidates",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out)
    cand_path = Path(args.cand)
    out_dir.mkdir(parents=True, exist_ok=True)
    global OUT
    OUT = out_dir

    pheno_path = PHENO if PHENO.exists() else PHENO_FALLBACK
    pheno = load_pheno(pheno_path)
    ab_set = {s for s, g in pheno.items() if g == "abnormal"}
    norm_set = {s for s, g in pheno.items() if g == "normal"}
    ctrl_set = {s for s, g in pheno.items() if g == "control"}
    n_ab_tot, n_norm_tot, n_ctrl_tot = len(ab_set), len(norm_set), len(ctrl_set)
    print(f"PHENOTYPE_FILE {pheno_path}", flush=True)
    print(
        f"PHENOTYPE abnormal={n_ab_tot} normal={n_norm_tot} control={n_ctrl_tot}",
        flush=True,
    )

    print("Loading DN union from candidate TSV ...", flush=True)
    sites = load_union_sites(cand_path, args.universe)
    print(
        f"Union sites ({args.universe}; ab or normal carrier): {len(sites)}",
        flush=True,
    )

    print("Scanning VCF for carriers ...", flush=True)
    carriers = scan_vcf_carriers(VCF, sites, pheno)

    # --- site-level table + patterns ---
    site_rows = []
    long_rows = []
    sample_n_sites: Counter[str] = Counter()
    sample_genes: dict[str, set[str]] = defaultdict(set)
    gene_samples: dict[str, set[str]] = defaultdict(set)
    gene_n_sites: Counter[str] = Counter()
    gene_site_keys: dict[str, set[tuple]] = defaultdict(set)

    for key, meta in sites.items():
        carr = carriers.get(key, {})
        ab_ids = [s for s in carr if s in ab_set]
        norm_ids = [s for s in carr if s in norm_set]
        ctrl_ids = [s for s in carr if s in ctrl_set]
        n_ab, n_norm, n_ctrl = len(ab_ids), len(norm_ids), len(ctrl_ids)
        if n_ab == 0 and n_norm == 0:
            continue
        # 2x2: ab carrier, ab non, norm carrier, norm non
        a, b = n_ab, n_ab_tot - n_ab
        c, d = n_norm, n_norm_tot - n_norm
        p_two = fisher_two_sided(a, b, c, d)
        p_ab = fisher_one_sided_row1_greater(a, b, c, d)
        p_norm = fisher_one_sided_row1_greater(c, d, a, b)
        or_ = fisher_or(a, b, c, d)
        rate_ab = n_ab / n_ab_tot if n_ab_tot else 0.0
        rate_norm = n_norm / n_norm_tot if n_norm_tot else 0.0
        genes = gene_tokens(meta["gene"])
        for s, gt in carr.items():
            if s not in ab_set and s not in norm_set:
                continue
            long_rows.append(
                {
                    "chrom": meta["chrom"],
                    "pos": meta["pos"],
                    "ref": meta["ref"],
                    "alt": meta["alt"],
                    "gene": meta["gene"],
                    "func_class": meta.get("func_class", ""),
                    "VCF_Sample_ID": s,
                    "Condition": pheno[s],
                    "GT": gt,
                    "clnsig": meta["clnsig"],
                    "gnomad_af": meta["gnomad_af"],
                    "revel": meta["revel"],
                    "cadd": meta["cadd"],
                    "sift": meta["sift"],
                    "polyphen_hdiv": meta["polyphen_hdiv"],
                }
            )
            sample_n_sites[s] += 1
            for g in genes:
                sample_genes[s].add(g)
                gene_samples[g].add(s)
                gene_site_keys[g].add(key)
        for g in genes:
            gene_n_sites[g] += 1

        site_rows.append(
            {
                "chrom": meta["chrom"],
                "pos": meta["pos"],
                "ref": meta["ref"],
                "alt": meta["alt"],
                "gene": meta["gene"],
                "clnsig": meta["clnsig"],
                "gnomad_af": meta["gnomad_af"],
                "func_class": meta.get("func_class", ""),
                "revel": meta["revel"],
                "cadd": meta["cadd"],
                "sift": meta["sift"],
                "polyphen_hdiv": meta["polyphen_hdiv"],
                "n_abnormal": n_ab,
                "n_normal": n_norm,
                "n_control": n_ctrl,
                "rate_abnormal": rate_ab,
                "rate_normal": rate_norm,
                "or_ab_vs_normal": or_,
                "p_fisher_two_sided": p_two,
                "p_ab_greater": p_ab,
                "p_norm_greater": p_norm,
            }
        )

    apply_full_union_fdr(site_rows, n_ab_key="n_abnormal", n_norm_key="n_normal")
    site_rows.sort(
        key=lambda r: (
            min(r["fdr_two_sided"], r["fdr_ab_greater"], r["fdr_norm_greater"]),
            r["p_fisher_two_sided"],
        )
    )

    # write site patterns
    site_path = OUT / "dn_site_patterns_ab_vs_normal.tsv"
    with open(site_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(site_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(site_rows)
    print(
        f"WROTE {site_path} n={len(site_rows)} FDR_m_full_union={len(site_rows)}",
        flush=True,
    )

    # long table can be large — write it
    long_path = OUT / "dn_carrier_long_ab_normal.tsv"
    with open(long_path, "w", encoding="utf-8", newline="") as f:
        if long_rows:
            w = csv.DictWriter(f, fieldnames=list(long_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(long_rows)
    print(f"WROTE {long_path} n={len(long_rows)}", flush=True)

    # --- gene burden ---
    gene_len: dict[str, int] = {}
    if GENE_LEN.exists():
        with open(GENE_LEN, encoding="utf-8") as f:
            r = csv.DictReader(f, delimiter="\t")
            for row in r:
                gene_len[row["gene"]] = int(float(row["cds_bp"]))

    gene_rows = []
    for g, samps in gene_samples.items():
        ab_c = sum(1 for s in samps if s in ab_set)
        norm_c = sum(1 for s in samps if s in norm_set)
        a, b = ab_c, n_ab_tot - ab_c
        c, d = norm_c, n_norm_tot - norm_c
        p_two = fisher_two_sided(a, b, c, d)
        p_ab = fisher_one_sided_row1_greater(a, b, c, d)
        p_norm = fisher_one_sided_row1_greater(c, d, a, b)
        or_ = fisher_or(a, b, c, d)
        cds = gene_len.get(g)
        # length-normalized carrier rates (per 100kb CDS)
        scale = (cds / 1e5) if cds and cds > 0 else float("nan")
        rate_ab = ab_c / n_ab_tot
        rate_norm = norm_c / n_norm_tot
        rate_ab_per_100kb = rate_ab / scale if scale == scale and scale > 0 else float("nan")
        rate_norm_per_100kb = rate_norm / scale if scale == scale and scale > 0 else float("nan")
        gene_rows.append(
            {
                "gene": g,
                "n_dn_sites": len(gene_site_keys[g]),
                "cds_bp": cds if cds is not None else "",
                "n_abnormal_carriers": ab_c,
                "n_normal_carriers": norm_c,
                "rate_abnormal": rate_ab,
                "rate_normal": rate_norm,
                "rate_ab_per_100kb": rate_ab_per_100kb if rate_ab_per_100kb == rate_ab_per_100kb else "",
                "rate_norm_per_100kb": rate_norm_per_100kb if rate_norm_per_100kb == rate_norm_per_100kb else "",
                "or_ab_vs_normal": or_,
                "p_fisher_two_sided": p_two,
                "p_ab_greater": p_ab,
                "p_norm_greater": p_norm,
            }
        )

    apply_full_union_fdr(
        gene_rows, n_ab_key="n_abnormal_carriers", n_norm_key="n_normal_carriers"
    )
    for r in gene_rows:
        ra = r["rate_ab_per_100kb"]
        rn = r["rate_norm_per_100kb"]
        if ra != "" and rn != "":
            if float(ra) > float(rn):
                r["length_adj_direction"] = "abnormal_higher"
            elif float(rn) > float(ra):
                r["length_adj_direction"] = "normal_higher"
            else:
                r["length_adj_direction"] = "equal"
        else:
            r["length_adj_direction"] = "no_cds_length"

    gene_rows.sort(
        key=lambda r: (
            min(r["fdr_two_sided"], r["fdr_ab_greater"], r["fdr_norm_greater"]),
            r["p_fisher_two_sided"],
        )
    )
    gene_path = OUT / "dn_gene_burden_ab_vs_normal.tsv"
    with open(gene_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(gene_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(gene_rows)
    print(
        f"WROTE {gene_path} n={len(gene_rows)} FDR_m_full_union={len(gene_rows)}",
        flush=True,
    )

    # gene x sample 0/1 (long sparse)
    gx_path = OUT / "dn_gene_sample_carriers.tsv"
    with open(gx_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene", "VCF_Sample_ID", "Condition"])
        for g, samps in sorted(gene_samples.items()):
            for s in sorted(samps):
                if s in ab_set or s in norm_set:
                    w.writerow([g, s, pheno[s]])
    print(f"WROTE {gx_path}", flush=True)

    # --- sample load ---
    sample_rows = []
    for s, g in pheno.items():
        if g not in {"abnormal", "normal"}:
            continue
        sample_rows.append(
            {
                "VCF_Sample_ID": s,
                "Condition": g,
                "n_dn_sites": sample_n_sites.get(s, 0),
                "n_dn_genes": len(sample_genes.get(s, ())),
            }
        )
    sample_path = OUT / "dn_sample_load_ab_normal.tsv"
    with open(sample_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(sample_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(sample_rows)

    ab_sites = [r["n_dn_sites"] for r in sample_rows if r["Condition"] == "abnormal"]
    norm_sites = [r["n_dn_sites"] for r in sample_rows if r["Condition"] == "normal"]
    ab_genes = [r["n_dn_genes"] for r in sample_rows if r["Condition"] == "abnormal"]
    norm_genes = [r["n_dn_genes"] for r in sample_rows if r["Condition"] == "normal"]
    u_s, p_s = mann_whitney_u(ab_sites, norm_sites)
    u_g, p_g = mann_whitney_u(ab_genes, norm_genes)

    def summarize(vals: list[float]) -> dict:
        if not vals:
            return {"n": 0}
        vs = sorted(vals)
        n = len(vs)
        def q(p):
            i = int(round((n - 1) * p))
            return vs[i]
        mean = sum(vs) / n
        var = sum((x - mean) ** 2 for x in vs) / n if n else 0
        return {
            "n": n,
            "mean": mean,
            "sd": math.sqrt(var),
            "median": q(0.5),
            "q25": q(0.25),
            "q75": q(0.75),
            "min": vs[0],
            "max": vs[-1],
        }

    sum_path = OUT / "dn_ab_vs_normal_summary.tsv"
    with open(sum_path, "w", encoding="utf-8") as f:
        f.write("metric\tvalue\n")
        f.write(f"phenotype_file\t{pheno_path}\n")
        f.write(f"candidate_file\t{cand_path}\n")
        f.write(f"universe\t{args.universe}\n")
        f.write(f"n_abnormal\t{n_ab_tot}\n")
        f.write(f"n_normal\t{n_norm_tot}\n")
        f.write("fdr_method\tfull_union_BH_two_sided_and_two_one_sided\n")
        f.write(
            "fdr_note\t"
            "For sites and genes separately: BH on full union for "
            "(1) two-sided difference p_fisher_two_sided -> fdr_two_sided; "
            "(2) one-sided ab higher p_ab_greater -> fdr_ab_greater; "
            "(3) one-sided normal higher p_norm_greater -> fdr_norm_greater. "
            "Direction is not used to shrink m.\n"
        )
        f.write(f"site_fdr_m\t{len(site_rows)}\n")
        f.write(f"gene_fdr_m\t{len(gene_rows)}\n")
        f.write(
            f"site_n_fdr05_two\t{sum(1 for r in site_rows if r['fdr_two_sided'] < 0.05)}\n"
        )
        f.write(
            f"site_n_fdr05_ab\t{sum(1 for r in site_rows if r['fdr_ab_greater'] < 0.05)}\n"
        )
        f.write(
            f"site_n_fdr05_norm\t{sum(1 for r in site_rows if r['fdr_norm_greater'] < 0.05)}\n"
        )
        f.write(
            f"gene_n_fdr05_two\t{sum(1 for r in gene_rows if r['fdr_two_sided'] < 0.05)}\n"
        )
        f.write(
            f"gene_n_fdr05_ab\t{sum(1 for r in gene_rows if r['fdr_ab_greater'] < 0.05)}\n"
        )
        f.write(
            f"gene_n_fdr05_norm\t{sum(1 for r in gene_rows if r['fdr_norm_greater'] < 0.05)}\n"
        )
        f.write(f"n_union_sites_scanned\t{len(sites)}\n")
        f.write(f"n_sites_with_ab_or_norm_carrier\t{len(site_rows)}\n")
        f.write(f"n_carrier_rows\t{len(long_rows)}\n")
        f.write(f"n_genes\t{len(gene_rows)}\n")
        site_pat = Counter(r["pattern"] for r in site_rows)
        for k, v in site_pat.most_common():
            f.write(f"site_pattern_{k}\t{v}\n")
        gene_pat = Counter(r["pattern"] for r in gene_rows)
        for k, v in gene_pat.most_common():
            f.write(f"gene_pattern_{k}\t{v}\n")
        f.write(f"mw_p_n_dn_sites\t{p_s}\n")
        f.write(f"mw_U_n_dn_sites\t{u_s}\n")
        f.write(f"mw_p_n_dn_genes\t{p_g}\n")
        for label, vals in [
            ("abnormal_sites", ab_sites),
            ("normal_sites", norm_sites),
            ("abnormal_genes", ab_genes),
            ("normal_genes", norm_genes),
        ]:
            sm = summarize(vals)
            for k, v in sm.items():
                f.write(f"{label}_{k}\t{v}\n")
        # normal heterogeneity: fraction above abnormal median
        ab_med = summarize(ab_sites).get("median", 0)
        if norm_sites:
            frac_hi = sum(1 for x in norm_sites if x > ab_med) / len(norm_sites)
            f.write(f"normal_frac_sites_gt_abnormal_median\t{frac_hi}\n")

    print(f"WROTE {sample_path}", flush=True)
    print(f"WROTE {sum_path}", flush=True)
    print(
        f"Sample load MW sites p={p_s:.4g} genes p={p_g:.4g} "
        f"ab_med_sites={summarize(ab_sites).get('median')} "
        f"norm_med_sites={summarize(norm_sites).get('median')}",
        flush=True,
    )
    print(f"Site patterns: {dict(site_pat)}", flush=True)
    print(f"Gene patterns: {dict(gene_pat)}", flush=True)

    # --- sensitivity / descriptive (raw columns) ---
    def revel_ok(v) -> bool | None:
        x = parse_af(v)
        if x is None:
            return None
        return x >= 0.5

    def af_strict(v) -> bool:
        x = parse_af(v)
        return x is None or x == 0.0

    def clnsig_blank(v: str) -> bool:
        return str(v).strip() in {"", ".", "NA", "na", "nan"}

    sens_rows = []
    for subset_name, pred in [
        ("all_union", lambda r: True),
        ("af_missing_or_zero", lambda r: af_strict(r["gnomad_af"])),
        ("func_lof_or_splicing", lambda r: r.get("func_class") in {"lof", "splicing"}),
        ("func_missense", lambda r: r.get("func_class") == "missense"),
        ("revel_ge_0.5", lambda r: revel_ok(r["revel"]) is True),
        ("clnsig_blank", lambda r: clnsig_blank(r["clnsig"])),
        ("clnsig_annotated", lambda r: not clnsig_blank(r["clnsig"])),
    ]:
        sub = [r for r in site_rows if pred(r)]
        pc = Counter(r["pattern"] for r in sub)
        sens_rows.append(
            {
                "subset": subset_name,
                "n_sites": len(sub),
                "abnormal_enriched": pc.get("abnormal_enriched", 0),
                "normal_enriched": pc.get("normal_enriched", 0),
                "shared": pc.get("shared", 0),
                "abnormal_only": pc.get("abnormal_only", 0),
                "normal_only": pc.get("normal_only", 0),
                "frac_clnsig_blank": (
                    sum(1 for r in sub if clnsig_blank(r["clnsig"])) / len(sub) if sub else ""
                ),
            }
        )
    sens_path = OUT / "dn_sensitivity_site_subsets.tsv"
    with open(sens_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(sens_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(sens_rows)
    print(f"WROTE {sens_path}", flush=True)

    # clnsig blank rate by pattern
    cln_path = OUT / "dn_clnsig_by_pattern.tsv"
    with open(cln_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["pattern", "n_sites", "n_clnsig_blank", "frac_blank"])
        by = defaultdict(list)
        for r in site_rows:
            by[r["pattern"]].append(r)
        for pat_name, rows in sorted(by.items()):
            nb = sum(1 for r in rows if clnsig_blank(r["clnsig"]))
            w.writerow([pat_name, len(rows), nb, nb / len(rows) if rows else 0])
    print(f"WROTE {cln_path}", flush=True)

    # talking points
    tp = OUT / "dn_ab_vs_normal_talking_points.txt"
    with open(tp, "w", encoding="utf-8") as f:
        f.write("DN-background SNV: abnormal vs normal (raw columns)\n")
        f.write("=" * 60 + "\n")
        if args.universe == "freq_gate":
            f.write(
                "DN background = loose dual filter (nonsyn/splice, drop Benign-only, "
                "gnomAD AF missing or <1%) plus cohort frequency gate: "
                "control carrier rate <5% and combined case (abnormal+normal) rate >5%. "
                "Not trio de novo.\n"
            )
        else:
            f.write(
                "DN background = loose dual filter (nonsyn/splice, drop Benign-only, "
                "gnomAD AF missing or <1%, own control carrier=0). Not trio de novo.\n"
            )
        f.write(f"Phenotype: {pheno_path}\n")
        f.write(
            f"Groups: abnormal={n_ab_tot} (47 label + CpG suspects), "
            f"normal={n_norm_tot}, control={n_ctrl_tot}\n"
        )
        f.write(
            "Valid FDR: full-union BH for two-sided difference and for each "
            "pre-specified one-sided test (ab higher; normal higher). "
            "m_site=all union sites; m_gene=all hit genes. No direction-based subsetting of m.\n"
        )
        f.write(
            f"FDR m: sites={len(site_rows)}, genes={len(gene_rows)}; "
            f"site FDR<0.05 two/ab/norm="
            f"{sum(1 for r in site_rows if r['fdr_two_sided'] < 0.05)}/"
            f"{sum(1 for r in site_rows if r['fdr_ab_greater'] < 0.05)}/"
            f"{sum(1 for r in site_rows if r['fdr_norm_greater'] < 0.05)}; "
            f"gene FDR<0.05 two/ab/norm="
            f"{sum(1 for r in gene_rows if r['fdr_two_sided'] < 0.05)}/"
            f"{sum(1 for r in gene_rows if r['fdr_ab_greater'] < 0.05)}/"
            f"{sum(1 for r in gene_rows if r['fdr_norm_greater'] < 0.05)}\n"
        )
        f.write(f"Sites with ab/norm carriers: {len(site_rows)}; genes: {len(gene_rows)}\n")
        f.write(f"Site patterns: {dict(site_pat)}\n")
        f.write(f"Gene patterns: {dict(gene_pat)}\n")
        f.write(
            f"Sample DN site load MW p={p_s:.4g}; "
            f"abnormal median={summarize(ab_sites).get('median')}, "
            f"normal median={summarize(norm_sites).get('median')}\n"
        )
        f.write("Main tests use: chrom/pos/ref/alt, gene, carrier counts, sample IDs.\n")
        f.write("Not used in discovery tests: clinvar_tiers, stratum, score_votes, score_labels.\n")
        f.write("clnsig / revel / cadd / sift used only in sensitivity or descriptive tables.\n")
    print(f"WROTE {tp}", flush=True)


if __name__ == "__main__":
    main()
