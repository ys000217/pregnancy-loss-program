#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compute unbiased SNV / SV / P-LP SV burden from VCF + annotation files."""

from __future__ import annotations

import csv
import gzip
import os
import re
from collections import Counter, defaultdict
from math import comb
from pathlib import Path

# Input data live under figure2 / ONT root (large files not in git).
# Outputs default to this module: analyses/burden/{tables,figures}.
MODULE = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get("BURDEN_ROOT", "D:/ONT/figure2"))
ONT = ROOT.parent

PHENO_PATH = ROOT / "sample_phenotype_648.tsv"
SNV_VCF = ROOT / "WGS_ONT_Intersection_648samples.vcf.gz"
SNV_ANNO = ROOT / "S956.snp.annovar.hg38_multianno.txt.gz"
SV_VCF = ONT / "clinical_649.GRCh38.correct.vcf"
SV_ANNO = ONT / "clinical_649.GRCh38.annotsv.tsv"

OUT_DIR = Path(os.environ.get("BURDEN_OUT", str(MODULE / "tables")))
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Rare P/LP SV: max gnomAD-SV AF from AnnotSV B_*_AFmax columns
POP_AF_RARE_THRESHOLD = 0.01
# Sensitivity: also require low frequency in this cohort
COHORT_AF_RARE_THRESHOLD = 0.05
POP_AF_COLS = ["B_gain_AFmax", "B_loss_AFmax", "B_ins_AFmax", "B_inv_AFmax"]


def norm_chr(chrom: str) -> str:
    return chrom.lower().replace("chr", "").strip()


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def load_phenotype(path: Path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def has_alt_allele(gt_field: str) -> bool:
    gt = gt_field.split(":")[0]
    return any(ch in gt for ch in "123456789")


def parse_max_pop_af(row: dict) -> float | None:
    vals = []
    for col in POP_AF_COLS:
        try:
            v = float(row.get(col, "") or "nan")
        except ValueError:
            v = float("nan")
        if v == v:
            vals.append(v)
    return max(vals) if vals else None


def is_pop_rare(max_pop_af: float | None) -> bool:
    return max_pop_af is None or max_pop_af < POP_AF_RARE_THRESHOLD


def is_cohort_rare(cohort_rate: float) -> bool:
    return cohort_rate < COHORT_AF_RARE_THRESHOLD


def parse_acmg_class(raw: str) -> int | None:
    if raw is None:
        return None
    raw = str(raw).strip()
    if not raw or raw.lower() in {"na", "nan", "."}:
        return None
    m = re.match(r"^full=(\d+)$", raw)
    if m:
        return int(m.group(1))
    try:
        return int(float(raw))
    except ValueError:
        return None


def fisher_exact_2x2(a: int, b: int, c: int, d: int) -> float:
    row1, row2, col1, n = a + b, c + d, a + c, a + b + c + d
    if n == 0 or min(a, b, c, d) < 0:
        return 1.0
    p_obs = 0.0
    for x in range(a, min(row1, col1) + 1):
        y, z, w = col1 - x, row1 - x, row2 - (col1 - x)
        if min(y, z, w) < 0:
            continue
        p_obs += comb(col1, x) * comb(n - col1, row1 - x) / comb(n, row1)
    return min(1.0, p_obs)


def bh_fdr(p_values: list[float]) -> list[float]:
    m = len(p_values)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: p_values[i])
    q = [0.0] * m
    prev = 1.0
    for rank in range(1, m + 1):
        idx = order[m - rank]
        val = p_values[idx] * m / (m - rank + 1)
        prev = min(prev, val)
        q[idx] = min(prev, 1.0)
    return q


def parse_info_field(info: str, key: str) -> str:
    for field in info.split(";"):
        if field.startswith(key + "="):
            return field.split("=", 1)[1]
    return ""


def build_functional_snv_keys(anno_path: Path) -> set[tuple[str, str, str, str]]:
    func_idx, exonic_idx = 10, 13
    keys: set[tuple[str, str, str, str]] = set()
    with open_text(anno_path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= exonic_idx:
                continue
            if parts[func_idx] != "exonic" or parts[exonic_idx] == "synonymous SNV":
                continue
            keys.add((norm_chr(parts[1]), parts[2], parts[4], parts[5]))
    return keys


def compute_snv_burden(samples: list[str], functional_keys: set, vcf_path: Path):
    counts = Counter({s: 0 for s in samples})
    matched = 0
    with open_text(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[6] != "PASS":
                continue
            key = (norm_chr(parts[0]), parts[1], parts[3], parts[4])
            if key not in functional_keys:
                continue
            matched += 1
            for sample, gt in zip(samples, parts[9:]):
                if has_alt_allele(gt):
                    counts[sample] += 1
    return counts, matched


def compute_sv_burden(samples: list[str], vcf_path: Path):
    total = Counter({s: 0 for s in samples})
    by_type: dict[str, Counter] = defaultdict(lambda: Counter({s: 0 for s in samples}))
    carriers_by_coord: dict[str, dict[str, int]] = {}
    carriers_by_id: dict[str, dict[str, int]] = {}
    n_variants = 0

    with open_text(vcf_path) as f:
        vcf_samples = None
        for line in f:
            if line.startswith("#CHROM"):
                vcf_samples = line.strip().split("\t")[9:]
                break
        if vcf_samples is None:
            raise RuntimeError("SV VCF header missing")

        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[6] != "PASS":
                continue
            svtype = parse_info_field(parts[7], "SVTYPE") or "UNK"
            sv_end = parse_info_field(parts[7], "END") or parts[1]
            sv_id = parts[2]
            coord_key = f"{parts[0]}:{parts[1]}:{sv_end}:{svtype}"
            gt_map = dict(zip(vcf_samples, parts[9:]))
            carried = {}
            for sample in samples:
                gt = gt_map.get(sample, "./.")
                g = gt.split(":")[0]
                if g in {".", "./.", ".|."}:
                    continue
                total[sample] += 1
                by_type[svtype][sample] += 1
                carried[sample] = 1
            if carried:
                n_variants += 1
                carriers_by_coord[coord_key] = {s: carried.get(s, 0) for s in samples}
                carriers_by_id[sv_id] = carriers_by_coord[coord_key]

    return total, by_type, n_variants, carriers_by_coord, carriers_by_id


def load_plp_sv_records(anno_path: Path):
    candidates: dict[str, dict] = {}
    with open(anno_path, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acmg = parse_acmg_class(row.get("ACMG_class"))
            if acmg is None or acmg < 4:
                continue
            coord_key = f"{row.get('SV_chrom')}:{row.get('SV_start')}:{row.get('SV_end')}:{row.get('SV_type')}"
            mode = row.get("Annotation_mode", "")
            score = float(row.get("AnnotSV_ranking_score") or 0)
            prev = candidates.get(coord_key)
            if prev is None or (mode == "full" and prev.get("mode") != "full") or (
                mode == prev.get("mode") and score > prev.get("score", 0)
            ):
                max_af = parse_max_pop_af(row)
                candidates[coord_key] = {
                    "coord_key": coord_key,
                    "chrom": row.get("SV_chrom"),
                    "start": row.get("SV_start"),
                    "end": row.get("SV_end"),
                    "svtype": row.get("SV_type"),
                    "sv_id": row.get("ID"),
                    "gene": row.get("Gene_name"),
                    "acmg_class": acmg,
                    "acmg_label": row.get("ACMG_class"),
                    "ranking_score": score,
                    "mode": mode,
                    "max_pop_af": max_af if max_af is not None else "",
                    "pop_rare": int(is_pop_rare(max_af)),
                }
    return list(candidates.values())


def match_carriers(rec: dict, by_coord: dict, by_id: dict, samples: list[str]):
    if rec["coord_key"] in by_coord:
        return by_coord[rec["coord_key"]]
    if rec.get("sv_id") and rec["sv_id"] in by_id:
        return by_id[rec["sv_id"]]
    return {s: 0 for s in samples}


def annotate_plp_cohort_stats(plp_records, samples, by_coord, by_id):
    n = len(samples)
    for rec in plp_records:
        carriers = match_carriers(rec, by_coord, by_id, samples)
        cohort_carrier = sum(carriers.values())
        cohort_rate = cohort_carrier / n if n else 0.0
        rec["cohort_carrier"] = cohort_carrier
        rec["cohort_rate"] = cohort_rate
        rec["cohort_rare"] = int(is_cohort_rare(cohort_rate))
        rec["strict_rare"] = int(rec.get("pop_rare") == 1 and is_cohort_rare(cohort_rate))
    return plp_records


def compute_plp_burden_and_enrichment(
    samples,
    pheno_map,
    plp_records,
    by_coord,
    by_id,
    record_filter=None,
):
    case_samples = [s for s in samples if pheno_map.get(s) in {"abnormal", "normal"}]
    control_samples = [s for s in samples if pheno_map.get(s) == "control"]
    plp_counts = Counter({s: 0 for s in samples})
    enrichment_rows = []

    for rec in plp_records:
        if record_filter is not None and not record_filter(rec):
            continue
        carriers = match_carriers(rec, by_coord, by_id, samples)
        for s, v in carriers.items():
            plp_counts[s] += v
        a = sum(carriers[s] for s in case_samples)
        b = len(case_samples) - a
        c = sum(carriers[s] for s in control_samples)
        d = len(control_samples) - c
        p = fisher_exact_2x2(a, b, c, d)
        enrichment_rows.append(
            {
                **rec,
                "case_carrier": a,
                "case_total": len(case_samples),
                "control_carrier": c,
                "control_total": len(control_samples),
                "case_rate": a / len(case_samples) if case_samples else 0,
                "control_rate": c / len(control_samples) if control_samples else 0,
                "p_value": p,
            }
        )

    qvals = bh_fdr([r["p_value"] for r in enrichment_rows]) if enrichment_rows else []
    for r, q in zip(enrichment_rows, qvals):
        r["fdr"] = q
    enrichment_rows.sort(key=lambda x: x["p_value"])
    return plp_counts, enrichment_rows


def write_sample_burden(
    path,
    pheno_rows,
    snv_counts,
    sv_total,
    sv_by_type,
    plp_counts,
    plp_rare_counts=None,
    plp_strict_counts=None,
):
    sv_types = sorted(sv_by_type.keys())
    fields = ["Sample_ID", "Condition", "Gestational_Week", "SNV_nonsyn_count", "SV_total"]
    fields += [f"SV_{t}" for t in sv_types]
    fields += [
        "SV_plp_count",
        "SV_plp_binary",
        "SV_plp_rare_count",
        "SV_plp_rare_binary",
        "SV_plp_strict_rare_count",
        "SV_plp_strict_rare_binary",
    ]
    plp_rare_counts = plp_rare_counts or Counter()
    plp_strict_counts = plp_strict_counts or Counter()
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in pheno_rows:
            sid = row["VCF_Sample_ID"]
            rare_n = plp_rare_counts.get(sid, 0)
            strict_n = plp_strict_counts.get(sid, 0)
            out = {
                "Sample_ID": sid,
                "Condition": row["Condition"],
                "Gestational_Week": row.get("Gestational_Week", ""),
                "SNV_nonsyn_count": snv_counts.get(sid, 0),
                "SV_total": sv_total.get(sid, 0),
                "SV_plp_count": plp_counts.get(sid, 0),
                "SV_plp_binary": int(plp_counts.get(sid, 0) > 0),
                "SV_plp_rare_count": rare_n,
                "SV_plp_rare_binary": int(rare_n > 0),
                "SV_plp_strict_rare_count": strict_n,
                "SV_plp_strict_rare_binary": int(strict_n > 0),
            }
            for t in sv_types:
                out[f"SV_{t}"] = sv_by_type[t].get(sid, 0)
            w.writerow(out)


def load_snv_counts_from_burden(path: Path, samples: list[str]) -> Counter:
    counts = Counter({s: 0 for s in samples})
    if not path.exists():
        return counts
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["Sample_ID"] in counts:
                counts[row["Sample_ID"]] = int(float(row["SNV_nonsyn_count"]))
    return counts


def main():
    print("[1/4] Loading phenotype...")
    pheno_rows = load_phenotype(PHENO_PATH)
    samples = [r["VCF_Sample_ID"] for r in pheno_rows]
    pheno_map = {r["VCF_Sample_ID"]: r["Condition"] for r in pheno_rows}
    print(f"  samples: {len(samples)}")

    out_sample = OUT_DIR / "sample_burden.tsv"
    skip_snv = os.environ.get("SKIP_SNV") == "1"
    func_keys = None
    if skip_snv and out_sample.exists():
        print("[2/4] SNV burden... SKIPPED (loading existing counts)")
        snv_counts = load_snv_counts_from_burden(out_sample, samples)
        snv_matched = -1
    else:
        print("[2/4] SNV burden...")
        func_keys = build_functional_snv_keys(SNV_ANNO)
        print(f"  functional nonsyn keys: {len(func_keys):,}")
        snv_counts, snv_matched = compute_snv_burden(samples, func_keys, SNV_VCF)
        print(f"  matched VCF variants: {snv_matched:,}")

    print("[3/4] SV burden...")
    sv_total, sv_by_type, n_sv, by_coord, by_id = compute_sv_burden(samples, SV_VCF)
    print(f"  PASS SV records: {n_sv:,}")

    print("[4/4] P/LP SV (all + rare)...")
    plp_records = load_plp_sv_records(SV_ANNO)
    plp_records = annotate_plp_cohort_stats(plp_records, samples, by_coord, by_id)
    print(f"  unique P/LP SVs: {len(plp_records)}")
    pop_rare_n = sum(1 for r in plp_records if r.get("pop_rare") == 1)
    strict_catalog = [r for r in plp_records if r.get("strict_rare") == 1]
    print(f"  population-rare P/LP (no gnomAD AF or AF<{POP_AF_RARE_THRESHOLD}): {pop_rare_n}")
    print(f"  strict rare P/LP (pop rare + cohort<{COHORT_AF_RARE_THRESHOLD:.0%}): {len(strict_catalog)}")

    plp_counts, enrichment_rows = compute_plp_burden_and_enrichment(
        samples, pheno_map, plp_records, by_coord, by_id, record_filter=None
    )
    plp_rare_counts, enrichment_rare = compute_plp_burden_and_enrichment(
        samples,
        pheno_map,
        plp_records,
        by_coord,
        by_id,
        record_filter=lambda r: r.get("pop_rare") == 1,
    )
    plp_strict_counts, enrichment_strict = compute_plp_burden_and_enrichment(
        samples,
        pheno_map,
        plp_records,
        by_coord,
        by_id,
        record_filter=lambda r: r.get("strict_rare") == 1,
    )

    write_sample_burden(
        out_sample,
        pheno_rows,
        snv_counts,
        sv_total,
        sv_by_type,
        plp_counts,
        plp_rare_counts,
        plp_strict_counts,
    )
    print(f"  wrote {out_sample}")

    def write_enrichment(rows, path):
        if not rows:
            return
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    write_enrichment(enrichment_rows, OUT_DIR / "sv_plp_enrichment_case_vs_control.tsv")
    write_enrichment(enrichment_rare, OUT_DIR / "sv_plp_rare_enrichment_case_vs_control.tsv")
    write_enrichment(enrichment_strict, OUT_DIR / "sv_plp_strict_rare_enrichment_case_vs_control.tsv")
    write_enrichment(plp_records, OUT_DIR / "sv_plp_catalog_with_rarity.tsv")
    print("  wrote P/LP enrichment tables (all / rare / strict-rare)")

    with open(OUT_DIR / "compute_burden_summary.txt", "w", encoding="utf-8") as f:
        if func_keys is not None:
            f.write(f"functional_snv_keys\t{len(func_keys)}\n")
        if snv_matched >= 0:
            f.write(f"snv_variants_in_vcf\t{snv_matched}\n")
        f.write(f"sv_variants_pass\t{n_sv}\n")
        f.write(f"unique_plp_sv\t{len(plp_records)}\n")
        f.write(f"unique_plp_pop_rare\t{pop_rare_n}\n")
        f.write(f"unique_plp_strict_rare\t{len(strict_catalog)}\n")
        f.write(f"pop_af_rare_threshold\t{POP_AF_RARE_THRESHOLD}\n")
        f.write(f"cohort_af_rare_threshold\t{COHORT_AF_RARE_THRESHOLD}\n")
    print("Done.")


if __name__ == "__main__":
    main()
