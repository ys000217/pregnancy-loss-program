#!/usr/bin/env python3
"""Quick check rare P/LP definitions."""
import csv
import re
from pathlib import Path
from collections import Counter

ROOT = Path("/mnt/d/ONT/figure2") if Path("/mnt/d/ONT/figure2").exists() else Path("D:/ONT/figure2")
ONT = ROOT.parent

def parse_acmg(raw):
    if not raw:
        return None
    raw = str(raw).strip()
    m = re.match(r"^full=(\d+)$", raw)
    if m:
        return int(m.group(1))
    try:
        return int(float(raw))
    except ValueError:
        return None

def max_pop_af(row):
    vals = []
    for c in ["B_gain_AFmax", "B_loss_AFmax", "B_ins_AFmax", "B_inv_AFmax"]:
        try:
            v = float(row.get(c, "") or "nan")
        except ValueError:
            v = float("nan")
        if v == v:
            vals.append(v)
    return max(vals) if vals else None

# load plp metadata
plp = {}
with open(ONT / "clinical_649.GRCh38.annotsv.tsv", newline="", encoding="utf-8", errors="replace") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        acmg = parse_acmg(row.get("ACMG_class"))
        if acmg is None or acmg < 4:
            continue
        key = f"{row['SV_chrom']}:{row['SV_start']}:{row['SV_end']}:{row['SV_type']}"
        rec = (max_pop_af(row), row.get("Annotation_mode"))
        prev = plp.get(key)
        if prev is None or (rec[1] == "full" and prev[1] != "full"):
            plp[key] = rec

# cohort carriers from enrichment file if exists
enrich = ROOT / "burden_analysis/tables/sv_plp_enrichment_case_vs_control.tsv"
cohort_n = 648
if enrich.exists():
    cohort_rate = {}
    with open(enrich, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = row["coord_key"]
            total = int(row["case_carrier"]) + int(row["control_carrier"])
            cohort_rate[key.replace("|", ":").split("|")[0]] = total / cohort_n
    # coord_key format chrom:start:end:type
    keys_rate = {}
    with open(enrich, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            ck = row["coord_key"]
            total = int(row["case_carrier"]) + int(row["control_carrier"])
            keys_rate[ck] = total / cohort_n

def pop_rare(af):
    return af is None or af < 0.01

def cohort_rare(rate):
    return rate < 0.05

# recompute cohort from sample burden - skip, use enrichment
print("unique plp", len(plp))
pop_only = sum(1 for af, _ in plp.values() if pop_rare(af))
print("pop rare (no AF or AF<1%):", pop_only)

# read enrichment for cohort rates
rates = {}
with open(enrich, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        rates[row["coord_key"]] = (
            int(row["case_carrier"]) + int(row["control_carrier"])
        ) / cohort_n

both = 0
pop = 0
cohort = 0
for key in plp:
    ck = key
    af = plp[key][0]
    rate = rates.get(ck, 1.0)
    pr = pop_rare(af)
    cr = cohort_rare(rate)
    if pr:
        pop += 1
    if cr:
        cohort += 1
    if pr and cr:
        both += 1
print("cohort rare <5%:", cohort)
print("both pop+cohort rare:", both)