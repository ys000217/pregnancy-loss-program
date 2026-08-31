#!/usr/bin/env python3
import csv
import re
from pathlib import Path

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

path = Path("/mnt/d/ONT/clinical_649.GRCh38.annotsv.tsv")
if not path.exists():
    path = Path("D:/ONT/clinical_649.GRCh38.annotsv.tsv")

af_cols = ["B_gain_AFmax", "B_loss_AFmax", "B_ins_AFmax", "B_inv_AFmax"]
rows = []
with open(path, newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        acmg = parse_acmg(row.get("ACMG_class"))
        if acmg is None or acmg < 4:
            continue
        key = f"{row['SV_chrom']}:{row['SV_start']}:{row['SV_end']}:{row['SV_type']}"
        vals = []
        for c in af_cols:
            try:
                v = float(row.get(c, "") or "nan")
            except ValueError:
                v = float("nan")
            if v == v:
                vals.append(v)
        max_af = max(vals) if vals else None
        rows.append((key, max_af, row.get("Annotation_mode"), acmg))

print("plp rows", len(rows))
print("unique sv", len(set(x[0] for x in rows)))
best = {}
for key, max_af, mode, acmg in rows:
    prev = best.get(key)
    if prev is None or (mode == "full" and prev[2] != "full"):
        best[key] = (max_af, mode, acmg)

common = [k for k, v in best.items() if v[0] is not None and v[0] >= 0.01]
rare = [k for k, v in best.items() if v[0] is not None and v[0] < 0.01]
novel = [k for k, v in best.items() if v[0] is None]
print("dedup unique", len(best), "common>=1%", len(common), "rare<1%", len(rare), "no_af", len(novel))
afs = [v[0] for v in best.values() if v[0] is not None]
if afs:
    print("AF stats min/median/max:", min(afs), sorted(afs)[len(afs)//2], max(afs))
print("sample AF values:", sorted(set(afs))[:20])