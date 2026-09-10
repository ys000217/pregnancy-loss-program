#!/usr/bin/env python3
"""Export Arnadottir 2025 Nature pathogenic SSV tables from supplementary xlsx."""
from __future__ import annotations

import csv
from pathlib import Path

from openpyxl import load_workbook

OUT = Path(r"D:/ONT/analyses/burden/tables/external_arnadottir2025")
XLSX = OUT / "41586_2025_9031_MOESM4_ESM.xlsx"


def dump_sheet(wb, name: str, dest: Path, skip_title: bool = True) -> list[tuple]:
    ws = wb[name]
    rows = list(ws.iter_rows(values_only=True))
    if skip_title:
        rows = rows[1:]
    header = [str(x) if x is not None else "" for x in rows[0]]
    body = rows[1:]
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for r in body:
            w.writerow(["" if x is None else x for x in r[: len(header)]])
    return [tuple(header)] + [tuple(r[: len(header)]) for r in body]


def main() -> None:
    wb = load_workbook(XLSX, data_only=True)
    dump_sheet(wb, "Supp.Table.2", OUT / "supp_table2_cnv.tsv")
    dump_sheet(wb, "Supp.Table.4", OUT / "supp_table4_pathogenic_ssv.tsv")
    dump_sheet(wb, "Supp.Table.8", OUT / "supp_table8_lethality_score.tsv")
    dump_sheet(wb, "Supp.Table.9", OUT / "supp_table9_euploid_likely_causal.tsv")

    # unique pathogenic SSVs (one row per gene+protein effect+genotype)
    with (OUT / "supp_table4_pathogenic_ssv.tsv").open(encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = list(r)
    seen = set()
    uniq = []
    for row in rows:
        gene = (row.get("Gene") or "").strip()
        if not gene or gene.startswith("SSV="):
            continue
        key = (
            gene,
            row.get("Effect on protein (MANE transcript)", ""),
            row.get("Genotype", ""),
        )
        if key in seen:
            continue
        seen.add(key)
        uniq.append(row)
    dest = OUT / "pathogenic_ssv_unique.tsv"
    with dest.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(uniq[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(uniq)
    print("unique_pathogenic_ssv", len(uniq))
    for u in uniq:
        print(
            f"  {u.get('Gene')}\t{u.get('Genotype')}\t{u.get('Variant type')}\t"
            f"{u.get('Effect on protein (MANE transcript)')}\t{u.get('Fetal karyotype')}"
        )
    print("WROTE", dest)


if __name__ == "__main__":
    main()
