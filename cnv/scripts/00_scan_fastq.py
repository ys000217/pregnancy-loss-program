#!/usr/bin/env python3
"""Scan Novogene NGS_Rawdata: 03.release/raw_data/{lib}_L{lane}_{1,2}.fq.gz

One manifest row per library. Multi-lane R1/R2 paths are comma-separated.
Column ont_id starts as the library ID (FDSW...); replace with ONT folder names (0002C).
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

LANE_R1 = re.compile(r"^(.+)_L(\d+)_1\.f(ast)?q\.gz$", re.I)


def find_grouped(root: Path) -> list[tuple[str, str, str, int]]:
    """Return (lib_id, r1_csv, r2_csv, n_lanes)."""
    libs: dict[str, list[tuple[int, Path, Path]]] = defaultdict(list)
    r1_files = sorted(root.glob("**/03.release/raw_data/*_1.fq.gz"))
    r1_files += sorted(root.glob("**/03.release/raw_data/*_1.fastq.gz"))
    # Some batches put each sample in a subfolder: raw_data/0002C/0002C_L1_1.fq.gz
    r1_files += sorted(root.glob("**/03.release/raw_data/*/*_1.fq.gz"))
    r1_files += sorted(root.glob("**/03.release/raw_data/*/*_1.fastq.gz"))
    seen: set[Path] = set()
    for r1 in r1_files:
        if r1 in seen:
            continue
        seen.add(r1)
        m = LANE_R1.match(r1.name)
        if m:
            lib, lane = m.group(1), int(m.group(2))
            r2 = r1.with_name(r1.name.replace("_1.fq.gz", "_2.fq.gz").replace("_1.fastq.gz", "_2.fastq.gz"))
        else:
            lib = r1.name.replace("_1.fq.gz", "").replace("_1.fastq.gz", "")
            lane = 0
            r2 = r1.with_name(r1.name.replace("_1.fq.gz", "_2.fq.gz").replace("_1.fastq.gz", "_2.fastq.gz"))
        if not r2.exists():
            print(f"WARN unpaired R1: {r1}")
            continue
        libs[lib].append((lane, r1, r2))

    rows = []
    for lib in sorted(libs):
        lanes = sorted(libs[lib], key=lambda x: x[0])
        r1s = ",".join(str(p) for _, p, _ in lanes)
        r2s = ",".join(str(p) for _, _, p in lanes)
        rows.append((lib, r1s, r2s, len(lanes)))
    return rows


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--ngs-root", required=True, type=Path)
    p.add_argument("-o", "--out", required=True, type=Path)
    args = p.parse_args()

    rows = find_grouped(args.ngs_root)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["ont_id", "wgs_r1", "wgs_r2", "sex", "n_lanes", "ngs_lib"])
        for lib, r1s, r2s, n in rows:
            w.writerow([lib, r1s, r2s, "", n, lib])
    print(f"wrote {len(rows)} libraries to {args.out}")
    if not rows:
        print("No **/03.release/raw_data/*_1.fq.gz found")
        return
    print("Replace column ont_id (currently FDSW library ID) with ONT folder names, e.g. 0002C.")
    print("Keep ngs_lib as-is. Multi-lane FASTQs are comma-separated; 01_wgs_align.sh concatenates them.")
    print("Mapping hint: head a batch ReadMe.txt, or NGS_release_structure_check.tsv")


if __name__ == "__main__":
    main()
