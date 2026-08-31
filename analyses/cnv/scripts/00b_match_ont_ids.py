#!/usr/bin/env python3
"""Set ont_id to the ONT BAM folder name when it matches ngs_lib (case-insensitive).

Keeps unmatched WGS libraries with ont_id still = ngs_lib.
Writes manifest.paired.tsv: only rows with an existing {ONT_BAM_ROOT}/{id}/{id}.merged.bam
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def ont_folders(root: Path) -> dict[str, str]:
    """lower -> actual folder name"""
    out: dict[str, str] = {}
    if not root.is_dir():
        raise SystemExit(f"ONT_BAM_ROOT not a directory: {root}")
    for p in root.iterdir():
        if p.is_dir():
            out[p.name.lower()] = p.name
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--ont-root", required=True, type=Path)
    p.add_argument("-o", "--out", required=True, type=Path)
    p.add_argument("--paired-only", action="store_true", help="keep only libs with ONT BAM")
    args = p.parse_args()

    folders = ont_folders(args.ont_root)
    rows = []
    with args.manifest.open(newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(r.fieldnames or [])
        for row in r:
            lib = (row.get("ngs_lib") or row.get("ont_id") or "").strip()
            key = lib.lower()
            if key in folders:
                row["ont_id"] = folders[key]
                row["ont_match"] = "folder"
            else:
                row["ont_id"] = lib
                row["ont_match"] = "none"
            bam = args.ont_root / row["ont_id"] / f"{row['ont_id']}.merged.bam"
            row["has_ont_bam"] = "yes" if bam.is_file() else "no"
            rows.append(row)

    if "ont_match" not in fieldnames:
        fieldnames.extend(["ont_match", "has_ont_bam"])

    paired = [x for x in rows if x["has_ont_bam"] == "yes"]
    write_rows = paired if args.paired_only else rows

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(write_rows)

    n_yes = sum(1 for x in rows if x["has_ont_bam"] == "yes")
    n_folder = sum(1 for x in rows if x["ont_match"] == "folder")
    print(f"WGS libraries: {len(rows)}")
    print(f"name matches an ONT folder: {n_folder}")
    print(f"has {{{{id}}}}.merged.bam: {n_yes}")
    print(f"wrote {len(write_rows)} rows -> {args.out}")
    if n_yes == 0:
        print("No CS/FDSW id equals an ONT folder. Need a CS<->0002C mapping table.")


if __name__ == "__main__":
    main()
