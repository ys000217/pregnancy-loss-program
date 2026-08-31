#!/usr/bin/env python3
"""ONT-centric pairing: unique sample IDs (any string: 0002C, CS…, 107, A0281, …).

Two BAM trees, then match WGS FASTQ libraries by the same unique id.

Batch A: --ont-root  (Dorado)  {id}/{id}.merged.bam
Batch B: --ont-extra (ONT_Rawdata)  **/ {id}_mods_merged.bam

Dedup by sample id (case-insensitive). Prefer Dorado {id}.merged.bam over _mods_merged.bam.

Match WGS:
  1) ngs_lib == sample id (case-insensitive)
  2) ReadMe.txt SAMPLE <-> FDSW aliases

Writes ont_bam_index.tsv, ont_wgs_coverage.tsv, manifest.from_ont.tsv (includes ont_bam).
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

def sample_id_from_bam(path: Path) -> str:
    n = path.name
    low = n.lower()
    if low.endswith("_mods_merged.bam"):
        return n[: -len("_mods_merged.bam")]
    if low.endswith(".merged.bam"):
        return n[: -len(".merged.bam")]
    return path.parent.name


def rank_bam(path: Path) -> int:
    n = path.name.lower()
    if n.endswith("_mods_merged.bam"):
        return 2
    if n.endswith(".merged.bam"):
        return 1
    return 3


def discover_bams(roots: list[Path]) -> list[Path]:
    found: list[Path] = []
    for root in roots:
        if not root.is_dir():
            print(f"WARN skip missing root: {root}")
            continue
        for p in root.rglob("*merged.bam"):
            if p.is_file() and p.name.lower().endswith("merged.bam"):
                found.append(p)
    return found


def dedupe(bams: list[Path]) -> dict[str, Path]:
    """lower(id) -> best bam path. Also map to canonical id from preferred file."""
    best: dict[str, tuple[int, Path, str]] = {}
    for bam in bams:
        sid = sample_id_from_bam(bam)
        key = sid.lower()
        r = rank_bam(bam)
        prev = best.get(key)
        if prev is None or r < prev[0]:
            best[key] = (r, bam, sid)
    return {k: v[1] for k, v in best.items()}


def canonical_ids(bams: dict[str, Path]) -> dict[str, str]:
    out = {}
    for k, bam in bams.items():
        out[k] = sample_id_from_bam(bam)
        # Prefer parent folder if it matches (Dorado 0002C/0002C.merged.bam)
        parent = bam.parent.name
        if parent.lower() == k:
            out[k] = parent
    return out


def load_wgs(manifest: Path) -> dict[str, dict]:
    by_lib: dict[str, dict] = {}
    with manifest.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            lib = (row.get("ngs_lib") or row.get("ont_id") or "").strip()
            if lib:
                by_lib[lib.lower()] = row
    return by_lib


def load_readme_aliases(ngs_root: Path) -> dict[str, str]:
    alias: dict[str, str] = {}
    for readme in ngs_root.glob("**/03.release/ReadMe.txt"):
        text = readme.read_text(encoding="utf-8", errors="replace")
        for line in text.splitlines():
            line = line.strip()
            if not line or line.startswith("=") or "项目" in line or line.startswith("raw_data"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            a, b = parts[0], parts[1]
            alias.setdefault(a.lower(), a)
            alias.setdefault(b.lower(), a)
    return alias


def find_wgs(key: str, wgs: dict, aliases: dict) -> tuple[dict | None, str]:
    row = wgs.get(key)
    if row:
        return row, "name"
    if key in aliases:
        alt = aliases[key]
        row = wgs.get(alt.lower())
        if row:
            return row, "readme"
    for a, b in aliases.items():
        if a == key or b.lower() == key:
            row = wgs.get(a) or wgs.get(b.lower())
            if row:
                return row, "readme"
    return None, ""


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--ont-root", required=True, type=Path, help="Dorado processed BAMs")
    p.add_argument(
        "--ont-extra",
        action="append",
        default=[],
        type=Path,
        help="Extra tree to rglob *merged.bam (e.g. ONT_Rawdata). Repeatable.",
    )
    p.add_argument("--ngs-root", required=True, type=Path)
    p.add_argument("--wgs-manifest", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    args = p.parse_args()

    roots = [args.ont_root] + list(args.ont_extra)
    raw = discover_bams(roots)
    by_key = dedupe(raw)
    ids = canonical_ids(by_key)

    wgs = load_wgs(args.wgs_manifest)
    aliases = load_readme_aliases(args.ngs_root)
    args.outdir.mkdir(parents=True, exist_ok=True)

    idx_path = args.outdir / "ont_bam_index.tsv"
    cov_path = args.outdir / "ont_wgs_coverage.tsv"
    man_path = args.outdir / "manifest.from_ont.tsv"
    fieldnames = [
        "ont_id", "ont_bam", "wgs_r1", "wgs_r2", "sex", "n_lanes", "ngs_lib", "match",
    ]

    n_ok = 0
    missing = []
    n_dorado = sum(1 for b in by_key.values() if rank_bam(b) == 1)
    n_mods = sum(1 for b in by_key.values() if rank_bam(b) == 2)

    with idx_path.open("w", newline="") as ifh, cov_path.open("w", newline="") as cfh, man_path.open("w", newline="") as mfh:
        iw = csv.writer(ifh, delimiter="\t")
        cw = csv.writer(cfh, delimiter="\t")
        mw = csv.DictWriter(mfh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        iw.writerow(["ont_id", "ont_bam", "rank", "n_raw_hits"])
        cw.writerow(["ont_id", "has_wgs", "match", "ngs_lib", "ont_bam"])
        mw.writeheader()

        raw_counts: dict[str, int] = {}
        for bam in raw:
            raw_counts[sample_id_from_bam(bam).lower()] = raw_counts.get(sample_id_from_bam(bam).lower(), 0) + 1

        for key in sorted(by_key, key=lambda k: ids[k]):
            bam = by_key[key]
            oid = ids[key]
            iw.writerow([oid, str(bam), rank_bam(bam), raw_counts.get(key, 1)])
            row, how = find_wgs(key, wgs, aliases)
            if row:
                n_ok += 1
                out = {
                    "ont_id": oid,
                    "ont_bam": str(bam),
                    "wgs_r1": row["wgs_r1"],
                    "wgs_r2": row["wgs_r2"],
                    "sex": row.get("sex") or "",
                    "n_lanes": row.get("n_lanes") or "",
                    "ngs_lib": row.get("ngs_lib") or oid,
                    "match": how,
                }
                mw.writerow(out)
                cw.writerow([oid, "yes", how, out["ngs_lib"], str(bam)])
            else:
                missing.append(oid)
                cw.writerow([oid, "no", "", "", str(bam)])

    print(f"raw *merged.bam files: {len(raw)}")
    print(f"unique ONT samples after dedup: {len(by_key)}")
    print(f"  preferred Dorado .merged.bam: {n_dorado}")
    print(f"  preferred _mods_merged.bam only: {n_mods}")
    print(f"with WGS: {n_ok}")
    print(f"without WGS: {len(missing)}")
    print(f"index: {idx_path}")
    print(f"coverage: {cov_path}")
    print(f"pipeline manifest: {man_path}")
    if missing:
        print("missing examples:", ", ".join(missing[:30]))
        if len(missing) > 30:
            print(f"... and {len(missing) - 30} more")


if __name__ == "__main__":
    main()
