#!/usr/bin/env python3
"""Merge ONT + WGS DEL/DUP callsets for one 10x paired sample.

Confidence tiers (cohort-facing):
  LARGE_HIGH  both platforms, same type, RO>=0.5, size>=100kb, outside hard mask,
              and size<=MAX unless both-depth (+ optional SV) passes the size cap
  SHARED_SV   both platforms have SV breakpoints, size < 100kb (small SV burden)
  MEDIUM      cross-platform same-type overlap that is neither LARGE_HIGH nor SHARED_SV
  ONT_SV      ONT Sniffles only
  ONT_DEPTH   ONT depth only, >=100kb, unmasked, within size cap
  WGS_DEPTH   WGS depth only, >=100kb, unmasked, within size cap
  MASKED      would be depth/LARGE candidate but hits hard mask or sex-Y rule
  DROP        not written (WGS-only small SV, tiny singleton depth, oversize uncapped)

Legacy name HIGH is no longer emitted. AnnotSV / cohort tables use LARGE_HIGH
({sample}.cnv.high.bed).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path

RO_DEFAULT = 0.5
MIN_SV = 50
MIN_DEPTH = 100_000
MAX_EVENT_DEFAULT = 10_000_000
MASK_FRAC = 0.50

CHRY = "NC_000024.10"
FEMALE_SEX = {"f", "female", "xx", "2", "woman"}


@dataclass
class Call:
    chrom: str
    start: int  # 0-based
    end: int
    svtype: str  # DEL or DUP
    source: str
    size: int = 0
    extra: str = ""

    def __post_init__(self) -> None:
        if self.size <= 0:
            self.size = max(0, self.end - self.start)


@dataclass
class Cluster:
    chrom: str
    start: int
    end: int
    svtype: str
    members: list[Call] = field(default_factory=list)

    @property
    def size(self) -> int:
        return self.end - self.start

    @property
    def sources(self) -> set[str]:
        return {m.source for m in self.members}


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        elif item:
            out[item] = "1"
    return out


def parse_sv_vcf(path: Path, source: str) -> list[Call]:
    calls: list[Call] = []
    if not path.exists():
        return calls
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            chrom, pos, _id, _ref, _alt, _qual, filt, info = cols[:8]
            if filt not in (".", "PASS"):
                continue
            inf = parse_info(info)
            svtype = inf.get("SVTYPE", "")
            if svtype not in ("DEL", "DUP"):
                continue
            start = int(pos) - 1
            if "END" in inf:
                end = int(inf["END"])
            elif "SVLEN" in inf:
                end = start + abs(int(float(inf["SVLEN"])))
            else:
                continue
            if end <= start:
                continue
            calls.append(Call(chrom, start, end, svtype, source, extra=filt))
    return calls


def parse_cnvpytor(path: Path, source: str, min_size: int) -> list[Call]:
    """CNVpytor -call TSV: type, region, size, cnv_level, pvals..."""
    calls: list[Call] = []
    if not path.exists():
        return calls
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) < 3:
                continue
            kind = parts[0].lower()
            if kind.startswith("del"):
                svtype = "DEL"
            elif kind.startswith("dup"):
                svtype = "DUP"
            else:
                continue
            region = parts[1]
            m = re.match(r"([^:]+):(\d+)-(\d+)", region)
            if not m:
                continue
            chrom, a, b = m.group(1), int(m.group(2)), int(m.group(3))
            start, end = min(a, b) - 1, max(a, b)
            if end - start < min_size:
                continue
            level = parts[3] if len(parts) > 3 else ""
            calls.append(Call(chrom, start, end, svtype, source, extra=str(level)))
    return calls


def parse_spectre_bed(path: Path, source: str, min_size: int) -> list[Call]:
    calls: list[Call] = []
    if not path.exists():
        return calls
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("chrom"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            blob = "\t".join(cols[3:]).upper()
            if "DEL" in blob or "LOSS" in blob or "DELETION" in blob:
                svtype = "DEL"
            elif "DUP" in blob or "GAIN" in blob or "DUPLICATION" in blob:
                svtype = "DUP"
            else:
                continue
            if end - start < min_size:
                continue
            calls.append(Call(chrom, start, end, svtype, source))
    return calls


def load_mask_bed(path: Path | None) -> dict[str, list[Interval]]:
    by_chrom: dict[str, list[Interval]] = {}
    if path is None or not path.exists():
        return by_chrom
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            by_chrom.setdefault(chrom, []).append(Interval(chrom, start, end))
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda x: (x.start, x.end))
    return by_chrom


def mask_coverage_fraction(chrom: str, start: int, end: int, masks: dict[str, list[Interval]]) -> float:
    size = end - start
    if size <= 0:
        return 0.0
    intervals = masks.get(chrom, [])
    covered = 0
    for iv in intervals:
        inter = max(0, min(end, iv.end) - max(start, iv.start))
        covered += inter
    return min(1.0, covered / size)


def is_female(sex: str) -> bool:
    return (sex or "").strip().lower() in FEMALE_SEX


def reciprocal_overlap(a: Call | Cluster, b: Call | Cluster) -> float:
    if a.chrom != b.chrom:
        return 0.0
    inter = max(0, min(a.end, b.end) - max(a.start, b.start))
    if inter == 0:
        return 0.0
    return min(inter / max(a.size, 1), inter / max(b.size, 1))


def cluster_calls(calls: list[Call], ro: float) -> list[Cluster]:
    remaining = sorted(calls, key=lambda c: (c.chrom, c.start, c.end))
    clusters: list[Cluster] = []
    while remaining:
        seed = remaining.pop(0)
        cl = Cluster(seed.chrom, seed.start, seed.end, seed.svtype, [seed])
        changed = True
        while changed:
            changed = False
            kept: list[Call] = []
            for c in remaining:
                if c.svtype == cl.svtype and reciprocal_overlap(cl, c) >= ro:
                    cl.members.append(c)
                    cl.start = min(cl.start, c.start)
                    cl.end = max(cl.end, c.end)
                    changed = True
                else:
                    kept.append(c)
            remaining = kept
        clusters.append(cl)
    return clusters


def hard_blocked(cl: Cluster, sex: str, masks: dict[str, list[Interval]], mask_frac: float) -> bool:
    if is_female(sex) and cl.chrom == CHRY:
        return True
    return mask_coverage_fraction(cl.chrom, cl.start, cl.end, masks) >= mask_frac


def passes_size_cap(cl: Cluster, max_event: int, ont_depth: bool, wgs_depth: bool, ont_sv: bool, wgs_sv: bool) -> bool:
    """Mb-scale events need both-depth; oversize without SV support is dropped."""
    if cl.size <= max_event:
        return True
    both_depth = ont_depth and wgs_depth
    any_sv = ont_sv or wgs_sv
    return both_depth and any_sv


def assign_confidence(
    cl: Cluster,
    min_depth: int,
    max_event: int,
    sex: str,
    masks: dict[str, list[Interval]],
    mask_frac: float,
) -> str | None:
    src = cl.sources
    ont_sv = "ont_sv" in src
    ont_depth = "ont_depth" in src
    wgs_sv = "wgs_sv" in src
    wgs_depth = "wgs_depth" in src
    ont = ont_sv or ont_depth
    wgs = wgs_sv or wgs_depth
    both_sv = ont_sv and wgs_sv
    size = cl.size
    blocked = hard_blocked(cl, sex, masks, mask_frac)

    if ont and wgs:
        if both_sv and size < min_depth:
            return "SHARED_SV"
        if size >= min_depth:
            if blocked:
                return "MASKED"
            if not passes_size_cap(cl, max_event, ont_depth, wgs_depth, ont_sv, wgs_sv):
                return "MASKED"
            return "LARGE_HIGH"
        # cross-platform <100 kb without both SV (rare)
        return "MEDIUM"

    if ont_sv and not wgs:
        return "ONT_SV"

    if ont_depth and not ont_sv and not wgs and size >= min_depth:
        if blocked or not passes_size_cap(cl, max_event, True, False, False, False):
            return "MASKED"
        return "ONT_DEPTH"

    if wgs_depth and not ont and size >= min_depth:
        if blocked or not passes_size_cap(cl, max_event, False, True, False, False):
            return "MASKED"
        return "WGS_DEPTH"

    return None


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--ont-sv", type=Path, required=True)
    p.add_argument("--wgs-sv", type=Path, required=True)
    p.add_argument("--ont-cnvpytor", type=Path, required=True)
    p.add_argument("--wgs-cnvpytor", type=Path, required=True)
    p.add_argument("--ont-spectre", type=Path, default=None)
    p.add_argument("--ro", type=float, default=RO_DEFAULT)
    p.add_argument("--min-depth", type=int, default=MIN_DEPTH)
    p.add_argument("--max-event", type=int, default=MAX_EVENT_DEFAULT)
    p.add_argument("--hard-mask", type=Path, default=None)
    p.add_argument("--mask-frac", type=float, default=MASK_FRAC)
    p.add_argument("--sex", default="", help="F/female/XX drops chrY; M/male/XY keeps Y under mask only")
    p.add_argument("-o", "--outdir", type=Path, required=True)
    args = p.parse_args()

    masks = load_mask_bed(args.hard_mask)
    calls: list[Call] = []
    calls += parse_sv_vcf(args.ont_sv, "ont_sv")
    calls += parse_sv_vcf(args.wgs_sv, "wgs_sv")
    calls += parse_cnvpytor(args.ont_cnvpytor, "ont_depth", args.min_depth)
    calls += parse_cnvpytor(args.wgs_cnvpytor, "wgs_depth", args.min_depth)
    if args.ont_spectre:
        calls += parse_spectre_bed(args.ont_spectre, "ont_depth", args.min_depth)

    clusters = cluster_calls(calls, args.ro)
    args.outdir.mkdir(parents=True, exist_ok=True)
    all_path = args.outdir / f"{args.sample}.cnv.bed"
    high_path = args.outdir / f"{args.sample}.cnv.high.bed"
    shared_path = args.outdir / f"{args.sample}.cnv.shared_sv.bed"
    counts = {
        "LARGE_HIGH": 0,
        "SHARED_SV": 0,
        "MEDIUM": 0,
        "ONT_SV": 0,
        "ONT_DEPTH": 0,
        "WGS_DEPTH": 0,
        "MASKED": 0,
        "DROP": 0,
    }

    header = [
        "chrom", "start", "end", "svtype", "size", "confidence",
        "sources", "n_ont_sv", "n_ont_depth", "n_wgs_sv", "n_wgs_depth",
    ]
    with (
        all_path.open("w", newline="") as fh_all,
        high_path.open("w", newline="") as fh_high,
        shared_path.open("w", newline="") as fh_shared,
    ):
        w_all = csv.writer(fh_all, delimiter="\t")
        w_high = csv.writer(fh_high, delimiter="\t")
        w_shared = csv.writer(fh_shared, delimiter="\t")
        w_all.writerow(header)
        w_high.writerow(header)
        w_shared.writerow(header)
        for cl in sorted(clusters, key=lambda x: (x.chrom, x.start)):
            conf = assign_confidence(
                cl, args.min_depth, args.max_event, args.sex, masks, args.mask_frac
            )
            if conf is None:
                counts["DROP"] += 1
                continue
            counts[conf] += 1
            srcs = sorted(cl.sources)
            row = [
                cl.chrom, cl.start, cl.end, cl.svtype, cl.size, conf,
                ",".join(srcs),
                sum(m.source == "ont_sv" for m in cl.members),
                sum(m.source == "ont_depth" for m in cl.members),
                sum(m.source == "wgs_sv" for m in cl.members),
                sum(m.source == "wgs_depth" for m in cl.members),
            ]
            w_all.writerow(row)
            if conf == "LARGE_HIGH":
                w_high.writerow(row)
            elif conf == "SHARED_SV":
                w_shared.writerow(row)

    summary = args.outdir / f"{args.sample}.cnv.summary.tsv"
    with summary.open("w") as fh:
        fh.write(
            "sample\tinput_calls\tclusters_kept\tLARGE_HIGH\tSHARED_SV\tMEDIUM\t"
            "ONT_SV\tONT_DEPTH\tWGS_DEPTH\tMASKED\tDROP\tsex\n"
        )
        kept = (
            counts["LARGE_HIGH"]
            + counts["SHARED_SV"]
            + counts["MEDIUM"]
            + counts["ONT_SV"]
            + counts["ONT_DEPTH"]
            + counts["WGS_DEPTH"]
            + counts["MASKED"]
        )
        fh.write(
            f"{args.sample}\t{len(calls)}\t{kept}\t{counts['LARGE_HIGH']}\t{counts['SHARED_SV']}\t"
            f"{counts['MEDIUM']}\t{counts['ONT_SV']}\t{counts['ONT_DEPTH']}\t"
            f"{counts['WGS_DEPTH']}\t{counts['MASKED']}\t{counts['DROP']}\t{args.sex or 'NA'}\n"
        )
    print(f"merged {args.sample}: {counts} -> {all_path}")
    print(f"LARGE_HIGH={counts['LARGE_HIGH']} SHARED_SV={counts['SHARED_SV']} MASKED={counts['MASKED']}")


if __name__ == "__main__":
    main()
