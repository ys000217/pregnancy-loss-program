#!/usr/bin/env python3
"""Merge ONT + WGS DEL/DUP callsets for one 10x paired sample.

Confidence tiers (cohort-facing):
  LARGE_HIGH  both platforms have depth (ont_depth+wgs_depth), same type, RO>=0.5,
              size>=100kb, outside hard mask, primary autosomes;
              >=1 Mb also needs both-platform 500 kb depth when those callsets exist;
              >10 Mb needs dual depth only (no SV required)
  SHARED_SV   both platforms have SV breakpoints, size < 100kb (small SV burden)
  MEDIUM      cross-platform same-type overlap that is not dual-depth LARGE_HIGH
              (e.g. ONT SV + WGS depth, or >=1 Mb missing one side's 500 kb)
  ONT_SV / ONT_DEPTH / WGS_DEPTH / MASKED / DROP
              ONT-only SV/depth are KEPT (not dropped) — long-read unique events

Sex chromosomes (X/Y) are dropped by default (DROP_SEX_CHROM).
Hard mask is required for LARGE_HIGH (--require-hard-mask).
Cohort breakpoint clustering of LARGE_HIGH is 09_cluster_breakpoints.py, not this script.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

RO_DEFAULT = 0.5
MIN_DEPTH = 100_000
MAX_EVENT_DEFAULT = 10_000_000
MASK_FRAC = 0.50
LARGE_BIN_MIN = 1_000_000  # ≥1 Mb: prefer 500 kb depth support when available

# RefSeq primary autosomes only (sex chr dropped unless --keep-sex-chrom)
PRIMARY_AUTOSOMES = {
    "NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000004.12",
    "NC_000005.10", "NC_000006.12", "NC_000007.14", "NC_000008.11",
    "NC_000009.12", "NC_000010.11", "NC_000011.10", "NC_000012.12",
    "NC_000013.11", "NC_000014.9", "NC_000015.10", "NC_000016.10",
    "NC_000017.11", "NC_000018.10", "NC_000019.10", "NC_000020.11",
    "NC_000021.9", "NC_000022.11",
}
SEX_CHROMS = {"NC_000023.11", "NC_000024.10"}

# CNVpytor -call columns (0-based): type, region, size, level, e1, e2, e3, e4, q0, pN, dG
DEFAULT_Q0_MAX = 0.5
DEFAULT_PN_MAX = 0.5
DEFAULT_EVAL_MAX = 1e-4


@dataclass
class Call:
    chrom: str
    start: int  # 0-based
    end: int
    svtype: str
    source: str
    size: int = 0
    extra: str = ""
    bin_size: int = 0  # 100000 / 500000 for depth calls

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

    @property
    def has_large_bin_depth(self) -> bool:
        """True if either platform has a 500 kb depth member (legacy helper)."""
        return self.has_ont_large_depth or self.has_wgs_large_depth

    @property
    def has_ont_large_depth(self) -> bool:
        return any(m.source == "ont_depth" and m.bin_size >= 500_000 for m in self.members)

    @property
    def has_wgs_large_depth(self) -> bool:
        return any(m.source == "wgs_depth" and m.bin_size >= 500_000 for m in self.members)


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


def _f(parts: list[str], i: int, default: float | None = None) -> float | None:
    if len(parts) <= i:
        return default
    try:
        return float(parts[i])
    except ValueError:
        return default


def parse_cnvpytor(
    path: Path,
    source: str,
    min_size: int,
    bin_size: int,
    q0_max: float,
    pn_max: float,
    eval_max: float,
    apply_qc: bool,
) -> list[Call]:
    """CNVpytor -call TSV with optional Q0 / pN / e-val1 filters."""
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
            if apply_qc:
                q0 = _f(parts, 8)
                pn = _f(parts, 9)
                e1 = _f(parts, 4)
                if q0 is not None and q0 > q0_max:
                    continue
                if pn is not None and pn > pn_max:
                    continue
                if e1 is not None and e1 > eval_max:
                    continue
            level = parts[3] if len(parts) > 3 else ""
            calls.append(
                Call(chrom, start, end, svtype, source, extra=str(level), bin_size=bin_size)
            )
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


def merge_intervals(intervals: list[Interval]) -> list[Interval]:
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x.start, x.end))
    out = [intervals[0]]
    for iv in intervals[1:]:
        last = out[-1]
        if iv.start <= last.end:
            out[-1] = Interval(last.chrom, last.start, max(last.end, iv.end))
        else:
            out.append(iv)
    return out


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
    for chrom in list(by_chrom):
        by_chrom[chrom] = merge_intervals(by_chrom[chrom])
    return by_chrom


def mask_coverage_fraction(chrom: str, start: int, end: int, masks: dict[str, list[Interval]]) -> float:
    size = end - start
    if size <= 0:
        return 0.0
    covered = 0
    for iv in masks.get(chrom, []):
        covered += max(0, min(end, iv.end) - max(start, iv.start))
    return min(1.0, covered / size)


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


def hard_blocked(cl: Cluster, masks: dict[str, list[Interval]], mask_frac: float) -> bool:
    return mask_coverage_fraction(cl.chrom, cl.start, cl.end, masks) >= mask_frac


def passes_size_cap(cl: Cluster, max_event: int, ont_depth: bool, wgs_depth: bool) -> bool:
    """Events > MAX (default 10 Mb): dual-platform depth is enough.

    Aneuploidy / chromosome-arm dosage usually has no Sniffles/DELLY breakpoint.
    Requiring SV here would drop the most important pregnancy-loss events.
    """
    if cl.size <= max_event:
        return True
    return ont_depth and wgs_depth


def assign_confidence(
    cl: Cluster,
    min_depth: int,
    max_event: int,
    masks: dict[str, list[Interval]],
    mask_frac: float,
    require_large_bin_for_mb: bool,
) -> str | None:
    src = cl.sources
    ont_sv = "ont_sv" in src
    ont_depth = "ont_depth" in src
    wgs_sv = "wgs_sv" in src
    wgs_depth = "wgs_depth" in src
    ont = ont_sv or ont_depth
    wgs = wgs_sv or wgs_depth
    both_sv = ont_sv and wgs_sv
    dual_depth = ont_depth and wgs_depth
    size = cl.size
    blocked = hard_blocked(cl, masks, mask_frac)

    if ont and wgs:
        if both_sv and size < min_depth:
            return "SHARED_SV"
        if size >= min_depth:
            if blocked:
                return "MASKED"
            # LARGE_HIGH requires both-platform depth (not SV+depth or SV+SV alone)
            if not dual_depth:
                return "MEDIUM"
            if not passes_size_cap(cl, max_event, ont_depth, wgs_depth):
                return "MEDIUM"
            # ≥1 Mb: both ONT and WGS 500 kb depth when large callsets were provided
            if require_large_bin_for_mb and size >= LARGE_BIN_MIN:
                if not (cl.has_ont_large_depth and cl.has_wgs_large_depth):
                    return "MEDIUM"
            return "LARGE_HIGH"
        return "MEDIUM"

    # Platform-unique calls are retained (ONT advantage for loci WGS misses)
    if not wgs:
        if ont_depth and size >= min_depth:
            if blocked or not passes_size_cap(cl, max_event, True, False):
                return "MASKED"
            return "ONT_DEPTH"
        if ont_sv:
            return "ONT_SV"

    if wgs_depth and not ont and size >= min_depth:
        if blocked or not passes_size_cap(cl, max_event, False, True):
            return "MASKED"
        return "WGS_DEPTH"

    return None


def keep_chrom(chrom: str, keep_sex: bool) -> bool:
    if chrom in PRIMARY_AUTOSOMES:
        return True
    if keep_sex and chrom in SEX_CHROMS:
        return True
    return False


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--ont-sv", type=Path, required=True)
    p.add_argument("--wgs-sv", type=Path, required=True)
    p.add_argument("--ont-cnvpytor", type=Path, required=True)
    p.add_argument("--wgs-cnvpytor", type=Path, required=True)
    p.add_argument("--ont-cnvpytor-large", type=Path, default=None, help="500 kb ONT CNVpytor TSV")
    p.add_argument("--wgs-cnvpytor-large", type=Path, default=None, help="500 kb WGS CNVpytor TSV")
    p.add_argument("--ont-spectre", type=Path, default=None)
    p.add_argument("--ro", type=float, default=RO_DEFAULT)
    p.add_argument("--min-depth", type=int, default=MIN_DEPTH)
    p.add_argument("--max-event", type=int, default=MAX_EVENT_DEFAULT)
    p.add_argument("--hard-mask", type=Path, default=None)
    p.add_argument("--mask-frac", type=float, default=MASK_FRAC)
    p.add_argument("--require-hard-mask", action="store_true",
                   help="Fail if hard-mask missing (needed for LARGE_HIGH)")
    p.add_argument("--keep-sex-chrom", action="store_true",
                   help="Keep chrX/Y; default drops sex chromosomes")
    p.add_argument("--cnvpytor-q0-max", type=float, default=DEFAULT_Q0_MAX)
    p.add_argument("--cnvpytor-pn-max", type=float, default=DEFAULT_PN_MAX)
    p.add_argument("--cnvpytor-eval-max", type=float, default=DEFAULT_EVAL_MAX)
    p.add_argument("--no-cnvpytor-qc", action="store_true", help="Disable Q0/pN/e-val filters")
    p.add_argument("--sex", default="", help="Retained for summary only; sex chr dropped by default")
    p.add_argument("-o", "--outdir", type=Path, required=True)
    args = p.parse_args()

    if args.require_hard_mask and (args.hard_mask is None or not args.hard_mask.exists() or not args.hard_mask.stat().st_size):
        print("ERROR hard mask required but missing; set HARD_MASK_BED", file=sys.stderr)
        sys.exit(2)

    masks = load_mask_bed(args.hard_mask)
    apply_qc = not args.no_cnvpytor_qc
    qc_kw = dict(
        q0_max=args.cnvpytor_q0_max,
        pn_max=args.cnvpytor_pn_max,
        eval_max=args.cnvpytor_eval_max,
        apply_qc=apply_qc,
    )

    calls: list[Call] = []
    calls += parse_sv_vcf(args.ont_sv, "ont_sv")
    calls += parse_sv_vcf(args.wgs_sv, "wgs_sv")
    calls += parse_cnvpytor(args.ont_cnvpytor, "ont_depth", args.min_depth, 100_000, **qc_kw)
    calls += parse_cnvpytor(args.wgs_cnvpytor, "wgs_depth", args.min_depth, 100_000, **qc_kw)
    if args.ont_cnvpytor_large:
        calls += parse_cnvpytor(args.ont_cnvpytor_large, "ont_depth", args.min_depth, 500_000, **qc_kw)
    if args.wgs_cnvpytor_large:
        calls += parse_cnvpytor(args.wgs_cnvpytor_large, "wgs_depth", args.min_depth, 500_000, **qc_kw)
    # Only enforce ≥1 Mb ↔ 500 kb concordance when both large TSVs actually exist
    require_large_bin = bool(
        args.ont_cnvpytor_large
        and args.wgs_cnvpytor_large
        and args.ont_cnvpytor_large.is_file()
        and args.wgs_cnvpytor_large.is_file()
        and args.ont_cnvpytor_large.stat().st_size > 0
        and args.wgs_cnvpytor_large.stat().st_size > 0
    )

    if args.ont_spectre:
        calls += parse_spectre_bed(args.ont_spectre, "ont_depth", args.min_depth)

    before = len(calls)
    calls = [c for c in calls if keep_chrom(c.chrom, args.keep_sex_chrom)]
    dropped_chr = before - len(calls)

    clusters = cluster_calls(calls, args.ro)
    args.outdir.mkdir(parents=True, exist_ok=True)
    all_path = args.outdir / f"{args.sample}.cnv.bed"
    high_path = args.outdir / f"{args.sample}.cnv.high.bed"
    shared_path = args.outdir / f"{args.sample}.cnv.shared_sv.bed"
    counts = {
        "LARGE_HIGH": 0, "SHARED_SV": 0, "MEDIUM": 0, "ONT_SV": 0,
        "ONT_DEPTH": 0, "WGS_DEPTH": 0, "MASKED": 0, "DROP": 0,
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
                cl, args.min_depth, args.max_event, masks, args.mask_frac, require_large_bin
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
            "ONT_SV\tONT_DEPTH\tWGS_DEPTH\tMASKED\tDROP\tsex\tdropped_non_autosome\n"
        )
        kept = sum(counts[k] for k in counts if k != "DROP")
        fh.write(
            f"{args.sample}\t{len(calls)}\t{kept}\t{counts['LARGE_HIGH']}\t{counts['SHARED_SV']}\t"
            f"{counts['MEDIUM']}\t{counts['ONT_SV']}\t{counts['ONT_DEPTH']}\t"
            f"{counts['WGS_DEPTH']}\t{counts['MASKED']}\t{counts['DROP']}\t"
            f"{args.sex or 'NA'}\t{dropped_chr}\n"
        )
    print(f"merged {args.sample}: {counts} -> {all_path}")
    print(
        f"LARGE_HIGH={counts['LARGE_HIGH']} SHARED_SV={counts['SHARED_SV']} "
        f"MASKED={counts['MASKED']} dropped_non_autosome={dropped_chr} "
        f"cnvpytor_qc={apply_qc} large_bin={require_large_bin}"
    )


if __name__ == "__main__":
    main()
