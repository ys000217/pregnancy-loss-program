#!/usr/bin/env python3
"""Cluster LARGE_HIGH calls across samples by similar breakpoints.

Within-sample merge (07) uses 50% RO and can union slightly offset bins.
Cohort burden needs the Nature-style rule: same chrom/type, both ends close,
and still reciprocal-overlap so a 200 kb DEL is not grouped with a multi-Mb DEL
that happens to start nearby.

Default pad = 100 kb (CNVpytor primary bin). Events ≥1 Mb on both sides use
the 500 kb pad. Cluster coordinates are median start/end of members (not union).
"""
from __future__ import annotations

import argparse
import csv
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path

BP_PAD = 100_000
BP_PAD_LARGE = 500_000
LARGE_MIN = 1_000_000
RO_DEFAULT = 0.50


@dataclass
class Event:
    sample: str
    chrom: str
    start: int
    end: int
    svtype: str
    size: int
    extra: dict[str, str]


def reciprocal_overlap(a: Event, b: Event) -> float:
    if a.chrom != b.chrom:
        return 0.0
    inter = max(0, min(a.end, b.end) - max(a.start, b.start))
    if inter == 0:
        return 0.0
    return min(inter / max(a.size, 1), inter / max(b.size, 1))


def similar_breakpoints(
    a: Event,
    b: Event,
    pad: int,
    pad_large: int,
    large_min: int,
    ro: float,
) -> bool:
    if a.chrom != b.chrom or a.svtype != b.svtype:
        return False
    window = pad_large if min(a.size, b.size) >= large_min else pad
    if abs(a.start - b.start) > window or abs(a.end - b.end) > window:
        return False
    return reciprocal_overlap(a, b) >= ro


def cluster_events(
    events: list[Event],
    pad: int,
    pad_large: int,
    large_min: int,
    ro: float,
) -> list[list[Event]]:
    n = len(events)
    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    for i in range(n):
        for j in range(i + 1, n):
            if similar_breakpoints(events[i], events[j], pad, pad_large, large_min, ro):
                union(i, j)

    buckets: dict[int, list[Event]] = {}
    for i, ev in enumerate(events):
        buckets.setdefault(find(i), []).append(ev)
    out = list(buckets.values())
    out.sort(key=lambda g: (g[0].chrom, min(e.start for e in g), g[0].svtype))
    return out


def parse_high_bed(path: Path) -> list[Event]:
    sample = path.name.replace(".cnv.high.bed", "")
    events: list[Event] = []
    with path.open() as fh:
        header = fh.readline()
        if not header:
            return events
        cols = header.rstrip("\n").split("\t")
        named = cols[0] == "chrom"
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            if named:
                row = dict(zip(cols, parts))
                chrom, start, end = row["chrom"], int(row["start"]), int(row["end"])
                svtype = row.get("svtype", "DEL")
                size = int(row["size"]) if row.get("size") else end - start
                extra = {k: row[k] for k in row if k not in ("chrom", "start", "end")}
            else:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                svtype = parts[3]
                size = int(parts[4]) if len(parts) > 4 else end - start
                extra = {"raw": "\t".join(parts[3:])}
            if end <= start:
                continue
            events.append(Event(sample, chrom, start, end, svtype, size, extra))
    return events


def collect_beds(indir: Path, glob_pat: str) -> list[Path]:
    return sorted(p for p in indir.glob(glob_pat) if p.is_file() and p.stat().st_size > 0)


def write_outputs(groups: list[list[Event]], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    cl_path = outdir / "cohort.large_high.clusters.tsv"
    mem_path = outdir / "cohort.large_high.members.tsv"
    sum_path = outdir / "cohort.large_high.cluster_summary.tsv"

    n_singleton = sum(1 for g in groups if len({e.sample for e in g}) == 1)
    n_recurrent = len(groups) - n_singleton
    samples_all = sorted({e.sample for g in groups for e in g})

    with cl_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "cluster_id", "chrom", "median_start", "median_end", "svtype",
            "median_size", "n_calls", "n_samples", "samples",
            "min_start", "max_end", "start_span", "end_span",
        ])
        for i, g in enumerate(groups, start=1):
            starts = [e.start for e in g]
            ends = [e.end for e in g]
            samples = sorted({e.sample for e in g})
            med_s = int(statistics.median(starts))
            med_e = int(statistics.median(ends))
            w.writerow([
                f"CL{i:05d}", g[0].chrom, med_s, med_e, g[0].svtype,
                med_e - med_s, len(g), len(samples), ",".join(samples),
                min(starts), max(ends), max(starts) - min(starts), max(ends) - min(ends),
            ])

    with mem_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "cluster_id", "sample", "chrom", "start", "end", "svtype", "size", "n_samples_in_cluster",
        ])
        for i, g in enumerate(groups, start=1):
            n_s = len({e.sample for e in g})
            cid = f"CL{i:05d}"
            for e in sorted(g, key=lambda x: (x.sample, x.start)):
                w.writerow([cid, e.sample, e.chrom, e.start, e.end, e.svtype, e.size, n_s])

    with sum_path.open("w") as fh:
        fh.write("n_input_samples\tn_calls\tn_clusters\tn_singleton\tn_recurrent\n")
        fh.write(
            f"{len(samples_all)}\t{sum(len(g) for g in groups)}\t{len(groups)}\t"
            f"{n_singleton}\t{n_recurrent}\n"
        )

    print(
        f"clustered {sum(len(g) for g in groups)} LARGE_HIGH calls from {len(samples_all)} samples "
        f"-> {len(groups)} loci ({n_recurrent} recurrent, {n_singleton} singleton) -> {cl_path}"
    )


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--indir", type=Path, required=True, help="Directory of per-sample cnv.high.bed")
    p.add_argument("--glob", default="*.cnv.high.bed")
    p.add_argument("-o", "--outdir", type=Path, required=True)
    p.add_argument("--bp-pad", type=int, default=BP_PAD, help="Max |Δstart| and |Δend| (bp) for <1 Mb")
    p.add_argument("--bp-pad-large", type=int, default=BP_PAD_LARGE, help="Pad when both events ≥1 Mb")
    p.add_argument("--large-min", type=int, default=LARGE_MIN)
    p.add_argument("--ro", type=float, default=RO_DEFAULT)
    args = p.parse_args()

    beds = collect_beds(args.indir, args.glob)
    if not beds:
        print(f"ERROR no files matching {args.indir}/{args.glob}", file=sys.stderr)
        sys.exit(2)

    events: list[Event] = []
    for path in beds:
        events.extend(parse_high_bed(path))
    if not events:
        args.outdir.mkdir(parents=True, exist_ok=True)
        print(f"WARN 0 LARGE_HIGH events in {len(beds)} beds; wrote empty tables to {args.outdir}")
        write_outputs([], args.outdir)
        return

    groups = cluster_events(events, args.bp_pad, args.bp_pad_large, args.large_min, args.ro)
    write_outputs(groups, args.outdir)


if __name__ == "__main__":
    main()
