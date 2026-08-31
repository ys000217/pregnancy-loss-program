#!/usr/bin/env python3
"""QC gates from mosdepth dist/summary. Invoked by 02_qc_coverage.sh."""
from __future__ import annotations

import sys
from pathlib import Path


def stats(qcdir: Path, prefix: str):
    dist = qcdir / f"{prefix}.mosdepth.global.dist.txt"
    breadth_5 = None
    rows = []
    with dist.open() as fh:
        for line in fh:
            chrom, depth, prop = line.split()
            if chrom == "total":
                rows.append((float(depth), float(prop)))
    rows.sort()
    median = 0.0
    for depth, prop in rows:
        if prop >= 0.5:
            median = depth
        if depth == 5:
            breadth_5 = prop
    if breadth_5 is None:
        for depth, prop in rows:
            if depth >= 5:
                breadth_5 = prop
                break
    summary = qcdir / f"{prefix}.mosdepth.summary.txt"
    mean = None
    with summary.open() as fh:
        fh.readline()
        for line in fh:
            if line.startswith("total"):
                mean = float(line.split()[3])
                break
    return mean, median, breadth_5


def main() -> None:
    sample, qcdir = sys.argv[1], Path(sys.argv[2])
    min_med, min_br = float(sys.argv[3]), float(sys.argv[4])
    do_abort = sys.argv[5] not in ("0", "false", "no")
    abort_med, abort_br = float(sys.argv[6]), float(sys.argv[7])

    out = qcdir / f"{sample}.qc.tsv"
    warn = False
    abort = False
    with out.open("w") as fh:
        fh.write("sample\tplatform\tmean\tmedian\tbreadth_ge5x\tpass\n")
        for plat in ("ont", "wgs"):
            mean, median, br = stats(qcdir, f"{sample}.{plat}")
            ok = median is not None and median >= min_med
            if plat == "wgs":
                ok = ok and br is not None and br >= min_br
            warn = warn or (not ok)
            too_thin = median is None or median < abort_med
            if plat == "wgs":
                too_thin = too_thin or br is None or br < abort_br
            abort = abort or too_thin
            fh.write(f"{sample}\t{plat}\t{mean}\t{median}\t{br}\t{int(ok)}\n")
            print(f"{sample} {plat}: mean={mean:.2f} median={median} breadth>=5x={br} pass={ok}")

    if warn:
        print(
            f"WARN {sample} below warning gates (median>={min_med}, WGS breadth>=5x>={min_br}); continuing",
            file=sys.stderr,
        )
    if abort and do_abort:
        print(
            f"ERROR {sample} coverage too thin to call (median>={abort_med}, WGS breadth>=5x>={abort_br})",
            file=sys.stderr,
        )
        sys.exit(2)
    print(f"QC table: {out}")


if __name__ == "__main__":
    main()
