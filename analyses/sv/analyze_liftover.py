#!/usr/bin/env python3
"""
Analyze liftover drop rate from CN1 -> GRCh38 breakpoint BED lifting.
Usage:
  python3 analyze_liftover.py CN1_breakpoints.bed GRCh38_breakpoints.bed [unmap_file]
"""
import sys, os
from collections import defaultdict, Counter

def parse_id(bp_id):
    m = bp_id.rsplit(':', 1)
    return (m[0], m[1]) if len(m) == 2 else (bp_id, 'S')

def main():
    src, lifted = sys.argv[1], sys.argv[2]
    unmap = sys.argv[3] if len(sys.argv) > 3 else None

    sv_sides = defaultdict(set)   # sv -> {side, ...}
    total_bp = 0
    svtype_of = {}
    for line in open(src):
        f = line.rstrip('\n').split('\t')
        if len(f) < 6:
            continue
        sv, side = parse_id(f[3])
        sv_sides[sv].add(side)
        svtype_of[sv] = f[4]
        total_bp += 1

    # lifted breakpoints: set of (sv, side)
    lifted_pairs = set()
    for line in open(lifted):
        f = line.rstrip('\n').split('\t')
        if len(f) < 6:
            continue
        lifted_pairs.add(parse_id(f[3]))

    # group lifted by sv (efficient)
    lifted_by_sv = defaultdict(set)
    for sv, side in lifted_pairs:
        lifted_by_sv[sv].add(side)

    # unmapped
    unmap_pairs = set()
    if unmap and os.path.exists(unmap):
        for line in open(unmap):
            f = line.rstrip('\n').split('\t')
            if len(f) < 4:
                continue
            unmap_pairs.add(parse_id(f[3]))

    n_sv = len(sv_sides)
    fully_lifted = partially = fully_dropped = 0
    for sv, want in sv_sides.items():
        got = lifted_by_sv.get(sv, set())
        if got == want:
            fully_lifted += 1
        elif got:
            partially += 1
        else:
            fully_dropped += 1

    lifted_bp = len(lifted_pairs)
    dropped_bp = total_bp - lifted_bp

    print("=== Breakpoint-level ===")
    print("total breakpoints      : %d" % total_bp)
    print("lifted breakpoints     : %d (%.3f%%)" % (lifted_bp, 100.0*lifted_bp/total_bp))
    print("dropped breakpoints    : %d (%.3f%%)" % (dropped_bp, 100.0*dropped_bp/total_bp))
    print("")
    print("=== SV-level ===")
    print("total SVs              : %d" % n_sv)
    print("fully lifted           : %d (%.3f%%)" % (fully_lifted, 100.0*fully_lifted/n_sv))
    print("partially lifted       : %d (%.3f%%)" % (partially, 100.0*partially/n_sv))
    print("fully dropped          : %d (%.3f%%)" % (fully_dropped, 100.0*fully_dropped/n_sv))

    print("")
    print("=== Breakpoint drop by SVTYPE ===")
    drop_by_type = Counter()
    tot_by_type = Counter()
    for sv, want in sv_sides.items():
        st = svtype_of.get(sv, '?')
        tot_by_type[st] += len(want)
        got = lifted_by_sv.get(sv, set())
        drop_by_type[st] += len(want) - len(got)
    for st in sorted(tot_by_type):
        print("%-5s total_bp=%6d dropped_bp=%6d (%.3f%%)" %
              (st, tot_by_type[st], drop_by_type[st], 100.0*drop_by_type[st]/tot_by_type[st]))

if __name__ == '__main__':
    main()
