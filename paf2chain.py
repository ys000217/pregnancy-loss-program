#!/usr/bin/env python3
"""
Convert minimap2 PAF (with cg tag, i.e. run with -c) to UCSC chain format.

Direction convention (matches CrossMap & UCSC liftOver):
  chain tName = SOURCE genome (input coords)   <- PAF tname (the reference)
  chain qName = TARGET genome (output coords)  <- PAF qname (the query)

So to lift CN1 -> GRCh38, run:
  minimap2 -cx asm5 -c -t8 CN1.fa GRCh38.fa > CN1_to_GRCh38.paf
then:
  python3 paf2chain.py CN1_to_GRCh38.paf CN1_to_GRCh38.chain

Usage:
  python3 paf2chain.py input.paf output.chain [--min-span N]
"""
import sys, re

def parse_cigar(cg):
    return [(int(m.group(1)), m.group(2))
            for m in re.finditer(r'(\d+)([MIDNSHP=X])', cg)]

def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    opts = [a for a in sys.argv[1:] if a.startswith('--')]
    min_span = 0
    for o in opts:
        if o.startswith('--min-span='):
            min_span = int(o.split('=', 1)[1])
    if len(args) != 2:
        sys.stderr.write("usage: paf2chain.py input.paf output.chain [--min-span=N]\n")
        sys.exit(1)
    paf_path, chain_path = args

    out = open(chain_path, 'w')
    cid = 0
    n_in = n_out = 0
    for line in open(paf_path):
        n_in += 1
        f = line.rstrip('\n').split('\t')
        if len(f) < 12:
            continue
        qname, qlen = f[0], int(f[1])
        qs, qe = int(f[2]), int(f[3])
        strand = f[4]
        tname, tlen = f[5], int(f[6])
        ts, te = int(f[7]), int(f[8])
        cg = None
        for tag in f[12:]:
            if tag.startswith('cg:Z:'):
                cg = tag[5:]
                break
        if cg is None:
            continue
        if te - ts < min_span:
            continue
        if strand not in ('+', '-'):
            continue

        ops = parse_cigar(cg)
        # build blocks (M-runs) and gaps
        blocks = []
        dts = []   # source(tName) gap after each block
        dqs = []   # target(qName) gap after each block
        i = 0
        n = len(ops)
        while i < n:
            ln, op = ops[i]
            if op in 'M=X':
                size = 0
                while i < n and ops[i][1] in 'M=X':
                    size += ops[i][0]
                    i += 1
                dt = dq = 0
                while i < n and ops[i][1] in 'IDN':
                    ln2, op2 = ops[i]
                    if op2 == 'I':
                        dq += ln2
                    elif op2 in 'DN':
                        dt += ln2
                    i += 1
                blocks.append(size)
                dts.append(dt)
                dqs.append(dq)
            else:
                i += 1  # leading clip/gap

        if not blocks:
            continue

        # header
        if strand == '+':
            q_start, q_end = qs, qe
        else:
            q_start, q_end = qe, qs

        out.write("chain 1000 %s %d + %d %d %s %d %s %d %d %d\n" %
                  (tname, tlen, ts, te, qname, qlen, strand, q_start, q_end, cid))
        nb = len(blocks)
        for j in range(nb):
            if j == nb - 1:
                out.write("%d\n" % blocks[j])
            else:
                out.write("%d %d %d\n" % (blocks[j], dts[j], dqs[j]))
        cid += 1
        n_out += 1

    out.close()
    sys.stderr.write("paf2chain: read %d PAF lines, wrote %d chains\n" % (n_in, n_out))

if __name__ == '__main__':
    main()
