#!/usr/bin/env python3
"""Write a CNVpytor -conf module for RefSeq-named GRCh38 (NC_*, not chr1)."""
from __future__ import annotations

import argparse
from pathlib import Path


def primary_chroms(fai: Path) -> list[tuple[str, int, str]]:
    out: list[tuple[str, int, str]] = []
    for line in fai.read_text().splitlines():
        if not line.strip():
            continue
        name, length_s = line.split("\t")[0], line.split("\t")[1]
        if not name.startswith("NC_0000"):
            continue
        num = int(name.split(".")[0].replace("NC_0000", ""))
        if num < 1 or num > 24:
            continue
        kind = "S" if num >= 23 else "A"
        out.append((name, int(length_s), kind))
    out.sort(key=lambda x: (0 if x[2] == "A" else 1, x[0]))
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--fai", type=Path, required=True)
    p.add_argument("--gc", type=Path, required=True, help="path to -make_gc_file pytor")
    p.add_argument("-o", type=Path, required=True)
    args = p.parse_args()
    chroms = primary_chroms(args.fai)
    if len(chroms) < 22:
        raise SystemExit(f"ERROR expected 22–24 NC_0000* chromosomes in {args.fai}, got {len(chroms)}")
    rows = ",\n            ".join(
        f'("{n}", ({L}, "{k}"))' for n, L, k in chroms
    )
    gc = args.gc.resolve().as_posix()
    args.o.parent.mkdir(parents=True, exist_ok=True)
    args.o.write_text(
        "# Generated for CNVpytor -conf. Chromosome names match this project's FASTA/BAM.\n"
        "from collections import OrderedDict\n\n"
        "import_reference_genomes = {\n"
        '    "GRCh38p14_refseq": {\n'
        '        "name": "GRCh38.p14 RefSeq",\n'
        '        "species": "human",\n'
        "        \"chromosomes\": OrderedDict([\n"
        f"            {rows}\n"
        "        ]),\n"
        f'        "gc_file": r"{gc}",\n'
        "    }\n"
        "}\n"
    )
    print(f"wrote {args.o} ({len(chroms)} chromosomes)")


if __name__ == "__main__":
    main()
