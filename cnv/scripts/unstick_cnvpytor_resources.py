#!/usr/bin/env python3
"""Satisfy CNVpytor Genome.check_resources() without `cnvpytor -download`.

Conda CNVpytor 1.3.x: -download crashes (PosixPath.split) and the package
often ships without gc_hg19/hg38 files. check_resources() then logs ERROR
and main() returns 0 before -gc/-rd. Empty placeholder files are enough
for that existence check. Calling still uses our RefSeq -conf + GC pytor.
"""
from __future__ import annotations

from pathlib import Path

from cnvpytor.genome import Genome


def main() -> None:
    n = 0
    for spec in Genome.reference_genomes.values():
        for key in ("gc_file", "mask_file"):
            if key not in spec:
                continue
            path = Path(spec[key])
            if path.exists() and path.stat().st_size > 0:
                continue
            path.parent.mkdir(parents=True, exist_ok=True)
            if not path.exists():
                path.touch()
                n += 1
                print(f"placeholder {path}")
    if not Genome.check_resources():
        raise SystemExit("ERROR check_resources() still false after placeholders")
    print(f"CNVpytor bundled resources OK (created {n} placeholders)")


if __name__ == "__main__":
    main()
