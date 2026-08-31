#!/usr/bin/env python3
"""Overwrite BAM-autodetected 'hg38' in a sample .pytor without Root().

Root.__init__ opens bundled gc_hg38/mask_hg38 when the pytor still says hg38.
Those conda placeholders are empty, IO aborts, and set_reference_genome never
runs. Write the signals through IO only, then later -his uses our -conf GC.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from cnvpytor.genome import Genome
from cnvpytor.io import IO


def main() -> None:
    conf, root, genome_id = sys.argv[1], sys.argv[2], sys.argv[3]
    if not Path(root).is_file():
        raise SystemExit(f"ERROR pytor missing: {root}")
    Genome.load_reference_genomes(conf)
    if genome_id not in Genome.reference_genomes:
        raise SystemExit(f"ERROR genome id {genome_id!r} not in {conf}")
    spec = Genome.reference_genomes[genome_id]
    gc = spec.get("gc_file")
    if not gc or not Path(gc).is_file() or Path(gc).stat().st_size < 1000:
        raise SystemExit(f"ERROR bad gc_file for {genome_id}: {gc!r}")

    io = IO(root)
    io.create_signal(None, None, "reference genome", np.array([np.bytes_(genome_id)]))
    # GC from conf; mask off (no RefSeq mask, and hg38 placeholders are empty)
    io.create_signal(None, None, "use reference", np.array([1, 0]).astype("uint8"))

    # Re-open read-only check
    io2 = IO(root, ro=True)
    name = np.array(io2.get_signal(None, None, "reference genome")).astype(str)[0]
    use = np.array(io2.get_signal(None, None, "use reference"))
    if name != genome_id:
        raise SystemExit(f"ERROR write failed: still {name!r}")
    print(
        f"forced genome={name} use_gc={int(use[0])} use_mask={int(use[1])} gc_file={gc}",
        flush=True,
    )


if __name__ == "__main__":
    main()
