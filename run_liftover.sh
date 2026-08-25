#!/bin/bash
set -e
export PATH="$HOME/miniconda3/envs/liftover/bin:$PATH"
cd /mnt/e/genotype_data/liftover

echo "=== [1/3] PAF -> chain ==="
python3 /mnt/d/ONT/paf2chain.py CN1_to_GRCh38.paf CN1_to_GRCh38.chain
echo "chain lines:"
grep -c "^chain" CN1_to_GRCh38.chain || true

echo "=== [2/3] CrossMap bed liftover of breakpoints ==="
CrossMap bed CN1_to_GRCh38.chain CN1_breakpoints.bed GRCh38_breakpoints.bed

echo "=== [3/3] summary ==="
echo "--- lifted breakpoints ---"
wc -l GRCh38_breakpoints.bed
echo "--- unmapped breakpoints ---"
wc -l GRCh38_breakpoints.bed.unmap 2>/dev/null || echo "0 (none)"
echo "=== DONE ==="
