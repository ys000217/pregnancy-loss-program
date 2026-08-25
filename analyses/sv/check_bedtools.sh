#!/bin/bash
echo "=== bedtools/bcftools 可用性 ==="
ls /home/administrator/miniconda3/envs/*/bin/bedtools 2>/dev/null
ls /home/administrator/miniconda3/envs/*/bin/bcftools 2>/dev/null
command -v bedtools || echo "bedtools 不在 PATH"
