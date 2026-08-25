#!/bin/bash
export PATH=/home/administrator/miniconda3/envs/annotsv/bin:$PATH
echo "=== AnnotSV -help ==="
AnnotSV -help 2>&1 | head -60
