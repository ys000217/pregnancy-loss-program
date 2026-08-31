#!/usr/bin/env bash
# Fix Windows CRLF -> LF if scripts were copied from Windows without git eol=lf.
# Run once on the cluster: bash scripts/fix_crlf.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
find "${ROOT}" -type f \( -name '*.sh' -o -name 'config.sh' \) -print0 |
  while IFS= read -r -d '' f; do
    sed -i 's/\r$//' "${f}"
    echo "fixed ${f}"
  done
echo "Done. Re-run: bash scripts/check_env.sh"
