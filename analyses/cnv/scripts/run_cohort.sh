#!/usr/bin/env bash
# Submit one sample per line. Wrap inner command in your scheduler if needed.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "${ROOT}/scripts/lib_source_config.sh"

if [[ ! -s "${MANIFEST}" ]]; then
  echo "ERROR ${MANIFEST} missing" >&2
  exit 1
fi

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r ont_id rest; do
  echo "sample ${ont_id}"
  bash "${ROOT}/scripts/run_sample.sh" "${ont_id}"
done
