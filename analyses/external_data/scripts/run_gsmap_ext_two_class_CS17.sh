#!/bin/bash
# Run gsMap on CS17 for external miscarriage vs pregnancy GWAS hit sets.
# Report step often fails; still copy Cauchy/spatial if present.
set -uo pipefail
source /home/administrator/gsmap_env/bin/activate

WORKDIR=/mnt/d/gsMap/CS17_HESTA
SAMPLE=CS17_E1S1_HESTA
RESOURCE=/mnt/d/gsMap/gsMap_resource
H5AD=/mnt/f/spatio-transcriptome-humanembryos/CS17_E1S1_HESTA.h5ad
OUT=/mnt/d/ONT/analyses/external_data/results/gsmap_external
LOGDIR=/mnt/d/gsMap/RPL_GWAS
mkdir -p "$OUT"

run_one () {
  local TRAIT="$1"
  local SUMSTATS="$2"
  local LOG="$LOGDIR/${TRAIT}_CS17_run.stdout.log"
  echo "[$(date)] START $TRAIT" | tee "$LOG"
  gsmap quick_mode \
    --workdir "$WORKDIR" \
    --sample_name "$SAMPLE" \
    --gsMap_resource_dir "$RESOURCE" \
    --hdf5_path "$H5AD" \
    --annotation celltype \
    --data_layer counts \
    --trait_name "$TRAIT" \
    --sumstats_file "$SUMSTATS" \
    --max_processes 8 \
    >> "$LOG" 2>&1 || echo "[$(date)] gsmap non-zero (report often fails)" | tee -a "$LOG"
  CAUCHY="$WORKDIR/$SAMPLE/cauchy_combination/${SAMPLE}_${TRAIT}.Cauchy.csv.gz"
  SPATIAL="$WORKDIR/$SAMPLE/spatial_ldsc/${SAMPLE}_${TRAIT}.csv.gz"
  if [[ ! -f "$CAUCHY" ]]; then
    echo "[$(date)] FAIL $TRAIT: missing Cauchy" | tee -a "$LOG"
    return 1
  fi
  cp -f "$CAUCHY" "$OUT/${TRAIT}_cauchy.csv.gz"
  cp -f "$SPATIAL" "$OUT/${TRAIT}_spatial_ldsc.csv.gz"
  echo "[$(date)] DONE $TRAIT" | tee -a "$LOG"
  zcat "$CAUCHY" | tee -a "$LOG"
}

run_one EXT_miscarriage /mnt/d/gsMap/RPL_GWAS/EXT_miscarriage_boost.sumstats.gz
run_one EXT_pregnancy /mnt/d/gsMap/RPL_GWAS/EXT_pregnancy_boost.sumstats.gz

echo "[$(date)] ALL DONE"
