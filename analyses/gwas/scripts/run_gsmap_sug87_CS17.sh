#!/bin/bash
# gsMap on human embryo CS17 HESTA ST + RPL suggestive-87 boosted sumstats
set -euo pipefail
source /home/administrator/gsmap_env/bin/activate

WORKDIR=/mnt/d/gsMap/CS17_HESTA
SAMPLE=CS17_E1S1_HESTA
TRAIT=RPL_sug87
H5AD=/mnt/f/spatio-transcriptome-humanembryos/CS17_E1S1_HESTA.h5ad
RESOURCE=/mnt/d/gsMap/gsMap_resource
SUMSTATS=/mnt/d/gsMap/RPL_GWAS/RPL_suggestive87_boost.sumstats.gz
OUTDIR=/mnt/d/ONT/figure3/gwas/gsmap_CS17_sug87_results
LOG=/mnt/d/gsMap/RPL_GWAS/rpl_sug87_CS17_run.stdout.log

mkdir -p "$WORKDIR" "$OUTDIR"
echo "[$(date)] start quick_mode sample=$SAMPLE trait=$TRAIT" | tee "$LOG"

# Human ST: no homolog_file. annotation=celltype, layer=counts.
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
  2>&1 | tee -a "$LOG"

SPATIAL="$WORKDIR/$SAMPLE/spatial_ldsc/${SAMPLE}_${TRAIT}.csv.gz"
CAUCHY="$WORKDIR/$SAMPLE/cauchy_combination/${SAMPLE}_${TRAIT}.Cauchy.csv.gz"
cp -f "$SPATIAL" "$OUTDIR/CS17_RPL_sug87_spatial_ldsc.csv.gz"
cp -f "$CAUCHY" "$OUTDIR/CS17_RPL_sug87_cauchy.csv.gz"
cp -f "$CAUCHY" /mnt/d/ONT/analyses/gwas/results/CS17_RPL_sug87_cauchy.csv.gz 2>/dev/null || true

echo "[$(date)] DONE" | tee -a "$LOG"
zcat "$CAUCHY" | tee -a "$LOG"
