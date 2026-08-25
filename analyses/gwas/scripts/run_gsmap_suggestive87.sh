#!/bin/bash
# Run gsMap quick_mode (spatial LDSC + Cauchy) for suggestive-87 boosted sumstats.
# Matches prior successful RPL_OLA1 / RPL runs (ldscore_save_format=quick_mode).
set -euo pipefail
source /home/administrator/gsmap_env/bin/activate

WORKDIR=/mnt/d/gsMap/example_quick_mode/Mouse_Embryo
SAMPLE=E16.5_E1S1.MOSTA
TRAIT=RPL_sug87
SUMSTATS=/mnt/d/gsMap/RPL_GWAS/RPL_suggestive87_boost.sumstats.gz
RESOURCE=/mnt/d/gsMap/gsMap_resource
H5AD=/mnt/d/gsMap/gsMap_example_data/ST/E16.5_E1S1.MOSTA.h5ad
HOMOLOG=/mnt/d/gsMap/gsMap_resource/homologs/mouse_human_homologs.txt
OUTDIR=/mnt/d/ONT/figure3/gwas/gsmap_RPL_results
LOG=/mnt/d/gsMap/RPL_GWAS/rpl_sug87_run.stdout.log

mkdir -p "$OUTDIR"
echo "[$(date)] start quick_mode trait=$TRAIT" | tee "$LOG"

gsmap quick_mode \
  --workdir "$WORKDIR" \
  --sample_name "$SAMPLE" \
  --gsMap_resource_dir "$RESOURCE" \
  --hdf5_path "$H5AD" \
  --annotation annotation \
  --data_layer count \
  --trait_name "$TRAIT" \
  --sumstats_file "$SUMSTATS" \
  --homolog_file "$HOMOLOG" \
  --max_processes 8 \
  2>&1 | tee -a "$LOG"

SPATIAL="$WORKDIR/$SAMPLE/spatial_ldsc/${SAMPLE}_${TRAIT}.csv.gz"
CAUCHY="$WORKDIR/$SAMPLE/cauchy_combination/${SAMPLE}_${TRAIT}.Cauchy.csv.gz"
cp -f "$SUMSTATS" "$OUTDIR/RPL_suggestive87_boost.sumstats.gz"
cp -f "$SPATIAL" "$OUTDIR/RPL_sug87_spatial_ldsc_spot_level.csv.gz"
cp -f "$CAUCHY" "$OUTDIR/RPL_sug87_cauchy_celltype_level.csv.gz"
# also copy small cauchy into analyses module results
cp -f "$CAUCHY" /mnt/d/ONT/analyses/gwas/results/RPL_sug87_cauchy_celltype_level.csv.gz

echo "[$(date)] DONE" | tee -a "$LOG"
zcat "$CAUCHY" | head -20 | tee -a "$LOG"
