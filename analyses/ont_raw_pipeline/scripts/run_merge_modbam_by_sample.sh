#!/bin/bash
#JSUB -q normal
#JSUB -n 8
#JSUB -J merge_modbam_by_sample
#JSUB -o /data/home/4250160/hpc_logs/merge_modbam_by_sample.out.%J
#JSUB -e /data/home/4250160/hpc_logs/merge_modbam_by_sample.err.%J
#JSUB -cwd /data/home/4250160/modbam
#JSUB -R span[hosts=1]

set -eo pipefail

SAMPLE="${SAMPLE:?}"

THREADS=8
MIN_AGE_MINUTES=60

ROOT_DIR="/data/home/4250160/modbam"
cd "$ROOT_DIR"

RUN_LOG_DIR="$ROOT_DIR/merge_logs"
mkdir -p "$RUN_LOG_DIR"
LOG_FILE="$RUN_LOG_DIR/${SAMPLE}.log"

MICROMAMBA="/data/home/4250160/bin/micromamba"
export MAMBA_ROOT_PREFIX="/data/home/4250160/.local/share/mamba"
export MAMBA_ENVS_DIRS="/data/home/4250160/.local/share/mamba/envs"
export MAMBA_PKGS_DIRS="/data/home/4250160/.local/share/mamba/pkgs:/data/home/4250160/.mamba/pkgs"

if [ ! -x "$MICROMAMBA" ]; then
  echo "ERROR: micromamba not executable: $MICROMAMBA" >&2
  exit 2
fi

if ! "$MICROMAMBA" env list | awk 'NR>1{print $1}' | grep -qx bam_merge; then
  echo "ERROR: micromamba env 'bam_merge' not found." >&2
  exit 3
fi

samtools() {
  "$MICROMAMBA" run -n bam_merge samtools "$@"
}

has_index() { [ -s "${1}.bai" ] || [ -s "${1}.csi" ]; }

OUT_BAM="$SAMPLE/${SAMPLE}.merged.bam"

{
  echo "============================================================"
  echo "RUN MODE (NO CLEANUP: will NOT delete or move originals, or any other files, after merge)"
  echo "Sample            : $SAMPLE"
  echo "Working directory : $ROOT_DIR"
  echo "Threads           : $THREADS"
  echo "Min batch age     : ${MIN_AGE_MINUTES} minutes"
  echo "Per-sample log    : $LOG_FILE"
  echo "Output            : $OUT_BAM"
  echo "Start             : $(date)"
  echo
} | tee -a "$LOG_FILE"

if [ ! -d "$SAMPLE" ]; then
  echo "ERROR: sample directory not found: $ROOT_DIR/$SAMPLE" | tee -a "$LOG_FILE"
  exit 1
fi

if [ -s "$OUT_BAM" ] && has_index "$OUT_BAM"; then
  echo "Skip $SAMPLE (already has merged.bam + index)" | tee -a "$LOG_FILE"
  exit 0
fi

mapfile -t BATCH_DIRS < <(find "$SAMPLE" -maxdepth 1 -type d -name "*_mods.bam" | sort)
NUM_BATCHES=${#BATCH_DIRS[@]}
if [ "$NUM_BATCHES" -eq 0 ]; then
  echo "ERROR: $SAMPLE has no *_mods.bam batch directories" | tee -a "$LOG_FILE"
  exit 1
fi

STABLE_BATCHES=$(find "$SAMPLE" -maxdepth 1 -type d -name "*_mods.bam" -mmin +"$MIN_AGE_MINUTES" | wc -l)
if [ "$STABLE_BATCHES" -ne "$NUM_BATCHES" ]; then
  echo "ERROR: $SAMPLE not stable yet ($STABLE_BATCHES / $NUM_BATCHES batches older than ${MIN_AGE_MINUTES} min)" | tee -a "$LOG_FILE"
  exit 1
fi

BAMS=()
for BATCH_DIR in "${BATCH_DIRS[@]}"; do
  mapfile -t BATCH_BAMS < <(find "$BATCH_DIR" -type f -path "*/bam_pass/*.bam" | sort)
  if [ "${#BATCH_BAMS[@]}" -eq 0 ]; then
    echo "ERROR: $SAMPLE batch $BATCH_DIR has no bam_pass/*.bam" | tee -a "$LOG_FILE"
    exit 1
  fi
  BAMS+=("${BATCH_BAMS[@]}")
done
NUM_BAMS=${#BAMS[@]}
if [ "$NUM_BAMS" -eq 0 ]; then
  echo "ERROR: $SAMPLE has batch dirs but no bam_pass/*.bam" | tee -a "$LOG_FILE"
  exit 1
fi

{
  echo "Batches: $NUM_BATCHES (all stable >${MIN_AGE_MINUTES} min)"
  echo "BAMs: $NUM_BAMS"
  printf '  %s\n' "${BAMS[@]}"
  echo
} | tee -a "$LOG_FILE"

# Only remove an incomplete output merged BAM/index. Never delete inputs or other files.
rm -f "$OUT_BAM" "$OUT_BAM.bai" "$OUT_BAM.csi"

if [ "$NUM_BAMS" -eq 1 ]; then
  echo "Mode: SINGLE → copy & index (no cleanup)" | tee -a "$LOG_FILE"
  cp -p "${BAMS[0]}" "$OUT_BAM"
  samtools index -@ "$THREADS" "$OUT_BAM" 2>&1 | tee -a "$LOG_FILE"
else
  echo "Mode: MULTI → samtools merge (no cleanup)" | tee -a "$LOG_FILE"
  samtools merge \
    -@ "$THREADS" \
    -f \
    -c -p \
    --write-index \
    -o "$OUT_BAM" \
    "${BAMS[@]}" \
    2>&1 | tee -a "$LOG_FILE"
fi

if [ ! -s "$OUT_BAM" ] || ! has_index "$OUT_BAM"; then
  echo "ERROR: output BAM/index missing or empty" | tee -a "$LOG_FILE"
  exit 1
fi

if ! samtools quickcheck "$OUT_BAM"; then
  echo "ERROR: samtools quickcheck failed for $OUT_BAM" | tee -a "$LOG_FILE"
  exit 1
fi

echo "Output verified OK. Original batch files and all other files kept." | tee -a "$LOG_FILE"
echo "Done : $(date)" | tee -a "$LOG_FILE"
