#!/bin/bash
#JSUB -q normal
#JSUB -n 56
#JSUB -J convert_pod5_by_sample
#JSUB -o /data/home/4250160/hpc_logs/convert_pod5_by_sample.out.%J
#JSUB -e /data/home/4250160/hpc_logs/convert_pod5_by_sample.err.%J
#JSUB -cwd /data/home/4250160
#JSUB -R span[hosts=1]

source ~/.bashrc
source /data/home/4250160/miniconda3/etc/profile.d/conda.sh

SAMPLE="${SAMPLE:?}"

BASE_DIR="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_Lab_shared_data/ONT_trio/02.ONT/01.RawData"
POD5_ROOT="/data/home/4250160/pod5/${SAMPLE}"

echo "[INFO] Start SAMPLE=${SAMPLE}"
echo "[INFO] POD5_ROOT=${POD5_ROOT}"
echo "[INFO] Start time: $(date)"

conda activate modkit

find "$BASE_DIR" -type d -path "*/Release-*/Data-*/${SAMPLE}" | sort | while read -r SAMPLE_DIR; do
    echo "[INFO] Found SAMPLE_DIR: $SAMPLE_DIR"

    find "$SAMPLE_DIR" -maxdepth 4 -type d -name "*fast5_pass" | sort | while read -r FAST5_DIR; do
        BATCH_NAME=$(basename "$(dirname "$FAST5_DIR")")
        OUT_DIR="${POD5_ROOT}/${BATCH_NAME}"

        echo "[INFO]   Processing batch: ${BATCH_NAME}"
        echo "[INFO]   FAST5_DIR: ${FAST5_DIR}"
        echo "[INFO]   OUT_DIR: ${OUT_DIR}"
        echo "[INFO]   Start pod5 convert: $(date)"

        mkdir -p "$OUT_DIR"
        pod5 convert fast5 "$FAST5_DIR" -o "$OUT_DIR" --threads 56

        echo "[INFO]   Finished batch: ${BATCH_NAME} at $(date)"
    done
done

conda deactivate

echo "[INFO] Finished SAMPLE=${SAMPLE} at $(date)"
