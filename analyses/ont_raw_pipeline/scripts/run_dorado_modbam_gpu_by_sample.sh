#!/bin/bash
#JSUB -q gpu
#JSUB -n 8
#JSUB -gpgpu "1"
#JSUB -J dorado_modbam_by_sample
#JSUB -o /data/home/4250160/hpc_logs/dorado_modbam_by_sample.out.%J
#JSUB -e /data/home/4250160/hpc_logs/dorado_modbam_by_sample.err.%J
#JSUB -cwd /data/home/4250160

source ~/.bashrc

SAMPLE="${SAMPLE:?}"

SAMPLE_DIR="/data/home/4250160/pod5/${SAMPLE}"
cd "$SAMPLE_DIR" || exit 1

REF="/data/home/4250160/refseq/grch38.p14/GRCh38.p14.fasta"
BASE_MODEL="/data/home/4250160/dna_r10.4.1_e8.2_400bps_hac@v5.2.0"
MOD_MODEL="/data/home/4250160/dna_r10.4.1_e8.2_400bps_hac@v5.2.0_5mC_5hmC@v2"
MODEL_DIR="/data/home/4250160"

OUT_ROOT="/data/home/4250160/modbam/${SAMPLE}"
mkdir -p "$OUT_ROOT"

conda activate dorado

find . -type f -name "*.pod5" | sort | while read -r POD5_FILE; do
    BATCH_DIR=$(dirname "$POD5_FILE")
    BATCH_NAME=$(basename "$BATCH_DIR")
    OUT_DIR="${OUT_ROOT}/${BATCH_NAME}_mods.bam"
    rm -rf "$OUT_DIR"
    mkdir -p "$OUT_DIR"
    /data/home/4250160/run_dorado.sh basecaller \
        --device cuda:0 \
        --reference "$REF" \
        --models-directory "$MODEL_DIR" \
        "$BASE_MODEL" \
        --modified-bases-models "$MOD_MODEL" \
        "$POD5_FILE" \
        -o "$OUT_DIR"
done

conda deactivate

if find "/data/home/4250160/modbam/${SAMPLE}" -type f -path "*/bam_pass/*.bam" | grep -q .; then
    rm -rf "/data/home/4250160/pod5/${SAMPLE}"
fi
