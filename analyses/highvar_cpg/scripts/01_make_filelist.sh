#!/bin/bash

# =========================================================
# STEP1: Build sample-file index for methylation data
# =========================================================

BASE_DIR="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/ONT_processed_data/methlybed_modkit"

OUT_DIR="/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/figure3/meQTL/prepare_methylation"

OUT_FILE="${OUT_DIR}/sample_file_map.tsv"

mkdir -p ${OUT_DIR}

echo -e "sample_id\tfile_path" > ${OUT_FILE}

echo "Searching .cpg.5mC.bed.gz files..."

find ${BASE_DIR} -type f -name "*.cpg.5mC.bed.gz" | while read file; do

    # extract sample name
    fname=$(basename ${file})
    sample_id=${fname%%.cpg.5mC.bed.gz}

    echo -e "${sample_id}\t${file}" >> ${OUT_FILE}

done

# summary
echo "Done!"
echo "Total files found:"
wc -l ${OUT_FILE}

echo "Output saved to:"
echo ${OUT_FILE}
