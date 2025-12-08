#!/usr/bin/env bash
# 03_trim_reads.sh
# Trim adapters or low-quality reads

set -euo pipefail

INPUT_DIR="data/fastq"
OUTPUT_DIR="data/trimmed_fastq"
LOG_DIR="results/trim/logs"
THREADS="${1:-8}"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

echo "[fastp] Input directory: ${INPUT_DIR}"
echo "[fastp] Output directory: ${OUTPUT_DIR}"
echo "[fastp] Using ${THREADS} threads per fastp run."

# Find all *.fastq.gz files and run fastp in parallel
for fq in "${INPUT_DIR}"/*.fastq.gz; do
    # If the glob does not match anything, bash will keep the literal pattern
    # so guard against that
    if [[ "${fq}" == "${INPUT_DIR}/*.fastq.gz" ]]; then
        echo "[fastp] No FASTQ files found in ${INPUT_DIR}."
        exit 1
    fi

    base=$(basename "$fq")
    sample="${base%.fastq.gz}"

    out_fq="${OUTPUT_DIR}/${sample}_trimmed.fastq.gz"
    out_html="${LOG_DIR}/${sample}_fastp.html"
    out_json="${LOG_DIR}/${sample}_fastp.json"

    echo "[fastp] Processing sample: ${sample}"
    echo "[fastp] Input file: ${fq}"
    echo "[fastp] Output file: ${out_fq}"

    fastp \
        -i "$fq" \
        -o "$out_fq" \
        -h "$out_html" \
        -j "$out_json" \
        -w "${THREADS}"

    echo "[fastp] Done processing sample: ${sample}"
    echo
done

echo "[fastp] All samples processed."