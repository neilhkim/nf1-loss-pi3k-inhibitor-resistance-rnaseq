#!/usr/bin/env bash
# 02_fastqc.sh
# Run FastQC on FASTQs

set -euo pipefail

RAW_DIR="data/fastq"
OUT_FIG_DIR="results/fastqc/figures"
OUT_LOG_DIR="results/fastqc/logs"

# Number of parallel FastQC processes to run at once
# Default is 8, can be overridden by first script argument 
# (e.g., ./02_fastqc.sh 8)
N_JOBS="${1:-8}"

mkdir -p "${OUT_FIG_DIR}" "${OUT_LOG_DIR}"

echo "[FastQC] Input directory: ${RAW_DIR}"
echo "[FastQC] Output figures: ${OUT_FIG_DIR}"
echo "[FastQC] Output logs: ${OUT_LOG_DIR}"
echo "[FastQC] Using ${N_JOBS} parallel jobs."
# echo "[FastQC] Starting FastQC analysis..."

# for fq in ${RAW_DIR}/*.fastq.gz; do 
#     sample=$(basename "${fq}" .fastq.gz)
#     log_files="${LOG_DIR}/${sample}.log"

#     echo "[FastQC] Processing: ${sample}"
#     fastqc "$fq" \
#         --outdir "${OUT_DIR}" \
#         > "${log_files}" 2>&1

#     echo "[FastQC] Done: ${sample}"
# done

# Find all *.fastq.gz files and run FastQC in parallel
shopt -s nullglob
files=("${RAW_DIR}"/*.fastq.gz)

if (( ${#files[@]} == 0 )); then
    echo "[FastQC] No FASTQ files found in ${RAW_DIR}."
    exit 1
fi

i=0
for file in "${files[@]}"; do
    mime=$(file --mime-type -b "$file")
    echo "$mime"
    base=$(basename "$file")
    fastqc -o "${OUT_FIG_DIR}" "$file" 2>"${OUT_LOG_DIR}/fastqc_${base}.log" &
    ((i++))
    if (( i % N_JOBS == 0 )); then
        wait
    fi
done
wait

echo "[FastQC] All samples processed."