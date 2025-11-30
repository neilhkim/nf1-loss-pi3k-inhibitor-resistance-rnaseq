#!/usr/bin/env bash
# Run FastQC on FASTQs

set -euo pipefail

RAW_DIR="data/fastq"
OUT_DIR="results/fastqc"
# LOG_DIR="results/logs/fastqc"

# Number of parallel FastQC processes to run at once
# Default is 8, can be overridden by first script argument 
# (e.g., ./02_fastqc.sh 8)
N_JOBS="${1:-8}"

mkdir -p "${OUT_DIR}"
# mkdir -p "${LOG_DIR}"

echo "[FastQC] Input directory: ${RAW_DIR}"
echo "[FastQC] Output directory: ${OUT_DIR}"
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
find "${RAW_DIR}" -maxdepth 1 -type f -name "*.fastq.gz" -print0 \
    | xargs -0 -n1 -P "$THREADS" -I{} bash -c '
        file="$1"
        mime=$(file --mime-type -b "$file")
        echo "$mime"
        fastqc -o "'"${OUT_DIR}"'" "$file"
    ' _ {}

echo "[FastQC] All samples processed."