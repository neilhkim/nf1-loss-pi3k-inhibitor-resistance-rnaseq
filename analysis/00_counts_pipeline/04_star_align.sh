#!/usr/bin/env bash
# Align with STAR

set -euo pipefail

# Config
RUNS_TSV="data/metadata/runs.tsv"
STAR_INDEX="data/star_index"
OUT_DIR="results/star"
LOG_DIR="analysis/00_counts_pipeline/logs/star"

# How many STAR processes to run in parallel
# With 14 logical cores and 30 GB RAM, 3 jobs * 4 threads is reasonable.
STAR_THREADS=4
MAX_JOBS=3

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

echo "[`date`] Starting STAR alignment"
echo "Using runs table: ${RUNS_TSV}"
echo "STAR index: ${STAR_INDEX}"
echo "Output directory: ${OUT_DIR}"
echo "Log directory: ${LOG_DIR}"
echo "Threads per STAR job: ${STAR_THREADS}, max concurrent jobs: ${MAX_JOBS}"
echo

# Function that runs STAR for a single run
run_star_for_run() {
    local run_id="$1"
    local sample_id="$2"
    local group_id="$3"
    local fastq="$4"

    # Output prefix and BAM path
    local prefix="${OUT_DIR}/${run_id}_"
    local bam="${prefix}Aligned.sortedByCoord.out.bam"
    local log_file="${LOG_DIR}/${run_id}.log"

    if [[ -f "${bam}" ]]; then
        echo "[`date`] [${run_id}] BAM file already exists. Skipping."
        echo " FASTQ: ${fastq}"
        echo " Prefix: ${prefix}"
        echo " Log: ${log_file}"

        STAR \
            --runThreadN "${STAR_THREADS}" \
            --genomeDir "${STAR_INDEX}" \
            --readFilesIn "${fastq}" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${prefix}" \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \ 
            --limitBAMsortRAM 12000000000 \
            > "${log_file}" 2>&1

        if [[ -f "${bam}" ]]; then
            echo "[`date`] [${run_id}] Finished STAR successfully."
        else
            echo "[`date`] [${run_id}] ERROR: BAM not found at ${bam}. Check log: ${log_file}"
            return 1
        fi
}

# Function to wait if we already have MAX_JOBS running STAR processes
wait_for_slot() {
    while true; do 
        # Count background jobs started by this script
        local njobs
        njobs=$(jobs -rp | wc -l)

        if (( njobs < MAX_JOBS )); then
            break
        fi

        sleep 5
    done
}

# Main loop over runs.tsv (skip header)
tail -n +2 "${RUNS_TSV}" | while IFS=$'\t' read -r run_id sample_id group_id fastq; do
    if [[ -z "${run_id}" ]]; then
        continue
    fi

    if [[ ! -f "${fastq}" ]]; then
        echo "[`date`] [${run_id}] WARNING: FASTQ file not found at ${fastq}. Skipping."
        continue
    fi

    wait_for_slot
    run_star_from_run "${run_id}" "${sample_id}" "${group_id}" "${fastq}" &
done

# Wait for all background jobs to finish
wait

echo 
echo "[`date`] All STAR alignments completed."