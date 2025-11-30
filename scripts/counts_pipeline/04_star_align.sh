#!/usr/bin/env bash
# 04_star_align.sh
# Align with STAR

set -euo pipefail

# Config
RUNS_TSV="data/metadata/runs.tsv"
STAR_INDEX="data/star_index"
TRIMMED_DIR="data/trimmed_fastq"
OUT_DIR="results/star"
LOG_DIR="analysis/00_counts_pipeline/logs/star"

# How many STAR processes to run in parallel
# My WSL kept failing so I am trying 1 jobs * 1 thread to limit resource usage. 
STAR_THREADS=1
MAX_JOBS=1
LIMIT_RAM=3000000000

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

echo "[$(date)] Starting STAR alignment"
echo "Using runs table: ${RUNS_TSV}"
echo "STAR index: ${STAR_INDEX}"
echo "Trimmed FASTQ dir: ${TRIMMED_DIR}"
echo "Output directory: ${OUT_DIR}"
echo "Log directory: ${LOG_DIR}"
echo "Threads per STAR job: ${STAR_THREADS}, max concurrent jobs: ${MAX_JOBS}"
echo


# Function that runs STAR for a single run
run_star_for_run() {
    local run_id="$1"
    local sample_id="$2"
    local group_id="$3"
    local fq_trimmed="$4"

    local prefix="${OUT_DIR}/${run_id}_"
    local bam="${prefix}Aligned.sortedByCoord.out.bam"
    local log_file="${LOG_DIR}/${run_id}.log"

    if [[ -f "${bam}" ]]; then
        echo "[$(date)] [${run_id}] BAM already exists at ${bam}. Skipping."
        return 0
    fi

    echo "[$(date)] [${run_id}] Starting STAR"
    echo "  Sample: ${sample_id} (group: ${group_id})"
    echo "  FASTQ:  ${fq_trimmed}"
    echo "  Prefix: ${prefix}"
    echo "  Log:    ${log_file}"

    STAR \
        --runThreadN "${STAR_THREADS}" \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "${fq_trimmed}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${prefix}" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --limitBAMsortRAM "${LIMIT_RAM}" \
        > "${log_file}" 2>&1

    if [[ -f "${bam}" ]]; then
        echo "[$(date)] [${run_id}] Finished STAR successfully."
    else
        echo "[$(date)] [${run_id}] ERROR: BAM not found at ${bam}. Check log: ${log_file}"
        return 1
    fi
}

# Function to wait until we have a free slot
wait_for_slot() {
    while true; do
        local njobs
        njobs=$(jobs -rp | wc -l)

        if (( njobs < MAX_JOBS )); then
            break
        fi

        sleep 5
    done
}

# Main loop over runs.tsv (skip header)
# Expected columns: run_id  sample_id  group_id  fastq_path
tail -n +2 "${RUNS_TSV}" | while IFS=$'\t' read -r run_id sample_id group_id fastq_path; do
    if [[ -z "${run_id}" ]]; then
        continue
    fi

    # Construct trimmed FASTQ path from run_id
    fq_trimmed="${TRIMMED_DIR}/${run_id}_trimmed.fastq.gz"

    if [[ ! -f "${fq_trimmed}" ]]; then
        echo "[$(date)] [${run_id}] WARNING: trimmed FASTQ not found at ${fq_trimmed}. Skipping."
        continue
    fi

    wait_for_slot
    run_star_for_run "${run_id}" "${sample_id}" "${group_id}" "${fq_trimmed}" &
done

# Wait for all background jobs to finish
wait

echo
echo "[$(date)] All STAR alignments completed."