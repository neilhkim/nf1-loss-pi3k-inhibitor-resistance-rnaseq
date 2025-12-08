#!/usr/bin/env bash
# 05_featurecounts.sh
# Run featureCounts on BAM files

set -euo pipefail

RUNS_TSV="data/metadata/runs.tsv"
BAM_DIR="data/aligned"
GTF="data/reference/gencode.v45.annotation.gtf"
OUT_DIR="data/counts/featurecounts"
LOG_DIR="results/featurecounts/logs"

# featureCounts options
THREADS="${1:-4}" # Number of threads to use (1st arg). Default: 4
STRAND=0 # 0 = unstranded, 1 = stranded, 2 = reversely stranded

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

echo "[featureCounts] Runs table: ${RUNS_TSV}"
echo "[featureCounts] BAM directory: ${BAM_DIR}"
echo "[featureCounts] GTF annotation: ${GTF}"
echo "[featureCounts] Output directory: ${OUT_DIR}"
echo "[featureCounts] Threads: ${THREADS}, Strandness: ${STRAND}"
echo

# Build list of BAM files in the same order as runs.tsv
bam_list=()

while IFS=$'\t' read -r run_id sample_id group_id fastq; do
    # skip empty lines
    [[ -z "${run_id}" ]] && continue
    # skip header if present
    if [[ "${run_id}" == "run_id" ]]; then
        continue
    fi

    bam="${BAM_DIR}/${run_id}_Aligned.sortedByCoord.out.bam"

    # # Test - skip if filename starts with SRR19987603
    # if [[ "${run_id}" == SRR19987603* ]]; then
    #     echo "[featureCounts] Skipping BAM file for ${run_id} as per test condition."
    #     continue
    # fi

    if [[ ! -f "${bam}" ]]; then
        echo "[featureCounts] WARNING: BAM file not found for ${run_id} at ${bam}. Skipping"
        continue
    fi

    bam_list+=("${bam}")
done < "${RUNS_TSV}"

if ((${#bam_list[@]} == 0)); then
    echo "[featureCounts] ERROR: No BAM files found. Exiting."
    exit 1
fi

echo "[featureCounts] Will count genes for ${#bam_list[@]} BAM files."
printf ' %s\n' "${bam_list[@]}"
echo

featureCounts \
    -T "${THREADS}" \
    -a "${GTF}" \
    -o "${OUT_DIR}/featurecounts_gene_counts.txt" \
    -g gene_id \
    -t exon \
    -s "${STRAND}" \
    "${bam_list[@]}" 2>"${LOG_DIR}/featurecounts.log"

echo
echo "[featureCounts] Done. Output written to ${OUT_DIR}/featurecounts_gene_counts.txt"
echo "[featureCounts] Log written to ${LOG_DIR}/featurecounts.log"