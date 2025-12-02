#!/usr/bin/env bash
# scripts/download_fastq_all.sh
# Download all FASTQ files for the project using prefetch and fasterq-dump

set -euo pipefail

# Directory to store FASTQ files
OUTDIR="data/fastq"
mkdir -p "${OUTDIR}"

# List of all SRR runs for GSE207514
SRR_LIST="
SRR19987570
SRR19987571
SRR19987572
SRR19987573
SRR19987574
SRR19987575
SRR19987576
SRR19987577
SRR19987578
SRR19987579
SRR19987580
SRR19987581
SRR19987582
SRR19987583
SRR19987584
SRR19987585
SRR19987586
SRR19987587
SRR19987588
SRR19987589
SRR19987590
SRR19987591
SRR19987592
SRR19987593
SRR19987594
SRR19987595
SRR19987596
SRR19987597
SRR19987598
SRR19987599
SRR19987600
SRR19987601
SRR19987602
SRR19987603
SRR19987604
SRR19987605
"

# Number of threads for fasterq-dump
THREADS=8

for SRR in ${SRR_LIST}; do 
    FASTQ_GZ="${OUTDIR}/${SRR}.fastq.gz"

    # If compressed FASTQ already exists, skip
    if [ -f "${FASTQ_GZ}" ]; then
        echo "File ${FASTQ_GZ} already exists. Skipping download."
        continue
    fi

    echo "=== Processing ${SRR} ==="

    # 1) Download .sra to local cache (if not already cached)
    echo "Running prefetch for ${SRR}..."
    prefetch "${SRR}"

    # 2) Convert .sra to FASTQ using fasterq-dump
    echo "Running fasterq-dump for ${SRR}..."
    fasterq-dump "${SRR}" --threads "${THREADS}" --outdir "${OUTDIR}"

    # 3) Compress the FASTQ file
    echo "Compressing ${SRR}.fastq..."
    gzip -f "${OUTDIR}/${SRR}.fastq"

    echo "Finished processing ${SRR}."
    echo
done

echo "All requested SRR files have been processed."