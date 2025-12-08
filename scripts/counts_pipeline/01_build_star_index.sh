#!/usr/bin/env bash
# 01_build_star_index.sh
# Build STAR genome index from downloaded FASTA + GTF

set -euo pipefail

REF_DIR="data/reference"
STAR_INDEX_DIR="data/star_index"

THREADS="${1:-8}"
SJDB_OVERHANG="${2:-100}"

if [[ ! -f "${REF_DIR}/GRCh38.primary_assembly.genome.fa" ]] || [[ ! -f "${REF_DIR}/gencode.v45.annotation.gtf" ]]; then
  echo "[star-index] ERROR: Missing reference files in ${REF_DIR}. Run 00_download_reference.sh first."
  exit 1
fi

mkdir -p "${STAR_INDEX_DIR}"

echo "[star-index] Building STAR index at ${STAR_INDEX_DIR} (threads=${THREADS}, sjdbOverhang=${SJDB_OVERHANG})..."
STAR \
  --runThreadN "${THREADS}" \
  --runMode genomeGenerate \
  --genomeDir "${STAR_INDEX_DIR}" \
  --genomeFastaFiles "${REF_DIR}/GRCh38.primary_assembly.genome.fa" \
  --sjdbGTFfile "${REF_DIR}/gencode.v45.annotation.gtf" \
  --sjdbOverhang "${SJDB_OVERHANG}"

echo "[star-index] STAR index build completed."
