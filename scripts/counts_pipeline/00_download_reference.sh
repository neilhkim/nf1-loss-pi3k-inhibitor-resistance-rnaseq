#!/usr/bin/env bash
# 00_download_reference.sh
# Fetch GENCODE GRCh38 reference (FASTA + GTF) and build STAR index

set -euo pipefail

REF_DIR="data/reference"

mkdir -p "${REF_DIR}" 

echo "[reference] Downloading GENCODE v45 annotation GTF..."
wget -q -O "${REF_DIR}/gencode.v45.annotation.gtf.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz"
gunzip -f "${REF_DIR}/gencode.v45.annotation.gtf.gz"

echo "[reference] Downloading GRCh38 primary assembly FASTA..."
wget -q -O "${REF_DIR}/GRCh38.primary_assembly.genome.fa.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz"
gunzip -f "${REF_DIR}/GRCh38.primary_assembly.genome.fa.gz"

echo "[reference] Reference download completed. To build STAR index, run 01_build_star_index.sh."
