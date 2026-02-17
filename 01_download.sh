#!/usr/bin/env bash
# =============================================================================
# 01_download.sh
# Workshop: WGS De Novo Assembly — Dr. Omics Lab
# Purpose:  Download SRR501130 paired-end reads from ENA FTP
#
# Why ENA (not SRA)?
#   NCBI SRA merges paired-end reads into a single file.
#   ENA provides them pre-split: _1 (forward) and _2 (reverse).
#
# Usage: bash scripts/01_download.sh
# =============================================================================

set -euo pipefail

OUTDIR="data/raw"
mkdir -p "${OUTDIR}" logs

echo "======================================================"
echo " Step 1: Downloading SRR501130 from ENA FTP"
echo " Organism: Staphylococcus aureus GR1"
echo "======================================================"

cd "${OUTDIR}"

# -c = continue partial downloads (safe to re-run if connection drops)
echo "[1/3] Downloading forward reads (SRR501130_1)..."
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR501/SRR501130/SRR501130_1.fastq.gz

echo ""
echo "[2/3] Downloading reverse reads (SRR501130_2)..."
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR501/SRR501130/SRR501130_2.fastq.gz

echo ""
echo "[3/3] Extracting compressed files..."
gunzip *.gz

echo ""
echo "======================================================"
echo " Download complete!"
echo " Files in ${OUTDIR}/:"
ls -lh *.fastq 2>/dev/null || echo "(check current directory)"
echo "======================================================"
