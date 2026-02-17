#!/usr/bin/env bash
# =============================================================================
# 02_fastqc.sh
# Workshop: WGS De Novo Assembly — Dr. Omics Lab
# Purpose:  Quality assessment of raw reads with FastQC
#
# Key decision: If FastQC shows green on:
#   1. Per base sequence quality  (Phred >= 30)
#   2. Overrepresented sequences  (no hits)
#   3. Adapter content            (absent)
# → NO TRIMMING NEEDED before de novo WGS assembly.
#   (Trimming risks discarding real genomic data without a reference.)
#
# Usage: bash scripts/02_fastqc.sh
# =============================================================================

set -euo pipefail

INDIR="data/raw"
OUTDIR="results/fastqc"
mkdir -p "${OUTDIR}" logs

echo "======================================================"
echo " Step 2: Quality Control with FastQC"
echo "======================================================"

# Install FastQC if not present
if ! command -v fastqc &>/dev/null; then
    echo "[0/2] Installing FastQC..."
    sudo apt-get update -qq
    sudo apt-get install -y fastqc
fi

echo ""
echo "[1/2] Running FastQC on raw reads..."
fastqc "${INDIR}"/*.fastq --outdir "${OUTDIR}" 2>&1 | tee logs/02_fastqc.log

echo ""
echo "[2/2] Done! Open these reports in your browser:"
ls "${OUTDIR}"/*.html 2>/dev/null

echo ""
echo "======================================================"
echo " Three parameters to evaluate:"
echo "   1. Per base sequence quality  → aim for Phred >= 30"
echo "   2. Overrepresented sequences  → 'No hit' is fine"
echo "   3. Adapter content            → should be absent"
echo ""
echo " Traffic light: Green = proceed | Yellow = caution"
echo "                Red = consider trimming"
echo "======================================================"
