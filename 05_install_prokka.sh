#!/usr/bin/env bash
# =============================================================================
# 05_install_prokka.sh
# Workshop: WGS De Novo Assembly — Dr. Omics Lab
# Purpose:  Install Prokka and all dependencies for bacterial genome annotation
#
# Prokka uses multiple tools internally:
#   Prodigal  → CDS / ORF prediction
#   BLAST+    → protein similarity search
#   HMMER     → protein domain / family detection
#   Barrnap   → rRNA gene detection (16S, 23S, 5S)
#   Aragorn   → tRNA gene detection
#   BioPerl   → Perl bioinformatics libraries
#
# Usage: bash scripts/05_install_prokka.sh
# =============================================================================

set -euo pipefail

TOOLS_DIR="tools"
mkdir -p "${TOOLS_DIR}" logs

echo "======================================================"
echo " Step 5: Installing Prokka + Dependencies"
echo "======================================================"

# ── 1. System dependencies ────────────────────────────────────────────────────
echo "[1/5] Installing system dependencies..."
sudo apt update
sudo apt install -y \
    git bioperl perl-doc \
    hmmer barrnap aragorn \
    ncbi-blast+ parallel prodigal \
    wget unzip cpanminus

# ── 2. Perl module ────────────────────────────────────────────────────────────
echo ""
echo "[2/5] Installing Bio::SearchIO::hmmer3 Perl module..."
sudo cpanm Bio::SearchIO::hmmer3 2>&1 | tee logs/05_cpanm.log

# ── 3. Clone Prokka ───────────────────────────────────────────────────────────
echo ""
echo "[3/5] Cloning Prokka from GitHub..."
cd "${TOOLS_DIR}"
git clone https://github.com/tseemann/prokka.git

# ── 4. Set up databases ───────────────────────────────────────────────────────
echo ""
echo "[4/5] Setting up Prokka internal databases..."
echo "      (This indexes protein/rRNA references — CRITICAL for speed!)"
cd prokka
./bin/prokka --setupdb 2>&1 | tee "../../logs/05_prokka_setupdb.log"

# ── 5. Fix tbl2asn if missing ─────────────────────────────────────────────────
echo ""
echo "[5/5] Checking for tbl2asn (required for GenBank output)..."
if ! command -v tbl2asn &>/dev/null; then
    echo "      tbl2asn not found — downloading..."
    wget -q https://ftp.ncbi.nlm.nih.gov/toolbox/ncbi_tools/converters/by_program/table2asn/linux64.table2asn.gz
    gunzip linux64.table2asn.gz
    mv linux64.table2asn tbl2asn
    chmod +x tbl2asn
    sudo mv tbl2asn /usr/local/bin/
    echo "      ✓ tbl2asn installed to /usr/local/bin/"
else
    echo "      ✓ tbl2asn already available"
fi

echo ""
echo "======================================================"
echo " Prokka installation complete!"
echo ""
echo " Executable: $(realpath bin/prokka)"
echo ""
echo " To annotate your assembly:"
echo "   1. Copy contigs.fasta to tools/prokka/bin/"
echo "   2. cd tools/prokka/bin"
echo "   3. ./prokka contigs.fasta --outdir prokka_result"
echo "======================================================"
