#!/usr/bin/env bash
# =============================================================================
# 07_blast_uniprot.sh
# Workshop: WGS De Novo Assembly — Dr. Omics Lab
# Purpose:  BLAST Prokka protein predictions against UniProt SwissProt
#           for deep functional annotation
#
# Why SwissProt?
#   It's the manually curated, high-confidence section of UniProt.
#   Every entry is reviewed by a biologist — much more reliable
#   than TrEMBL (automated, unreviewed).
#
# Usage: bash scripts/07_blast_uniprot.sh
# Note:  SwissProt FASTA download = ~90 MB compressed
# =============================================================================

set -euo pipefail

BLAST_DIR="results/blast"
PROKKA_DIR="results/annotation/prokka_result"
mkdir -p "${BLAST_DIR}" logs

# Find protein predictions from Prokka (.faa file)
QUERY=$(find "${PROKKA_DIR}" -name "*.faa" | head -1)

echo "======================================================"
echo " Step 7: Functional Annotation via BLAST"
echo " Query:  Prokka protein predictions (.faa)"
echo " DB:     UniProt SwissProt (curated proteins)"
echo "======================================================"

# ── Install BLAST+ ────────────────────────────────────────────────────────────
if ! command -v blastp &>/dev/null; then
    echo "[0/4] Installing NCBI BLAST+..."
    sudo apt-get install -y ncbi-blast+
fi

cd "${BLAST_DIR}"

# ── 1. Download SwissProt FASTA ───────────────────────────────────────────────
echo ""
echo "[1/4] Downloading UniProt SwissProt database (~90 MB)..."
wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip -f uniprot_sprot.fasta.gz

# ── 2. Build BLAST database ───────────────────────────────────────────────────
echo ""
echo "[2/4] Building BLAST protein database from SwissProt..."
makeblastdb \
    -in uniprot_sprot.fasta \
    -dbtype prot \
    -out uniprodb \
    2>&1 | tee "../../logs/07_makeblastdb.log"

# ── 3. Run BLASTP ─────────────────────────────────────────────────────────────
echo ""
echo "[3/4] Running BLASTP (query vs SwissProt)..."
echo "      This may take 15–60 minutes..."

if [[ -z "${QUERY}" ]]; then
    echo "WARNING: No .faa file found in ${PROKKA_DIR}"
    echo "         Using 'longest_orfs.pep' as in the workshop."
    QUERY="longest_orfs.pep"
fi

blastp \
    -query "${QUERY}" \
    -db uniprodb \
    -num_alignments 1 \
    -num_threads 4 \
    -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen" \
    -out result.txt \
    2>&1 | tee "../../logs/07_blastp.log"

# ── 4. Quick summary ──────────────────────────────────────────────────────────
echo ""
echo "[4/4] Results summary:"
echo "      Total hits: $(wc -l < result.txt)"
echo "      Output: ${BLAST_DIR}/result.txt"

echo ""
echo "======================================================"
echo " BLAST complete!"
echo ""
echo " Output columns (outfmt 6):"
echo "   qseqid  sseqid  pident  length  mismatch  gaps"
echo "   qstart  qend  sstart  send  evalue  bitscore  qlen"
echo ""
echo " % Query Coverage formula:"
echo "   coverage = (length - mismatch - gaps) / qlen * 100"
echo ""
echo " To look up a SwissProt hit:"
echo "   https://www.uniprot.org/uniprot/<sseqid>"
echo "======================================================"
