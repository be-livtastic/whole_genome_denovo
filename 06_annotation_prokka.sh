#!/usr/bin/env bash
# =============================================================================
# 06_annotation_prokka.sh
# Workshop: WGS De Novo Assembly — Dr. Omics Lab
# Purpose:  Annotate the SPAdes assembly using Prokka
#
# Input:   results/assembly/spades_result/contigs.fasta
# Output:  results/annotation/prokka_result/
#
# Usage: bash scripts/06_annotation_prokka.sh
# Runtime: 2–4 minutes (after --setupdb; hours without it!)
# =============================================================================

set -euo pipefail

# ── EDIT THIS PATH to your Prokka executable ─────────────────────────────────
PROKKA_BIN="tools/prokka/bin/prokka"
# Workshop approach: cd prokka/bin && ./prokka contigs.fasta --outdir prokka_result
# ─────────────────────────────────────────────────────────────────────────────

CONTIGS="results/assembly/spades_result/contigs.fasta"
OUTDIR="results/annotation/prokka_result"
mkdir -p logs

echo "======================================================"
echo " Step 6: Genome Annotation with Prokka"
echo " Organism: Staphylococcus aureus GR1"
echo "======================================================"
echo ""
echo " Prokka annotation pipeline (internal):"
echo "   1. Prodigal  → predict protein-coding ORFs"
echo "   2. BLAST+    → match proteins to known database"
echo "   3. HMMER     → identify protein domains/families"
echo "   4. Barrnap   → detect ribosomal RNA genes"
echo "   5. Aragorn   → detect tRNA genes"
echo "   6. Combine all → write output files"
echo ""

# Check inputs
if [[ ! -f "${CONTIGS}" ]]; then
    echo "ERROR: Assembly not found: ${CONTIGS}"
    echo "       Run 04_assembly_spades.sh first."
    exit 1
fi

if [[ ! -f "${PROKKA_BIN}" ]]; then
    echo "ERROR: Prokka not found at: ${PROKKA_BIN}"
    echo "       Run 05_install_prokka.sh first."
    exit 1
fi

echo "[1/1] Running Prokka annotation..."
"${PROKKA_BIN}" \
    "${CONTIGS}" \
    --outdir "${PROKKA_RESULT_DIR:-${OUTDIR}}" \
    2>&1 | tee logs/06_prokka.log

# ── Results summary ───────────────────────────────────────────────────────────
TXT_REPORT=$(find "${OUTDIR}" -name "*.txt" | head -1)
if [[ -f "${TXT_REPORT}" ]]; then
    echo ""
    echo "======================================================"
    echo " Annotation complete! Summary:"
    echo ""
    cat "${TXT_REPORT}"
    echo ""
    echo " Output files in ${OUTDIR}/:"
    echo "   .gff   Gene features — open in IGV or Artemis"
    echo "   .gbk   GenBank format — open in Geneious"
    echo "   .faa   Protein FASTA — use for BLAST or phylogenetics"
    echo "   .ffn   CDS nucleotide FASTA"
    echo "   .txt   This summary"
    echo "   .tbl   Feature table for NCBI GenBank submission"
    echo "======================================================"
else
    echo "WARNING: Could not find .txt summary. Check logs/06_prokka.log"
fi
