#!/usr/bin/env bash
# =============================================================================
# 04_assembly_spades.sh
# Workshop: WGS De Novo Assembly — Dr. Omics Lab
# Purpose:  Assemble SRR501130 reads into a draft genome using SPAdes
#
# Assembly strategy:
#   --isolate : For pure bacterial cultures (single organism)
#   No trimming performed — raw reads used directly (QC passed)
#
# Usage: bash scripts/04_assembly_spades.sh
# Runtime: 30–90 minutes depending on hardware
# =============================================================================

set -euo pipefail

# ── EDIT THIS PATH to your SPAdes bin location ──────────────────────────────
SPADES_BIN="tools/spades-4.2.0/bin/spades.py"
# Example if running from spades bin dir: ./spades.py
# Workshop example: /mnt/f/Dr_OmicsLab/NGS/Whole-genome-Assembly/spades-4.2.0/bin/spades.py
# ─────────────────────────────────────────────────────────────────────────────

READ1="data/raw/SRR501130_1.fastq"
READ2="data/raw/SRR501130_2.fastq"
OUTDIR="results/assembly/spades_result"
mkdir -p logs

echo "======================================================"
echo " Step 4: De Novo Assembly with SPAdes"
echo " Organism: Staphylococcus aureus GR1"
echo "======================================================"
echo ""
echo " Parameters:"
echo "   --isolate   Pure culture (single bacteria, not metagenome)"
echo "   -m 250      Max memory: 250 GB"
echo "   -1 / -2     Forward and reverse paired-end reads"
echo "   -o          Output directory"
echo ""

# Check inputs
for f in "${READ1}" "${READ2}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: File not found: ${f}"
        echo "       Run 01_download.sh first."
        exit 1
    fi
done

if [[ ! -f "${SPADES_BIN}" ]]; then
    echo "ERROR: SPAdes not found at: ${SPADES_BIN}"
    echo "       Edit SPADES_BIN at the top of this script."
    exit 1
fi

echo "[1/1] Running SPAdes assembly (this may take 30–90 minutes)..."
echo ""

"${SPADES_BIN}" \
    --isolate \
    -m 250 \
    -1 "${READ1}" \
    -2 "${READ2}" \
    -o "${OUTDIR}" \
    2>&1 | tee logs/04_spades_assembly.log

# ── Summary ───────────────────────────────────────────────────────────────────
CONTIGS="${OUTDIR}/contigs.fasta"
if [[ -f "${CONTIGS}" ]]; then
    echo ""
    echo "======================================================"
    echo " Assembly complete!"
    echo ""
    echo " Key output files in ${OUTDIR}/:"
    echo "   contigs.fasta          ← PRIMARY RESULT (use for Prokka)"
    echo "   scaffolds.fasta        ← Linked contigs with N-gaps"
    echo "   assembly_graph.gfa     ← Visualise in Bandage"
    echo "   spades.log             ← Detailed process log"
    echo "   K21/ K33/ K55/         ← k-mer intermediates (ignore)"
    echo ""
    echo " Contig header format example:"
    echo "   >NODE_1_length_163380_cov_139.527782"
    echo "    ^node  ^length in bp    ^read coverage (higher = better)"
    echo ""
    echo " Quick stats:"
    grep -c "^>" "${CONTIGS}" | xargs -I{} echo "   Total contigs: {}"
    echo "======================================================"
    echo ""
    echo " NEXT STEP: Copy contigs.fasta to prokka/bin/ and annotate."
else
    echo "ERROR: contigs.fasta not found. Check logs/04_spades_assembly.log"
    exit 1
fi

# ── HAMMER note ───────────────────────────────────────────────────────────────
echo " TIP: If you see HAMMER errors, add --only-assembler to skip"
echo "      the read correction step."
