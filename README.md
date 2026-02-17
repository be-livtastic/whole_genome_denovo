# 🧬 Whole Genome De Novo Assembly — *Staphylococcus aureus* GR1

A complete bioinformatics pipeline for de novo assembly and annotation of a bacterial genome from raw Illumina paired-end reads. Developed during a 5-day live code-along workshop with [Dr. Omics Lab](https://www.ncbi.nlm.nih.gov/sra?term=SRX146853).

> **De novo** means we reconstruct the genome **without a reference** — purely from overlapping reads.

---

## 🧫 Dataset

| Field | Value |
|---|---|
| **Organism** | *Staphylococcus aureus* GR1 |
| **SRA Accession** | [SRR501130]((https://r.search.yahoo.com/_ylt=Awrji.wi7ZNp8QEAq2VXNyoA;_ylu=Y29sbwNiZjEEcG9zAzEEdnRpZAMEc2VjA3Ny/RV=2/RE=1772511779/RO=10/RU=https%3a%2f%2fdromicslabs.com%2f/RK=2/RS=Z_2xkO.BInIX46CAv7AuJ_rU00Q-)) |
| **Experiment** | SRX146853 |
| **Sequencing** | Illumina paired-end (76 bp reads) |
| **Total reads** | ~10,900,000 spots |
| **GC content** | ~33% |
| **Download source** | ENA (European Nucleotide Archive) |

---

## 🗺️ Pipeline Overview

```
  Raw Paired-End Reads (FASTQ)
           │
           ▼
  ① Quality Control — FastQC
     (assess base quality, adapter content, overrepresented seqs)
           │
           ▼
  ② De Novo Assembly — SPAdes v4.2.0
     (--isolate mode, paired-end, outputs contigs + scaffolds)
           │
           ▼
  ③ Genome Annotation — Prokka
     (gene prediction via Prodigal, rRNA via Barrnap, tRNA via Aragorn)
           │
           ▼
  ④ Functional Annotation — BLAST vs UniProt/SwissProt
     (protein-level function assignment, %query coverage calculation)
           │
           ▼
  Annotated Genome (.gff / .gbk / .faa) + BLAST results table
```

> **Note:** Unlike RNA-seq or variant calling workflows, **no trimming** was performed before assembly. For de novo WGS without a reference, trimming risks discarding real genomic information. The raw reads passed FastQC QC thresholds (median Phred ≥ 30, no adapter contamination flagged).

---

## 📁 Repository Structure

```
wgs-denovo-assembly/
│
├── README.md                        ← This file
│
├── scripts/
│   ├── 01_download.sh               ← wget from ENA FTP
│   ├── 02_fastqc.sh                 ← Quality assessment
│   ├── 03_install_spades.sh         ← SPAdes source install (WSL/Ubuntu)
│   ├── 04_assembly_spades.sh        ← De novo genome assembly
│   ├── 05_install_prokka.sh         ← Prokka + dependency install
│   ├── 06_annotation_prokka.sh      ← Genome annotation
│   └── 07_blast_uniprot.sh          ← Functional BLAST annotation
│
├── data/
│   ├── raw/                         ← SRR501130_1.fastq, SRR501130_2.fastq
│   └── trimmed/                     ← (not used — see note above)
│
├── results/
│   ├── fastqc/                      ← HTML quality reports
│   ├── assembly/spades_result/      ← contigs.fasta, scaffolds.fasta, etc.
│   ├── annotation/prokka_result/    ← .gff, .gbk, .faa, .ffn, .txt, .tbl
│   └── blast/                       ← result.txt (tabular BLAST output)
│
└── logs/                            ← Tool stdout/stderr logs
```

---

## 🔧 Step-by-Step Workflow

### Step 1 — Download Raw Data

Data is downloaded directly from the **ENA (European Nucleotide Archive)** FTP server using `wget`. ENA provides paired-end files already separated into `_1` (forward) and `_2` (reverse), unlike NCBI SRA which merges them.

```bash
bash scripts/01_download.sh
```

```bash
# What the script runs:
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR501/SRR501130/SRR501130_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR501/SRR501130/SRR501130_2.fastq.gz

# Extract compressed files
gunzip *.gz
```

**Why `-c`?** The `-c` (continue) flag lets `wget` resume a partial download if your connection drops — essential for large sequencing files.

**Output:** `data/raw/SRR501130_1.fastq` and `SRR501130_2.fastq`

---

### Step 2 — Quality Control with FastQC

```bash
bash scripts/02_fastqc.sh
```

```bash
# Install
sudo apt-get update
sudo apt-get install fastqc

# Run on all .fastq files
fastqc *.fastq
```

Open the `.html` reports in your browser. The three key parameters to check:

| Parameter | What it means | Our result |
|---|---|---|
| **Per base sequence quality** | Phred score per position across all reads | ✅ Median ~30 (green zone) |
| **Overrepresented sequences** | Repeated sequences > 0.1% of total reads | ✅ No hits |
| **Adapter content** | Residual sequencing adapters | ✅ Absent |

> **FastQC traffic light:** 🟢 green = proceed, 🟡 yellow = caution, 🔴 red = trimming needed.
> Our data showed all green — **no trimming was required before assembly**.

> 💡 **Phred score Q30** = 99.9% base call accuracy (1 error per 1,000 bases). This is the standard threshold for high-quality NGS data.

---

### Step 3 — Install SPAdes (from source)

SPAdes is not available in standard Ubuntu repositories and must be compiled from source.

```bash
bash scripts/03_install_spades.sh
```

```bash
# Download SPAdes v4.2.0 source tarball
wget -c https://github.com/ablab/spades/archive/refs/tags/v4.2.0.tar.gz
tar -xzf v4.2.0.tar.gz
cd spades-4.2.0/

# Install build dependencies
sudo apt-get install g++ cmake zlib1g-dev
sudo apt install -y build-essential cmake libbz2-dev libz-dev \
    libcurl4-openssl-dev libssl-dev

# Compile SPAdes
./spades_compile.sh
```

**Verify:** `cd bin && ./spades.py` — the help manual should print.

---

### Step 4 — De Novo Assembly with SPAdes

```bash
bash scripts/04_assembly_spades.sh
```

```bash
# Run SPAdes in isolate mode
/path/to/spades-4.2.0/bin/spades.py \
    --isolate \
    -m 250 \
    -1 SRR501130_1.fastq \
    -2 SRR501130_2.fastq \
    -o spades_result
```

**Parameter explanations:**

| Flag | Value | Meaning |
|---|---|---|
| `--isolate` | — | Optimised for pure single-organism cultures (not metagenomes) |
| `-m` | `250` | Max memory limit in GB |
| `-1` / `-2` | FASTQ files | Forward and reverse paired-end reads |
| `-o` | `spades_result` | Output directory |

**Key output files:**

| File | Description |
|---|---|
| `contigs.fasta` | ⭐ Primary result — assembled contiguous sequences |
| `scaffolds.fasta` | Contigs linked using paired-end distance info; gaps filled with `N`s |
| `assembly_graph.gfa` | Assembly graph — visualise in [Bandage](https://github.com/rrwick/Bandage) to check circularity |
| `spades.log` | Full log with timestamps — check here if assembly fails |
| `K21/, K33/, K55/` | Intermediate k-mer assemblies (can be ignored) |

> 💡 **How SPAdes works:** Reads are broken into **k-mers** (short substrings of length k). SPAdes builds a **de Bruijn graph** from overlapping k-mers, then traverses the graph to reconstruct longer sequences. It tries multiple k-mer sizes automatically and merges the best results.

> 💡 **contigs.fasta header format:** `>NODE_1_length_163380_cov_139.527782` — the node number, contig length in bp, and average read coverage. Higher coverage = more read support = higher confidence.

---

### Step 5 — Install Prokka

```bash
bash scripts/05_install_prokka.sh
```

```bash
# System dependencies
sudo apt update
sudo apt install -y git bioperl perl-doc hmmer barrnap aragorn \
    ncbi-blast+ parallel prodigal wget unzip
sudo apt install -y cpanminus
sudo cpanm Bio::SearchIO::hmmer3

# Clone Prokka from GitHub
git clone https://github.com/tseemann/prokka.git
cd prokka

# Critical: set up internal databases (indexes proteins/rRNA refs)
./bin/prokka --setupdb
```

> ⚠️ **`--setupdb` is critical.** Without it, Prokka runs database searches unindexed — annotation that should take 2–4 minutes can take hours.

> ⚠️ **tbl2asn error?** If you see `tbl2asn not found`, download it manually:
> ```bash
> wget https://ftp.ncbi.nlm.nih.gov/toolbox/ncbi_tools/converters/by_program/table2asn/linux64.table2asn.gz
> gunzip linux64.table2asn.gz && mv linux64.table2asn tbl2asn && chmod +x tbl2asn
> sudo mv tbl2asn /usr/local/bin/
> ```

**What each Prokka dependency does:**

| Tool | Role |
|---|---|
| **Prodigal** | Predicts protein-coding sequences (CDS) and ORFs |
| **BLAST+** | Sequence similarity searches against protein databases |
| **HMMER** | Profile HMM searches for protein domain/family identification |
| **Barrnap** | Detects ribosomal RNA genes (16S, 23S, 5S rRNA) |
| **Aragorn** | Identifies tRNA genes and anticodon sequences |
| **BioPerl** | Perl libraries for biological sequence parsing |

---

### Step 6 — Genome Annotation with Prokka

```bash
bash scripts/06_annotation_prokka.sh
```

```bash
# Copy contigs to prokka/bin (or use full path)
cp spades_result/contigs.fasta prokka/bin/
cd prokka/bin

# Annotate
./prokka contigs.fasta --outdir prokka_result
```

**Output files:**

| File | Format | Contents |
|---|---|---|
| `.gff` | GFF3 | Gene features with coordinates — use in genome browsers |
| `.gbk` | GenBank | Sequence + full annotation in one file (Geneious, Artemis) |
| `.faa` | FASTA | Protein sequences of all predicted CDSs |
| `.ffn` | FASTA | Nucleotide sequences of all predicted CDSs |
| `.txt` | Text | Summary: gene count, tRNA count, rRNA count |
| `.tbl` | Table | Feature table for NCBI GenBank submission |

> 💡 **Prokka annotation pipeline (internal):** Prodigal → CDS prediction → BLAST against protein DB → HMMER domain search → Barrnap rRNA detection → Aragorn tRNA detection → combine all → output.

---

### Step 7 — Functional Annotation with BLAST

For deeper functional insight, we BLAST predicted proteins against **UniProt/SwissProt** — a curated, high-quality protein database.

```bash
bash scripts/07_blast_uniprot.sh
```

```bash
# Install BLAST+
sudo apt-get install ncbi-blast+
mkdir blast && cd blast

# Download UniProt SwissProt database (~90 MB compressed)
wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

# a) Create searchable BLAST database from SwissProt FASTA
makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprodb

# b) Run BLAST: query = Prokka protein predictions, db = uniprodb
blastp \
    -query longest_orfs.pep \
    -db uniprodb \
    -num_alignments 1 \
    -num_threads 4 \
    -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen" \
    -out result.txt
```

**Output format columns (outfmt 6):**

| Column | Field | Meaning |
|---|---|---|
| 1 | qseqid | Query protein ID |
| 2 | sseqid | Subject (SwissProt) protein ID |
| 3 | pident | % identical residues |
| 4 | length | Alignment length |
| 5 | mismatch | Number of mismatches |
| 6 | gaps | Number of gaps |
| 7–8 | qstart/qend | Query alignment coordinates |
| 9–10 | sstart/send | Subject alignment coordinates |
| 11 | evalue | Statistical significance (lower = better) |
| 12 | bitscore | Alignment quality score |
| 13 | qlen | Total query length |

**Calculating % Query Coverage:**
```
% Query Coverage = (alignment_length - mismatches - gaps) / qlen × 100
```

---

## 📊 Results Summary

> Fill in after running the pipeline.

| Metric | Value |
|---|---|
| Organism | *Staphylococcus aureus* GR1 |
| SRA Accession | SRR501130 |
| Raw reads (total) | ~10,900,000 |
| Read length | 76 bp |
| Trimming performed | No (FastQC passed all thresholds) |
| Assembler | SPAdes v4.2.0 (`--isolate`) |
| Total contigs | _(fill in)_ |
| Total assembly size | _(fill in)_ bp |
| Largest contig | _(fill in)_ bp |
| N50 | _(fill in)_ bp |
| Annotator | Prokka |
| Predicted CDS | _(fill in)_ |
| tRNA genes | _(fill in)_ |
| rRNA genes | _(fill in)_ |
| BLAST DB | UniProt SwissProt |

> 💡 **N50 explained:** Sort all contigs by length (longest first). Walk down the list adding lengths. N50 is the length of the contig where your running total crosses 50% of the total assembly size. A higher N50 means a more contiguous assembly.

---

## 🧠 Concepts Glossary

| Term | Definition |
|---|---|
| **De novo assembly** | Reconstructing a genome from reads alone, without a reference genome |
| **FASTQ** | File format with DNA sequence + per-base Phred quality scores |
| **Phred score (Q)** | Q30 = 99.9% accuracy; Q20 = 99% accuracy |
| **Paired-end reads** | Each DNA fragment sequenced from both ends; improves assembly |
| **Contig** | Contiguous assembled sequence with no gaps |
| **Scaffold** | Contigs ordered and oriented using paired-end distances; gaps = `N`s |
| **k-mer** | Substring of length k; building block of SPAdes de Bruijn graph |
| **de Bruijn graph** | Graph-based data structure SPAdes uses to assemble reads |
| **Coverage** | How many reads cover each base position; higher = more confidence |
| **N50** | Assembly contiguity metric — see Results section above |
| **ORF** | Open Reading Frame — a sequence that could encode a protein |
| **CDS** | Coding Sequence — a confirmed protein-coding region |
| **GFF3** | Standard annotation file format for genome browsers |
| **GenBank (.gbk)** | Combined sequence + annotation format; used for DB submission |
| **HMM** | Hidden Markov Model; used by HMMER to detect protein families |
| **SwissProt** | Manually curated, high-quality protein sequence database (UniProt) |

---

## ⚙️ Tool Versions

| Tool | Version | Purpose |
|---|---|---|
| FastQC | ≥ 0.12 | Quality assessment |
| SPAdes | 4.2.0 | De novo assembly |
| Prokka | ≥ 1.14 | Genome annotation |
| BLAST+ | ≥ 2.12 | Protein similarity search |
| Prodigal | ≥ 2.6 | CDS/ORF prediction (Prokka dep.) |
| HMMER | ≥ 3.3 | Protein domain search (Prokka dep.) |
| Barrnap | ≥ 0.9 | rRNA detection (Prokka dep.) |
| Aragorn | ≥ 1.2 | tRNA detection (Prokka dep.) |

---

## 🔗 Resources

- **NCBI SRA record:** https://www.ncbi.nlm.nih.gov/sra?term=SRX146853
- **ENA record:** https://www.ebi.ac.uk/ena/browser/view/SRR501130
- **SPAdes GitHub:** https://github.com/ablab/spades
- **Prokka GitHub:** https://github.com/tseemann/prokka
- **FastQC docs:** https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **UniProt SwissProt:** https://www.uniprot.org/
- **Bandage (assembly graph viewer):** https://github.com/rrwick/Bandage
- **tbl2asn (if needed):** https://ftp.ncbi.nlm.nih.gov/toolbox/ncbi_tools/converters/by_program/table2asn/

---

*5-day whole genome de novo assembly workshop — Dr. Omics Lab*
