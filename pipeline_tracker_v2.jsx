import { useState } from "react";

const STEPS = [
  {
    id: "01",
    phase: "DATA",
    icon: "⬇",
    title: "Download Raw Reads",
    script: "bash scripts/01_download.sh",
    tool: "wget + gunzip",
    color: "#00b4d8",
    accent: "#0077b6",
    organism: "S. aureus GR1 · SRR501130",
    commands: [
      "wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR501/SRR501130/SRR501130_1.fastq.gz",
      "wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR501/SRR501130/SRR501130_2.fastq.gz",
      "gunzip *.gz",
    ],
    checks: [
      "SRR501130_1.fastq exists in data/raw/",
      "SRR501130_2.fastq exists in data/raw/",
      "Files are non-empty (ls -lh data/raw/)",
      "Used ENA (not SRA) — files already split into _1 and _2",
    ],
    note: "Use -c flag so wget can resume if your connection drops. ENA pre-splits paired-end reads; NCBI SRA merges them into one file.",
    output: "data/raw/SRR501130_{1,2}.fastq",
  },
  {
    id: "02",
    phase: "QC",
    icon: "🔬",
    title: "Quality Control",
    script: "bash scripts/02_fastqc.sh",
    tool: "FastQC",
    color: "#f4a261",
    accent: "#e76f51",
    organism: "76 bp reads · ~33% GC · 10.9M spots",
    commands: [
      "sudo apt-get update && sudo apt-get install fastqc",
      "fastqc *.fastq",
    ],
    checks: [
      "HTML reports open in browser",
      "Per base sequence quality: median Phred ≥ 30 (green zone) ✅",
      "Overrepresented sequences: 'No hit' ✅",
      "Adapter content: absent ✅",
      "Decision: NO TRIMMING needed (all green = proceed to assembly)",
    ],
    note: "For de novo WGS assembly without a reference, skipping trimming is standard — trimming risks discarding real genomic sequence. Our data passed all 3 thresholds.",
    output: "results/fastqc/SRR501130_{1,2}_fastqc.html",
  },
  {
    id: "03",
    phase: "SETUP",
    icon: "⚙",
    title: "Install SPAdes",
    script: "bash scripts/03_install_spades.sh",
    tool: "SPAdes v4.2.0 (from source)",
    color: "#9b5de5",
    accent: "#7209b7",
    organism: "Source compile — not in apt repos",
    commands: [
      "wget -c https://github.com/ablab/spades/archive/refs/tags/v4.2.0.tar.gz",
      "tar -xzf v4.2.0.tar.gz && cd spades-4.2.0/",
      "sudo apt-get install g++ cmake zlib1g-dev build-essential",
      "./spades_compile.sh",
    ],
    checks: [
      "spades_compile.sh runs without errors",
      "./spades.py prints the help manual (verify with: cd bin && ./spades.py)",
      "Note the full path to spades.py for next step",
      "HAMMER error? → add --only-assembler flag in assembly step",
    ],
    note: "SPAdes must be compiled from source. Note your full executable path — you'll need it in the assembly command (e.g. /mnt/f/.../spades-4.2.0/bin/spades.py).",
    output: "tools/spades-4.2.0/bin/spades.py",
  },
  {
    id: "04",
    phase: "ASSEMBLY",
    icon: "🧩",
    title: "De Novo Assembly",
    script: "bash scripts/04_assembly_spades.sh",
    tool: "SPAdes --isolate",
    color: "#2dc653",
    accent: "#007f5f",
    organism: "S. aureus GR1 · pure culture",
    commands: [
      "/path/to/spades.py --isolate -m 250 -1 SRR501130_1.fastq -2 SRR501130_2.fastq -o spades_result",
    ],
    checks: [
      "spades_result/contigs.fasta exists and is non-empty",
      "spades_result/scaffolds.fasta exists",
      "spades_result/assembly_graph.gfa exists (view in Bandage)",
      "No ERROR lines in spades_result/spades.log",
      "Record: number of contigs, largest contig length, N50",
    ],
    note: "Header format: >NODE_1_length_163380_cov_139.527 — length = bp, cov = read coverage. Higher coverage = more confident sequence. contigs.fasta is your primary output — use this for Prokka, not scaffolds.",
    output: "results/assembly/spades_result/contigs.fasta",
  },
  {
    id: "05",
    phase: "SETUP",
    icon: "⚙",
    title: "Install Prokka",
    script: "bash scripts/05_install_prokka.sh",
    tool: "Prokka + deps",
    color: "#fb5607",
    accent: "#c1121f",
    organism: "git clone + --setupdb",
    commands: [
      "sudo apt install -y git bioperl hmmer barrnap aragorn ncbi-blast+ prodigal",
      "sudo cpanm Bio::SearchIO::hmmer3",
      "git clone https://github.com/tseemann/prokka.git",
      "cd prokka && ./bin/prokka --setupdb",
    ],
    checks: [
      "--setupdb completes without errors (indexes databases!)",
      "If tbl2asn error: download manually from NCBI FTP",
      "cd prokka/bin && ./prokka --version prints a version number",
      "Copy (not move) contigs.fasta into prokka/bin/",
    ],
    note: "--setupdb is critical! Without it Prokka runs unindexed searches that take hours. With it, annotation takes 2–4 minutes. tbl2asn error is common — see script for the fix.",
    output: "tools/prokka/bin/prokka (executable)",
  },
  {
    id: "06",
    phase: "ANNOTATION",
    icon: "🏷",
    title: "Genome Annotation",
    script: "bash scripts/06_annotation_prokka.sh",
    tool: "Prokka",
    color: "#ff006e",
    accent: "#8338ec",
    organism: "Prodigal · BLAST · HMMER · Barrnap · Aragorn",
    commands: [
      "cp spades_result/contigs.fasta prokka/bin/",
      "cd prokka/bin",
      "./prokka contigs.fasta --outdir prokka_result",
    ],
    checks: [
      "prokka_result/ directory contains output files",
      ".gff file is non-empty (has gene features)",
      ".txt summary shows CDS count > 0",
      "Record: # CDSs, tRNAs, rRNAs from the .txt summary",
      "Update Results Summary table in README",
    ],
    note: "Prokka runs: Prodigal (CDS) → BLAST (function) → HMMER (domains) → Barrnap (rRNA) → Aragorn (tRNA). .gff for genome browsers, .gbk for Geneious/Artemis, .faa for BLAST/phylogenetics.",
    output: "results/annotation/prokka_result/*.{gff,gbk,faa,ffn,txt,tbl}",
  },
  {
    id: "07",
    phase: "BLAST",
    icon: "💥",
    title: "BLAST vs UniProt",
    script: "bash scripts/07_blast_uniprot.sh",
    tool: "BLASTP + SwissProt",
    color: "#ffbe0b",
    accent: "#fb5607",
    organism: "UniProt SwissProt (manually curated)",
    commands: [
      "wget -c https://ftp.uniprot.org/.../uniprot_sprot.fasta.gz && gunzip",
      "makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprodb",
      'blastp -query longest_orfs.pep -db uniprodb -num_alignments 1 -num_threads 4 -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen" -out result.txt',
    ],
    checks: [
      "uniprodb.p* database files created by makeblastdb",
      "result.txt is non-empty (has BLAST hits)",
      "Calculate % query coverage: (length - mismatch - gaps) / qlen × 100",
      "Look up top hits on https://www.uniprot.org/",
    ],
    note: "SwissProt = manually reviewed, high-confidence proteins. outfmt 6 = tabular output. % Query Coverage = (alignment_length − mismatches − gaps) ÷ qlen × 100. Higher coverage + low e-value = strong functional assignment.",
    output: "results/blast/result.txt",
  },
];

const phaseColors = { DATA: "#00b4d8", QC: "#f4a261", SETUP: "#9b5de5", ASSEMBLY: "#2dc653", ANNOTATION: "#ff006e", BLAST: "#ffbe0b" };

export default function App() {
  const [open, setOpen] = useState(null);
  const [done, setDone] = useState({});
  const [notes, setNotes] = useState({});

  const toggle = (stepId, ci) => {
    const k = `${stepId}-${ci}`;
    setDone(p => ({ ...p, [k]: !p[k] }));
  };

  const stepDone = sid => STEPS.find(s => s.id === sid).checks.filter((_, i) => done[`${sid}-${i}`]).length;
  const totalDone = Object.values(done).filter(Boolean).length;
  const totalAll = STEPS.reduce((a, s) => a + s.checks.length, 0);
  const pct = Math.round((totalDone / totalAll) * 100);

  return (
    <div style={{ minHeight: "100vh", background: "#0d0d0d", fontFamily: "'Courier New', Courier, monospace", color: "#e0e0e0" }}>

      {/* ── Header ── */}
      <div style={{ borderBottom: "1px solid #222", padding: "28px 32px 24px" }}>
        <div style={{ display: "flex", alignItems: "flex-start", justifyContent: "space-between", flexWrap: "wrap", gap: 16 }}>
          <div>
            <div style={{ fontSize: 10, letterSpacing: 4, color: "#555", textTransform: "uppercase", marginBottom: 6 }}>Dr. Omics Lab · 5-Day Workshop</div>
            <h1 style={{ margin: 0, fontSize: 22, fontWeight: 700, color: "#f0f0f0", letterSpacing: -0.5 }}>
              WGS De Novo Assembly
            </h1>
            <div style={{ marginTop: 6, fontSize: 12, color: "#666" }}>
              <span style={{ color: "#2dc653" }}>▶</span> <em>Staphylococcus aureus</em> GR1 · SRR501130 · Illumina PE 76 bp
            </div>
          </div>

          {/* Progress ring */}
          <div style={{ textAlign: "right" }}>
            <div style={{ fontSize: 28, fontWeight: 800, color: pct === 100 ? "#2dc653" : "#f0f0f0", lineHeight: 1 }}>{pct}<span style={{ fontSize: 14 }}>%</span></div>
            <div style={{ fontSize: 10, color: "#555", letterSpacing: 2 }}>{totalDone}/{totalAll} CHECKS</div>
          </div>
        </div>

        {/* Progress bar */}
        <div style={{ marginTop: 20, height: 3, background: "#1a1a1a", borderRadius: 2, overflow: "hidden" }}>
          <div style={{ height: "100%", width: `${pct}%`, background: "linear-gradient(90deg, #00b4d8, #2dc653)", borderRadius: 2, transition: "width 0.5s ease" }} />
        </div>

        {/* Phase legend */}
        <div style={{ marginTop: 14, display: "flex", gap: 16, flexWrap: "wrap" }}>
          {Object.entries(phaseColors).map(([ph, col]) => (
            <span key={ph} style={{ fontSize: 10, color: col, letterSpacing: 2 }}>■ {ph}</span>
          ))}
        </div>
      </div>

      {/* ── Steps ── */}
      <div style={{ padding: "20px 32px 40px" }}>
        {STEPS.map((step, idx) => {
          const isOpen = open === step.id;
          const doneCount = stepDone(step.id);
          const total = step.checks.length;
          const complete = doneCount === total;

          return (
            <div key={step.id} style={{ marginBottom: 8, borderRadius: 4, border: `1px solid ${complete ? step.color + "44" : "#1e1e1e"}`, overflow: "hidden", transition: "border-color 0.3s" }}>

              {/* Step header */}
              <div
                onClick={() => setOpen(isOpen ? null : step.id)}
                style={{ cursor: "pointer", padding: "14px 20px", display: "flex", alignItems: "center", gap: 14, background: isOpen ? "#111" : "#0a0a0a", transition: "background 0.2s" }}
              >
                {/* Step number */}
                <div style={{ width: 36, height: 36, borderRadius: 3, background: complete ? step.color : "#1a1a1a", border: `1px solid ${complete ? step.color : "#2a2a2a"}`, display: "flex", alignItems: "center", justifyContent: "center", fontSize: complete ? 16 : 11, color: complete ? "#000" : "#444", fontWeight: 800, flexShrink: 0, transition: "all 0.3s" }}>
                  {complete ? "✓" : step.id}
                </div>

                {/* Phase badge + title */}
                <div style={{ flex: 1, minWidth: 0 }}>
                  <div style={{ display: "flex", alignItems: "center", gap: 8, marginBottom: 2 }}>
                    <span style={{ fontSize: 9, letterSpacing: 2, color: phaseColors[step.phase], background: phaseColors[step.phase] + "18", padding: "2px 6px", borderRadius: 2 }}>{step.phase}</span>
                    <span style={{ fontSize: 14, fontWeight: 700, color: "#e0e0e0" }}>{step.icon} {step.title}</span>
                  </div>
                  <div style={{ fontSize: 10, color: "#444" }}>{step.organism}</div>
                </div>

                {/* Progress + toggle */}
                <div style={{ display: "flex", alignItems: "center", gap: 12, flexShrink: 0 }}>
                  <div style={{ fontSize: 11, color: doneCount > 0 ? step.color : "#333" }}>{doneCount}/{total}</div>
                  <div style={{ width: 48, height: 3, background: "#1a1a1a", borderRadius: 2 }}>
                    <div style={{ width: `${(doneCount / total) * 100}%`, height: "100%", background: step.color, borderRadius: 2, transition: "width 0.3s" }} />
                  </div>
                  <span style={{ color: "#333", fontSize: 10, transform: isOpen ? "rotate(180deg)" : "none", display: "inline-block", transition: "transform 0.2s" }}>▼</span>
                </div>
              </div>

              {/* Expanded */}
              {isOpen && (
                <div style={{ padding: "0 20px 20px", background: "#0a0a0a", borderTop: "1px solid #1a1a1a" }}>

                  {/* Command block */}
                  <div style={{ marginTop: 16 }}>
                    <div style={{ fontSize: 9, letterSpacing: 3, color: "#444", marginBottom: 8 }}>COMMANDS</div>
                    <div style={{ background: "#050505", border: "1px solid #1a1a1a", borderRadius: 3, padding: "12px 14px" }}>
                      {step.commands.map((cmd, i) => (
                        <div key={i} style={{ marginBottom: i < step.commands.length - 1 ? 6 : 0 }}>
                          <span style={{ color: step.color, marginRight: 6 }}>$</span>
                          <span style={{ fontSize: 11, color: "#aaa", wordBreak: "break-all" }}>{cmd}</span>
                        </div>
                      ))}
                    </div>
                  </div>

                  {/* Checklist */}
                  <div style={{ marginTop: 16 }}>
                    <div style={{ fontSize: 9, letterSpacing: 3, color: "#444", marginBottom: 8 }}>VERIFICATION CHECKLIST</div>
                    {step.checks.map((check, ci) => {
                      const k = `${step.id}-${ci}`;
                      const isDone = done[k];
                      return (
                        <div
                          key={ci}
                          onClick={() => toggle(step.id, ci)}
                          style={{ display: "flex", alignItems: "flex-start", gap: 10, padding: "8px 10px", marginBottom: 4, borderRadius: 3, cursor: "pointer", background: isDone ? step.color + "12" : "#111", border: `1px solid ${isDone ? step.color + "33" : "#1a1a1a"}`, transition: "all 0.2s" }}
                        >
                          <div style={{ width: 16, height: 16, borderRadius: 2, border: `1px solid ${isDone ? step.color : "#2a2a2a"}`, background: isDone ? step.color : "transparent", flexShrink: 0, marginTop: 1, display: "flex", alignItems: "center", justifyContent: "center", fontSize: 10, color: "#000", fontWeight: 800, transition: "all 0.2s" }}>
                            {isDone ? "✓" : ""}
                          </div>
                          <span style={{ fontSize: 12, color: isDone ? "#888" : "#bbb", textDecoration: isDone ? "line-through" : "none", lineHeight: 1.4 }}>{check}</span>
                        </div>
                      );
                    })}
                  </div>

                  {/* Note */}
                  <div style={{ marginTop: 14, background: "#111", border: `1px solid ${step.color}22`, borderLeft: `3px solid ${step.color}`, borderRadius: "0 3px 3px 0", padding: "10px 14px", fontSize: 11, color: "#777", lineHeight: 1.6 }}>
                    <span style={{ color: step.color, marginRight: 6, fontWeight: 700 }}>NOTE</span>{step.note}
                  </div>

                  {/* Output */}
                  <div style={{ marginTop: 12, display: "flex", gap: 8, alignItems: "center", flexWrap: "wrap" }}>
                    <span style={{ fontSize: 9, letterSpacing: 3, color: "#333" }}>OUTPUT</span>
                    <code style={{ fontSize: 10, color: step.color, background: step.color + "11", padding: "3px 8px", borderRadius: 2 }}>{step.output}</code>
                  </div>

                  {/* Notes textarea */}
                  <div style={{ marginTop: 14 }}>
                    <div style={{ fontSize: 9, letterSpacing: 3, color: "#333", marginBottom: 6 }}>MY NOTES</div>
                    <textarea
                      placeholder="e.g. N50=84,231 bp · 52 contigs · 2,821 CDSs predicted..."
                      value={notes[step.id] || ""}
                      onChange={e => setNotes(p => ({ ...p, [step.id]: e.target.value }))}
                      style={{ width: "100%", minHeight: 52, background: "#050505", border: "1px solid #1a1a1a", borderRadius: 3, padding: "8px 10px", fontSize: 11, color: "#888", fontFamily: "inherit", resize: "vertical", outline: "none", boxSizing: "border-box" }}
                    />
                  </div>
                </div>
              )}
            </div>
          );
        })}
      </div>

      {/* Footer */}
      <div style={{ borderTop: "1px solid #111", padding: "16px 32px", display: "flex", justifyContent: "space-between", alignItems: "center", flexWrap: "wrap", gap: 8 }}>
        <span style={{ fontSize: 10, color: "#333", letterSpacing: 2 }}>github.com/your-username/wgs-denovo-assembly</span>
        <span style={{ fontSize: 10, color: "#333" }}>⚠ Update Results Summary in README after completing pipeline</span>
      </div>
    </div>
  );
}
