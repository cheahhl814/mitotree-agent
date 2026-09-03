#!/usr/bin/env python3
"""
core_mito_phylo.py — Core mitogenome phylogenetics pipeline

Extracts PCG (and optionally rRNA) sequences from MitoAnnotator archives,
aligns each gene with MAFFT, trims with trimAl, concatenates with AMAS,
and reports alignment statistics. Produces a supermatrix FASTA + partition
file for IQ-TREE.

Modes
-----
  nt (default)  Nucleotide supermatrix: 13 PCGs + 12S rRNA + 16S rRNA.
                For intragenera / species-level analyses.
  aa            Amino acid supermatrix: 13 PCGs translated with the
                vertebrate mitochondrial code (NCBI table 2). No rRNA.
                For family-level or broader analyses.

Input sources
-------------
  * MitoAnnotator .tar.gz archives (contain seq_genes.fa with annotated genes)
  * Bare per-species gene FASTAs via --fasta LABEL:FILE

Usage
-----
    python3 core_mito_phylo.py <annotation_dir> [<annotation_dir> ...]
        [--fasta LABEL:FILE ...]
        [--mode nt|aa]
        [-o OUTPUT_DIR]
        [--genes GENE [GENE ...]]
        [--mafft PATH] [--trimal PATH] [--amas PATH]

Example
-------
    python3 scripts/core_mito_phylo.py \\
        analyses/MLTree_Betta_waseri_complex/annotation_mitoannotator/ \\
        --fasta Betta_hipposideros:analyses/MLTree_Betta_waseri_complex/Betta_hipposideros_annotation/re_extracted_CDS.fasta \\
        --mode nt \\
        -o analyses/MLTree_Betta_waseri_complex/mitogenome_phylo/
"""

import argparse
import os
import sys
import tarfile
import subprocess
import shutil
from pathlib import Path
from collections import defaultdict

# ─── Default configuration ────────────────────────────────────────────────────

DEFAULT_PCG_ORDER = [
    "ND1", "ND2", "COXI", "COXII", "ATPase8", "ATPase6",
    "COXIII", "ND3", "ND4L", "ND4", "ND5", "ND6", "Cytb",
]

# rRNA gene names as parsed from MitoAnnotator seq_genes.fa headers
# (header ">12S RNA ..." → first word "12S"; same for "16S")
DEFAULT_RRNA_ORDER = ["12S", "16S"]

# nt mode: rRNAs first (mitogenome conventional order), then PCGs
DEFAULT_NT_GENE_ORDER = DEFAULT_RRNA_ORDER + DEFAULT_PCG_ORDER

# aa mode: PCGs only (rRNAs are not translated)
DEFAULT_AA_GENE_ORDER = DEFAULT_PCG_ORDER

# AMAS.py search order (first found is used; AMAS is bundled by mitotree-agent's pixi env
# so this is normally auto-resolved via the env's bin path).
AMAS_SEARCH_PATHS = [
    os.environ.get("AMAS_PY_PATH", ""),  # override via env var if AMAS isn't on PATH
    "/usr/local/bin/AMAS.py",            # common system path
    "AMAS.py",                           # fallback: rely on PATH (resolved by shutil.which below)
]

# ─── Vertebrate mitochondrial genetic code (NCBI table 2) ────────────────────
# Differences from standard code: ATA→M, TGA→W, AGA→*, AGG→*

MITO_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "*", "AGG": "*",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_mito(dna_seq):
    """
    Translate a DNA sequence using the vertebrate mitochondrial code.
    Trims to a multiple of 3, translates, removes any trailing stop codon.
    Ambiguous codons → X.
    """
    dna = dna_seq.upper().replace("-", "")
    remainder = len(dna) % 3
    if remainder:
        dna = dna[:-remainder]
    aa = [MITO_CODON_TABLE.get(dna[i:i+3], "X") for i in range(0, len(dna), 3)]
    if aa and aa[-1] == "*":
        aa.pop()
    return "".join(aa)


# ─── FASTA I/O ────────────────────────────────────────────────────────────────

def parse_fasta(source, is_text=False):
    """
    Parse FASTA from a file path or text string.
    Returns {first_word_of_header: sequence}.
    """
    seqs = {}
    current_id = None
    current_seq = []

    if is_text:
        lines = source.splitlines()
    else:
        with open(source) as f:
            lines = f.readlines()

    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                seqs[current_id] = "".join(current_seq)
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(current_seq)
    return seqs


def write_fasta(path, records):
    with open(path, "w") as f:
        for seq_id, seq in records.items():
            f.write(f">{seq_id}\n{seq}\n")


# ─── Gene extraction ─────────────────────────────────────────────────────────

def genes_from_archive(archive_path, gene_set):
    """Extract PCG sequences from a MitoAnnotator .tar.gz archive."""
    with tarfile.open(archive_path) as tar:
        for member in tar.getmembers():
            if member.name.endswith("seq_genes.fa"):
                f = tar.extractfile(member)
                raw = parse_fasta(f.read().decode(), is_text=True)
                return {k: v for k, v in raw.items() if k in gene_set}
    raise FileNotFoundError(f"seq_genes.fa not found inside {archive_path.name}")


def genes_from_fasta(fasta_path, gene_set):
    """Extract PCG sequences from a bare per-species gene FASTA."""
    raw = parse_fasta(fasta_path)
    return {k: v for k, v in raw.items() if k in gene_set}


# ─── Source discovery ─────────────────────────────────────────────────────────

def discover_archives(dirs):
    """
    Find all .tar.gz files in dirs.
    Returns list of ("archive", label, path) tuples.
    Label = filename without .tar.gz extension.
    """
    sources = []
    for d in dirs:
        d = Path(d)
        if not d.exists():
            print(f"  Warning: {d} not found, skipping", file=sys.stderr)
            continue
        for archive in sorted(d.glob("*.tar.gz")):
            label = archive.name[:-7]  # strip ".tar.gz"
            sources.append(("archive", label, archive))
    return sources


def parse_fasta_specs(specs):
    """
    Parse LABEL:PATH strings from --fasta arguments.
    Returns list of ("fasta", label, path) tuples.
    """
    sources = []
    for spec in (specs or []):
        if ":" not in spec:
            sys.exit(f"--fasta argument must be LABEL:FILE, got: {spec!r}")
        label, path = spec.split(":", 1)
        sources.append(("fasta", label, Path(path)))
    return sources


# ─── Tool helpers ─────────────────────────────────────────────────────────────

def find_amas(user_path):
    if user_path:
        return user_path
    for p in AMAS_SEARCH_PATHS:
        if Path(p).exists():
            return p
    found = shutil.which("AMAS.py")
    if found:
        return found
    sys.exit(
        "AMAS.py not found. Add its path to AMAS_SEARCH_PATHS in this script "
        "or pass --amas PATH."
    )


def amas_python(amas_py):
    """Return the Python interpreter co-located with AMAS.py (pixi env)."""
    py = Path(amas_py).parent / "python3"
    return str(py) if py.exists() else "python3"


def run(cmd, label):
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"\n  ERROR in {label}:", file=sys.stderr)
        print(result.stderr[-1000:], file=sys.stderr)
        sys.exit(1)
    return result


# ─── Pipeline steps ───────────────────────────────────────────────────────────

def step_extract(all_sources, gene_order, out_dir):
    """Step 1: extract per-gene sequences from all sources."""
    gene_set  = set(gene_order)
    gene_seqs = defaultdict(dict)   # gene → {label: sequence}

    for kind, label, path in all_sources:
        print(f"    {label}")
        if kind == "archive":
            genes = genes_from_archive(path, gene_set)
        else:
            genes = genes_from_fasta(path, gene_set)
        for gene, seq in genes.items():
            gene_seqs[gene][label] = seq

    n = len(all_sources)
    print(f"\n  Coverage across {n} samples:")
    for gene in gene_order:
        found = len(gene_seqs.get(gene, {}))
        flag  = "" if found == n else f"  ← MISSING IN {n - found} SAMPLE(S)"
        print(f"    {gene:10s}  {found}/{n}{flag}")

    out_dir.mkdir(parents=True, exist_ok=True)
    gene_fastas = {}
    for gene in gene_order:
        if gene not in gene_seqs:
            print(f"  Warning: {gene} absent from all samples — skipping")
            continue
        fp = out_dir / f"{gene}.fasta"
        write_fasta(fp, gene_seqs[gene])
        gene_fastas[gene] = fp
    return gene_fastas


def step_translate(gene_fastas, out_dir):
    """Translate PCG nucleotide FASTAs to amino acid using vertebrate mito code."""
    out_dir.mkdir(parents=True, exist_ok=True)
    translated = {}
    for gene, fasta in gene_fastas.items():
        out     = out_dir / f"{gene}.aa.fasta"
        nt_seqs = parse_fasta(fasta)
        aa_seqs = {label: translate_mito(seq) for label, seq in nt_seqs.items()}
        write_fasta(out, aa_seqs)
        translated[gene] = out
        print(f"    {gene}")
    return translated


def step_align(gene_fastas, out_dir, mafft_argv):
    """Step 2: align each gene with MAFFT --auto."""
    out_dir.mkdir(parents=True, exist_ok=True)
    aligned = {}
    for gene, fasta in gene_fastas.items():
        out    = out_dir / f"{gene}.aln.fasta"
        result = run(mafft_argv + ["--auto", "--quiet", str(fasta)], f"mafft {gene}")
        out.write_text(result.stdout)
        aligned[gene] = out
        print(f"    {gene}")
    return aligned


def step_trim(aligned, out_dir, trimal_argv):
    """Step 3: trim each alignment with trimAl -automated1."""
    out_dir.mkdir(parents=True, exist_ok=True)
    trimmed = {}
    for gene, aln in aligned.items():
        out = out_dir / f"{gene}.trim.fasta"
        run(trimal_argv + ["-in", str(aln), "-out", str(out), "-automated1"],
            f"trimal {gene}")
        trimmed[gene] = out
        print(f"    {gene}")
    return trimmed


def step_concatenate(trimmed, gene_order, out_dir, amas_py, data_type="dna"):
    """Step 4: concatenate trimmed alignments into a supermatrix with AMAS."""
    out_dir.mkdir(parents=True, exist_ok=True)
    infiles     = [str(trimmed[g]) for g in gene_order if g in trimmed]
    supermatrix = out_dir / "supermatrix.fasta"
    partitions  = out_dir / "partitions.nex"
    run(
        [amas_python(amas_py), amas_py, "concat",
         "-f", "fasta", "-d", data_type,
         "-i"] + infiles +
        ["-t", str(supermatrix), "-p", str(partitions)],
        "AMAS concat",
    )
    
    # Convert AMAS RAxML-style partitions to NEXUS format for IQ-TREE compatibility
    if partitions.exists():
        lines = partitions.read_text().splitlines()
        nexus_lines = ["#NEXUS", "begin sets;"]
        for line in lines:
            if '=' in line:
                # Remove spaces and split by '='
                parts = line.replace(" ", "").split('=')
                name = parts[0].replace("DNA,", "").replace("AA,", "")
                coords = parts[1]
                nexus_lines.append(f"    charset {name} = {coords};")
        nexus_lines.append("end;")
        partitions.write_text("\n".join(nexus_lines))
        
    return supermatrix, partitions


def step_stats(trimmed, gene_order, supermatrix, out_dir, amas_py, data_type="dna"):
    """Step 5: per-gene and supermatrix statistics."""
    out_dir.mkdir(parents=True, exist_ok=True)

    infiles = [str(trimmed[g]) for g in gene_order if g in trimmed]

    py = amas_python(amas_py)
    run(
        [py, amas_py, "summary",
         "-f", "fasta", "-d", data_type,
         "-i"] + infiles +
        ["-o", str(out_dir / "per_gene_stats.txt")],
        "AMAS per-gene summary",
    )
    run(
        [py, amas_py, "summary",
         "-f", "fasta", "-d", data_type,
         "-i", str(supermatrix),
         "-o", str(out_dir / "supermatrix_stats.txt")],
        "AMAS supermatrix summary",
    )

    seqs = parse_fasta(supermatrix)
    n_taxa = len(seqs)
    sm_len = len(next(iter(seqs.values()))) if seqs else 0
    missing = sum(
        c in "-?N" for s in seqs.values() for c in s.upper()
    ) / max(n_taxa * sm_len, 1) * 100

    print(f"\n  Supermatrix : {n_taxa} taxa × {sm_len} bp")
    print(f"  Missing data: {missing:.1f}%")
    print(f"  Per-gene stats : {out_dir / 'per_gene_stats.txt'}")
    print(f"  Supermatrix stats: {out_dir / 'supermatrix_stats.txt'}")


# ─── Entry point ─────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "source_dirs", nargs="*", default=[],
        help="Directory/ies containing MitoAnnotator .tar.gz archives (optional if --fasta covers all samples)",
    )
    ap.add_argument(
        "--fasta", metavar="LABEL:FILE", action="append",
        help="Additional per-species gene FASTA (may be repeated)",
    )
    ap.add_argument(
        "--mode", choices=["nt", "aa"], default="nt",
        help=(
            "nt  (default) nucleotide supermatrix: 13 PCGs + 12S + 16S rRNA. "
            "aa  amino acid supermatrix: 13 PCGs translated with vertebrate "
            "mitochondrial code (NCBI table 2), no rRNA."
        ),
    )
    ap.add_argument(
        "-o", "--outdir", default="mitogenome_phylo",
        help="Output directory (default: mitogenome_phylo)",
    )
    ap.add_argument(
        "--genes", nargs="+", default=None,
        help=(
            "Override gene list and order. "
            "Default for --mode nt: 12S 16S + 13 PCGs. "
            "Default for --mode aa: 13 PCGs only."
        ),
    )
    ap.add_argument("--mafft",  default="mafft",  help='mafft command, e.g. "pixi run mafft" (default: mafft)')
    ap.add_argument("--trimal", default="trimal", help='trimal command, e.g. "pixi run trimal" (default: trimal)')
    ap.add_argument("--amas",   default=None,     help="Path to AMAS.py")
    args = ap.parse_args()

    # Resolve gene order and AMAS data type from mode
    if args.genes:
        gene_order = args.genes
    elif args.mode == "aa":
        gene_order = DEFAULT_AA_GENE_ORDER
    else:
        gene_order = DEFAULT_NT_GENE_ORDER
    data_type = "aa" if args.mode == "aa" else "dna"

    amas_py = find_amas(args.amas)
    out     = Path(args.outdir)

    print(f"=== Mitogenome Phylogenetics Pipeline  [mode: {args.mode}] ===\n")

    print("Step 1  Discovering annotation sources")
    all_sources = discover_archives(args.source_dirs) + parse_fasta_specs(args.fasta)
    if not all_sources:
        sys.exit("No annotation sources found. Check input directories.")
    print(f"  {len(all_sources)} sample(s)\n")

    print("Step 2  Extracting sequences")
    gene_fastas = step_extract(all_sources, gene_order, out / "01_per_gene")

    if args.mode == "aa":
        print("\nStep 2b  Translating PCGs (vertebrate mito code)")
        gene_fastas = step_translate(gene_fastas, out / "01b_translated")

    print("\nStep 3  Aligning with MAFFT")
    aligned = step_align(gene_fastas, out / "02_aligned", args.mafft.split())

    print("\nStep 4  Trimming with trimAl")
    trimmed = step_trim(aligned, out / "03_trimmed", args.trimal.split())

    print("\nStep 5  Concatenating with AMAS")
    supermatrix, partitions = step_concatenate(
        trimmed, gene_order, out / "04_concatenated", amas_py, data_type
    )

    print("\nStep 6  Computing statistics")
    step_stats(trimmed, gene_order, supermatrix, out / "05_stats", amas_py, data_type)

    model_flag = "-m LG+G -B 1000" if args.mode == "aa" else "-m GTR+G -B 1000"
    print(f"\n=== Done ===")
    print(f"  Supermatrix : {supermatrix}")
    print(f"  Partitions  : {partitions}")
    print(
        f"\nNext step:\n"
        f"  iqtree -s {supermatrix} -p {partitions} {model_flag} -T AUTO"
    )


if __name__ == "__main__":
    main()
