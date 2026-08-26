#!/usr/bin/env python3
"""
blast_extract_genes.py -- BLASTN-extract the 13 PCGs + 12S/16S rRNA from
whole mitogenome FASTA(s) fetched by fetch_mitogenomes.py, using a reference
gene set from the user's OWN already-annotated data.

Trimmed down from PHYLOGENETIC/scripts/blast_extract_cds.py: reuses only its
core BLASTN-extraction logic (makeblastdb on the target, BLASTN each
reference gene, best-hit strand-correct + in-frame ATG-trim for PCGs).
Dropped: tRNAscan-SE annotation, D-loop detection, GFF3/BED output -- none
of those are consumed anywhere in mitotree-agent, which only needs the 13
PCGs + 2 rRNAs as build_supermatrix.py --fasta input.

The reference gene set does NOT come from a GFF3/GenBank annotation file --
it comes directly from build_supermatrix.py's own genes_from_archive /
genes_from_fasta loaders (a MitoAnnotator archive, or a bare per-gene FASTA
the user already has), or from a prior mitotree-agent run's
build/01_per_gene/<gene>.fasta directory. No annotation-format parsing
needed since the reference is already split per gene.

Usage:
  python3 blast_extract_genes.py \\
      --targets-dir taxon_sampling/raw/ \\
      --reference-archive annotation_mitoannotator/Sample.tar.gz \\
      -o taxon_sampling/genes/

  python3 blast_extract_genes.py \\
      --targets-dir taxon_sampling/raw/ \\
      --reference-per-gene-dir build/01_per_gene/ \\
      -o taxon_sampling/genes/
"""
import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "build" / "mitotree-builder" / "scripts"))
from build_supermatrix import (  # noqa: E402
    DEFAULT_PCG_ORDER, DEFAULT_RRNA_ORDER, genes_from_archive, genes_from_fasta, parse_fasta,
)

PCG_NAMES = set(DEFAULT_PCG_ORDER)
RRNA_NAMES = set(DEFAULT_RRNA_ORDER)
GENE_NAMES = PCG_NAMES | RRNA_NAMES
EVALUE = "1e-20"


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


def read_ref_genes_from_dir(gene_dir: Path) -> dict[str, str]:
    """One <gene>.fasta per file (e.g. build/01_per_gene/), single record each."""
    genes = {}
    for gene in GENE_NAMES:
        fp = gene_dir / f"{gene}.fasta"
        if fp.exists():
            seqs = parse_fasta(str(fp))
            if seqs:
                genes[gene] = next(iter(seqs.values()))
    return genes


def read_target(path: Path) -> tuple[str, str]:
    """Read a single-record whole-genome FASTA. Returns (label, sequence)."""
    seqs = parse_fasta(str(path))
    if not seqs:
        return path.stem, ""
    label = next(iter(seqs.keys()))
    return label, next(iter(seqs.values()))


def blast_extract(target_seq: str, ref_genes: dict[str, str], blastn_cmd: list[str],
                   makeblastdb_cmd: list[str], tmpdir: Path) -> dict[str, str]:
    """For each reference gene, BLASTN against the target genome and keep the
    best hit, strand-corrected and (for PCGs) trimmed to the first in-frame ATG."""
    target_fasta = tmpdir / "target.fasta"
    target_fasta.write_text(">target\n" + target_seq + "\n")
    db = str(tmpdir / "blastdb")
    result = subprocess.run(makeblastdb_cmd + ["-in", str(target_fasta), "-dbtype", "nucl", "-out", db],
                             capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: makeblastdb failed:\n{result.stderr}", file=sys.stderr)
        return {}

    extracted = {}
    for gene, ref_seq in sorted(ref_genes.items()):
        qfile = tmpdir / f"{gene}.fa"
        qfile.write_text(f">{gene}\n{ref_seq}\n")
        result = subprocess.run(
            blastn_cmd + ["-task", "blastn", "-query", str(qfile), "-db", db,
                          "-outfmt", "6 sseqid sstart send pident length qcovs bitscore",
                          "-evalue", EVALUE, "-max_target_seqs", "1"],
            capture_output=True, text=True,
        )
        hits = [h for h in result.stdout.strip().split("\n") if h]
        if not hits:
            print(f"  {gene:10s}  WARNING: no BLAST hit")
            continue

        fields = hits[0].split("\t")
        sstart, send = int(fields[1]), int(fields[2])
        rev = sstart > send
        lo, hi = (send, sstart) if rev else (sstart, send)

        seq = target_seq[lo - 1:hi]
        if rev:
            seq = reverse_complement(seq)

        if gene in PCG_NAMES:
            atg = seq.lower().find("atg")
            if 0 <= atg < 9:
                seq = seq[atg:]

        diff = len(seq) - len(ref_seq)
        flag = "OK" if abs(diff) <= 9 else f"{diff:+d} bp"
        print(f"  {gene:10s}  {lo:7d}-{hi:7d} {'−' if rev else '+'}  {len(seq):5d} bp  [{flag}]")
        extracted[gene] = seq

    return extracted


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--targets-dir", required=True, help="directory of whole-genome FASTA(s), one sample per file (e.g. fetch_mitogenomes.py's -o dir)")
    ref_group = ap.add_mutually_exclusive_group(required=True)
    ref_group.add_argument("--reference-archive", help="MitoAnnotator .tar.gz to source reference gene sequences from")
    ref_group.add_argument("--reference-fasta", help="a single sample's bare per-gene FASTA (multiple gene records in one file)")
    ref_group.add_argument("--reference-per-gene-dir", help="a directory with one <gene>.fasta per gene, e.g. build/01_per_gene/")
    ap.add_argument("-o", "--outdir", required=True, help="output directory for one <label>.genes.fasta per target")
    ap.add_argument("--blastn-cmd", default="blastn")
    ap.add_argument("--makeblastdb-cmd", default="makeblastdb")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.reference_archive:
        ref_genes = genes_from_archive(Path(args.reference_archive), GENE_NAMES)
    elif args.reference_fasta:
        ref_genes = genes_from_fasta(Path(args.reference_fasta), GENE_NAMES)
    else:
        ref_genes = read_ref_genes_from_dir(Path(args.reference_per_gene_dir))

    missing_ref = GENE_NAMES - set(ref_genes)
    print(f"Reference: {len(ref_genes)}/{len(GENE_NAMES)} genes loaded"
          + (f" (missing: {', '.join(sorted(missing_ref))})" if missing_ref else ""))
    if not ref_genes:
        print("ERROR: reference has zero usable genes -- check the reference source.", file=sys.stderr)
        sys.exit(1)

    blastn_cmd = args.blastn_cmd.split()
    makeblastdb_cmd = args.makeblastdb_cmd.split()

    targets = sorted(Path(args.targets_dir).glob("*.fasta"))
    if not targets:
        print(f"ERROR: no *.fasta files in {args.targets_dir}", file=sys.stderr)
        sys.exit(1)

    no_hit_report = []
    for target_path in targets:
        label, target_seq = read_target(target_path)
        if not target_seq:
            print(f"WARNING: {target_path} is empty or unparseable, skipping")
            continue
        print(f"\n{label} ({len(target_seq):,} bp)")
        with tempfile.TemporaryDirectory(prefix="mitotree_blast_") as tmp:
            extracted = blast_extract(target_seq, ref_genes, blastn_cmd, makeblastdb_cmd, Path(tmp))

        missing = GENE_NAMES - set(extracted)
        if missing:
            no_hit_report.append(f"{label}: {', '.join(sorted(missing))}")
        if not extracted:
            print(f"  ALL genes missed for {label} -- likely too divergent from the reference for BLASTN, not an annotation gap.")
            continue

        out_fasta = outdir / f"{label}.genes.fasta"
        with open(out_fasta, "w") as f:
            for gene, seq in extracted.items():
                f.write(f">{gene}\n{seq}\n")
        print(f"  -> {out_fasta}")

    if no_hit_report:
        (outdir / "no_hit_genes.txt").write_text("\n".join(no_hit_report) + "\n")
        print(f"\n{len(no_hit_report)} sample(s) had at least one gene with no BLAST hit -- see {outdir / 'no_hit_genes.txt'}")


if __name__ == "__main__":
    main()
