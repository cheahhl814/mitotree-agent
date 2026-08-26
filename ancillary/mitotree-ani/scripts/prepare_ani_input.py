#!/usr/bin/env python3
"""
prepare_ani_input.py — Build a genome list file for fastANI from --genome specs.

Mirrors the --fasta LABEL:FILE argparse pattern from build_supermatrix.py for consistency.
Outputs a plain text file with one FASTA path per line, suitable for fastANI --ql/--rl.

Usage:
  python3 prepare_ani_input.py --genome Sample1:mito1.fasta --genome Sample2:mito2.fasta -o genome_list.txt
"""
import argparse
import sys
from pathlib import Path


def parse_genome_spec(spec: str) -> tuple[str, str]:
    """Parse a LABEL:FILE spec into (label, path)."""
    if ":" not in spec:
        raise ValueError(f"Invalid genome spec '{spec}': expected LABEL:FILE format")
    label, path = spec.split(":", 1)
    if not label or not path:
        raise ValueError(f"Invalid genome spec '{spec}': LABEL and FILE must be non-empty")
    return label, path


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--genome",
        action="append",
        required=True,
        metavar="LABEL:FILE",
        help="Whole mitogenome FASTA spec (repeatable), e.g. 'Sample1:mito1.fasta'",
    )
    ap.add_argument(
        "-o", "--output",
        required=True,
        help="Output genome list file (one FASTA path per line)",
    )
    args = ap.parse_args()

    genome_paths = []
    for spec in args.genome:
        try:
            label, path = parse_genome_spec(spec)
        except ValueError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            sys.exit(1)
        
        if not Path(path).exists():
            print(f"WARNING: genome file '{path}' not found (label: {label})", file=sys.stderr)
        
        genome_paths.append(path)

    if len(genome_paths) < 2:
        print("ERROR: at least 2 --genome specs are required for pairwise ANI", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(genome_paths) + "\n")

    print(f"Wrote {len(genome_paths)} genome path(s) to {output_path}")


if __name__ == "__main__":
    main()
