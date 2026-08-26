#!/usr/bin/env python3
"""
aggregate_codon_usage.py — Combine per-gene codon usage TSVs into one file.

Optional helper script — not required for the core workflow, but useful for
downstream analysis or visualization of codon usage patterns across genes.

Reads all *.codon_usage.tsv files in the input directory and writes a single
aggregated TSV with an extra 'Gene' column.

Usage:
  python3 aggregate_codon_usage.py -i codon_usage/ -o all_genes_codon_usage.tsv
"""
import argparse
import sys
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "-i", "--input-dir",
        required=True,
        type=Path,
        help="Directory containing <gene>.codon_usage.tsv files",
    )
    ap.add_argument(
        "-o", "--output",
        required=True,
        type=Path,
        help="Output aggregated TSV file",
    )
    args = ap.parse_args()

    if not args.input_dir.exists():
        print(f"ERROR: input directory not found: {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    tsv_files = list(args.input_dir.glob("*.codon_usage.tsv"))
    if not tsv_files:
        print(f"ERROR: no *.codon_usage.tsv files found in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    dfs = []
    for tsv_path in sorted(tsv_files):
        gene_name = tsv_path.stem.replace(".codon_usage", "")
        df = pd.read_csv(tsv_path, sep="\t")
        df.insert(0, "Gene", gene_name)
        dfs.append(df)

    combined = pd.concat(dfs, ignore_index=True)
    
    output_path = args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_path, sep="\t", index=False)

    print(f"Aggregated {len(dfs)} gene(s) ({len(combined)} rows) into {output_path}")


if __name__ == "__main__":
    main()
