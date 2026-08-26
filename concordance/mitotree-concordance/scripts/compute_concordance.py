#!/usr/bin/env python3
"""
compute_concordance.py — Gene-tree concordance for mitotree-agent.

Computes Robinson-Foulds distances between each gene tree and the species tree,
normalizing per gene (normalized RF = raw_rf / (2*(n_tips-3))), and writes
rf_distances.tsv in the exact format qc/mitotree-qc/scripts/qc_report.py expects.

Usage:
  python3 compute_concordance.py \
      --species-tree <path> \
      --gene-trees <path> \
      --output-rf <path> \
      [--gotree-cmd "pixi run gotree"]

The gene-trees file should contain one Newick tree per line (as produced by
concatenating IQ-TREE .treefile outputs). Gene names are derived from a
parallel og_names.txt file if present, otherwise inferred from tree order.
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path


def run(argv: list[str]) -> subprocess.CompletedProcess:
    """Run a subprocess and return the result."""
    return subprocess.run(argv, capture_output=True, text=True)


def list_tips(gotree_argv: list[str], tree_path: str) -> set[str]:
    """Get the set of tip names from a tree using gotree stats tips."""
    result = run(gotree_argv + ["stats", "tips", "-i", tree_path])
    if result.returncode != 0:
        print(f"ERROR: gotree stats tips failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)
    lines = result.stdout.strip().splitlines()[1:]  # drop header
    # Format: tree<TAB>id<TAB>nneigh<TAB>name<TAB>ExternalBranch<TAB>RootToTip
    return {line.split("\t")[3] for line in lines if line.strip()}


def prune_species_tree(
    gotree_argv: list[str],
    species_tree: str,
    gene_tree: str,
    out_path: str,
) -> None:
    """Prune the species tree to match the gene tree's tip set."""
    # gotree prune -c <compared_tree> -i <reference> -o <output>
    # This removes tips from the reference that are not in the compared tree
    result = run(
        gotree_argv
        + ["prune", "-c", gene_tree, "-i", species_tree, "-o", out_path]
    )
    if result.returncode != 0:
        print(f"ERROR: gotree prune failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)


def compute_rf(
    gotree_argv: list[str],
    species_tree: str,
    gene_tree: str,
) -> int:
    """Compute raw RF distance between two trees using gotree compare trees --rf."""
    # gotree compare trees -i <ref> -c <compared> --rf
    # Output: one line per compared tree with just the RF distance value
    result = run(
        gotree_argv
        + ["compare", "trees", "-i", species_tree, "-c", gene_tree, "--rf"]
    )
    if result.returncode != 0:
        print(f"ERROR: gotree compare trees failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)
    # Parse output: just the RF value on each line (e.g., "0" or "4")
    lines = result.stdout.strip().splitlines()
    if not lines:
        print("ERROR: gotree compare trees produced no output", file=sys.stderr)
        sys.exit(1)
    # Take the RF value from the first (and only) comparison line
    return int(lines[0].strip())


def read_gene_trees(gene_trees_path: str) -> list[str]:
    """Read gene trees from a file (one Newick per line)."""
    with open(gene_trees_path) as f:
        return [line.strip() for line in f if line.strip()]


def read_og_names(og_names_path: str) -> list[str]:
    """Read orthogroup/gene names from a file (one per line)."""
    with open(og_names_path) as f:
        return [line.strip() for line in f if line.strip()]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--species-tree",
        required=True,
        help="Species tree Newick file (from build/mitotree-builder)",
    )
    ap.add_argument(
        "--gene-trees",
        required=True,
        help="Gene trees file (one Newick per line, from IQ-TREE)",
    )
    ap.add_argument(
        "--output-rf",
        required=True,
        help="Output TSV path for RF distances",
    )
    ap.add_argument(
        "--gotree-cmd",
        default="gotree",
        help='e.g. "pixi run gotree" (default: gotree)',
    )
    ap.add_argument(
        "--og-names",
        default=None,
        help="Optional: file with gene names (one per line, same order as gene trees)",
    )
    args = ap.parse_args()

    gotree = args.gotree_cmd.split()
    output_rf = Path(args.output_rf)
    output_rf.parent.mkdir(parents=True, exist_ok=True)

    # Read gene trees
    gene_trees = read_gene_trees(args.gene_trees)
    if not gene_trees:
        print(f"ERROR: no gene trees found in {args.gene_trees}", file=sys.stderr)
        sys.exit(1)

    # Read gene names if provided
    gene_names = None
    if args.og_names and Path(args.og_names).exists():
        gene_names = read_og_names(args.og_names)
        if len(gene_names) != len(gene_trees):
            print(
                f"WARNING: og_names ({len(gene_names)}) and gene_trees ({len(gene_trees)}) "
                "have different lengths — using index-based names",
                file=sys.stderr,
            )
            gene_names = None

    # Write temporary gene trees (one file per gene for gotree comparison)
    # We need to write each gene tree to a temp file for gotree to read
    results = []
    for i, gt_newick in enumerate(gene_trees):
        gene_name = gene_names[i] if gene_names else f"gene_{i}"

        # Write gene tree to temp file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".nwk", delete=False
        ) as gt_f:
            gt_f.write(gt_newick + "\n")
            gt_path = gt_f.name

        try:
            # Get gene tree tips
            gt_tips = list_tips(gotree, gt_path)
            n_tips = len(gt_tips)

            if n_tips < 4:
                # Can't compute meaningful RF with <4 tips
                results.append((gene_name, n_tips, None, None))
                continue

            # Prune species tree to match gene tree tips
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".pruned.nwk", delete=False
            ) as sp_f:
                sp_pruned_path = sp_f.name

            prune_species_tree(gotree, args.species_tree, gt_path, sp_pruned_path)

            # Get pruned species tree tips (should match gene tree tips)
            sp_pruned_tips = list_tips(gotree, sp_pruned_path)
            n_shared = len(sp_pruned_tips)

            # Compute raw RF
            raw_rf = compute_rf(gotree, sp_pruned_path, gt_path)

            # Normalize: normalized_rf = raw_rf / (2 * (n_tips - 3))
            # for unrooted binary trees
            max_rf = 2 * (n_shared - 3)
            normalized_rf = raw_rf / max_rf if max_rf > 0 else None

            results.append((gene_name, n_shared, raw_rf, normalized_rf))

            # Clean up temp files
            Path(gt_path).unlink(missing_ok=True)
            Path(sp_pruned_path).unlink(missing_ok=True)

        except Exception as e:
            print(
                f"WARNING: failed to process gene {gene_name}: {e}",
                file=sys.stderr,
            )
            results.append((gene_name, None, None, None))

    # Write output TSV
    with open(output_rf, "w") as f:
        f.write("orthogroup\tn_shared_taxa\trf_distance\tnormalised_rf\n")
        for gene_name, n_shared, raw_rf, norm_rf in results:
            n_shared_str = str(n_shared) if n_shared is not None else "NA"
            raw_rf_str = str(raw_rf) if raw_rf is not None else "NA"
            norm_rf_str = f"{norm_rf:.6f}" if norm_rf is not None else "NA"
            f.write(f"{gene_name}\t{n_shared_str}\t{raw_rf_str}\t{norm_rf_str}\n")

    print(f"RF distances written to {output_rf}")
    print(f"Processed {len(results)} gene trees")


if __name__ == "__main__":
    main()
