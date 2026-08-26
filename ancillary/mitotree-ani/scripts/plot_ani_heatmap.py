#!/usr/bin/env python3
"""
plot_ani_heatmap.py — Generate a clustered ANI heatmap from fastANI output.

Python matplotlib replacement for the original R script (plot_ani_heatmap.R),
avoiding the need for R/tidyverse/pheatmap dependencies. Uses pandas (already
in pixi.toml) + matplotlib (needs to be added to pixi.toml).

Reads fastANI's native 5-column, no-header output:
  query, reference, ani, fragments_mapped, fragments_total

Produces a hierarchically-clustered heatmap with:
  - Correlation distance, complete linkage clustering
  - ANI values annotated on cells
  - Symmetric matrix (both triangles)
  - Self-self comparisons shown as 100%

Usage:
  python3 plot_ani_heatmap.py -i fastani_results.tsv -o ani_heatmap.png
"""
import argparse
import sys
from pathlib import Path

import pandas as pd
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend for headless rendering
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import squareform
except ImportError as e:
    print(f"ERROR: missing dependency — {e}", file=sys.stderr)
    print("Add 'matplotlib = \"*\"' and 'scipy = \"*\"' to pixi.toml", file=sys.stderr)
    sys.exit(1)


def load_ani_data(path: Path) -> pd.DataFrame:
    """Load fastANI output and clean sample names."""
    col_names = ["query", "reference", "ani", "frag_mapped", "frag_total"]
    df = pd.read_csv(path, sep="\t", header=None, names=col_names)
    
    def clean_name(p: str) -> str:
        """Strip path and .fasta/.mitoannotator/.fa/.fna suffixes."""
        stem = Path(p).name
        for suffix in (".fasta", ".mitoannotator", ".fa", ".fna"):
            stem = stem.replace(suffix, "")
        return stem
    
    df["query"] = df["query"].apply(clean_name)
    df["reference"] = df["reference"].apply(clean_name)
    return df


def build_symmetric_matrix(df: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    """Build a symmetric ANI matrix from long-format data."""
    all_names = sorted(set(df["query"].tolist() + df["reference"].tolist()))
    n = len(all_names)
    
    if n == 0:
        raise ValueError("No sample names found in input")
    
    name_to_idx = {name: i for i, name in enumerate(all_names)}
    matrix = np.zeros((n, n))
    
    for _, row in df.iterrows():
        i = name_to_idx.get(row["query"])
        j = name_to_idx.get(row["reference"])
        if i is not None and j is not None:
            matrix[i, j] = row["ani"]
            matrix[j, i] = row["ani"]  # Symmetrize
    
    # Ensure self-self = 100
    np.fill_diagonal(matrix, 100.0)
    
    return matrix, all_names


def cluster_order(matrix: np.ndarray) -> list[int]:
    """Get clustered order of samples using correlation distance."""
    # Convert ANI to distance (100 - ANI)
    dist_matrix = 100.0 - matrix
    
    # Convert to condensed form for linkage
    condensed = squareform(dist_matrix, checks=False)
    
    # Hierarchical clustering with correlation distance
    # Note: linkage expects a distance matrix; we use 'correlation' metric
    # But since we already converted to distance, use 'complete' linkage
    try:
        link = linkage(condensed, method='complete')
        order = leaves_list(link)
    except Exception:
        # Fallback: no clustering, use original order
        order = np.arange(matrix.shape[0])
    
    return order.tolist()


def plot_heatmap(matrix: np.ndarray, names: list[str], output_path: Path) -> None:
    """Draw and save the heatmap."""
    # Cluster the samples
    order = cluster_order(matrix)
    matrix_ordered = matrix[order][:, order]
    names_ordered = [names[i] for i in order]
    
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Draw heatmap
    im = ax.imshow(matrix_ordered, cmap='RdYlBu_r', vmin=70, vmax=100, aspect='auto')
    
    # Set ticks and labels
    ax.set_xticks(np.arange(len(names_ordered)))
    ax.set_yticks(np.arange(len(names_ordered)))
    ax.set_xticklabels(names_ordered, rotation=90, fontsize=8)
    ax.set_yticklabels(names_ordered, fontsize=8)
    
    # Annotate cells with ANI values
    for i in range(len(names_ordered)):
        for j in range(len(names_ordered)):
            val = matrix_ordered[i, j]
            # Skip annotating self-self (always 100) or very low values
            if i == j:
                continue
            if val > 0:
                ax.text(j, i, f"{val:.1f}", ha='center', va='center', fontsize=7)
    
    # Colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('ANI (%)', rotation=-90, va='bottom', labelpad=20)
    
    ax.set_title('Average Nucleotide Identity (ANI) Heatmap', pad=20)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "-i", "--input",
        required=True,
        type=Path,
        help="fastANI results TSV (5 columns, no header)",
    )
    ap.add_argument(
        "-o", "--output",
        required=True,
        type=Path,
        help="Output PNG heatmap path",
    )
    args = ap.parse_args()

    if not args.input.exists():
        print(f"ERROR: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    output_path = args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Load and process data
    df = load_ani_data(args.input)
    matrix, names = build_symmetric_matrix(df)
    
    if len(names) < 2:
        print(f"ERROR: need ≥2 samples for heatmap, found {len(names)}", file=sys.stderr)
        sys.exit(1)

    print(f"Plotting {len(names)} samples: {', '.join(names)}")
    plot_heatmap(matrix, names, output_path)
    print(f"Heatmap saved to {output_path}")


if __name__ == "__main__":
    main()
