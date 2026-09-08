#!/usr/bin/env python3
"""
qc_report.py — Post-run QC summary for mitotree-agent mitogenome phylogenomics.

Reads $RUN_DIR/{build,concordance,ani}/ outputs and writes:
  <results_dir>/run_summary.json    machine-readable
  <results_dir>/run_summary.txt     human-readable plain text

Concordance and ANI sections are optional (mitotree-concordance / mitotree-ani
sub-skills) — if their output dirs are absent, that section reports "not found"
rather than failing the whole QC run.

Usage:
  python3 qc_report.py [results_dir]   # default: analysis_results

Exit codes:
  0  all QC checks passed
  1  one or more QC flags raised (review run_summary.txt)
  2  argument / path error
"""

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ── helpers ───────────────────────────────────────────────────────────────────

def _parse_newick_support(nwk: str) -> list:
    """Extract numeric internal-node support values from a Newick string.

    Matches patterns like )100: or )95.2: — skips tips and root labels.
    """
    return [float(m) for m in re.findall(r'\)(\d+(?:\.\d+)?)[:\)]', nwk)]


def render_ascii_tree(tree_path: Path, width: int = 120) -> str:
    """Return a gotree ASCII rendering of a Newick file, or a fallback message."""
    if not shutil.which("gotree"):
        return "  (gotree not found on PATH — install via pixi to enable ASCII tree rendering)"
    if not tree_path.exists():
        return f"  (tree file not found: {tree_path})"
    try:
        result = subprocess.run(
            [
                "gotree", "draw", "text",
                "--with-node-labels",
                "--with-branch-support",
                "--no-branch-lengths",
                "-w", str(width),
                "-i", str(tree_path),
            ],
            capture_output=True, text=True, timeout=30,
        )
        return result.stdout.rstrip()
    except Exception as exc:
        return f"  (gotree error: {exc})"


def _jsonify(obj):
    """Convert numpy scalars to Python native types for JSON serialisation."""
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.ndarray,)):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj)} is not JSON serialisable")


# ── section parsers ───────────────────────────────────────────────────────────

def parse_supermatrix(results_dir: Path) -> dict:
    stats_dir = results_dir / "build" / "05_stats"
    out = {}

    sm_file = stats_dir / "supermatrix_stats.txt"
    if sm_file.exists():
        df = pd.read_csv(sm_file, sep="\t")
        r = df.iloc[0]
        out["n_taxa"] = int(r["No_of_taxa"])
        out["alignment_length_bp"] = int(r["Alignment_length"])
        out["missing_percent"] = float(r["Missing_percent"])
        out["variable_sites"] = int(r["No_variable_sites"])
        out["parsimony_informative_sites"] = int(r["Parsimony_informative_sites"])
    else:
        out["error"] = f"not found: {sm_file}"

    pg_file = stats_dir / "per_gene_stats.txt"
    if pg_file.exists():
        df = pd.read_csv(pg_file, sep="\t")
        df["gene"] = df["Alignment_name"].str.replace(r"\.trim\.fasta$", "", regex=True)
        out["n_genes"] = len(df)
        out["genes"] = [
            {
                "gene": row["gene"],
                "length_bp": int(row["Alignment_length"]),
                "missing_percent": float(row["Missing_percent"]),
                "parsimony_informative_sites": int(row["Parsimony_informative_sites"]),
            }
            for _, row in df.iterrows()
        ]
    return out


def parse_species_tree(results_dir: Path) -> dict:
    tree_file = results_dir / "build" / "species_tree.nwk"
    if not tree_file.exists():
        return {"error": f"not found: {tree_file}"}

    nwk = tree_file.read_text().strip()
    vals = _parse_newick_support(nwk)
    if not vals:
        return {"error": "no support values found in species_tree.nwk"}

    # Normalize support scale: fraction-encoded trees carry 0–1 support
    # fractions, percent-encoded trees (IQ-TREE .treefile/.contree — this
    # pipeline's normal input) carry 0–100. Detect the scale (max <= 1 =>
    # fractions) and rescale to percent before thresholding/formatting.
    if max(vals) <= 1.0:
        vals = [v * 100.0 for v in vals]

    below_95 = [v for v in vals if v < 95]
    below_70 = [v for v in vals if v < 70]
    return {
        "n_internal_nodes_with_support": len(vals),
        "support_min": float(min(vals)),
        "support_median": float(np.median(vals)),
        "support_mean": float(np.mean(vals)),
        "nodes_below_95": sorted(below_95),
        "nodes_below_70": sorted(below_70),
        "flag_low_support": len(below_95) > 0,
    }


def parse_rf_distances(results_dir: Path) -> dict:
    rf_file = results_dir / "concordance" / "rf_distances.tsv"
    if not rf_file.exists():
        return {"error": f"not found: {rf_file}"}

    df = pd.read_csv(rf_file, sep="\t")
    gene_col = df.columns[0]   # 'orthogroup'
    nrf_col = "normalised_rf"

    vals = df[nrf_col].values.astype(float)
    q25, q75 = np.percentile(vals, [25, 75])
    threshold = float(np.median(vals) + 2 * (q75 - q25))

    outlier_rows = df[df[nrf_col] > threshold]
    outliers = [
        {"gene": str(r[gene_col]), "normalised_rf": float(r[nrf_col])}
        for _, r in outlier_rows.iterrows()
    ]
    return {
        "n_genes": len(df),
        "rf_mean": float(np.mean(vals)),
        "rf_max": float(np.max(vals)),
        "outlier_threshold": threshold,
        "outlier_genes": outliers,
        "flag_rf_outliers": len(outliers) > 0,
    }


def parse_concordance(results_dir: Path) -> dict:
    cf_file = results_dir / "concordance" / "concordance.cf.stat"
    if not cf_file.exists():
        return {"error": f"not found: {cf_file}"}

    df = pd.read_csv(cf_file, sep="\t", comment="#")
    for col in ("gCF", "sCF", "Label"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["gCF"])   # root row has all NA

    # misleading: UFBoot ≥95 but gCF <50 or sCF <33
    misleading = df[(df["Label"] >= 95) & ((df["gCF"] < 50) | (df["sCF"] < 33))]
    low_gcf = df[df["gCF"] < 50]

    def _node_row(r):
        return {
            "branch_id": int(r["ID"]),
            "ufboot": None if pd.isna(r["Label"]) else float(r["Label"]),
            "gcf": float(r["gCF"]),
            "scf": float(r["sCF"]),
        }

    return {
        "n_internal_branches": len(df),
        "gcf_mean": float(df["gCF"].mean()),
        "gcf_min": float(df["gCF"].min()),
        "scf_mean": float(df["sCF"].mean()),
        "scf_min": float(df["sCF"].min()),
        "misleading_nodes": [_node_row(r) for _, r in misleading.iterrows()],
        "low_gcf_nodes": [_node_row(r) for _, r in low_gcf.iterrows()],
        "flag_misleading_concordance": len(misleading) > 0,
    }


def parse_ani(results_dir: Path) -> dict:
    ani_file = results_dir / "ani" / "fastani_results.tsv"
    if not ani_file.exists():
        return {"error": f"not found: {ani_file}"}

    df = pd.read_csv(
        ani_file, sep="\t", header=None,
        names=["query", "reference", "ani", "frag_mapped", "frag_total"],
    )

    def _clean(p: str) -> str:
        stem = Path(p).name
        for suffix in (".fasta", ".mitoannotator", ".fa", ".fna"):
            stem = stem.replace(suffix, "")
        return stem

    df["q"] = df["query"].apply(_clean)
    df["r"] = df["reference"].apply(_clean)
    cross = df[df["q"] != df["r"]].copy()

    near_identical = cross[cross["ani"] > 99][["q", "r", "ani"]]
    return {
        "n_pairwise_comparisons": len(cross),
        "ani_min": float(cross["ani"].min()),
        "ani_max": float(cross["ani"].max()),
        "ani_mean": float(cross["ani"].mean()),
        "near_identical_pairs": [
            {"taxon_a": r["q"], "taxon_b": r["r"], "ani": float(r["ani"])}
            for _, r in near_identical.iterrows()
        ],
        "flag_near_identical": len(near_identical) > 0,
    }


# ── text report ───────────────────────────────────────────────────────────────

def _flag(val) -> str:
    return "WARN" if val else "OK"


def write_text_report(summary: dict, out_path: Path) -> None:
    lines = []
    a = lines.append

    a("=" * 62)
    a("MitoTree Run Summary")
    a("=" * 62)

    # Supermatrix
    a("\n── Supermatrix ─────────────────────────────────────────────")
    sm = summary.get("supermatrix", {})
    if "error" in sm:
        a(f"  ERROR: {sm['error']}")
    else:
        a(f"  Taxa              : {sm.get('n_taxa', '?')}")
        a(f"  Genes             : {sm.get('n_genes', '?')}")
        a(f"  Total length      : {sm.get('alignment_length_bp', '?'):,} bp")
        a(f"  Missing data      : {sm.get('missing_percent', '?'):.2f}%")
        a(f"  Variable sites    : {sm.get('variable_sites', '?'):,}")
        a(f"  Parsimony-inform. : {sm.get('parsimony_informative_sites', '?'):,}")
        if "genes" in sm:
            a(f"\n  Per-gene summary (trimmed alignments):")
            for g in sm["genes"]:
                a(f"    {g['gene']:<12}  {g['length_bp']:>5} bp  "
                  f"PI={g['parsimony_informative_sites']:>4}  "
                  f"missing={g['missing_percent']:.1f}%")

    # Species tree
    a("\n── Species Tree Bootstrap ──────────────────────────────────")
    st = summary.get("species_tree", {})
    if "error" in st:
        a(f"  ERROR: {st['error']}")
    else:
        a(f"  Internal nodes    : {st['n_internal_nodes_with_support']}")
        a(f"  UFBoot min/median : {st['support_min']:.0f} / {st['support_median']:.0f}")
        n_below = len(st['nodes_below_95'])
        a(f"  Nodes < 95 UFBoot : {n_below}  [{_flag(st['flag_low_support'])}]")
        if st["nodes_below_95"]:
            a(f"    Values          : {sorted(st['nodes_below_95'])}")
        if st["nodes_below_70"]:
            a(f"  Nodes < 70 UFBoot : {len(st['nodes_below_70'])}  [WARN]")
            a(f"    Values          : {sorted(st['nodes_below_70'])}")

    # RF distances
    a("\n── Gene Tree RF Distances ──────────────────────────────────")
    rf = summary.get("rf_distances", {})
    if "error" in rf:
        a(f"  ERROR: {rf['error']}")
    else:
        a(f"  Genes             : {rf['n_genes']}")
        a(f"  Mean norm. RF     : {rf['rf_mean']:.4f}")
        a(f"  Max norm. RF      : {rf['rf_max']:.4f}")
        a(f"  Outlier threshold : {rf['outlier_threshold']:.4f}  "
          f"(median + 2×IQR)")
        outliers = rf["outlier_genes"]
        if outliers:
            genes_str = ", ".join(
                f"{g['gene']} ({g['normalised_rf']:.3f})" for g in outliers
            )
            a(f"  Outlier genes     : {genes_str}  [WARN]")
        else:
            a(f"  Outlier genes     : none  [OK]")

    # Concordance factors
    a("\n── Concordance Factors ─────────────────────────────────────")
    cf = summary.get("concordance", {})
    if "error" in cf:
        a(f"  ERROR: {cf['error']}")
    else:
        a(f"  Internal branches : {cf['n_internal_branches']}")
        a(f"  Mean gCF / sCF    : {cf['gcf_mean']:.1f}% / {cf['scf_mean']:.1f}%")
        a(f"  Min  gCF / sCF    : {cf['gcf_min']:.1f}% / {cf['scf_min']:.1f}%")
        mis = cf["misleading_nodes"]
        a(f"  Misleading nodes  : {len(mis)}  "
          f"(UFBoot≥95 but gCF<50 or sCF<33)  [{_flag(cf['flag_misleading_concordance'])}]")
        for n in mis:
            ub = f"{n['ufboot']:.0f}" if n["ufboot"] is not None else "N/A"
            a(f"    Branch {n['branch_id']:>3}: UFBoot={ub}  "
              f"gCF={n['gcf']:.1f}%  sCF={n['scf']:.1f}%")
        low = cf["low_gcf_nodes"]
        if low:
            a(f"  Low-gCF nodes (<50%): {len(low)}")
            for n in low:
                ub = f"{n['ufboot']:.0f}" if n["ufboot"] is not None else "N/A"
                a(f"    Branch {n['branch_id']:>3}: UFBoot={ub}  "
                  f"gCF={n['gcf']:.1f}%  sCF={n['scf']:.1f}%")

    # ANI
    a("\n── Average Nucleotide Identity ─────────────────────────────")
    ani = summary.get("ani", {})
    if "error" in ani:
        a(f"  ERROR: {ani['error']}")
    else:
        a(f"  Pairwise comparisons : {ani['n_pairwise_comparisons']}")
        a(f"  ANI range            : {ani['ani_min']:.2f}% – {ani['ani_max']:.2f}%")
        a(f"  ANI mean             : {ani['ani_mean']:.2f}%")
        pairs = ani["near_identical_pairs"]
        a(f"  Near-identical >99%  : {len(pairs)} pairs  [{_flag(ani['flag_near_identical'])}]")
        for p in pairs:
            a(f"    {p['taxon_a']}  ↔  {p['taxon_b']}  ({p['ani']:.2f}%)")

    # ASCII species tree — topology + support (* = UFBoot ≥ 70)
    a("\n── Species Tree (ASCII) ────────────────────────────────────")
    tree_file = Path(summary["results_dir"]) / "build" / "species_tree.nwk"
    a(render_ascii_tree(tree_file))

    a("\n" + "=" * 62)

    out_path.write_text("\n".join(lines) + "\n")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Summarise an MitoTree pipeline run into run_summary.json/.txt",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "results_dir",
        nargs="?",
        default="analysis_results",
        help="MitoTree output directory to read results from (default: analysis_results)",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        metavar="DIR",
        help="Directory to write run_summary.json/.txt (default: same as results_dir)",
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir).resolve()
    if not results_dir.exists():
        print(f"Error: results directory not found: {results_dir}", file=sys.stderr)
        sys.exit(2)

    output_dir = Path(args.output_dir).resolve() if args.output_dir else results_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Summarising: {results_dir}")

    summary = {
        "results_dir": str(results_dir),
        "supermatrix": parse_supermatrix(results_dir),
        "species_tree": parse_species_tree(results_dir),
        "rf_distances": parse_rf_distances(results_dir),
        "concordance": parse_concordance(results_dir),
        "ani": parse_ani(results_dir),
    }

    json_out = output_dir / "run_summary.json"
    txt_out = output_dir / "run_summary.txt"

    json_out.write_text(json.dumps(summary, indent=2, default=_jsonify))
    write_text_report(summary, txt_out)

    print(f"Written: {json_out}")
    print(f"Written: {txt_out}")

    active_flags = [
        summary["species_tree"].get("flag_low_support"),
        summary["rf_distances"].get("flag_rf_outliers"),
        summary["concordance"].get("flag_misleading_concordance"),
        summary["ani"].get("flag_near_identical"),
    ]
    if any(active_flags):
        print("\nQC flags raised — review run_summary.txt for details.", file=sys.stderr)
        sys.exit(1)
    else:
        print("\nAll QC checks passed.")


if __name__ == "__main__":
    main()
