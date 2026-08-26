#!/usr/bin/env python3
"""
preflight_check.py -- Validate mitotree-agent inputs before the expensive
align/trim/tree build runs.

Reuses build/mitotree-builder/scripts/build_supermatrix.py's own source
discovery (discover_archives, parse_fasta_specs) and extraction (step_extract)
functions as a dry run -- no mafft/trimal/AMAS invoked, just the fast
FASTA/archive parsing + per-gene coverage table. This dry run is not just a
speed optimization: it is the main way to catch the single most common
mitotree-agent input problem early. core_mito_phylo.py extracts genes by
matching known gene-name headers (from a MitoAnnotator archive's seq_genes.fa,
or a bare FASTA whose records are already split per gene) -- it does not
annotate raw genomes. A whole, unannotated mitogenome straight off NCBI (one
FASTA record per sample, no per-gene boundaries) will match zero genes for
every sample, which this script reports as a NO-GO with a specific, actionable
explanation instead of a confusing downstream near-empty-supermatrix failure.

Usage:
  python3 preflight_check.py [source_dirs ...] [--fasta LABEL:FILE ...]
      [--mode nt|aa] [--outgroup LABEL] [--genes GENE ...]
      [--mafft-cmd CMD] [--trimal-cmd CMD] [--amas AMAS_PY]
      [--iqtree-cmd CMD] [--gotree-cmd CMD] [--fastani-cmd CMD]
      -o RUN_DIR

Writes RUN_DIR/params.json and RUN_DIR/preflight.md.
Exit codes: 0 = GO or GO-WITH-WARNINGS, 1 = NO-GO.
"""
import argparse
import json
import shutil
import sys
from pathlib import Path

# Reuse the build phase's own discovery/extraction logic instead of
# reimplementing it -- keeps preflight's coverage numbers identical to what
# build/mitotree-builder will actually see.
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "build" / "mitotree-builder" / "scripts"))
from build_supermatrix import (  # noqa: E402
    DEFAULT_NT_GENE_ORDER, DEFAULT_AA_GENE_ORDER,
    discover_archives, parse_fasta_specs, step_extract, find_amas,
)

TOOL_CHECKS = {
    "mafft": "mafft",
    "trimal": "trimal",
    "iqtree": "iqtree",
    "gotree": "gotree",
    "fastANI": "fastANI",
    "Rscript": "Rscript",
}


def audit_tools(amas_arg):
    found = {}
    for name, exe in TOOL_CHECKS.items():
        found[name] = shutil.which(exe)
    try:
        found["AMAS.py"] = find_amas(amas_arg)
    except SystemExit:
        found["AMAS.py"] = None
    return found


def check_gene_coverage(all_sources, gene_order, tmp_out):
    gene_fastas = step_extract(all_sources, gene_order, tmp_out)
    n = len(all_sources)
    coverage = {}
    for gene in gene_order:
        found_dir = tmp_out / f"{gene}.fasta"
        if found_dir.exists():
            n_found = sum(1 for line in found_dir.read_text().splitlines() if line.startswith(">"))
        else:
            n_found = 0
        coverage[gene] = (n_found, n)
    return coverage


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source_dirs", nargs="*", default=[])
    ap.add_argument("--fasta", metavar="LABEL:FILE", action="append")
    ap.add_argument("--mode", choices=["nt", "aa"], default="nt")
    ap.add_argument("--outgroup", default=None)
    ap.add_argument("--genes", nargs="+", default=None)
    ap.add_argument("--mafft-cmd", default="pixi run mafft")
    ap.add_argument("--trimal-cmd", default="pixi run trimal")
    ap.add_argument("--amas", default=None)
    ap.add_argument("--iqtree-cmd", default="pixi run iqtree")
    ap.add_argument("--gotree-cmd", default="pixi run gotree")
    ap.add_argument("--fastani-cmd", default="pixi run fastANI")
    ap.add_argument("-o", "--outdir", required=True, help="RUN_DIR to write params.json + preflight.md into")
    args = ap.parse_args()

    run_dir = Path(args.outdir)
    run_dir.mkdir(parents=True, exist_ok=True)
    gene_order = args.genes or (DEFAULT_AA_GENE_ORDER if args.mode == "aa" else DEFAULT_NT_GENE_ORDER)

    warnings, blockers = [], []

    # --- tool audit ---
    tools = audit_tools(args.amas)
    for name, path in tools.items():
        if path is None:
            blockers.append(f"tool not found: {name}")

    # --- source discovery ---
    all_sources = discover_archives(args.source_dirs) + parse_fasta_specs(args.fasta)
    if not all_sources:
        blockers.append("no sample sources discovered (check source_dirs / --fasta)")
        sample_labels = []
        coverage = {}
    else:
        sample_labels = [label for _, label, _ in all_sources]

        # --- gene-coverage dry run (SP1 evidence + the raw-genome detector) ---
        tmp_out = run_dir / "_preflight_extract_dryrun"
        coverage = check_gene_coverage(all_sources, gene_order, tmp_out)
        shutil.rmtree(tmp_out, ignore_errors=True)

        total_found = sum(f for f, _ in coverage.values())
        if total_found == 0:
            blockers.append(
                "zero gene coverage across ALL samples and ALL genes -- this pipeline extracts "
                "genes by matching known gene-name headers (12S, 16S, ND1, ... from a MitoAnnotator "
                "archive's seq_genes.fa, or a bare FASTA already split per gene). It does NOT "
                "annotate raw genomes. This looks like whole, unannotated mitogenome FASTA(s) "
                "straight from NCBI -- run MitoAnnotator (or another gene-annotation tool) first "
                "to produce per-gene records, then re-run preflight."
            )
        else:
            for gene, (found, n) in coverage.items():
                if found < n:
                    warnings.append(f"gene '{gene}' found in {found}/{n} samples")

        # --- outgroup check ---
        if args.outgroup and args.outgroup not in sample_labels:
            blockers.append(
                f"--outgroup '{args.outgroup}' not among discovered sample labels: {sample_labels}"
            )

    # --- aa-mode rRNA reminder (not a blocker) ---
    if args.mode == "aa" and any(g in gene_order for g in ("12S", "16S")):
        warnings.append("mode=aa: rRNA genes (12S/16S) are excluded from the amino-acid supermatrix (not translated) -- expected, not an error.")

    verdict = "NO-GO" if blockers else ("GO-WITH-WARNINGS" if warnings else "GO")

    params = {
        "mode": args.mode,
        "source_dirs": args.source_dirs,
        "fasta_specs": args.fasta or [],
        "gene_order": gene_order,
        "outgroup": args.outgroup,
        "sample_labels": sample_labels,
        "mafft_cmd": args.mafft_cmd,
        "trimal_cmd": args.trimal_cmd,
        "amas_cmd": tools.get("AMAS.py"),
        "iqtree_cmd": args.iqtree_cmd,
        "gotree_cmd": args.gotree_cmd,
        "fastani_cmd": args.fastani_cmd,
    }
    (run_dir / "params.json").write_text(json.dumps(params, indent=2))

    lines = [
        "# mitotree-agent Preflight Report", "",
        f"**Verdict: {verdict}**", "",
        "## Tool audit",
    ]
    for name, path in tools.items():
        lines.append(f"- {name}: {'FOUND ' + path if path else 'MISSING'}")
    lines += ["", "## Sample sources", f"{len(all_sources)} discovered: {sample_labels}", "", "## Gene coverage"]
    for gene, (found, n) in coverage.items():
        lines.append(f"- {gene}: {found}/{n}")
    if blockers:
        lines += ["", "## Blockers (NO-GO)"] + [f"- {b}" for b in blockers]
    if warnings:
        lines += ["", "## Warnings"] + [f"- {w}" for w in warnings]
    (run_dir / "preflight.md").write_text("\n".join(lines) + "\n")

    print(f"Verdict: {verdict}")
    print(f"Wrote {run_dir / 'params.json'} and {run_dir / 'preflight.md'}")
    sys.exit(1 if verdict == "NO-GO" else 0)


if __name__ == "__main__":
    main()
