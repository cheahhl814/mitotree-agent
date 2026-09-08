#!/usr/bin/env python3
"""Retrieve whole complete-mitogenome sequences from NCBI to broaden a
mitotree-agent dataset (same genus/family/etc.) and optionally add an
outgroup, with quality filtering, near-duplicate collapse, and a per-taxon
representative cap. Adapted from phylogenetics-agent's taxon_sampling.py --
same pipeline shape, but for whole mitogenomes rather than a single marker,
and with an NCBI-side length filter (SLEN) instead of a marker-title search.

This deliberately fetches raw whole-genome FASTA only -- gene-level
extraction (13 PCGs + 12S/16S) happens downstream in blast_extract_genes.py
against a reference gene set, not here. Most complete-mitogenome submissions
on NCBI are annotated in GenBank format but NOT split into per-gene FASTA
records, which is exactly the gap blast_extract_genes.py closes.

Primary-submission-only: RefSeq records (NC_/NG_-style accessions) are
excluded at the esearch level (`NOT srcdb_refseq[PROP]`) and again as a
post-fetch guard; excluded records are still logged in the metadata CSV with
drop_reason=refseq_excluded for traceability. When a primary submission and
its RefSeq counterpart describe the same sequence, the PRIMARY accession is
the one kept.

Pipeline: NCBI esearch (SLEN-filtered, RefSeq-excluded) / efetch (GenBank
          flat, for metadata + sequence) -> RefSeq guard -> name-quality
          filter (drop cf./sp./hybrid/uncultured) -> sequence-quality ranking
          (length band, ambiguity, type material, recency) -> >=99.9%
          redundancy collapse (vsearch --cluster_fast) -> per-taxon cap ->
          outgroup fetch.

Outputs: raw whole-genome FASTA(s) (one file per kept accession, ready for
blast_extract_genes.py), a per-taxon summary CSV, and a per-accession
metadata CSV (kept AND dropped, for full traceability).

Requires: esearch/efetch (entrez-direct), vsearch on PATH (or pass
--esearch-cmd/--efetch-cmd/--vsearch-cmd, e.g. "pixi run vsearch").
"""
import argparse
import csv
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

REFSEQ_PREFIX_RE = re.compile(r"^[A-Z]{2}_")
SP_RE = re.compile(r"\bspp?\.(\s|$)")


def run(argv: list[str]) -> str:
    result = subprocess.run(argv, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"WARNING: command failed: {argv}\n{result.stderr}", file=sys.stderr)
        return ""
    return result.stdout


def run_piped(argv_list: list[list[str]]) -> str:
    """Run a list of argv commands piped left-to-right, no shell involved."""
    procs = []
    prev_stdout = None
    for i, argv in enumerate(argv_list):
        is_last = i == len(argv_list) - 1
        p = subprocess.Popen(argv, stdin=prev_stdout, stdout=subprocess.PIPE if is_last else subprocess.PIPE,
                              stderr=subprocess.PIPE, text=True)
        if prev_stdout is not None:
            prev_stdout.close()
        prev_stdout = p.stdout
        procs.append(p)
    out, err = procs[-1].communicate()
    for p in procs[:-1]:
        p.wait()
    if procs[-1].returncode != 0:
        print(f"WARNING: pipeline failed: {argv_list}\n{err}", file=sys.stderr)
        return ""
    return out


def fetch_gb_records(taxon: str, min_length: int, max_length: int, esearch_cmd: list[str],
                      efetch_cmd: list[str], retmax: int) -> list[str]:
    # NOT srcdb_refseq[PROP] keeps primary GenBank submissions only — RefSeq
    # mirrors are excluded at the source (live-verified: Betta 12 -> 9 records).
    query = (
        f'"{taxon}"[Organism] AND mitochondrion[Filter] AND "complete genome"[Title] '
        f'AND {min_length}:{max_length}[SLEN] NOT srcdb_refseq[PROP]'
    )
    esearch_argv = [*esearch_cmd, "-db", "nucleotide", "-query", query]
    efetch_argv = [*efetch_cmd, "-format", "gb", "-start", "1", "-stop", str(retmax)]
    text = run_piped([esearch_argv, efetch_argv])
    if not text.strip():
        return []
    return [rec for rec in re.split(r"\n//\n", text) if rec.strip()]


def parse_gb_record(text: str) -> dict | None:
    locus_m = re.search(r"^LOCUS\s+(\S+)\s+(\d+)\s+bp.*?(\d{2}-[A-Z]{3}-\d{4})?\s*$", text, re.MULTILINE)
    if not locus_m:
        return None
    accession, length, sub_date = locus_m.group(1), int(locus_m.group(2)), locus_m.group(3) or ""

    org_m = re.search(r"^\s*ORGANISM\s+(.+)$", text, re.MULTILINE)
    organism = org_m.group(1).strip() if org_m else ""

    country_m = re.search(r'/country="([^"]+)"', text)
    coldate_m = re.search(r'/collection_date="([^"]+)"', text)
    latlon_m = re.search(r'/lat_lon="([^"]+)"', text)
    type_m = re.search(r"/type_material=", text) or re.search(r"type strain", text, re.IGNORECASE)

    # Records are split on "\n//\n" upstream, which strips the trailing "//"
    # terminator from every record -- so ORIGIN's sequence block just runs to
    # the end of `text` here, not up to a "^//" anchor.
    origin_m = re.search(r"^ORIGIN.*?\n(.*)", text, re.DOTALL | re.MULTILINE)
    seq = ""
    if origin_m:
        seq = re.sub(r"[\d\s]", "", origin_m.group(1)).upper()

    n_ambig = sum(1 for c in seq if c not in "ACGT")

    return {
        "accession": accession,
        "organism": organism,
        "length": length,
        "sequence": seq,
        "n_ambig": n_ambig,
        "country": country_m.group(1) if country_m else "NA",
        "collection_date": coldate_m.group(1) if coldate_m else "NA",
        "lat_lon": latlon_m.group(1) if latlon_m else "NA",
        "is_refseq": bool(REFSEQ_PREFIX_RE.match(accession)),
        "is_type_material": bool(type_m),
        "submission_date": sub_date,
    }


def name_quality_ok(organism: str, keep_cf: bool, keep_sp: bool) -> bool:
    if not organism:
        return False
    if re.search(r"uncultured|environmental sample|metagenome", organism, re.IGNORECASE):
        return False
    if "×" in organism or " hybrid" in organism.lower():
        return False
    if not keep_cf and " cf. " in organism:
        return False
    if not keep_sp and SP_RE.search(organism):
        return False
    return True


def priority_key(rec: dict, expected_length: int):
    """Lower is better. RefSeq never appears here (excluded upstream); type
    material still outranks ordinary primary submissions."""
    ambig_frac = (rec["n_ambig"] / rec["length"]) if rec["length"] else 1.0
    length_dev = abs(rec["length"] - expected_length) / expected_length if expected_length else 0.0
    return (
        0 if rec["is_type_material"] else 1,
        round(ambig_frac, 4),
        round(length_dev, 4),
        rec["submission_date"],
    )


def collapse_redundant(records: list[dict], vsearch_cmd: list[str], threshold: float, tmp_prefix: str) -> tuple[list[dict], dict]:
    """>= `threshold` identical sequences within this taxon's candidate pool are
    collapsed to the single highest-priority (already-sorted) representative."""
    if len(records) < 2:
        return records, {}

    fasta_path = f"{tmp_prefix}.fasta"
    uc_path = f"{tmp_prefix}.uc"
    with open(fasta_path, "w") as f:
        for r in records:
            f.write(f">{r['accession']}\n{r['sequence']}\n")

    run([*vsearch_cmd, "--cluster_fast", fasta_path, "--id", str(threshold), "--uc", uc_path, "--quiet"])

    cluster_of = {}
    try:
        with open(uc_path) as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols[0] in ("H", "S"):
                    cluster_of[cols[8]] = cols[1]
    except FileNotFoundError:
        return records, {}

    by_cluster = defaultdict(list)
    for r in records:
        by_cluster[cluster_of.get(r["accession"], r["accession"])].append(r)

    kept, dropped_of = [], {}
    for cluster_records in by_cluster.values():
        cluster_records.sort(key=lambda r: r["_priority"])
        best = cluster_records[0]
        kept.append(best)
        for r in cluster_records[1:]:
            dropped_of[r["accession"]] = best["accession"]
    return kept, dropped_of


def fetch_outgroup(outgroup_taxon: str, min_length: int, max_length: int, esearch_cmd: list[str], efetch_cmd: list[str]) -> dict | None:
    raw = fetch_gb_records(outgroup_taxon, min_length, max_length, esearch_cmd, efetch_cmd, retmax=5)
    candidates = [r for rec in raw if (r := parse_gb_record(rec))]
    candidates = [r for r in candidates if r["sequence"] and not r["is_refseq"]]
    if not candidates:
        return None
    for r in candidates:
        r["_priority"] = priority_key(r, expected_length=candidates[0]["length"])
    candidates.sort(key=lambda r: r["_priority"])
    return candidates[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--expand-taxon", required=True, help='Taxon name to search within, e.g. "Betta" or "Osphronemidae"')
    ap.add_argument("--min-length", type=int, default=14000, help="minimum mitogenome length in bp (default: 14000)")
    ap.add_argument("--max-length", type=int, default=20000, help="maximum mitogenome length in bp (default: 20000)")
    ap.add_argument("--n-per-taxon", type=int, default=2)
    ap.add_argument("--redundancy-threshold", type=float, default=0.999)
    ap.add_argument("--keep-cf", action="store_true")
    ap.add_argument("--keep-sp", action="store_true")
    ap.add_argument("--outgroup-taxon", help="explicit outgroup taxon name")
    ap.add_argument("--retmax", type=int, default=200, help="max NCBI hits to consider before filtering")
    ap.add_argument("-o", "--outdir", required=True, help="output directory for one raw FASTA per kept accession")
    ap.add_argument("--summary-out", required=True, help="per-taxon summary CSV")
    ap.add_argument("--metadata-out", required=True, help="per-accession metadata CSV (kept + dropped)")
    ap.add_argument("--tmp-prefix", default="/tmp/mitotree_taxon_sampling")
    ap.add_argument("--esearch-cmd", default="esearch")
    ap.add_argument("--efetch-cmd", default="efetch")
    ap.add_argument("--vsearch-cmd", default="vsearch")
    args = ap.parse_args()
    args.esearch_cmd = args.esearch_cmd.split()
    args.efetch_cmd = args.efetch_cmd.split()
    args.vsearch_cmd = args.vsearch_cmd.split()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw = fetch_gb_records(args.expand_taxon, args.min_length, args.max_length, args.esearch_cmd, args.efetch_cmd, args.retmax)
    all_records = [r for rec in raw if (r := parse_gb_record(rec))]
    all_records = [r for r in all_records if r["sequence"]]
    # Post-fetch RefSeq guard (belt-and-braces behind the NOT srcdb_refseq[PROP]
    # query filter; also catches records that slip through a custom --esearch-cmd).
    # Excluded records stay in all_records so they appear in metadata.csv.
    refseq_dropped = {r["accession"]: "refseq_excluded" for r in all_records if r["is_refseq"]}
    working_records = [r for r in all_records if not r["is_refseq"]]
    print(f"Fetched: {len(all_records)} candidate complete-mitogenome records for '{args.expand_taxon}' "
          f"({args.min_length}-{args.max_length} bp)")
    print(f"Primary-submission-only: {len(refseq_dropped)} RefSeq record(s) excluded")

    filtered = [r for r in working_records if name_quality_ok(r["organism"], args.keep_cf, args.keep_sp)]
    name_dropped = {r["accession"]: "name_quality" for r in working_records if r not in filtered}
    print(f"After name-quality filter: {len(filtered)} (dropped {len(name_dropped)})")

    by_taxon = defaultdict(list)
    for r in filtered:
        by_taxon[r["organism"]].append(r)

    kept_all, dup_dropped, cap_dropped = [], {}, {}
    summary_rows = []
    for taxon, recs in by_taxon.items():
        expected_length = sorted(x["length"] for x in recs)[len(recs) // 2]
        for r in recs:
            r["_priority"] = priority_key(r, expected_length)
        recs.sort(key=lambda r: r["_priority"])

        collapsed, dropped_of = collapse_redundant(recs, args.vsearch_cmd, args.redundancy_threshold, f"{args.tmp_prefix}_{re.sub(r'[^A-Za-z0-9]', '_', taxon)}")
        dup_dropped.update(dropped_of)
        collapsed.sort(key=lambda r: r["_priority"])

        reps = collapsed[:args.n_per_taxon]
        for r in collapsed[args.n_per_taxon:]:
            cap_dropped[r["accession"]] = "over_cap"
        kept_all.extend(reps)

        best_q = "type_material" if any(r["is_type_material"] for r in reps) else "none"
        n_dropped_taxon = len(recs) - len(reps)
        summary_rows.append({
            "taxon": taxon, "n_kept": len(reps), "best_quality": best_q, "n_dropped": n_dropped_taxon,
        })

    outgroup_rec = None
    if args.outgroup_taxon:
        outgroup_rec = fetch_outgroup(args.outgroup_taxon, args.min_length, args.max_length, args.esearch_cmd, args.efetch_cmd)
    if outgroup_rec:
        kept_all.append(outgroup_rec)
        summary_rows.append({"taxon": outgroup_rec["organism"], "n_kept": 1, "best_quality": "outgroup", "n_dropped": 0})

    for r in kept_all:
        label = f"{r['organism'].replace(' ', '_')}_{r['accession']}"
        (outdir / f"{label}.fasta").write_text(f">{label}\n{r['sequence']}\n")
    print(f"\nWritten {len(kept_all)} whole-mitogenome FASTA(s) across {len(summary_rows)} taxa -> {outdir}/")

    with open(args.summary_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["taxon", "n_kept", "best_quality", "n_dropped"])
        w.writeheader()
        w.writerows(sorted(summary_rows, key=lambda r: -r["n_kept"]))
    print(f"Summary table -> {args.summary_out}")

    kept_accessions = {r["accession"] for r in kept_all}
    with open(args.metadata_out, "w", newline="") as f:
        fieldnames = ["accession", "organism", "length", "n_ambig", "country", "collection_date",
                      "lat_lon", "is_refseq", "is_type_material", "submission_date", "kept", "drop_reason"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in all_records:
            drop_reason = ("refseq_excluded" if r["accession"] in refseq_dropped else
                name_dropped.get(r["accession"]) or
                (f"duplicate_of_{dup_dropped[r['accession']]}" if r["accession"] in dup_dropped else "") or
                cap_dropped.get(r["accession"], ""))
            w.writerow({
                "accession": r["accession"], "organism": r["organism"], "length": r["length"],
                "n_ambig": r["n_ambig"], "country": r["country"], "collection_date": r["collection_date"],
                "lat_lon": r["lat_lon"], "is_refseq": r["is_refseq"], "is_type_material": r["is_type_material"],
                "submission_date": r["submission_date"],
                "kept": r["accession"] in kept_accessions, "drop_reason": drop_reason,
            })
        if outgroup_rec:
            w.writerow({
                "accession": outgroup_rec["accession"], "organism": outgroup_rec["organism"],
                "length": outgroup_rec["length"], "n_ambig": outgroup_rec["n_ambig"],
                "country": outgroup_rec["country"], "collection_date": outgroup_rec["collection_date"],
                "lat_lon": outgroup_rec["lat_lon"], "is_refseq": outgroup_rec["is_refseq"],
                "is_type_material": outgroup_rec["is_type_material"], "submission_date": outgroup_rec["submission_date"],
                "kept": True, "drop_reason": "",
            })
    print(f"Per-accession metadata -> {args.metadata_out}")


if __name__ == "__main__":
    main()
