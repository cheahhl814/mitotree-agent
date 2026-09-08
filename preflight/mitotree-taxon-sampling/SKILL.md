---
name: mitotree-taxon-sampling
description: 'Broaden a mitotree-agent dataset by fetching whole complete-mitogenome sequences from NCBI (same genus/family/etc., primary GenBank submissions only — RefSeq records excluded) plus an optional outgroup, then BLASTN-extract the 13 PCGs + 12S/16S rRNA against a reference gene set from the user''s own already-annotated data. Triggers: "add more mitogenomes", "pull mitogenomes from NCBI", "expand the mitotree dataset", "mitotree-taxon-sampling".'
version: 1.1.0
updated: "2026-09-08"
triggers:
  - "add more mitogenomes"
  - "pull mitogenomes from NCBI"
  - "expand the mitotree dataset"
  - "fetch mitogenomes"
  - "mitotree-taxon-sampling"
---

# mitotree-taxon-sampling

> **v1.1.0.** Optional stage of `mitotree-agent`, run before `preflight/mitotree-preflight`. Two scripts: `fetch_mitogenomes.py` (NCBI retrieval, adapted from `phylogenetics-agent`'s `taxon_sampling.py`) and `blast_extract_genes.py` (BLASTN gene extraction, trimmed down from `PHYLOGENETIC/scripts/blast_extract_cds.py` to just its core homology-extraction logic, reference-loading reused directly from `build_supermatrix.py`).
>
> **Primary-submission-only (v1.1.0):** RefSeq records (`NC_`/`NG_`-style accessions) are **excluded**, not preferred. The esearch query carries `NOT srcdb_refseq[PROP]` and a post-fetch guard drops any straggler, logging it in `metadata.csv` with `drop_reason=refseq_excluded`. When a primary submission and its byte-identical RefSeq mirror both exist, the **primary accession is the one kept**. Type material still outranks ordinary primary submissions in the quality ranking.

## When to Use This Skill

Use before `preflight/mitotree-preflight`, when the user wants to broaden a mitogenome dataset with more taxa from NCBI. Most NCBI mitogenome submissions are **not** split into per-gene FASTA records the way `build/mitotree-builder` needs — this skill fetches the raw whole-genome sequence, then extracts the 13 PCGs + 12S/16S via BLASTN against a reference gene set the user already has (from a prior `mitotree-agent` run or a MitoAnnotator archive), sidestepping the need to run MitoAnnotator on every newly-fetched genome.

Requires the user already have **at least one** already-annotated sample to serve as the BLAST reference — this skill broadens an existing dataset, it doesn't bootstrap one from nothing.

## 0. Inputs / Outputs contract

### Inputs (consumed)

| Path / Arg | Required? | Notes |
| --- | --- | --- |
| `--expand-taxon TAXON` | yes | e.g. `"Betta"` or `"Osphronemidae"` — esearch query, `mitochondrion[Filter] AND "complete genome"[Title]`. |
| `--min-length`/`--max-length` | no (default 14000/20000) | Pushed into the esearch query itself (`SLEN` field) — NCBI does the length filtering, not a local post-filter. |
| `--outgroup-taxon TAXON` | no | Same retrieval path, best single candidate kept. |
| One of `--reference-archive` (MitoAnnotator `.tar.gz`) / `--reference-fasta LABEL:FILE` / `--reference-per-gene-dir` (e.g. a prior run's `build/01_per_gene/`) | yes, for the extraction step | Reference gene sequences for BLASTN — loaded via `build_supermatrix.py`'s own `genes_from_archive`/`genes_from_fasta`, or a small per-gene-file loader for the `build/01_per_gene/` layout. |

### Outputs (produced)

| Path | Format |
| --- | --- |
| `$RUN_DIR/taxon_sampling/raw/<label>.fasta` | Whole-genome FASTA, one per kept accession (post name-quality filter + ≥99.9% dedup + per-taxon cap) |
| `$RUN_DIR/taxon_sampling/summary.csv` | Per-taxon: n_kept, best_quality, n_dropped |
| `$RUN_DIR/taxon_sampling/metadata.csv` | Per-accession, kept AND dropped, full provenance |
| `$RUN_DIR/taxon_sampling/genes/<label>.genes.fasta` | BLASTN-extracted per-gene FASTA — directly usable as `build_supermatrix.py --fasta LABEL:FILE` |
| `$RUN_DIR/taxon_sampling/genes/no_hit_genes.txt` | Per-label list of genes with zero BLAST hit, if any |

## 0.5 Ask-User Stop Points

### SP1 — No reference provided, or reference has an incomplete gene set

| Trigger | Evidence | Action |
| --- | --- | --- |
| `blast_extract_genes.py` reports `<15 genes loaded` from the reference (live-verified: a real MitoAnnotator archive loaded all 15/15) | Reference gene count printed before extraction starts | **Recommend**: use the most completely-annotated sample on hand as the reference. Ask: "(A) proceed with the partial reference (missing genes will be `0/n` for every fetched sample, not just this one), (B) point me at a more complete reference" |

**Auto-pick when**: reference has all 15 genes (13 PCG + 2 rRNA). No ask.

### SP2 — A fetched sample gets zero BLAST hits across every gene

| Trigger | Evidence | Action |
| --- | --- | --- |
| `blast_extract_genes.py` prints `ALL genes missed for <label>` (this is a distinct signature from a few individual gene misses) | The target is too evolutionarily divergent from the reference for BLASTN to find homology at all — not an annotation gap | Ask: "`<label>` matched none of the reference genes — likely too distant from the reference taxon for BLASTN. (A) drop this sample, (B) supply a closer reference and re-run just this sample" |

**Auto-pick when**: at least one gene extracted successfully. No ask — route individual gene misses through `mitotree-preflight` SP1 downstream instead.

### SP3 — Near-duplicate collapse threshold

| Trigger | Evidence | Action |
| --- | --- | --- |
| User wants a different redundancy threshold than the default ≥99.9% | n/a | Ask only if the user pushes back on the default; otherwise proceed with `--redundancy-threshold 0.999` (mirrors `phylo-taxon-sampling`'s convention). |

**Auto-pick when**: no objection to the default. No ask.

### SP4 — Length-range default doesn't fit the taxon

| Trigger | Evidence | Action |
| --- | --- | --- |
| Zero (or very few) hits returned at the default 14000–20000 bp range | esearch `Count` in the fetch log | **Recommend**: some taxa have unusually large/small mitogenomes (e.g. some invertebrates); widen `--min-length`/`--max-length` and retry rather than assuming zero NCBI coverage. |

**Auto-pick when**: hits found at the default range. No ask.

## Procedure

### 1. Fetch whole mitogenomes

```bash
python3 preflight/mitotree-taxon-sampling/scripts/fetch_mitogenomes.py \
    --expand-taxon "Betta" \
    --outgroup-taxon "Macropodus opercularis" \
    -o "$RUN_DIR/taxon_sampling/raw" \
    --summary-out "$RUN_DIR/taxon_sampling/summary.csv" \
    --metadata-out "$RUN_DIR/taxon_sampling/metadata.csv" \
    --esearch-cmd "pixi run esearch" --efetch-cmd "pixi run efetch" --vsearch-cmd "pixi run vsearch"
```

Live-verified against real NCBI data (`--expand-taxon "Betta"`, no network mocking): 12 candidate records fetched → **9 kept after the `NOT srcdb_refseq[PROP]` exclusion (3 RefSeq records dropped)**; all 9 kept accessions are primary GenBank submissions (AB/KR/LC/MT/MW/OQ/PX prefixes, zero `XX_`-style RefSeq accessions). Previously (v1.0.0) the behavior was the opposite — RefSeq was preferred over a byte-identical non-RefSeq duplicate.

### 2. Extract genes via BLASTN against a reference

```bash
python3 preflight/mitotree-taxon-sampling/scripts/blast_extract_genes.py \
    --targets-dir "$RUN_DIR/taxon_sampling/raw" \
    --reference-archive /path/to/a/MitoAnnotator/archive.tar.gz \
    -o "$RUN_DIR/taxon_sampling/genes" \
    --blastn-cmd "pixi run blastn" --makeblastdb-cmd "pixi run makeblastdb"
```

Live-verified against the 8 real fetched *Betta* genomes from step 1, using a real MitoAnnotator archive (`Betta_pi_BPI001.tar.gz`) as reference: 7/8 samples extracted all 15 genes cleanly; the 8th missed only `ATPase8` (a naturally short, fast-evolving gene) — a realistic, expected outcome at real cross-species divergence, not a bug.

### 3. Hand off to preflight

`$RUN_DIR/taxon_sampling/genes/*.genes.fasta` are ready to pass straight to `preflight/mitotree-preflight` as `--fasta LABEL:FILE` specs, alongside any existing MitoAnnotator archives.

## Interpretation Guidelines

- **`no_hit_genes.txt` lists a handful of individual genes across a few samples**: normal — route through `mitotree-preflight` SP1 like any other partial-coverage gap.
- **`ALL genes missed for <label>`**: SP2 — reference too divergent, not a coverage gap.
- **Reference loaded < 15/15 genes**: SP1 — the gap propagates to every fetched sample equally, fix the reference first.

## Troubleshooting — Signature library

| Signature | Likely cause | Fix |
| --- | --- | --- |
| `curl: (56) OpenSSL SSL_read: ... unexpected eof` in `esearch` output, but still followed by a valid `<ENTREZ_DIRECT>` block | Transient NCBI/network hiccup; entrez-direct retries automatically (verified live) | No action needed if the XML block with a real `Count` still appears after it; only investigate if the command exits with no XML at all. |
| `Warning: [blastn] Examining 5 or more matches is recommended` | Cosmetic — `-max_target_seqs 1` is intentional (best hit only), this warning is expected and harmless | Ignore. |
| Many genes `WARNING: no BLAST hit` across most samples, but a few genes (e.g. rRNAs, COXI) do hit | **Real bug found and fixed live**: `blastn`'s default `-task` is `megablast` (tuned for near-identical sequences) — it misses real cross-species homologs that classic `-task blastn` finds. Confirmed live: a 78%-identity Cytb hit was invisible to megablast but found immediately with `-task blastn`. | Already fixed in `blast_extract_genes.py` (always passes `-task blastn` explicitly). If you see this pattern anyway, check the script wasn't reverted or a custom `--blastn-cmd` bypasses the `-task` flag. |
| `ERROR: reference has zero usable genes` | Reference archive/FASTA path wrong, or its gene headers don't match the canonical PCG/rRNA name set | Check with `tar -tzf <archive>` + inspect `seq_genes.fa` headers, same as `mitotree-preflight`'s equivalent troubleshooting row. |
| `no *.fasta files in <targets-dir>` | `fetch_mitogenomes.py` wasn't run first, or `-o`/`--targets-dir` point at different directories | Run step 1 before step 2; confirm the paths match. |

## Verification

- [ ] `taxon_sampling/raw/` has ≥1 whole-genome FASTA, `summary.csv` and `metadata.csv` are non-empty.
- [ ] `taxon_sampling/genes/` has one `.genes.fasta` per successfully-extracted sample.
- [ ] Reference gene count was 15/15, or SP1 was surfaced to the user if not.

## What NOT to do

- Do **not** fall back to fetching GenBank-format annotations and parsing their own CDS/rRNA features — that was the original design and was deliberately dropped in favor of BLASTN extraction against the user's own reference, per explicit user direction, since it reuses this project's own already-working `blast_extract_cds.py`/`core_mito_phylo.py` logic instead of a new GenBank-feature parser.
- Do **not** omit `-task blastn` from the BLASTN call — the default `megablast` task silently misses real cross-species gene homologs (verified live).
- Do **not** re-add a RefSeq preference or remove the `NOT srcdb_refseq[PROP]` filter / the post-fetch `is_refseq` guard — primary-submission-only retrieval is the deliberate v1.1.0 contract (a byte-identical RefSeq mirror must never win over the primary accession it duplicates).

## Handoff

> Fetched `<N>` whole mitogenomes, extracted genes for `<M>`/`<N>`. `taxon_sampling/genes/*.genes.fasta` ready for `preflight/mitotree-preflight --fasta LABEL:FILE`.
