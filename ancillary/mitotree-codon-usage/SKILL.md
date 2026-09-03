---
name: mitotree-codon-usage
description: 'Calculate Relative Synonymous Codon Usage (RSCU) for nucleotide CDS sequences via calculate_rscu.py. Informational only — no QC flags. Nucleotide mode only (meaningless on amino-acid input). Triggers: "codon usage", "RSCU", "codon bias", "mitotree-codon-usage".'
version: 1.0.0
updated: "2026-08-26"
triggers:
  - "codon usage"
  - "RSCU"
  - "codon bias"
  - "mitotree-codon-usage"
---

# mitotree-codon-usage

> **v1.0.0.** Optional ancillary branch of `mitotree-agent`, run any time after `build/mitotree-builder` has produced per-gene nucleotide FASTAs (`build/01_per_gene/<gene>.fasta`). Wraps `calculate_rscu.py` (direct port from the original MitoTree Nextflow pipeline's `bin/calculate_rscu.py`). Outputs are **informational only** — `qc/mitotree-qc/scripts/qc_report.py` does NOT currently read codon-usage output (no QC flags), so this branch is purely for user interest.

## When to Use This Skill

Use this when you want to:
- Examine codon usage bias across mitochondrial protein-coding genes.
- Compare RSCU patterns between samples or genes.
- Investigate mitochondrial genetic code usage (default: vertebrate mitochondrial, NCBI table 2).

**Do NOT use this skill** if:
- You're working in **amino-acid (`aa`) mode** — RSCU is a codon-usage metric computed from nucleotide CDS sequences; it is meaningless on translated amino-acid alignments.
- Your input is rRNA (12S/16S) — RSCU applies only to protein-coding genes with codons.
- You expect QC flags — codon usage is informational only, not a QC checkpoint.

## 0. Inputs / Outputs contract

### Inputs

| Path / Arg | Required? | Notes |
| --- | --- | --- |
| `<gene>.fasta` (nucleotide CDS) | yes, one file per invocation | Typically one of `build/mitotree-builder`'s `build/01_per_gene/<gene>.fasta` outputs (PCGs only: ND1, ND2, COXI, etc. — not rRNA). |
| `--genetic-code N` | no (default: 2) | NCBI translation table: 1=Standard (nuclear), 2=Vertebrate Mitochondrial (default), 5=Invertebrate Mitochondrial. |

### Outputs

| Path | Format | Notes |
| --- | --- | --- |
| `$RUN_DIR/codon_usage/<gene>.codon_usage.tsv` | TSV (4 columns + header) | Columns: `Amino_Acid`, `Codon`, `Count`, `RSCU` — one row per codon (61 sense codons, excluding 3 stop codons). |

## 0.5 Ask-User Stop Points

### SP1 — Amino-acid mode detected

| Trigger | Evidence | Action |
| --- | --- | --- |
| User invokes this skill but `params.json.mode` = `aa` (or input FASTA headers look like amino-acid sequences: no T/U, or contains B/Z/X) | RSCU is meaningless on amino-acid input | **Refuse and explain**: "RSCU (Relative Synonymous Codon Usage) is computed from nucleotide codons — it cannot be calculated from amino-acid sequences. Re-run `build/mitotree-builder` in `nt` mode first, or provide nucleotide CDS FASTAs." |

**Auto-pick when**: input is confirmed nucleotide → proceed.

### SP2 — Genetic code mismatch

| Trigger | Evidence | Action |
| --- | --- | --- |
| User's clade uses a different mitochondrial genetic code (e.g. invertebrates use table 5, not table 2) | User mentions invertebrate/non-vertebrate clade | **Recommend**: "Default genetic code is 2 (Vertebrate Mitochondrial). For invertebrates, use `--genetic-code 5`. For nuclear genes, use `--genetic-code 1`. Proceed with default, or specify?" |

**Auto-pick when**: user explicitly provides `--genetic-code` → use their value.

### SP3 — Overwrite existing codon-usage output

| Trigger | Evidence | Action |
| --- | --- | --- |
| `$RUN_DIR/codon_usage/<gene>.codon_usage.tsv` already exists | A previous run is present | Ask before overwriting: "(A) overwrite, (B) write to a new subdirectory (e.g. `codon_usage_v2/`)" |

**Auto-pick when**: no existing file at that path.

## Procedure

### 1. Copy the script from the Nextflow reference

The script `bin/calculate_rscu.py` from the original MitoTree repository is copied into `scripts/calculate_rscu.py` as-is — it has no MitoTree-specific hardcoded assumptions and is directly reusable.

### 2. Run RSCU calculation per gene

```bash
run_dir="$RUN_DIR"
mkdir -p "$run_dir/codon_usage"

# For each gene FASTA from build/01_per_gene/
for gene_fasta in "$run_dir/build/01_per_gene"/*.fasta; do
    gene_name=$(basename "$gene_fasta" .fasta)
    
    # Skip rRNA genes (RSCU is meaningless for non-coding RNA)
    if [[ "$gene_name" == "12S"* ]] || [[ "$gene_name" == "16S"* ]]; then
        echo "Skipping rRNA gene: $gene_name"
        continue
    fi
    
    pixi run python3 ancillary/mitotree-codon-usage/scripts/calculate_rscu.py \
        --genetic_code 2 \
        "$gene_fasta" \
        > "$run_dir/codon_usage/${gene_name}.codon_usage.tsv"
done
```

### 3. (Optional) Aggregate per-gene outputs into a summary

```bash
# Combine all per-gene TSVs into one file (optional, for downstream analysis)
python3 ancillary/mitotree-codon-usage/scripts/aggregate_codon_usage.py \
    -i "$run_dir/codon_usage" \
    -o "$run_dir/codon_usage/all_genes_codon_usage.tsv"
```

## Interpretation Guidelines

- **RSCU = 1.0**: codon used exactly as expected under equal usage among synonymous codons.
- **RSCU > 1.0**: codon used more frequently than expected (preferred codon).
- **RSCU < 1.0**: codon used less frequently than expected (disfavored codon).
- **RSCU = 0**: codon never observed in this gene's CDS.

Typical vertebrate mitochondrial patterns:
- Strong bias toward A/T-ending codons (reflecting AT-rich mitogenomes).
- Some codons absent entirely (e.g. AGA/AGG are stop codons in vertebrate mitochondrial code).

## Troubleshooting — Signature library

| Signature | Likely cause | Fix |
| --- | --- | --- |
| `Unsupported genetic code N` | `--genetic_code` value not in {1, 2, 5} | Use 1 (Standard), 2 (Vertebrate Mitochondrial), or 5 (Invertebrate Mitochondrial). |
| Input FASTA contains amino-acid sequences | User ran this in `aa` mode or on translated FASTAs | Route through SP1 — re-run `build/mitotree-builder` in `nt` mode first. |
| `calculate_rscu.py` output is empty | Input FASTA is not a valid CDS (not multiple of 3, or contains invalid bases) | Check the input FASTA — it should be the raw nucleotide CDS from MitoAnnotator, not a trimmed/translated version. |
| RSCU values all = 0 or NaN | Input sequence too short or all codons are stops | Verify the gene is a protein-coding gene (not rRNA), and the sequence is in-frame. |

## Verification

- [ ] `$RUN_DIR/codon_usage/<gene>.codon_usage.tsv` exists for each PCG gene.
- [ ] Each TSV has 61 rows (one per sense codon) + 1 header row.
- [ ] RSCU values sum to the number of synonymous codons per amino acid (e.g. Leu has 6 codons, so their RSCU values sum to 6).
- [ ] Stop codons (TAA, TAG, and in vertebrate mito: AGA, AGG) have Count = 0 and RSCU = 0.

## What NOT to do

- Do **not** run this on amino-acid sequences — RSCU is inherently a nucleotide codon metric.
- Do **not** run this on rRNA genes (12S, 16S) — they don't have codons.
- Do **not** expect QC flags from this branch — `qc_report.py` does not read codon-usage output (confirmed via grep: no "codon" references in qc_report.py).

## Handoff

> Codon usage analysis complete. `<N>` gene(s) processed: `codon_usage/<gene>.codon_usage.tsv` written. This branch is informational only — no QC flags are raised. To fold into the QC report, note that `qc/mitotree-qc` does not currently read codon-usage output.
