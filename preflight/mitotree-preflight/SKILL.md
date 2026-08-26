---
name: mitotree-preflight
description: 'Validate inputs for mitotree-agent: MitoAnnotator .tar.gz archives and/or per-species gene FASTAs, mode (nt/aa), outgroup, gene coverage, and tool audit (mafft/trimal/amas/iqtree/gotree). Writes preflight.md + params.json. Triggers: internal — invoked by mitotree-agent orchestrator only.'
version: 1.0.0
updated: "2026-08-26"
triggers:
  - "mitotree preflight"
---

# mitotree-preflight

> **v1.0.0.** Phase 1 of `mitotree-agent`. Validates the mitogenome gene-set inputs and tool availability before the expensive multi-gene alignment/tree build runs.

## 0. Inputs / Outputs contract

### Inputs

| Path / Arg                                                                                        | Required?                | Notes                                                                                                                                                                                                                                                                                  |
| ------------------------------------------------------------------------------------------------- | ------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| One or more `annotation_dir` (MitoAnnotator `.tar.gz` archives) and/or `--fasta LABEL:FILE` specs | at least 1 sample source | Same discovery rules as `core_mito_phylo.py`: archives matched by `*.tar.gz` in each dir, bare FASTAs via explicit `LABEL:FILE`.                                                                                                                                                       |
| `--mode nt\|aa`                                                                                   | no (default `nt`)        | `aa` excludes rRNA (12S/16S) — they have no translation — and uses `LG+G`; `nt` uses `GTR+G`.                                                                                                                                                                                          |
| `--outgroup LABEL`                                                                                | no                       | Must match a sample label exactly if given; unrooted output otherwise (species tree left unrooted — no re-root step runs).                                                                                                                                                             |
| `--calibrations FILE`                                                                             | **not supported**        | Divergence-time dating (MCMCtree/PAML) is deliberately out of scope — mitogenome-only data does not reliably calibrate node ages (rate heterogeneity, saturation at deeper splits). Refuse and explain if the user asks for it; point to nuclear/multi-locus dating pipelines instead. |

### Outputs

| Path                    | Format                                                                                                        |
| ----------------------- | ------------------------------------------------------------------------------------------------------------- |
| `$RUN_DIR/preflight.md` | Markdown, verdict `GO`/`GO-WITH-WARNINGS`/`NO-GO`                                                             |
| `$RUN_DIR/params.json`  | JSON: `mode`, `source_dirs`, `fasta_specs`, `gene_order`, `outgroup`, `sample_labels`, resolved tool commands |

## 0.5 Ask-User Stop Points

### SP1 — Gene coverage gaps

| Trigger                                                                                                                   | Evidence                                                                                 | Action                                                                                                                                                                                                                                            |
| ------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| A PCG or rRNA gene (from the default 13 PCG + 12S/16S set) is missing in ≥1 sample, but at least some genes have coverage | `scripts/preflight_check.py`'s gene-coverage dry-run (`found/n` per gene, live-verified) | **Recommend**: drop the incomplete sample, or drop the gene from `--genes` if it's missing in many samples. Ask: "(A) drop sample(s) X, (B) drop gene Y from the analysis, (C) proceed anyway (AMAS will pad with gaps, inflating missing-data%)" |

**Auto-pick when**: gene coverage is 100% across all samples. No ask.

### SP2 — Outgroup not found

| Trigger                                                            | Evidence                                | Action                                                                                                    |
| ------------------------------------------------------------------ | --------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| `--outgroup LABEL` given but LABEL not in discovered sample labels | Sample label list from source discovery | Ask: "Outgroup 'LABEL' isn't among the discovered samples: [list]. Pick one, or proceed without rooting?" |

**Auto-pick when**: no `--outgroup` given → proceed unrooted, `build/mitotree-builder` skips its reroot stage.

### SP3 — Divergence-time dating requested

| Trigger                                                            | Evidence | Action                                                                                                                                                                                                                                                                                                                                                                                             |
| ------------------------------------------------------------------ | -------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| User asks for node ages, MCMCtree, PAML, or a time-calibrated tree | n/a      | **Refuse and explain**: mitogenome-only dating is unreliable (rate heterogeneity across the 13 PCGs + rRNAs, saturation at deeper splits, no fossil-informed nuclear rate to cross-check against). Recommend a nuclear multi-locus or genome-scale dataset for MCMCtree instead. Do not scaffold a dating stage even if asked twice — this is a deliberate scope exclusion, not a missing feature. |

**Auto-pick when**: n/a — always explain, never silently build a dating path.

### SP4 — aa mode with rRNA present in input

| Trigger                                                                              | Evidence                                        | Action                                                                                                                                                             |
| ------------------------------------------------------------------------------------ | ----------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--mode aa` and a `.tar.gz`/`--fasta` source's gene set includes `12S`/`16S` headers | grep the source for rRNA headers before running | **Auto-pick**: rRNA is silently excluded (never translated) — this is expected `aa` mode behavior, not an error. Mention it once in `preflight.md`, no ask needed. |

### SP5 — Zero gene coverage across all samples (likely unannotated raw mitogenomes)

| Trigger                                                                          | Evidence                                                                                                                                                                                                                                                                                                                                  | Action                                                                                                                                                                                                                                                                                                                                                                      |
| -------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Every gene shows `0/n` in the coverage dry-run (live-verified `NO-GO` signature) | This pipeline extracts genes by matching known gene-name headers (from a MitoAnnotator archive's `seq_genes.fa`, or a bare FASTA already split per gene). It does not annotate raw genomes — a very common situation for mitogenomes pulled straight from NCBI (a single whole-sequence FASTA record per sample, no per-gene boundaries). | **Refuse and explain, don't guess**: tell the user their input looks unannotated, and that they need to run MitoAnnotator (or another gene-annotation tool, e.g. MITOS) first to produce per-gene records before `mitotree-agent` can extract PCGs/rRNA. Always `NO-GO`, never a warning — there is nothing useful `build/mitotree-builder` can do with zero matched genes. |

**Auto-pick when**: n/a — always `NO-GO` and explain; not a case to proceed past even with a `--force`-style override, since the downstream supermatrix would be entirely empty.

## Procedure

### 1. Run the preflight check

`scripts/preflight_check.py` reuses `build/mitotree-builder/scripts/build_supermatrix.py`'s own `discover_archives`/`parse_fasta_specs`/`step_extract` functions directly (imported, not reimplemented) so the coverage numbers preflight reports are guaranteed identical to what the build phase will actually see. It runs the tool audit, the gene-coverage dry-run (SP1/SP5 evidence — no mafft/trimal/AMAS invoked, just FASTA/archive parsing), the outgroup check (SP2), and writes `params.json` + `preflight.md` in one pass.

```bash
python3 preflight/mitotree-preflight/scripts/preflight_check.py \
    annotation_mitoannotator/ \
    --fasta Betta_hipposideros:re_extracted_CDS.fasta \
    --mode nt \
    --outgroup Salmonella_enterica \
    --amas "$(pixi run which AMAS.py)" \
    -o "$RUN_DIR"
```

Live-verified against two synthetic 4-sample datasets: per-gene FASTAs (12S/16S/Cytb, 100% coverage) → `GO`; and whole-sequence "raw mitogenome" FASTAs with non-matching headers → `NO-GO` via SP5, with the exact unannotated-genome explanation written into `preflight.md`.

## Interpretation Guidelines

- **All tools found, 100% gene coverage, outgroup resolved or absent**: `GO`.
- **1-2 non-blocking warnings (e.g. one sample missing a single gene, no outgroup given)**: `GO-WITH-WARNINGS`.
- **Any tool missing, zero sources discovered, outgroup mismatch unresolved, or SP5's zero-coverage-across-everything signature**: `NO-GO`.

## Troubleshooting — Signature library

| Signature                                                  | Likely cause                                                                                             | Fix                                                                                                                         |
| ---------------------------------------------------------- | -------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| `no sample sources discovered`                             | Wrong `source_dirs` path, or archives don't end in `.tar.gz`                                             | Check `ls <dir>/*.tar.gz`; use `--fasta LABEL:FILE` for bare FASTAs instead.                                                |
| `AMAS.py not found`                                        | Not on PATH and not in the pixi env's bin                                                                | Pass `--amas /path/to/AMAS.py`, or confirm `amas` pixi dependency installed (`pixi run which AMAS.py`).                     |
| `seq_genes.fa not found inside <archive>`                  | Archive isn't a MitoAnnotator output (wrong tar layout)                                                  | Verify the archive was produced by MitoAnnotator, not repacked by hand.                                                     |
| One gene found in 0/n samples (but others have coverage)   | Gene name mismatch for that one gene — MitoAnnotator's header naming vs `--genes` override               | `tar -tzf <archive>` + inspect `seq_genes.fa` headers directly; match exact first-word tokens (e.g. `12S`, not `12S rRNA`). |
| **Every** gene 0/n across **every** sample (SP5, `NO-GO`)  | Input is whole, unannotated mitogenome FASTA(s), not per-gene extracted sequences (verified live)        | Run MitoAnnotator (or another gene-annotation tool) first, then re-run preflight against the annotated output.              |
| `ModuleNotFoundError: No module named 'build_supermatrix'` | `preflight_check.py` invoked from outside the skill tree, or `build/mitotree-builder/scripts/` was moved | Run from the skill root, or fix the `sys.path.insert` in `preflight_check.py` if the directory layout changed.              |

## Verification

- [ ] All 6 tools resolve (mafft, trimal, AMAS.py, iqtree, gotree, fastANI — fastANI only required if `ancillary/mitotree-ani` will run).
- [ ] `params.json` written with a non-empty `sample_labels` list.
- [ ] `preflight.md` verdict ≥ `GO-WITH-WARNINGS`, and if `NO-GO`, the report names the specific blocker (missing tool, zero sources, or SP5's unannotated-input explanation) rather than a generic failure.

## What NOT to do

- Do **not** scaffold a dating/MCMCtree stage even on repeated request (SP3) — this is a deliberate, explained scope exclusion.
- Do **not** silently drop rRNA in `aa` mode without mentioning it in `preflight.md` — it's correct behavior but should be visible to the user.
- Do **not** let a zero-coverage input (SP5) fall through as a generic `NO-GO` — name the likely cause (unannotated raw genome) explicitly, since the fix (run MitoAnnotator first) is not obvious from a bare coverage table alone.

## Handoff

> Preflight complete. Verdict: `<verdict>`. `params.json` + `preflight.md` written. Proceed to `build/mitotree-builder`.
