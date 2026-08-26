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

| Path / Arg | Required? | Notes |
| --- | --- | --- |
| One or more `annotation_dir` (MitoAnnotator `.tar.gz` archives) and/or `--fasta LABEL:FILE` specs | at least 1 sample source | Same discovery rules as `core_mito_phylo.py`: archives matched by `*.tar.gz` in each dir, bare FASTAs via explicit `LABEL:FILE`. |
| `--mode nt\|aa` | no (default `nt`) | `aa` excludes rRNA (12S/16S) — they have no translation — and uses `LG+G`; `nt` uses `GTR+G`. |
| `--outgroup LABEL` | no | Must match a sample label exactly if given; unrooted output otherwise (species tree left as OrthoFinder-style unrooted topology — no re-root step runs). |
| `--calibrations FILE` | **not supported** | Divergence-time dating (MCMCtree/PAML) is deliberately out of scope — mitogenome-only data does not reliably calibrate node ages (rate heterogeneity, saturation at deeper splits). Refuse and explain if the user asks for it; point to nuclear/multi-locus dating pipelines instead. |

### Outputs

| Path | Format |
| --- | --- |
| `$RUN_DIR/preflight.md` | Markdown, verdict `GO`/`GO-WITH-WARNINGS`/`NO-GO` |
| `$RUN_DIR/params.json` | JSON: `mode`, `source_dirs`, `fasta_specs`, `gene_order`, `outgroup`, `sample_labels`, resolved tool commands |

## 0.5 Ask-User Stop Points

### SP1 — Gene coverage gaps

| Trigger | Evidence | Action |
| --- | --- | --- |
| A PCG or rRNA gene (from the default 13 PCG + 12S/16S set) is missing in ≥1 sample | `core_mito_phylo.py`'s own coverage table (`found/n` per gene) run in dry-run mode | **Recommend**: drop the incomplete sample, or drop the gene from `--genes` if it's missing in many samples. Ask: "(A) drop sample(s) X, (B) drop gene Y from the analysis, (C) proceed anyway (AMAS will pad with gaps, inflating missing-data%)" |

**Auto-pick when**: gene coverage is 100% across all samples. No ask.

### SP2 — Outgroup not found

| Trigger | Evidence | Action |
| --- | --- | --- |
| `--outgroup LABEL` given but LABEL not in discovered sample labels | Sample label list from source discovery | Ask: "Outgroup 'LABEL' isn't among the discovered samples: [list]. Pick one, or proceed without rooting?" |

**Auto-pick when**: no `--outgroup` given → proceed unrooted, `build/mitotree-builder` skips its reroot stage.

### SP3 — Divergence-time dating requested

| Trigger | Evidence | Action |
| --- | --- | --- |
| User asks for node ages, MCMCtree, PAML, or a time-calibrated tree | n/a | **Refuse and explain**: mitogenome-only dating is unreliable (rate heterogeneity across the 13 PCGs + rRNAs, saturation at deeper splits, no fossil-informed nuclear rate to cross-check against). Recommend a nuclear multi-locus or genome-scale dataset for MCMCtree instead. Do not scaffold a dating stage even if asked twice — this is a deliberate scope exclusion, not a missing feature. |

**Auto-pick when**: n/a — always explain, never silently build a dating path.

### SP4 — aa mode with rRNA present in input

| Trigger | Evidence | Action |
| --- | --- | --- |
| `--mode aa` and a `.tar.gz`/`--fasta` source's gene set includes `12S`/`16S` headers | grep the source for rRNA headers before running | **Auto-pick**: rRNA is silently excluded (never translated) — this is expected `aa` mode behavior, not an error. Mention it once in `preflight.md`, no ask needed. |

## Procedure

### 1. Discover sources and tools

```bash
python3 -c "
import shutil, subprocess, sys
tools = {
    'mafft': 'mafft --version',
    'trimal': 'trimal --version',
    'iqtree': 'iqtree --version',
    'gotree': 'gotree version',
    'fastani': 'fastANI --version',
    'Rscript': 'Rscript --version',
}
for name, check in tools.items():
    exe = check.split()[0]
    found = shutil.which(exe)
    print(f'{name}: {\"FOUND \" + found if found else \"MISSING\"}')"
```

AMAS.py is not on `PATH` by name — resolve via `shutil.which('AMAS.py')` or the pixi env's `bin/AMAS.py` (see `build/mitotree-builder/scripts/build_supermatrix.py --amas` default search, mirrors `core_mito_phylo.py`'s `AMAS_SEARCH_PATHS`).

### 2. Gene coverage dry-run

Run `core_mito_phylo.py`'s Step 1+2 (discovery + extraction) without alignment, to get the coverage table without paying for mafft/trimal on incomplete data. If a bespoke dry-run flag doesn't exist yet, run the sources through the extraction step only and inspect stdout's coverage block before proceeding to SP1.

### 3. Write params.json

```json
{
  "mode": "nt",
  "source_dirs": ["annotation_mitoannotator/"],
  "fasta_specs": ["Betta_hipposideros:re_extracted_CDS.fasta"],
  "gene_order": ["12S", "16S", "ND1", "ND2", "COXI", "COXII", "ATPase8", "ATPase6", "COXIII", "ND3", "ND4L", "ND4", "ND5", "ND6", "Cytb"],
  "outgroup": null,
  "sample_labels": ["Betta_hipposideros", "..."],
  "mafft_cmd": "pixi run mafft",
  "trimal_cmd": "pixi run trimal",
  "amas_cmd": "/path/to/AMAS.py",
  "iqtree_cmd": "pixi run iqtree",
  "gotree_cmd": "pixi run gotree"
}
```

## Interpretation Guidelines

- **All tools found, 100% gene coverage, outgroup resolved or absent**: `GO`.
- **1-2 non-blocking warnings (e.g. one sample missing a single gene, no outgroup given)**: `GO-WITH-WARNINGS`.
- **Any tool missing, zero sources discovered, or outgroup mismatch unresolved**: `NO-GO`.

## Troubleshooting — Signature library

| Signature | Likely cause | Fix |
| --- | --- | --- |
| `No annotation sources found` | Wrong `source_dirs` path, or archives don't end in `.tar.gz` | Check `ls <dir>/*.tar.gz`; use `--fasta LABEL:FILE` for bare FASTAs instead. |
| `AMAS.py not found` | Not on PATH and not in the pixi env's bin | Pass `--amas /path/to/AMAS.py`, or confirm `amas` pixi dependency installed (`pixi run which AMAS.py`). |
| `seq_genes.fa not found inside <archive>` | Archive isn't a MitoAnnotator output (wrong tar layout) | Verify the archive was produced by MitoAnnotator, not repacked by hand. |
| Gene found in 0/n samples | Gene name mismatch — MitoAnnotator's header naming vs `--genes` override | `tar -tzf <archive>` + inspect `seq_genes.fa` headers directly; match exact first-word tokens (e.g. `12S`, not `12S rRNA`). |

## Verification

- [ ] All 6 tools resolve (mafft, trimal, AMAS.py, iqtree, gotree, fastANI — fastANI only required if `ancillary/mitotree-ani` will run).
- [ ] `params.json` written with a non-empty `sample_labels` list.
- [ ] `preflight.md` verdict ≥ `GO-WITH-WARNINGS`.

## What NOT to do

- Do **not** scaffold a dating/MCMCtree stage even on repeated request (SP3) — this is a deliberate, explained scope exclusion.
- Do **not** silently drop rRNA in `aa` mode without mentioning it in `preflight.md` — it's correct behavior but should be visible to the user.

## Handoff

> Preflight complete. Verdict: `<verdict>`. `params.json` + `preflight.md` written. Proceed to `build/mitotree-builder`.
