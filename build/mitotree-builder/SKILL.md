---
name: mitotree-builder
description: 'Build the mitogenome supermatrix (extract PCGs+rRNA, align with mafft, trim with trimAl, concatenate with AMAS) and species tree (IQ-TREE, GTR+G or LG+G, 1000 UFBoot), then optionally re-root on an outgroup via ape. Checkpointed across 3 stages for resume. Triggers: internal — invoked by mitotree-agent orchestrator only.'
version: 1.2.0
updated: "2026-11-16"
triggers:
  - "mitotree build"
disable-model-invocation: true
---

# mitotree-builder

> Phase 2 of `mitotree-agent`. Wraps `build_supermatrix.py` (extract → align → trim → concat → stats, live-verified against a synthetic 4-taxon/3-gene dataset), `iqtree` (species tree), and `reroot_tree.R` (ape-based re-rooting, live-verified).

## 0. Inputs / Outputs contract

### Inputs

| Path | Source | Required? |
| --- | --- | --- |
| `$RUN_DIR/params.json` | `preflight/mitotree-preflight` | yes |
| `$RUN_DIR/preflight.md` verdict ≥ `GO-WITH-WARNINGS` | `preflight/mitotree-preflight` | yes |

### Outputs

| Path | Format | Notes |
| --- | --- | --- |
| `$RUN_DIR/build/04_concatenated/supermatrix.fasta` + `partitions.nex` | FASTA + NEXUS | from `build_supermatrix.py` |
| `$RUN_DIR/build/05_stats/{per_gene_stats.txt,supermatrix_stats.txt}` | Text (AMAS summary) | consumed by `qc/mitotree-qc` |
| `$RUN_DIR/build/species_tree.raw.treefile` | Newick | raw IQ-TREE output, UFBoot support |
| `$RUN_DIR/build/species_tree.nwk` | Newick | final tree: re-rooted if `params.json.outgroup` set, else a copy of the raw tree |
| `$RUN_DIR/checkpoints/{supermatrix,tree,reroot}.done.json` | JSON | resume markers |

## 0.5 Ask-User Stop Points

### SP1 — Resume vs. rebuild when checkpoints already exist

| Trigger | Evidence | Action |
| --- | --- | --- |
| `checkpoints/tree.done.json` exists but `params.json` has changed since (different `mode`/`gene_order`/`outgroup`) | Compare `params.json` mtime/hash against the checkpoint's recorded params hash | Ask: "(A) resume from the existing checkpoint (ignore params.json changes), (B) rebuild from the changed stage onward, (C) rebuild everything" |

**Auto-pick when**: `params.json` unchanged since the newest checkpoint → resume automatically, skip completed stages.

### SP2 — IQ-TREE bootstrap replicate count for small datasets

| Trigger | Evidence | Action |
| --- | --- | --- |
| Sample count < 5 | `params.json.sample_labels` length | **Recommend**: UFBoot support values on trees this small are not very informative (verified live: a 4-taxon synthetic run produced a single internal node at 32% support — expected, not a bug, since there is exactly one possible unrooted topology's worth of signal to test). Ask: "(A) proceed with default `-B 1000`, (B) skip bootstrap (`-B` omitted) since the topology is essentially fixed at this taxon count" |

**Auto-pick when**: ≥5 samples → proceed with `-B 1000`, no ask.

## Procedure

### 1. Build the supermatrix (checkpoint: `supermatrix`)

```bash
run_dir="$RUN_DIR"
mode=$(python3 -c "import json;print(json.load(open('$run_dir/params.json'))['mode'])")

python3 build/mitotree-builder/scripts/build_supermatrix.py \
    $(python3 -c "import json;print(' '.join(json.load(open('$run_dir/params.json'))['source_dirs']))") \
    $(python3 -c "import json;print(' '.join('--fasta ' + s for s in json.load(open('$run_dir/params.json'))['fasta_specs']))") \
    --mode "$mode" \
    -o "$run_dir/build" \
    --mafft "pixi run mafft" --trimal "pixi run trimal" \
    --amas "$(pixi run which AMAS.py)"
```

`build_supermatrix.py` takes `source_dirs` as `nargs="*"` — an empty list is fine when every sample comes via `--fasta`.

### 2. Infer the species tree (checkpoint: `tree`)

```bash
model_flag=$([ "$mode" = aa ] && echo "-m LG+G -B 1000" || echo "-m GTR+G -B 1000")
pixi run iqtree -s "$run_dir/build/04_concatenated/supermatrix.fasta" \
    -p "$run_dir/build/04_concatenated/partitions.nex" \
    $model_flag -T AUTO --quiet
cp "$run_dir/build/04_concatenated/partitions.nex.treefile" "$run_dir/build/species_tree.raw.treefile"
```

### 3. Re-root on the outgroup, if given (checkpoint: `reroot`)

```bash
outgroup=$(python3 -c "import json;print(json.load(open('$run_dir/params.json')).get('outgroup') or '')")
if [ -n "$outgroup" ]; then
    Rscript build/mitotree-builder/scripts/reroot_tree.R \
        -t "$run_dir/build/species_tree.raw.treefile" \
        -o "$outgroup" \
        -w "$run_dir/build/species_tree.nwk"
else
    cp "$run_dir/build/species_tree.raw.treefile" "$run_dir/build/species_tree.nwk"
fi
```

## Interpretation Guidelines

- **AMAS "Missing_percent" > 20%**: check `05_stats/per_gene_stats.txt` for which gene is dragging it up — usually one incompletely-annotated sample, not a systemic problem.
- **UFBoot support all low (<50) on a small dataset**: expected at low taxon counts (verified live), not necessarily a bad tree — cross-check against `concordance/mitotree-concordance` (gCF/sCF) once genes are more numerous than 3-4.
- **Re-root fails with "outgroup not found in tree"**: the outgroup label in `params.json` doesn't match a tip label in the supermatrix — tip labels come from `--fasta LABEL:FILE`'s LABEL or the archive filename stem, not any FASTA header inside.

## Troubleshooting — Signature library

| Signature | Likely cause | Fix |
| --- | --- | --- |
| `AMAS.py not found` | `--amas` path stale or pixi env not installed | Re-run `pixi install`; re-resolve via `pixi run which AMAS.py`. |
| `ERROR in mafft <gene>` | `--mafft` was a bare string like `"pixi run mafft"` passed to a wrapper expecting a single binary path (pre-fix behavior) | Fixed as of v1.0.0 — `step_align`/`step_trim` split the command into argv-list form. If seen, confirm the script wasn't reverted to string-concatenation. |
| `outgroup '<name>' not found in tree` | Outgroup label typo, or the label refers to a FASTA header rather than the `--fasta LABEL:FILE` label / archive stem | Check tip labels with `gotree stats tips -i species_tree.raw.treefile`. |
| IQ-TREE exits nonzero on a re-run of the same output dir | Stale `.ckp.gz` checkpoint from a previous partial/failed run conflicts with a fresh `-s`/`-p` invocation (verified live) | `rm 04_concatenated/*.ckp.gz` before re-running Step 2 if the params (model, partitions) changed. |
| Supermatrix has 0 parsimony-informative sites | Genes too short / too conserved for the taxon set (common on tiny test data — not a real-data concern at typical mitogenome scales) | Expected on synthetic/smoke-test data; on real mitogenomes this would indicate too-shallow a taxon sample for the marker. |

## Verification

- [ ] `checkpoints/supermatrix.done.json`, `checkpoints/tree.done.json`, `checkpoints/reroot.done.json` all exist.
- [ ] `build/species_tree.nwk` parses (`gotree stats tips -i <file>` succeeds).
- [ ] `build/05_stats/supermatrix_stats.txt` shows the expected taxon count.

## What NOT to do

- Do **not** add a divergence-time dating stage here or anywhere in this skill — see `preflight/mitotree-preflight` SP3. This is the tree-building phase only.
- Do **not** pass `--mafft`/`--trimal` as bare command strings without splitting them into argv form first — `subprocess.run` with a space-containing single string as `argv[0]` fails with `FileNotFoundError` (verified live during this skill's own build).

## Handoff

> Build complete. `build/species_tree.nwk` written (`<N>` taxa). Proceed to `qc/mitotree-qc`, or to the optional `concordance/mitotree-concordance` / `ancillary/mitotree-ani` / `ancillary/mitotree-codon-usage` branches first if the user wants those folded into the QC report.
