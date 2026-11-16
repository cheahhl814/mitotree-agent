---
name: mitotree-prune
description: 'Prune the final tree from build/mitotree-builder down to a curated subset of taxa -- keep list, exclude list, or diff against a reference pruned tree -- via gotree prune, preserving branch lengths and support values. Ported unchanged from phylogenetics-agent''s prune/phylo-prune (a generic Newick-tree wrapper). Triggers: "prune the mitotree", "remove these taxa", "drop misidentified mitogenomes", "mitotree-prune".'
version: 1.2.0
updated: "2026-11-16"
triggers:
  - "prune the mitotree"
  - "remove these taxa"
  - "drop misidentified mitogenomes"
  - "exclude these sequences from the mitotree"
  - "mitotree-prune"
disable-model-invocation: true
---

# mitotree-prune

> Optional stage of `mitotree-agent`, run after `build/mitotree-builder` (before or alongside `visualize/mitotree-itol` / `visualize/mitotree-render`, both of which can consume the pruned tree instead of `build/species_tree.nwk`). `scripts/prune_tree.py` is copied unchanged from `phylogenetics-agent/prune/phylo-prune` — it's a generic `gotree prune` wrapper with no single-marker-specific logic, so it works on `build/species_tree.nwk` exactly as it did on `phylo-builder`'s output.

## When to Use This Skill

Use after `build/mitotree-builder` has produced `build/species_tree.nwk`, when the user wants to:
- Drop misidentified or low-quality mitogenome samples before publishing.
- Exclude specific taxa.
- Cut the tree down to a curated subset, matching an existing reference pruned tree from a related analysis.

Do NOT use this to change which sequences go *into* the supermatrix build — that's a `build/mitotree-builder` re-run with a filtered input set (fewer `--fasta`/archive sources). This skill only trims tips from an already-built tree.

## 0. Inputs / Outputs contract

### Inputs

| Path | Source | Required? |
| --- | --- | --- |
| `$RUN_DIR/checkpoints/reroot.done.json` | `build/mitotree-builder` | yes |
| `$RUN_DIR/build/species_tree.nwk` | `build/mitotree-builder` | yes |
| Exactly one of: `--keep NAME...`, `--keep-file FILE`, `--exclude-file FILE`, `--reference-tree FILE` | User | yes |

### Outputs

| Path | Format | Notes |
| --- | --- | --- |
| `$RUN_DIR/pruned/species_tree_pruned.nwk` | Newick | branch lengths + support values preserved |
| `$RUN_DIR/pruned/species_tree_pruned.nwk.report.txt` | Text | tip counts before/after, removed names, any requested names not found |

## 0.5 Ask-User Stop Points

### SP1 — Determine the keep/exclude mode

| Trigger | Evidence | Action |
| --- | --- | --- |
| User said "prune the tree" without specifying how | No default mode makes sense — keep vs. exclude are opposite operations | Ask: "How should I select taxa? (A) give me a list of names to **keep**, (B) give me a list of names to **exclude** (e.g. misidentified mitogenomes), (C) I have a reference pruned tree from a related analysis — match its tip set" |

**Auto-pick when**: the user's request already names the taxa and the verb ("keep only X, Y, Z" vs "remove X, Y, Z") — infer `--keep`/`--exclude-file` directly, no need to ask.

### SP2 — Requested name not found in the tree

| Trigger | Evidence | Action |
| --- | --- | --- |
| `prune_tree.py` exits with `N requested name(s) not found` (gotree itself silently ignores unknown names rather than erroring, so this check runs in the wrapper before calling gotree) | The name doesn't match any tip in `gotree stats tips` | **Recommend**: mitogenome tip names come from `fetch_mitogenomes.py`'s `Genus_species_ACCESSION` convention or the sample label passed at preflight — check `gotree stats tips -i build/species_tree.nwk` for the real name first. Ask: "(A) fix the name and re-run, (B) proceed anyway with `--force` (missing names are simply ignored, not an error)" |

**Auto-pick when**: zero missing names. No ask.

### SP3 — Result would leave too few tips

| Trigger | Evidence | Action |
| --- | --- | --- |
| gotree stderr contains `only made of two tips` / `only made of one tip` (gotree hard-refuses to prune below 3 tips) | The keep/exclude set is too aggressive for a meaningful tree | Ask: "Pruning to this set leaves fewer than 3 tips, which gotree refuses to produce. Which taxa should stay in, at minimum?" |

**Auto-pick when**: n/a — this always requires user input.

### SP4 — Overwrite existing pruned output

| Trigger | Evidence | Action |
| --- | --- | --- |
| `$RUN_DIR/pruned/species_tree_pruned.nwk` already exists | A previous prune run is present | Ask before overwriting: "(A) overwrite, (B) write to a new name (e.g. `_pruned_v2.nwk`)" |

**Auto-pick when**: no existing file at that path.

## Procedure

```bash
cd /path/to/mitotree-agent   # skill root, so `pixi run gotree` resolves
mkdir -p "$RUN_DIR/pruned"

# Keep mode
python3 prune/mitotree-prune/scripts/prune_tree.py "$RUN_DIR/build/species_tree.nwk" \
    --keep Taxon_a Taxon_b Taxon_c \
    -o "$RUN_DIR/pruned/species_tree_pruned.nwk" \
    --gotree-cmd "pixi run gotree"

# Exclude mode (misidentified mitogenomes, duplicates flagged by ancillary/mitotree-ani, etc.)
python3 prune/mitotree-prune/scripts/prune_tree.py "$RUN_DIR/build/species_tree.nwk" \
    --exclude-file /path/to/exclude_list.txt \
    -o "$RUN_DIR/pruned/species_tree_pruned.nwk" \
    --gotree-cmd "pixi run gotree"

# Match a reference pruned tree from a related analysis
python3 prune/mitotree-prune/scripts/prune_tree.py "$RUN_DIR/build/species_tree.nwk" \
    --reference-tree /path/to/reference_pruned.nwk \
    -o "$RUN_DIR/pruned/species_tree_pruned.nwk" \
    --gotree-cmd "pixi run gotree"
```

## Interpretation Guidelines

- **Tip count drops exactly to the expected size, 0 missing names**: clean prune, no follow-up needed.
- **Some requested names missing**: route through SP2 before proceeding — likely a typo or a mismatch against the actual sample label used at preflight/taxon-sampling.
- **gotree refuses (<3 tips)**: route through SP3.
- **Excluding ANI-flagged near-identical pairs**: `ancillary/mitotree-ani`'s `near_identical_pairs` list is a natural `--exclude-file` source — drop the lower-quality member of each pair.

## Troubleshooting — Signature library

| Signature | Likely cause | Fix |
| --- | --- | --- |
| `N requested name(s) not found in <tree>` | Typo, or the name doesn't match the exact label used when the sample was added (`--fasta LABEL:...`, or `fetch_mitogenomes.py`'s `Genus_species_ACCESSION` naming) | Run `gotree stats tips -i build/species_tree.nwk` first and match the exact name; re-run with the corrected name, or `--force` to ignore. |
| `the requested keep/exclude set leaves too few tips` | gotree refuses trees with < 3 tips | Keep at least 3 taxa. |
| `gotree: command not found` | Script invoked without `--gotree-cmd "pixi run gotree"` outside an activated pixi shell | Always pass `--gotree-cmd`, or run inside `pixi shell`. |
| Pruned tree's branch lengths look longer than expected on a formerly-single-child node | Expected, not a bug — gotree sums the collapsed node's branch length into the remaining edge | No fix needed; standard phylogenetic practice. |
| `ERROR: gotree stats tips failed` | Input tree path wrong, or `build/mitotree-builder` hasn't produced it yet | Check `checkpoints/reroot.done.json` exists first. |

## Verification

- [ ] `$RUN_DIR/pruned/species_tree_pruned.nwk` exists and parses (`gotree stats tips -i <file>` succeeds).
- [ ] Report's "tips after" count matches what the user asked for.
- [ ] No unexpected names in "removed" list.

## What NOT to do

- Do **not** hand-roll a Newick parser for this — `gotree` already handles it, same as `phylo-prune`.
- Do **not** silently proceed past SP2 without asking.

## Handoff

> Prune complete. `pruned/species_tree_pruned.nwk` written (`<N>` → `<M>` tips). To visualize the pruned tree instead of the full one, point `visualize/mitotree-itol` or `visualize/mitotree-render` at this file instead of `build/species_tree.nwk`.
