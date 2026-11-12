---
name: mitotree-itol
description: 'Generate iTOL (itol.embl.de) annotation datasets from the final tree produced by build/mitotree-builder: italic tip-label styling (always) and a colorstrip grouping taxa by genus/clade/any user-defined grouping (when a group source is given). Unlike phylogenetics-agent''s phylo-itol, no separate reroot step is needed here -- build/mitotree-builder already writes a rerooted species_tree.nwk. Tip names come from gotree stats tips, not regex Newick parsing. Triggers: "visualize the mitotree", "iTOL annotation", "colorstrip", "mitotree-itol".'
version: 1.0.0
updated: "2026-08-26"
triggers:
  - "visualize the mitotree"
  - "iTOL annotation"
  - "colorstrip"
  - "style the tree labels"
  - "mitotree-itol"
disable-model-invocation: true
---

# mitotree-itol

> **v1.0.0.** Optional visualization phase of `mitotree-agent`, run after `qc/mitotree-qc` (or any time after `build/mitotree-builder` completes). Ported from `phylogenetics-agent/visualize/phylo-itol` — the wrapper scripts (`itol_tip_labels.py`, `itol_colorstrip.py`) are generic Newick-tree tools, copied unchanged. **Simplified vs. the original**: `phylo-itol` reroots on the outgroup itself (`gotree reroot outgroup`) because `phylo-builder` doesn't root its tree. `mitotree-builder` already does that rooting at build time (its own `reroot_tree.R` step) and writes the result straight to `build/species_tree.nwk` — so this sub-skill just points at that file directly, no separate reroot step.

## When to Use This Skill

Use after `qc/mitotree-qc` (or right after `build/mitotree-builder`), when the user wants:
- Italicized tip labels for a publication-style tree (always available, no extra input needed).
- A colored strip grouping taxa by genus, clade, sample source, or any other tip-name-derived or user-supplied category.

Do NOT use if the user just wants the plain-text ASCII tree from `qc/mitotree-qc` — that's already sufficient for a quick look.

## 0. Inputs / Outputs contract

### Inputs

| Path | Source | Required? |
| --- | --- | --- |
| `$RUN_DIR/checkpoints/reroot.done.json` | `build/mitotree-builder` | yes |
| `$RUN_DIR/build/species_tree.nwk` | `build/mitotree-builder` | yes — already rooted if an outgroup was set, no separate rerooting needed here |
| A group source: `--groups-tsv <tip\tgroup file>` OR `--regex <pattern>` | User | only for the colorstrip; labels-only doesn't need one |

### Outputs

| Path | Format | Notes |
| --- | --- | --- |
| `$RUN_DIR/itol/labels.txt` | iTOL `TREE_COLORS` dataset | Always generated — italicizes every tip label |
| `$RUN_DIR/itol/colorstrip.txt` | iTOL `DATASET_COLORSTRIP` dataset | Only if a group source was given |
| `$RUN_DIR/itol/legend.md` | Markdown | group → color → tip-count table |

Upload `build/species_tree.nwk` and these `.txt` datasets together at https://itol.embl.de/upload.cgi (drag the tree first, then drag each dataset file onto the rendered tree).

## 0.5 Ask-User Stop Points

### SP1 — No group source given for the colorstrip

| Trigger | Evidence | Action |
| --- | --- | --- |
| User asked for a colorstrip/grouped visualization but gave neither `--groups-tsv` nor `--regex` | Colorstrip needs a group source, there's no sensible default | Ask: "How should tips be grouped for the colorstrip? Pick: (A) I'll give a regex (e.g. `^([A-Za-z]+)_` to group by genus, matching this skill's `Genus_species_ACCESSION` tip-naming convention), (B) I'll supply a `tip_name<TAB>group` TSV file, (C) skip the colorstrip — labels-only" |

**Auto-pick when**: user only asked for label styling (no grouping mentioned). Default: labels-only, skip colorstrip.

### SP2 — Regex leaves tips unmatched

| Trigger | Evidence | Action |
| --- | --- | --- |
| `itol_colorstrip.py`'s stderr reports `N/M tips have no group` | Some tip names don't match the given pattern | Ask: "`<N>` of `<M>` tips didn't match the regex and were grouped as `unclassified` (grey). Pick: (A) accept as-is, (B) give a corrected regex and re-run, (C) switch to a `--groups-tsv` file for full control" |

**Auto-pick when**: 0 unmatched tips. No ask.

### SP3 — Color override conflict

| Trigger | Evidence | Action |
| --- | --- | --- |
| `--colors-tsv` assigns a color to a group name that doesn't appear in the actual tip groups | Likely a typo in the group name | Ask: "`--colors-tsv` names group `<X>` but no tip was assigned to it. Pick: (A) fix the group name, (B) ignore (color unused, no harm), (C) abort" |

**Auto-pick when**: every override group name is actually used. No ask.

## Procedure

### 1. Always: generate the italic-label dataset

```bash
cd /path/to/mitotree-agent   # skill root, so `pixi run gotree` resolves
mkdir -p "$RUN_DIR/itol"
python3 visualize/mitotree-itol/scripts/itol_tip_labels.py "$RUN_DIR/build/species_tree.nwk" \
    -o "$RUN_DIR/itol/labels.txt" \
    --gotree-cmd "pixi run gotree"
```

### 2. If a group source was given: generate the colorstrip

```bash
# Regex mode (SP1 option A):
python3 visualize/mitotree-itol/scripts/itol_colorstrip.py "$RUN_DIR/build/species_tree.nwk" \
    -o "$RUN_DIR/itol/colorstrip.txt" \
    --regex '^([A-Za-z]+)_' --label "Genus" \
    --legend-out "$RUN_DIR/itol/legend.md" \
    --gotree-cmd "pixi run gotree"

# TSV mode (SP1 option B):
python3 visualize/mitotree-itol/scripts/itol_colorstrip.py "$RUN_DIR/build/species_tree.nwk" \
    -o "$RUN_DIR/itol/colorstrip.txt" \
    --groups-tsv /path/to/groups.tsv --label "<group label>" \
    --legend-out "$RUN_DIR/itol/legend.md" \
    --gotree-cmd "pixi run gotree"
```

## Interpretation Guidelines

- **Labels dataset**: always trivially "correct" — if it's empty, the tree file itself has 0 tips (upstream bug, not this sub-skill's).
- **Colorstrip with 0 unmatched tips**: every tip got a real group — good regex/TSV coverage.
- **Colorstrip with some `unclassified`**: still usable (rendered grey in iTOL), but the legend will show it — route through SP2 if the user cares about full coverage.

## Troubleshooting — Signature library

| Signature | Likely cause | Fix |
| --- | --- | --- |
| `ERROR: no tips found in <tree>` | Tree file path wrong, or `build/mitotree-builder` hasn't produced it yet | Check `checkpoints/reroot.done.json` exists first. |
| `gotree: command not found` | Script invoked without `--gotree-cmd "pixi run gotree"` outside an activated pixi shell | Always pass `--gotree-cmd`, or run inside `pixi shell`. |
| All tips end up `unclassified` | Regex doesn't match this tip-naming convention at all | Print a few tip names first (`gotree stats tips -i build/species_tree.nwk`) and adjust the regex against real examples before running the full colorstrip. |
| iTOL upload shows the colorstrip but all one grey color | Group names in `DATA` rows don't match any key iTOL recognizes — actually a cosmetic non-issue, iTOL derives the color from the 2nd column directly, not the group name; check the hex codes in the file are valid `#RRGGBB`. |
| Legend group counts don't sum to tip count | Should be impossible (every tip gets exactly one group) — if seen, the TSV had duplicate tip-name rows; the last one silently wins. | Deduplicate the `--groups-tsv` file. |
| Tree looks unrooted in iTOL even though an outgroup was given at preflight | `build/mitotree-builder`'s reroot stage didn't run — check `checkpoints/reroot.done.json` and `params.json`'s `outgroup` field | This is a `build/mitotree-builder` issue, not `mitotree-itol`'s — see that sub-skill's troubleshooting table, not this one. |

## Verification

- [ ] `$RUN_DIR/itol/labels.txt` exists, starts with `TREE_COLORS`.
- [ ] If a colorstrip was requested: `$RUN_DIR/itol/colorstrip.txt` exists, starts with `DATASET_COLORSTRIP`, and `legend.md` group counts sum to the tree's tip count.

## What NOT to do

- Do **not** hardcode any project-specific taxonomy into this sub-skill — that belongs in the user's own `--groups-tsv` file, kept outside the skill repo.
- Do **not** parse the Newick file with a hand-rolled regex — always shell out to `gotree stats tips`.
- Do **not** re-implement outgroup rerooting here — `build/mitotree-builder` already does it once, upstream. Rerooting it again here would be redundant and risks disagreeing with the tree `qc/mitotree-qc` already reported on.

## Handoff

> Visualization complete. `itol/labels.txt` written (plus `itol/colorstrip.txt` and `itol/legend.md` if a colorstrip was requested). Upload `build/species_tree.nwk` and these dataset(s) together at https://itol.embl.de/upload.cgi.
