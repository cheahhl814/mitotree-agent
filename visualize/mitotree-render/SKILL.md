---
name: mitotree-render
description: 'Render an offline, scriptable publication-format tree figure (SVG/PDF/PNG) from build/mitotree-builder''s species_tree.nwk, using toytree. The tree is already rooted by build/mitotree-builder -- --root-outgroup is available only as an optional override for exploring a different rooting, not required. Complements visualize/mitotree-itol (interactive, true italics). Triggers: "render the mitotree figure", "publication tree image", "export mitotree as PDF/SVG/PNG", "mitotree-render".'
version: 1.2.0
updated: "2026-11-16"
triggers:
  - "render the mitotree figure"
  - "publication tree image"
  - "export mitotree as PDF"
  - "export mitotree as SVG"
  - "mitotree-render"
disable-model-invocation: true
---

# mitotree-render

> Optional visualization phase of `mitotree-agent`, run after `qc/mitotree-qc` (or any time after `build/mitotree-builder` completes), alongside `visualize/mitotree-itol`. Ported from `phylogenetics-agent/visualize/phylo-render` — `render_tree.py` is copied unchanged (generic toytree wrapper). **Simplified vs. the original**: `phylo-render` roots on the outgroup itself because `phylo-builder` doesn't. `build/mitotree-builder` already writes a rerooted `species_tree.nwk`, so `--root-outgroup` here is only useful as an optional override (e.g. exploring an alternative rooting), not a required step.

## When to use this vs. `mitotree-itol`

| Need | Use |
| --- | --- |
| Interactive exploration, true italic species names, custom colorstrips | `visualize/mitotree-itol` (uploads to itol.embl.de) |
| Offline/scriptable figure, CI pipeline, no internet/browser dependency | `mitotree-render` (this sub-skill) |
| Both | Run both — they don't conflict, write to different subdirectories |

**Known limitation** (inherited from `phylo-render`, same toytree/toyplot version): tip labels are **not italicized**. toyplot's CSS-style validator rejects `font-style`. If italic species names are required, use `mitotree-itol` instead.

## 0. Inputs / Outputs contract

### Inputs

| Path | Source | Required? |
| --- | --- | --- |
| `$RUN_DIR/checkpoints/reroot.done.json` | `build/mitotree-builder` | yes |
| `$RUN_DIR/build/species_tree.nwk` | `build/mitotree-builder` | yes — already rooted |

### Outputs

| Path | Format | Notes |
| --- | --- | --- |
| `$RUN_DIR/figures/species_tree.svg` | SVG | vector, always included by default |
| `$RUN_DIR/figures/species_tree.pdf` | PDF | vector, always included by default |
| `$RUN_DIR/figures/species_tree.png` | PNG | raster, always included by default |

## 0.5 Ask-User Stop Points

### SP1 — No support values found on the tree

| Trigger | Evidence | Action |
| --- | --- | --- |
| Every internal node's `support` field is `NaN` (script prints a warning, not an error) | Small taxon count producing an uninformative/near-empty support tree (see `build/mitotree-builder` SP2), or the wrong tree file was passed | Ask: "No bootstrap support values found on this tree — figure will render without them. Pick: (A) proceed anyway (topology-only figure), (B) check `checkpoints/tree.done.json` — was `-B 1000` skipped per `build/mitotree-builder` SP2?, (C) abort" |

**Auto-pick when**: support values are present. No ask.

### SP2 — `--root-outgroup` override disagrees with the build's own rooting

| Trigger | Evidence | Action |
| --- | --- | --- |
| User explicitly passes `--root-outgroup` with a different taxon than `params.json`'s `outgroup` | The figure would then be rooted differently than `qc/mitotree-qc`'s report and `mitotree-itol`'s upload | Ask: "This will root the figure differently from the tree `build/mitotree-builder` produced. Pick: (A) proceed (exploratory alternative rooting, figure-only), (B) use the build's own rooting instead (omit `--root-outgroup`)" |

**Auto-pick when**: `--root-outgroup` omitted, or matches `params.json`'s `outgroup`. No ask.

## Procedure

```bash
mkdir -p "$RUN_DIR/figures"
pixi run python3 visualize/mitotree-render/scripts/render_tree.py "$RUN_DIR/build/species_tree.nwk" \
    -o "$RUN_DIR/figures/species_tree" \
    --formats svg,pdf,png
```

Useful flags: `--no-support-labels` (topology-only), `--no-ladderize`, `--width`/`--height`, `--font-size`, `--root-outgroup NAME[,NAME...]` (optional override only — see SP2, the tree is already rooted).

## Interpretation Guidelines

- **Support values shown, tree ladderized, all 3 formats written**: normal success.
- **Warning "no support values found"**: route to SP1.

## Troubleshooting — Signature library

| Signature | Likely cause | Fix |
| --- | --- | --- |
| `ModuleNotFoundError: toytree` | pixi env not reinstalled after `toytree` was added to `pixi.toml` | `pixi install`. |
| `ModuleNotFoundError: packaging` (during PNG export specifically) | toyplot's reportlab-based PNG backend needs `packaging`, not pulled in transitively by `toytree` alone | Already pinned in `pixi.toml`; if seen, re-run `pixi install`. |
| `ValueError: Not an allowed CSS style: font-style` | Something tried to pass `tip_labels_style={"font-style": ...}` | Don't — toyplot rejects it. Use `mitotree-itol` for italics instead. |
| `WARNING: no support values found` | Route to SP1 | See SP1. |
| Figure is very cramped / tip labels overlap | `--height` too small for the tip count | Default height is `28px * ntips` (min 300) — for very large trees, raise `--height` further or drop `--width`. |
| `--root-outgroup '<X>' matched no tip label` | Substring match against tip labels found nothing — typo | Check `gotree stats tips -i build/species_tree.nwk` for the real name; figure still renders, just with the build's original rooting unchanged. |

## Verification

- [ ] `$RUN_DIR/figures/species_tree.{svg,pdf,png}` all exist.
- [ ] Support labels present unless the user explicitly asked for `--no-support-labels`.

## What NOT to do

- Do **not** treat `--root-outgroup` as required — `build/mitotree-builder` already rooted the tree; only pass it to explore a *different* rooting (route through SP2 first).
- Do **not** try to force italic tip labels via `font-style` — use `mitotree-itol` instead.

## Handoff

> Render complete. `figures/species_tree.{svg,pdf,png}` written. For interactive exploration or true italic tip labels, also run `visualize/mitotree-itol`.
