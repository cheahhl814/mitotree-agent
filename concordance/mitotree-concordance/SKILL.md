---
name: mitotree-concordance
description: 'Assess gene-tree/species-tree concordance: infer per-gene trees (IQ-TREE, -m MFP -B 1000), compute Robinson-Foulds distances (gotree compare trees --rf, normalized per gene), and calculate gCF/sCF concordance factors (IQ-TREE --gcf/--scf). Writes concordance/rf_distances.tsv + concordance.cf.stat optionally plot a gCF-vs-sCF scatter figure, and write concordance_scatter.pdf/.png. Outputs feed qc/mitotree-qc. Triggers: "gene tree concordance", "check concordance factors", "gCF/sCF", "concordance scatter", "RF distances", "mitotree-concordance".'
version: 1.1.0
updated: "2026-09-08"
triggers:
  - "gene tree concordance"
  - "check concordance factors"
  - "gCF/sCF"
  - "RF distances"
  - "mitotree-concordance"
---

# mitotree-concordance

> **v1.0.0.** Optional post-build stage of `mitotree-agent`, run after `build/mitotree-builder` completes (`checkpoints/tree.done.json` exists). Wraps IQ-TREE for gene-tree inference and concordance-factor calculation, and gotree for Robinson-Foulds distance computation — no hand-rolled Newick parsing needed, since gotree already handles tree comparison correctly and IQ-TREE natively outputs gCF/sCF in the exact format `qc/mitotree-qc/scripts/qc_report.py` expects.

## When to Use This Skill

Use after `build/mitotree-builder` has produced `build/species_tree.nwk` and per-gene trimmed alignments (`build/03_trimmed/*.trim.fasta`), when the user wants to:

- Assess how much individual gene trees disagree with the species tree (RF distances).
- Calculate gene concordance factors (gCF) and site concordance factors (sCF) for each internal branch.
- Identify outlier genes whose topologies are strongly discordant (potential paralogs, misannotations, or alignment artifacts).

Do NOT use this to **fix** discordance — it's a diagnostic stage. High RF or low gCF/sCF flags genes worth inspecting in `qc/mitotree-qc`, but the fix (drop the gene, re-align, re-annotate) is a `build/mitotree-builder` re-run with filtered inputs.

## 0. Inputs / Outputs contract

### Inputs

| Path                                               | Source                   | Required?                 |
| -------------------------------------------------- | ------------------------ | ------------------------- |
| `$RUN_DIR/build/species_tree.nwk`                  | `build/mitotree-builder` | yes                       |
| `$RUN_DIR/build/03_trimmed/<gene>.trim.fasta`      | `build/mitotree-builder` | yes (≥1 gene)             |
| `$RUN_DIR/build/04_concatenated/supermatrix.fasta` | `build/mitotree-builder` | yes (for sCF calculation) |

### Outputs

| Path                                         | Format                     | Notes                                                                                |
| -------------------------------------------- | -------------------------- | ------------------------------------------------------------------------------------ |
| `$RUN_DIR/concordance/gene_trees.nwk`        | Newick (one line per gene) | IQ-TREE gene trees (`-m MFP -B 1000 -T 1 --quiet`)                                   |
| `$RUN_DIR/concordance/rf_distances.tsv`      | TSV, header row            | Column 0 = gene name, column `normalised_rf` (exact name, hardcoded in qc_report.py) |
| `$RUN_DIR/concordance/concordance.cf.tree`   | Newick                     | Species tree with gCF/sCF node labels                                                |
| `$RUN_DIR/concordance/concordance.cf.branch` | TSV                        | Per-branch concordance statistics                                                    |
| `$RUN_DIR/concordance/concordance.cf.stat`   | TSV, `#`-commented         | Overall run summary; columns include `ID`, `Label` (UFBoot), `gCF`, `sCF`            |
| `$RUN_DIR/concordance/concordance_scatter.pdf` + `.png` | Figure (optional, user-triggered) | gCF-vs-sCF scatter — one point per internal branch, UFBoot-category colour, 33% random-expectation reference lines, node labels |

## 0.5 Ask-User Stop Points

### SP1 — Concordance on a tiny taxon set

| Trigger                                | Evidence                   | Action                                                                                                                                                                                                                                                                                                                         |
| -------------------------------------- | -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `params.json.sample_labels` length < 5 | Taxon count from preflight | **Recommend**: gCF/sCF/RF values on <5 taxa are statistically thin (verified live: a 4-taxon synthetic run produced a single internal node's worth of signal — expected, not a bug). Ask: "(A) proceed anyway (output will be well-formed, just not biologically informative), (B) skip concordance until more taxa are added" |

**Auto-pick when**: ≥5 samples → proceed, no ask.

### SP2 — Gene tree inference fails for one gene

| Trigger                                          | Evidence                                                      | Action                                                                                                                                                                                                                           |
| ------------------------------------------------ | ------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| IQ-TREE exits nonzero on one `<gene>.trim.fasta` | IQ-TREE stderr (e.g. "ERROR: no parsimony-informative sites") | That gene's alignment is too short/too conserved for tree inference. Ask: "(A) skip this gene and continue with the rest, (B) drop it from the analysis entirely (re-run build without it), (C) investigate the alignment first" |

**Auto-pick when**: all genes infer successfully. No ask.

### SP3 — Overwrite existing concordance output

| Trigger                                                                         | Evidence                              | Action                                                                                      |
| ------------------------------------------------------------------------------- | ------------------------------------- | ------------------------------------------------------------------------------------------- |
| `$RUN_DIR/concordance/rf_distances.tsv` or `concordance.cf.stat` already exists | A previous concordance run is present | Ask before overwriting: "(A) overwrite, (B) write to a new subdir (e.g. `concordance_v2/`)" |

**Auto-pick when**: no existing concordance outputs. No ask.

## Procedure

### 1. Infer gene trees (one IQ-TREE process per gene)

```bash
run_dir="$RUN_DIR"
mkdir -p "$run_dir/concordance"

# Clear any stale gene_trees.nwk from a previous run
rm -f "$run_dir/concordance/gene_trees.nwk"
touch "$run_dir/concordance/gene_trees.nwk"

for aln in "$run_dir"/build/03_trimmed/*.trim.fasta; do
    gene=$(basename "$aln" .trim.fasta)
    pixi run iqtree -s "$aln" -m MFP -B 1000 -T 1 --quiet
    # Append the treefile content (one line per gene)
    cat "$aln.treefile" >> "$run_dir/concordance/gene_trees.nwk"
done
```

### 2. Compute Robinson-Foulds distances (gotree compare trees --rf)

```bash
python3 concordance/mitotree-concordance/scripts/compute_concordance.py \
    --species-tree "$run_dir/build/species_tree.nwk" \
    --gene-trees "$run_dir/concordance/gene_trees.nwk" \
    --output-rf "$run_dir/concordance/rf_distances.tsv" \
    --gotree-cmd "pixi run gotree"
```

This wrapper:

- Parses the multi-line `gene_trees.nwk` file (one Newick per gene).
- For each gene tree: prunes the species tree to the gene tree's tip set via `gotree prune -c <gene_tree>`.
- Runs `gotree compare trees -i <pruned_species> -c <gene_tree> --rf` to get raw RF.
- Normalizes: `normalised_rf = raw_rf / (2 * (n_tips - 3))` where `n_tips` is the number of tips actually compared.
- Writes `rf_distances.tsv` with header row and `normalised_rf` column (exact name qc_report.py expects).

### 3. Calculate gCF/sCF concordance factors (IQ-TREE --gcf/--scf)

```bash
mode=$(python3 -c "import json;print(json.load(open('$run_dir/params.json'))['mode'])")
seqtype_flag=$([ "$mode" = aa ] && echo "--seqtype AA" || echo "")

cd "$run_dir/concordance"
pixi run iqtree -t "$run_dir/build/species_tree.nwk" \
    --gcf "$run_dir/concordance/gene_trees.nwk" \
    -s "$run_dir/build/04_concatenated/supermatrix.fasta" \
    --scf 100 \
    $seqtype_flag \
    -T AUTO \
    --prefix concordance \
    --quiet
```

IQ-TREE writes:

- `concordance.cf.tree` — species tree with gCF/sCF node labels.
- `concordance.cf.branch` — per-branch statistics (TSV).
- `concordance.cf.stat` — overall summary; tab-separated with `#`-prefixed comment lines, columns including `ID`, `Label` (UFBoot), `gCF`, `sCF`.

### 4. Optional: gCF-vs-sCF scatter plot (user-triggered)

Matplotlib port of the manuscript `Fig4_concordance_scatter.R` (the pixi env has no ggplot2; same Python-replacement convention as `ancillary/mitotree-ani`'s heatmap):

```bash
# Auto-label the N most discordant nodes (lowest gCF, then sCF)
pixi run python3 concordance/mitotree-concordance/scripts/plot_concordance_scatter.py \
    "$RUN_DIR/concordance/concordance.cf.stat" \
    -o "$RUN_DIR/concordance/concordance_scatter" --label-nodes 3

# Or manual clade labels: generate the template, fill it in, re-run
pixi run python3 concordance/mitotree-concordance/scripts/plot_concordance_scatter.py \
    "$RUN_DIR/concordance/concordance.cf.stat" \
    -o "$RUN_DIR/concordance/concordance_scatter" --labels-template-out "$RUN_DIR/concordance/scatter_labels.tsv"
# ... edit scatter_labels.tsv (node_id<TAB>clade_label) ...
pixi run python3 concordance/mitotree-concordance/scripts/plot_concordance_scatter.py \
    "$RUN_DIR/concordance/concordance.cf.stat" \
    -o "$RUN_DIR/concordance/concordance_scatter" --labels-tsv "$RUN_DIR/concordance/scatter_labels.tsv"
```

Writes `concordance_scatter.pdf` + `.png`. Live-verified against the real Betta waseri-complex dataset (`Figure_4_data/concordance.cf.stat`, 14 nodes): node values match the R original exactly; greedy label placement keeps every label inside the axes without overlap.

## Interpretation Guidelines

- **normalised_rf ≈ 0**: gene tree topology matches the species tree exactly (on the shared tip set).
- **normalised_rf > 0.5**: strong topological discordance — inspect the gene's alignment (`build/03_trimmed/<gene>.trim.fasta`) for alignment errors, or check if it's a paralog misidentified as an ortholog.
- **gCF < 50%**: fewer than half of the gene trees support this branch — could be incomplete lineage sorting (ILS) or genuine biological discordance.
- **sCF < 33%**: fewer than a third of the alignment sites support this branch — the branch may be an artifact of model misspecification or saturation.
- **UFBoot ≥95% but gCF < 50% or sCF < 33%**: "misleading" node — high bootstrap support does not mean the underlying gene trees or sites agree; this is the single most important reliability check IQ-TREE offers beyond bootstrap alone (see `qc/mitotree-qc` SP2).

## Troubleshooting — Signature library

| Signature                                                     | Likely cause                                                                                         | Fix                                                                                                         |
| ------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| `ERROR: no parsimony-informative sites` (IQ-TREE on one gene) | That gene's trimmed alignment is too short or too conserved for tree inference                       | Skip that gene (SP2), or drop it from the build and re-run with fewer genes.                                |
| `gotree: command not found`                                   | Script invoked without `--gotree-cmd "pixi run gotree"` outside an activated pixi shell              | Always pass `--gotree-cmd`, or run inside `pixi shell`.                                                     |
| `normalised_rf` column missing or named differently           | Wrapper wrote a different column name                                                                | The wrapper must write `normalised_rf` exactly — qc_report.py's `parse_rf_distances()` hardcodes this name. |
| IQ-TREE `--gcf` fails with "gene trees file not found"        | `gene_trees.nwk` path wrong or empty                                                                 | Check Step 1 completed; `gene_trees.nwk` should have one Newick line per gene.                              |
| `concordance.cf.stat` has no `gCF`/`sCF` columns              | IQ-TREE invoked without `--scf 100` or `--gcf`                                                       | Verify Step 3 command matches the Procedure exactly.                                                        |
| RF distances all zero                                         | Gene trees are identical to the species tree (rare on real data, common on tiny synthetic test data) | Expected on 4-taxon synthetic datasets with one possible unrooted topology — not a bug.                     |

## Verification

- [ ] `$RUN_DIR/concordance/gene_trees.nwk` exists and has one Newick line per gene.
- [ ] `$RUN_DIR/concordance/rf_distances.tsv` parses (`head` shows header + rows, `normalised_rf` column present).
- [ ] `$RUN_DIR/concordance/concordance.cf.stat` parses (comment lines starting with `#`, data rows with `ID`, `Label`, `gCF`, `sCF` columns).
- [ ] Re-run `qc/mitotree-qc` to fold the concordance section into `qc-report.md`.
- [ ] If the scatter was requested: `concordance_scatter.pdf`/`.png` exist and the console table's node values match `concordance.cf.stat`.

## What NOT to do

- Do **not** port the original Nextflow pipeline's R script (`plot_rf_distances.R`) — gotree already computes RF distances correctly, and the Python wrapper normalizes them per the same formula the R script used (verified against the R source: `norm_rf = rf / (2*(n-3))`).
- Do **not** skip the species-tree pruning step before comparing — gene trees may have fewer tips than the species tree if a gene was missing in some samples, so prune the species tree down to each gene tree's tip set first via `gotree prune -c <gene_tree>` (same pattern `prune/phylo-prune` uses for its `--reference-tree` mode).
- Do **not** run concordance before `build/mitotree-builder` completes — the species tree and per-gene alignments must exist first.

## Handoff

> Concordance complete. `concordance/rf_distances.tsv` + `concordance.cf.stat` written. `<N>` genes evaluated, `<M>` outlier genes (normalised_rf > median + 2×IQR), `<K>` misleading nodes (UFBoot≥95 but gCF<50 or sCF<33). Re-run `qc/mitotree-qc` to fold these sections into `qc-report.md`. Offer the gCF-vs-sCF scatter (step 4) if a manuscript-style figure would help.
