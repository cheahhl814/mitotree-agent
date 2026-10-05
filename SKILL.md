---
name: mitotree-agent
description: Agentic mitogenome phylogenomics pipeline (multi-gene supermatrix -- extract PCGs+rRNA, align, trim, concatenate, species tree) with optional gene-tree concordance, ANI, and codon-usage branches. Rebuilt from the MitoTree Nextflow pipeline as a checkpointed BettaMt-pattern skill. Deliberately excludes divergence-time dating (MCMCtree/PAML) -- unreliable from mitogenome-only data.
version: 1.1.0
updated: "2026-10-05"
triggers:
  - "mitogenome phylogenomics"
  - "mitotree"
  - "mitogenome supermatrix"
  - "build a mitogenome tree"
  - "gene concordance factors"
  - "mitotree agent"
---

# MitoTree Agent

> **v1.0.0.** Checkpointed, agent-driven mitogenome phylogenomics pipeline (extract, align, trim, concatenate, tree, re-root) with optional gene-tree concordance, ANI, and codon-usage branches.

## Audience

This skill serves two purposes:

- **AI Agents**: Triggered by the phrases above. The agent must follow the strict evidence chain and run the Go/No-Go gates.
- **Human Users**: Provides a transparent audit trail.

## When to Use This Skill

Use this skill when you need to:

- Build a multi-gene mitochondrial-genome species tree from MitoAnnotator archives or per-species gene FASTAs (13 PCGs + 12S/16S rRNA).
- Assess gene-tree/species-tree concordance (RF distances, gCF/sCF) or check for near-identical samples via ANI.
- Swap between nucleotide (`nt`, intragenus/species-level) and amino-acid (`aa`, family-level+) supermatrix modes.

**Do NOT use this skill** if:

- You need a single-marker tree from one gene (16S, COI, Cytb alone) — use `phylogenetics-agent` instead; this skill is multi-gene supermatrix only.
- You need divergence-time estimates (node ages, a time-calibrated tree). This is a **deliberate scope exclusion**: mitogenome-only data does not reliably calibrate node ages (rate heterogeneity across the 13 PCGs + rRNAs, saturation at deeper splits, no independent nuclear-rate cross-check). Use a nuclear multi-locus or genome-scale dataset with MCMCtree instead.
- You want the original rigid Nextflow DAG for HPC cohort batch runs with no agent involvement — use `MitoTree/main.nf` directly instead.

## 0. Orchestrator — detect stage, route to the right sub-skill

This meta-skill is a **router**, not a doer. It does not duplicate logic from the sub-skills. Its job is to ask: **"what stage is the user at, and which sub-skill should they invoke next?"**

### 0.1 Locate the run directory

By convention the agent writes handoff files (see §B) to a run directory. Default: the current working directory. Override with `RUN_DIR` env var.

### 0.2 Detect the user's stage

```bash
# Stage detection ladder — first match wins
test -f "$RUN_DIR/qc-report.md" && STAGE="qc-done"
test -f "$RUN_DIR/checkpoints/reroot.done.json" && STAGE="qc"
test -f "$RUN_DIR/checkpoints/tree.done.json" && STAGE="build"
test -f "$RUN_DIR/checkpoints/supermatrix.done.json" && STAGE="build"
test -f "$RUN_DIR/params.json" && STAGE="build"
: "${STAGE:=preflight}"
```

If auto-detection is ambiguous, ask the user one short question (see **SP0** below).

### 0.5 Master ask-user stop points

This orchestrator has **one** user-facing stop point (SP0). All other stop points live in the sub-skills.

#### SP0 — Entry-stage ambiguity

| Trigger | Ask |
|---|---|
| Stage detection returns `preflight` AND no input sources given | "I don't see MitoAnnotator archives or `--fasta LABEL:FILE` specs. What are the sample inputs?" |

**Auto-pick when**: input sources are already given → default to `preflight` stage.

### 0.3 Route to the right sub-skill

| Stage | Action |
| --- | --- |
| `preflight` | Invoke `preflight/mitotree-preflight`. Validates inputs, gene coverage, tool audit. Writes `preflight.md` + `params.json`. |
| `build` | Invoke `build/mitotree-builder`. Extract → align → trim → concat → species tree → re-root, checkpointed (3 stages). |
| `qc` | Invoke `qc/mitotree-qc`. Builds `qc-report.md` from `build/` + any optional branches already run. |
| `qc-done` | Show the user `$RUN_DIR/qc-report.md`. If `concordance/mitotree-concordance`, `ancillary/mitotree-ani`, or `ancillary/mitotree-codon-usage` weren't run, offer them — then re-run `qc/mitotree-qc` to fold the new section in. Also offer `prune/mitotree-prune`, `visualize/mitotree-itol`, or `visualize/mitotree-render`. |

**`concordance/mitotree-concordance`, `ancillary/mitotree-ani`, and `ancillary/mitotree-codon-usage` are optional and user-triggered**, not part of the Preflight→Build→QC gate chain. Each can run any time after `build/mitotree-builder` completes (`checkpoints/tree.done.json` exists) — before or after `qc/mitotree-qc`. If run after `qc-report.md` already exists, re-run `qc/mitotree-qc` so its section is no longer "not found." They don't conflict with each other or gate one another. See their respective `SKILL.md` files.

**`prune/mitotree-prune`, `visualize/mitotree-itol`, and `visualize/mitotree-render` are optional and user-triggered**, also not part of the gate chain. All three operate on `build/species_tree.nwk` — which `build/mitotree-builder` already roots on the outgroup, if one was set, so (unlike `phylogenetics-agent`'s equivalents) neither visualize sub-skill needs its own rerooting step. Run `mitotree-prune` first if the user wants the visualization built on a curated taxon subset instead of the full tree. Ported unchanged/lightly-adapted from `phylogenetics-agent`'s `prune/phylo-prune`, `visualize/phylo-itol`, `visualize/phylo-render` — see their respective `SKILL.md` files.

**Do not skip preflight.** The preflight sub-skill computes evidence and writes `params.json` that the build sub-skill consumes.

**Do not add a dating stage.** See "Do NOT use this skill" above and `preflight/mitotree-preflight` SP3 — this is a deliberate, explained scope exclusion, not a missing feature to fill in later.

**`preflight/mitotree-taxon-sampling` is optional and user-triggered**, invoked *before* `preflight/mitotree-preflight` when the user wants the dataset broadened with more taxa from NCBI. Fetches whole complete-mitogenome FASTA (most NCBI submissions aren't split into per-gene records), then BLASTN-extracts the 13 PCGs + 12S/16S against a reference gene set from the user's own already-annotated data (a prior run's `build/01_per_gene/`, or a MitoAnnotator archive) — it does not bypass `mitotree-preflight`'s own validation. Requires at least one already-annotated sample on hand to serve as the reference. See `preflight/mitotree-taxon-sampling/SKILL.md`.

### 0.4 The run command

```bash
# Run the preflight
pixi run preflight

# Run the main build (after preflight passes)
pixi run run

# Run the QC (after the build finishes)
pixi run qc
```

## A. Pipeline Architecture

The analysis is divided into **three sequential phases**: Preflight → Build → QC, plus three optional post-build branches (Concordance, ANI, Codon Usage) that fan out from `build/mitotree-builder` and fan back in at `qc/mitotree-qc` — mirroring the original Nextflow DAG's parallel `done_ch.mix(...)` branches, but as independently invocable sub-skills instead of a fixed DAG.

### Phase 1: Preflight (The Audit)

**Goal**: Validate inputs, gene coverage, tool availability, and write `params.json` + `preflight.md`.

- **Sub-skill**: `preflight/mitotree-preflight`
- **Must-Verify**: mafft/trimal/AMAS.py/iqtree/gotree present; ≥1 sample source discovered; gene coverage checked.
- **Exit Gate**: `$RUN_DIR/preflight.md` exists with overall verdict ≥ `GO-WITH-WARNINGS`.

### Phase 2: Build (The Workflow)

**Goal**: Build the supermatrix, species tree, and (if an outgroup was given) re-root it.

- **Sub-skill**: `build/mitotree-builder`
- **Must-Verify**: Preflight verdict ≥ `GO-WITH-WARNINGS`.
- **Exit Gate**: `build/species_tree.nwk` exists.

### Phase 3: QC (The Report)

**Goal**: Build the final report, folding in any optional branches that were run.

- **Sub-skill**: `qc/mitotree-qc`
- **Must-Verify**: `build/species_tree.nwk` present.
- **Exit Gate**: `$RUN_DIR/qc-report.md` exists.

## B. The Evidence Chain & Handoff Contract

| From → To | Artifact | Owner | Consumer |
| --- | --- | --- | --- |
| Upstream → Preflight | Input archives/FASTAs | upstream | `preflight/mitotree-preflight` |
| Preflight → Build | `preflight.md` + `params.json` | `preflight/mitotree-preflight` | `build/mitotree-builder` |
| Build → QC | `build/species_tree.nwk` + `05_stats/` | `build/mitotree-builder` | `qc/mitotree-qc` |
| Build → Concordance/ANI/Codon-usage | `build/` outputs (gene alignments, species tree) | `build/mitotree-builder` | optional branches |
| Concordance/ANI → QC | `concordance/`, `ani/` dirs | optional branches | `qc/mitotree-qc` |
| QC → User | `qc-report.md` | `qc/mitotree-qc` | User |

## C. Output contract — no artifact

This orchestrator skill produces **no file of its own**. All artifacts are produced by the sub-skills.

## D. Go/No-Go Gates

The agent must stop and warn the user if:

- **Preflight `NO-GO`**: `preflight.md` shows `NO-GO` verdict.
- **User asks for divergence-time dating**: refuse and explain (see `preflight/mitotree-preflight` SP3) rather than scaffolding a dating stage.

## E. Common follow-ups

| User says | What to do |
| --- | --- |
| "Show me my last report" | `cat $RUN_DIR/qc-report.md` |
| "What did the run produce?" | Read `$RUN_DIR/qc-report.md` |
| "Can I date this tree?" / "divergence times" / "MCMCtree" | Explain the scope exclusion (mitogenome-only dating is unreliable) — see `preflight/mitotree-preflight` SP3. Do not build a dating stage. |
| "Check gene-tree concordance" / "gCF/sCF" / "RF distances" | Invoke `concordance/mitotree-concordance` — writes `concordance/rf_distances.tsv` + `concordance.cf.stat`, then re-run `qc/mitotree-qc` to fold it in. Offer the gCF-vs-sCF scatter figure (see next row) afterwards. |
| "Concordance scatter" / "gCF vs sCF plot" / "manuscript figure of concordance" | Invoke `concordance/mitotree-concordance` step 4 (`scripts/plot_concordance_scatter.py`) on an existing `concordance.cf.stat` — writes `concordance/concordance_scatter.pdf` + `.png` (UFBoot-category colours, 33% random-expectation lines; manual labels via `--labels-tsv`, or auto-label discordant nodes via `--label-nodes`). Requires `concordance/mitotree-concordance` to have run first. |
| "Check for duplicate/near-identical samples" / "ANI" | Invoke `ancillary/mitotree-ani` — writes `ani/fastani_results.tsv` + a heatmap, then re-run `qc/mitotree-qc`. |
| "Codon usage" / "RSCU" | Invoke `ancillary/mitotree-codon-usage` (nucleotide mode only — meaningless on translated `aa`-mode input). |
| "Add more taxa" / "pull mitogenomes from NCBI" / "expand the dataset" | Invoke `preflight/mitotree-taxon-sampling` *before* `preflight/mitotree-preflight` — writes `taxon_sampling/genes/*.genes.fasta` (ready as `--fasta LABEL:FILE` input), then hand off to the normal preflight gate. Needs an already-annotated reference sample on hand; ask for one (MitoAnnotator archive or a prior run's `build/01_per_gene/`) if not already clear from the request. |
| "Prune the tree" / "remove these taxa" / "drop misidentified mitogenomes" | Invoke `prune/mitotree-prune` — writes `pruned/species_tree_pruned.nwk` + a report. Ask which mode (keep list / exclude list / match a reference tree) if not already clear from the request. |
| "Visualize the tree" / "iTOL annotation" / "colorstrip" | Invoke `visualize/mitotree-itol` — writes `itol/labels.txt` (+ `colorstrip.txt`, `legend.md` if a grouping was given). Upload alongside `build/species_tree.nwk` at itol.embl.de. |
| "Render a publication figure" / "export tree as PDF/SVG/PNG" | Invoke `visualize/mitotree-render` — writes `figures/species_tree.{svg,pdf,png}`. No italics (toyplot limitation) — use `mitotree-itol` instead if italic species names are required. |

## Verification

- [ ] Preflight was performed; `preflight.md` overall verdict ≥ `GO-WITH-WARNINGS`.
- [ ] Build completed; `build/species_tree.nwk` exists.
- [ ] QC was performed; `qc-report.md` exists.

## Handoff pointers

After this orchestrator routes the user to a phase, the agent should explicitly state which sub-skill was invoked and what file it produces.

## Provenance

Rebuilt with the [bioinfo-skill-creator](https://github.com/cheahhl814/bioinfo-skill-creator) meta-skill from the `MitoTree` Nextflow pipeline (`PHYLOGENETIC/MitoTree/main.nf` + `bin/`), following the same pattern as `phylogenetics-agent`. Divergence-time dating (MCMCtree/PAML, OrthoFinder orthology inference was investigated and confirmed absent from the source pipeline — genes are extracted by fixed name, not inferred orthology) was deliberately dropped during the port.


## Dependencies → biodb-fetch (shared data lake)

Consumes [biodb-fetch](https://github.com/cheahhl814/biodb-fetch) v1.0.0 (deployed at `~/.agents/skills/biodb-fetch/`) for mitogenome phylogenomics — mitogenome/marker record retrieval from NCBI:

- `scripts/ncbi_fetch.py fetch --db nucleotide --uids <ids> --rettype fasta` (+ `search` for query discovery) — rate-limit owned, store dedup across runs

- **Shared store + contract v2**: payloads land in biodb-fetch's content-addressed `store/` (md5 `.sha256` sidecars); repeat requests return `cache=hit` in the `manifest.tsv` (accession/source/url/path/md5/retrieval_date/tool/format/contract_version/cache) — zero re-download across skills/batches; `--refresh` re-pulls; never parse biodb-fetch internals.
- **Degradation**: biodb-fetch unavailable (not deployed / offline / preflight NO-GO) → proceed exactly as before with user-supplied files; note the deviation in the run report.
