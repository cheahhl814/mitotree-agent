# MitoTree Agent

> **v1.0.0.** Checkpointed, agent-driven mitogenome phylogenomics pipeline (extract, align, trim, concatenate, tree, re-root) with optional gene-tree concordance, ANI, and codon-usage branches. No divergence-time dating — mitogenome-only data does not reliably calibrate node ages.

This is the **mitotree-agent** skill, built with the [bioinfo-skill-creator](https://github.com/cheahhl814/bioinfo-skill-creator) meta-skill, ported from the `MitoTree` Nextflow pipeline (`main.nf` + `bin/`).

## Pipeline Overview

| # | Phase | Goal | Sub-Skill | Artifact produced |
|:---|:---|:---|:---|:---|
| *opt.* | **Taxon Sampling** (broaden dataset) | Fetch whole complete-mitogenome FASTA from NCBI, then BLASTN-extract PCGs+rRNA against a reference gene set the user already has. Not gated — user-triggered, before preflight. | `preflight/mitotree-taxon-sampling` | `$RUN_DIR/taxon_sampling/genes/*.genes.fasta` |
| **1** | **Preflight** (The Audit) | Validate inputs, gene coverage, tool availability. | `preflight/mitotree-preflight` | `$RUN_DIR/preflight.md` + `$RUN_DIR/params.json` |
| **2** | **Build** (The Workflow) | Extract PCGs+rRNA → align (mafft) → trim (trimAl) → concatenate (AMAS) → species tree (IQ-TREE) → re-root (ape). | `build/mitotree-builder` | `$RUN_DIR/build/species_tree.nwk` |
| **3** | **QC** (The Report) | Build the final report: supermatrix stats, UFBoot support, and any optional-branch results. | `qc/mitotree-qc` | `$RUN_DIR/qc-report.md` |
| *opt.* | **Concordance** (gene-tree agreement) | Per-gene trees, RF distances vs. species tree, gCF/sCF concordance factors. Not gated — user-triggered after build. | `concordance/mitotree-concordance` | `$RUN_DIR/concordance/*` |
| *opt.* | **ANI** (sample identity check) | Pairwise Average Nucleotide Identity via fastANI — flags near-identical/duplicate samples. Independent of the phylogenomics branch, can run any time after preflight. | `ancillary/mitotree-ani` | `$RUN_DIR/ani/*` |
| *opt.* | **Codon Usage** (RSCU) | Relative synonymous codon usage per gene (nucleotide mode only). | `ancillary/mitotree-codon-usage` | `$RUN_DIR/codon_usage/*` |
| *opt.* | **Prune** (curate taxa) | Drop misidentified/excluded taxa or match a reference pruned tree, via gotree. Not gated — user-triggered after build. | `prune/mitotree-prune` | `$RUN_DIR/pruned/species_tree_pruned.nwk` |
| *opt.* | **Visualize — iTOL** (interactive annotation) | iTOL tip-label + colorstrip datasets for `build/species_tree.nwk` (already rooted — no separate reroot step needed). | `visualize/mitotree-itol` | `$RUN_DIR/itol/*` |
| *opt.* | **Visualize — Render** (static figure) | Offline SVG/PDF/PNG tree figure via toytree, with support values. | `visualize/mitotree-render` | `$RUN_DIR/render/*` |

## When to Use This Skill

✅ **Use this skill** when you need to:

- Build a multi-gene mitochondrial-genome species tree from MitoAnnotator archives or per-species gene FASTAs (13 PCGs + 12S/16S rRNA).
- Assess gene-tree/species-tree concordance (RF distances, gCF/sCF) or check for near-identical samples via ANI.
- Swap between nucleotide (`nt`, intragenus/species-level) and amino-acid (`aa`, family-level+) supermatrix modes.

❌ **Do NOT use this skill** for:

- Single-marker trees (one gene only) — use `phylogenetics-agent` instead.
- Divergence-time estimation (node ages, MCMCtree/PAML). **Deliberately out of scope**: mitogenome-only data does not reliably calibrate node ages (rate heterogeneity across the 13 PCGs + rRNAs, saturation at deeper splits, no independent nuclear-rate cross-check). Use a nuclear multi-locus or genome-scale dataset instead.

## Prerequisites

- A pixi environment with the required tools: `mafft, trimal, amas, iqtree, gotree, fastani, r-ape, r-getopt, pandas, numpy, matplotlib, scipy, biopython` (+ `entrez-direct, vsearch, blast` for taxon-sampling; `toytree, packaging` for `visualize/mitotree-render`), all optional stages.
- MitoAnnotator `.tar.gz` archives and/or per-species gene FASTAs (`--fasta LABEL:FILE`), at least 1 sample source (≥4 samples recommended for a meaningful tree).

## Installation

```bash
pixi init mitotree-agent-env
cd mitotree-agent-env
pixi workspace channel add conda-forge bioconda
pixi add --manifest-path /path/to/mitotree-agent/pixi.toml
```

## Usage

```bash
# Optional: broaden the dataset with more mitogenomes from NCBI (before Phase 1)
python3 preflight/mitotree-taxon-sampling/scripts/fetch_mitogenomes.py \
    --expand-taxon "Betta" --outgroup-taxon "Macropodus opercularis" \
    -o "$RUN_DIR/taxon_sampling/raw" \
    --summary-out "$RUN_DIR/taxon_sampling/summary.csv" \
    --metadata-out "$RUN_DIR/taxon_sampling/metadata.csv" \
    --esearch-cmd "pixi run esearch" --efetch-cmd "pixi run efetch" --vsearch-cmd "pixi run vsearch"
python3 preflight/mitotree-taxon-sampling/scripts/blast_extract_genes.py \
    --targets-dir "$RUN_DIR/taxon_sampling/raw" \
    --reference-archive /path/to/a/MitoAnnotator/archive.tar.gz \
    -o "$RUN_DIR/taxon_sampling/genes" \
    --blastn-cmd "pixi run blastn" --makeblastdb-cmd "pixi run makeblastdb"

# Phase 1: Preflight
RUN_DIR=/path/to/run
pixi run preflight

# Phase 2: Build (after Phase 1 passes)
python3 build/mitotree-builder/scripts/build_supermatrix.py \
    annotation_dir/ --fasta extra_sample:extra_sample.fasta \
    --mode nt -o "$RUN_DIR/build" \
    --mafft "pixi run mafft" --trimal "pixi run trimal" --amas "$(pixi run which AMAS.py)"

# Phase 3: QC (after the build finishes)
python3 qc/mitotree-qc/scripts/qc_report.py "$RUN_DIR" --output-dir "$RUN_DIR"

# Optional: gene-tree concordance (after the build finishes, user-triggered)
# See concordance/mitotree-concordance/SKILL.md

# Optional: ANI near-identical-sample check (after preflight, independent of the build)
# See ancillary/mitotree-ani/SKILL.md

# Optional: codon usage / RSCU (nucleotide mode only)
# See ancillary/mitotree-codon-usage/SKILL.md

# Optional: prune the final tree to a curated taxon subset (after the build finishes)
python3 prune/mitotree-prune/scripts/prune_tree.py "$RUN_DIR/build/species_tree.nwk" \
    --exclude-file exclude_list.txt \
    -o "$RUN_DIR/pruned/species_tree_pruned.nwk" \
    --gotree-cmd "pixi run gotree"

# Optional: iTOL interactive annotation datasets (already-rooted tree, no reroot step)
python3 visualize/mitotree-itol/scripts/itol_tip_labels.py "$RUN_DIR/build/species_tree.nwk" \
    -o "$RUN_DIR/itol/labels.txt" --gotree-cmd "pixi run gotree"

# Optional: offline static figure (SVG/PDF/PNG) via toytree
python3 visualize/mitotree-render/scripts/render_tree.py "$RUN_DIR/build/species_tree.nwk" \
    -o "$RUN_DIR/render/species_tree"
```

Each phase has explicit **Inputs/Outputs contract** at the top of its sub-skill SKILL.md. If the upstream artifact is missing, the sub-skill **must refuse to proceed** and tell the user which upstream step to run.

## Repository Layout

```
mitotree-agent/
├── SKILL.md                              # Master orchestrator (router only)
├── README.md                             # This file
├── pixi.toml                             # Pixi env: mafft, trimal, amas, iqtree, gotree, fastani, r-ape, r-getopt, pandas, numpy
├── preflight/
│   ├── mitotree-taxon-sampling/           # Optional — broaden dataset via NCBI + BLASTN gene extraction
│   │   ├── SKILL.md
│   │   └── scripts/
│   │       ├── fetch_mitogenomes.py
│   │       └── blast_extract_genes.py
│   └── mitotree-preflight/SKILL.md       # Phase 1 — input + gene-coverage validation
├── build/
│   └── mitotree-builder/
│       ├── SKILL.md                      # Phase 2 — supermatrix + species tree + re-root
│       └── scripts/
│           ├── build_supermatrix.py
│           └── reroot_tree.R
├── qc/
│   └── mitotree-qc/
│       ├── SKILL.md                      # Phase 3 — final report
│       └── scripts/
│           └── qc_report.py
├── concordance/
│   └── mitotree-concordance/             # Optional — gene-tree RF distances + gCF/sCF
│       ├── SKILL.md
│       └── scripts/
├── ancillary/
│   ├── mitotree-ani/                     # Optional — fastANI near-identical-sample check
│   │   ├── SKILL.md
│   │   └── scripts/
│   └── mitotree-codon-usage/             # Optional — RSCU codon usage (nt mode only)
│       ├── SKILL.md
│       └── scripts/
├── prune/
│   └── mitotree-prune/                   # Optional — curate taxa (gotree prune)
│       ├── SKILL.md
│       └── scripts/
├── visualize/
│   ├── mitotree-itol/                    # Optional — iTOL interactive annotation datasets
│   │   ├── SKILL.md
│   │   └── scripts/
│   └── mitotree-render/                  # Optional — offline SVG/PDF/PNG figure (toytree)
│       ├── SKILL.md
│       └── scripts/
├── debug/
│   └── mitotree-debug/SKILL.md           # Failure interpreter (signature library)
├── battle-test/
│   └── mitotree-battle-test/SKILL.md     # Structural smoke tests
├── docs-corpus/                          # Offline Markdown snapshot of upstream tool docs
│   ├── mafft/ trimal/ iqtree/ gotree/    # reused from phylogenetics-agent (same tools)
│   ├── entrez-direct/ vsearch/ blast/    # reused from phylogenetics-agent (same tools)
│   ├── toytree/                          # reused from phylogenetics-agent (same tool)
│   └── amas/ fastani/                    # new for this skill
└── test_smoke.py                         # Skill-level structural smoke test (NOT shipped; see battle-test/ for QA)
```

_Note: `test_smoke.py` is referenced above for historical reasons — the v1.0.0 mitotree-agent bundle does not actually ship a `test_smoke.py`. Use the battle-test sub-skill (`battle-test/mitotree-battle-test/SKILL.md`) for structural validation._

## Provenance

Ported with the [bioinfo-skill-creator](https://github.com/cheahhl814/bioinfo-skill-creator) meta-skill from the `MitoTree` Nextflow pipeline, following the same pattern as [phylogenetics-agent](https://github.com/cheahhl814/phylogenetics-agent). OrthoFinder orthology inference was investigated and confirmed absent from the source pipeline (genes are extracted by fixed name from MitoAnnotator archives, not inferred orthology) — no orthology-inference dependency was needed in the port. Divergence-time dating (MCMCtree/PAML) was deliberately dropped.

## License

Same as the rest of the AiX-BIO skill tree (private GitHub repo, MIT-style — see repository owner).
