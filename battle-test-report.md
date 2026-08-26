# mitotree-agent — Battle-Test Report

**Date**: 2026-08-26
**Verdict**: `PASS`

## 1. Frontmatter coherence

All 9 `SKILL.md` files parse as valid YAML frontmatter with `name`/`description`/`version`/`updated`/`triggers`, all at `version: 1.0.0`:

```
SKILL.md OK: mitotree-agent 1.0.0
preflight/mitotree-preflight/SKILL.md OK: mitotree-preflight 1.0.0
build/mitotree-builder/SKILL.md OK: mitotree-builder 1.0.0
qc/mitotree-qc/SKILL.md OK: mitotree-qc 1.0.0
debug/mitotree-debug/SKILL.md OK: mitotree-debug 1.0.0
battle-test/mitotree-battle-test/SKILL.md OK: mitotree-battle-test 1.0.0
concordance/mitotree-concordance/SKILL.md OK: mitotree-concordance 1.0.0
ancillary/mitotree-ani/SKILL.md OK: mitotree-ani 1.0.0
ancillary/mitotree-codon-usage/SKILL.md OK: mitotree-codon-usage 1.0.0
```

## 2. Contract wiring

- `preflight/mitotree-preflight` → `params.json` + `preflight.md` → consumed by `build/mitotree-builder`. ✅
- `build/mitotree-builder` → `build/species_tree.nwk` + `build/05_stats/*.txt` → consumed by `qc/mitotree-qc`. ✅
- `qc/mitotree-qc` → `qc-report.md` → orchestrator's `qc-done` stage. ✅
- `concordance/mitotree-concordance` → `concordance/rf_distances.tsv` + `concordance.cf.stat` → consumed (optionally) by `qc/mitotree-qc`. ✅ (live-verified: `qc_report.py`'s `parse_rf_distances`/`parse_concordance` read these correctly)
- `ancillary/mitotree-ani` → `ani/fastani_results.tsv` → consumed (optionally) by `qc/mitotree-qc`. ✅ (live-verified: `parse_ani` reads it correctly)
- `ancillary/mitotree-codon-usage` → `codon_usage/*.tsv` → informational only, not read by `qc_report.py` (confirmed no reference).

## 3. Signature-library completeness

All 7 non-orchestrator sub-skills have well over the 3-row minimum (6–14 rows each, counted via the same `grep -c '^| `' pattern as `phylo-battle-test`).

## 4. Docs-corpus freshness

```
mafft: OK (cli-reference.md github-readme.md README.md)
trimal: OK (cli-reference.md github-readme.md README.md)
amas: OK (cli-reference.md README.md)
iqtree: OK (cli-reference.md github-readme.md README.md)
gotree: OK (cli-reference.md github-readme.md README.md)
fastani: OK (cli-reference.md README.md)
```

`mafft`/`trimal`/`iqtree`/`gotree` entries were reused directly from `phylogenetics-agent/docs-corpus/` (same tools, same versions via the shared pixi channel pins). `amas`/`fastani` are new — verbatim `--help` captures from this skill's own pixi env.

## 5. Pixi parse

```
i Dry-run:lock file would not change
```

## 6. Git hygiene

Not yet a git repo — `git init` + first commit pending user confirmation before any push.

## Build process — investigated live

- **OrthoFinder not present**: read `core_mito_phylo.py` in full; genes are extracted by fixed name (13 PCGs + 12S/16S) from MitoAnnotator archive headers, not inferred orthology. No orthology-inference dependency was needed in the port (this corrected an earlier assumption made mid-design-discussion, before this build).
- **Dating scope exclusion**: `bin/annotate_calibrations.R`, `bin/bifurcate_tree.py`, and the `mcmctree_dating` Nextflow process were read but deliberately not ported, per explicit user instruction — mitogenome-only divergence-time estimation is unreliable. `preflight/mitotree-preflight` SP3 refuses and explains this if requested.
- **argv-list bug found and fixed live**: `build_supermatrix.py`'s `step_align`/`step_trim` originally built `subprocess.run([mafft_bin, ...])` assuming `mafft_bin` was a single executable path; passing `"pixi run mafft"` as one CLI arg raised `FileNotFoundError`. Fixed by accepting the command as a string and `.split()`-ing it into argv form before use — same convention as `prune_tree.py`'s `--gotree-cmd`. Verified fixed via a live end-to-end run.
- **`source_dirs` argparse bug found and fixed live**: was `nargs="+"` (required), which made bare-`--fasta`-only invocations (no MitoAnnotator archives) impossible despite the script's own docstring describing that as a valid mode. Changed to `nargs="*", default=[]`.

## Delegated work (herdr + Pi agents)

Two sub-skills were delegated to Pi agents running in parallel herdr panes (`pi` integration, model `qwen3.5:397b-cloud`), each given a self-contained brief covering: the `phylo-prune` SKILL.md format template, the exact upstream/downstream artifact contract (including the already-written, already-tested `qc_report.py` parser functions they had to match byte-for-byte on column names/paths), the relevant original Nextflow process definitions, and a live-testing requirement against a synthetic dataset.

- **`concordance/mitotree-concordance`** (agent `concordance-agent`): built `scripts/compute_concordance.py` wrapping `gotree compare trees --rf` (chosen over porting the original `plot_rf_distances.R`, avoiding a new R-plotting dependency — same "prefer an already-installed tool" principle used by `phylo-prune`) plus an `iqtree --gcf --scf` bash recipe. Live-tested end-to-end against a 4-taxon/3-gene synthetic build, including running the real `qc/mitotree-qc/scripts/qc_report.py` against its output and confirming correct parsing. Found and worked around a real discrepancy: `gotree compare trees --rf`'s actual output is a bare RF value per line, not the tab-separated index+RF format its own `--help` text implies.
- **`ancillary/mitotree-ani`** + **`ancillary/mitotree-codon-usage`** (agent `ani-codon-agent`): built `scripts/prepare_ani_input.py` + `scripts/plot_ani_heatmap.py` (a Python/matplotlib+scipy replacement for the original `plot_ani_heatmap.R`, since that R script needs `tidyverse`+`pheatmap` which weren't in `pixi.toml` — same dependency-minimization reasoning as the concordance agent's gotree choice) and copied `calculate_rscu.py` in as-is. First live-test attempt hit transient tool-stream errors from an oversized inline heredoc; re-prompted with a smaller-synthetic-data instruction and it completed successfully on retry.

**Post-delegation verification I ran myself** (not just trusting the agents' self-reports): re-ran `calculate_rscu.py` and `plot_ani_heatmap.py` against this skill's real pixi env (the agent had tested `calculate_rscu.py` in an ad-hoc `venv` since `biopython` wasn't yet in `pixi.toml`) — both produced correct output once `matplotlib`, `scipy`, and `biopython` were added to `pixi.toml` and `pixi install` re-run. `fastANI` produced real (non-trivial) pairwise ANI values (98.3–100%) against synthetic near-identical sequences, confirming the wrapper's argument-passing is correct end to end.

## Completed

- [x] `preflight/mitotree-preflight`, `build/mitotree-builder`, `qc/mitotree-qc`, `debug/mitotree-debug`, orchestrator `SKILL.md`, `battle-test/mitotree-battle-test` — written and (preflight/build/qc) live-tested against a synthetic 4-taxon/3-gene dataset.
- [x] `concordance/mitotree-concordance`, `ancillary/mitotree-ani`, `ancillary/mitotree-codon-usage` — delegated, built, live-tested, and independently re-verified.
- [x] `docs-corpus/` populated for all 6 tools (4 reused, 2 new).
- [x] `README.md`, `pixi.toml`, `.gitignore` written.
- [x] Structural battle-test: `PASS`.

## Addendum (2026-08-26, later same day): mitotree-preflight given a real script

The user noticed `preflight/mitotree-preflight/scripts/` was empty — a real gap: unlike every other phase, `mitotree-preflight`'s Procedure section only had inline one-liners and a hand-wavy "if a bespoke dry-run flag doesn't exist yet..." note for gene-coverage checking, with no actual script ever written.

**Fixed**: `scripts/preflight_check.py` now does the tool audit, source discovery, and gene-coverage dry-run for real — by importing `discover_archives`/`parse_fasta_specs`/`step_extract` directly from `build_supermatrix.py` rather than reimplementing them, so preflight's coverage numbers are guaranteed identical to what the build phase will see.

**Investigated live, per the user's follow-up**: the user pointed out this dry-run is valuable specifically because most mitogenomes retrieved from NCBI lack gene annotation (a single whole-sequence FASTA record per sample, not per-gene records). Built and added **SP5**: when every gene shows `0/n` across every sample, this is a distinct `NO-GO` signature from ordinary partial coverage gaps (SP1) — the pipeline extracts genes by header-name matching and does not annotate raw genomes, so a fully-zero coverage table means the input needs MitoAnnotator (or similar) run first, not a partial-sample/partial-gene fix.

Live-tested against two synthetic 4-sample datasets:
- Per-gene FASTAs (12S/16S/Cytb headers) → `GO`, 100% coverage.
- Whole-sequence "raw mitogenome" FASTAs (headers like `taxA_complete_mitogenome`, no gene-name match) → `NO-GO` via SP5, with the exact "this looks like whole, unannotated mitogenome FASTA(s)... run MitoAnnotator first" explanation written into `preflight.md`.

## Completed (addendum)

- [x] `preflight/mitotree-preflight/scripts/preflight_check.py` written, live-tested (GO and NO-GO/SP5 paths), frontmatter/signature-library re-verified.

## Addendum (2026-08-26, later same day): preflight/mitotree-taxon-sampling added

Optional NCBI-retrieval stage, per user direction: fetch whole complete-mitogenome FASTA (not GenBank-annotation parsing, which was the original design sketch) and extract the 13 PCGs + 12S/16S via BLASTN against a reference gene set the user already has, inspired by `PHYLOGENETIC/scripts/blast_extract_cds.py` and `core_mito_phylo.py`.

**Investigated live**:
- `fetch_mitogenomes.py` (adapted from `phylogenetics-agent/preflight/phylo-taxon-sampling/scripts/taxon_sampling.py`) run against **real NCBI data** (`--expand-taxon "Betta"`, no mocking): 12 candidate complete-mitogenome records fetched, correctly deduplicated to 8 via `vsearch --cluster_fast` (RefSeq `NC_026581` preferred over byte-identical duplicate `AB571120`, both real *Betta splendens* accessions).
- `blast_extract_genes.py` (trimmed from `blast_extract_cds.py` to just its BLASTN-extraction core; reference loading reuses `build_supermatrix.py`'s `genes_from_archive`/`genes_from_fasta` directly, no GFF3/GenBank-feature parsing needed) — first synthetic-data test passed cleanly (strand correction + ATG-trim verified on a fabricated 993bp genome with a reverse-complemented PCG).
- **Real bug found and fixed via the real-data test**: running the full extraction against the 8 real fetched *Betta* genomes, using a real MitoAnnotator archive (`Betta_pi_BPI001.tar.gz`) as reference, initially showed most samples missing 6-10+ genes each. Root cause: `blastn`'s default `-task` is `megablast` (tuned for near-identical sequences) — a real Cytb hit at 78% identity was invisible to megablast but found immediately with `-task blastn`. Fixed by always passing `-task blastn` explicitly. Re-run after the fix: 7/8 samples extracted all 15/15 genes cleanly; the 8th missed only `ATPase8` (a short, fast-evolving gene) — a realistic biological result at real cross-species divergence, not a bug.

**Completed**:
- [x] `preflight/mitotree-taxon-sampling/scripts/fetch_mitogenomes.py` — live-tested against real NCBI data.
- [x] `preflight/mitotree-taxon-sampling/scripts/blast_extract_genes.py` — live-tested against synthetic data and real fetched genomes + a real MitoAnnotator reference; `-task blastn` bug found and fixed.
- [x] `preflight/mitotree-taxon-sampling/SKILL.md` (v1.0.0, 4 stop points, 5-row troubleshooting table including the `-task blastn` finding).
- [x] Orchestrator `SKILL.md`, `README.md`, `battle-test/mitotree-battle-test/SKILL.md` updated; `pixi.toml` gained `entrez-direct`/`vsearch`/`blast`; `docs-corpus/{entrez-direct,vsearch,blast}` reused from `phylogenetics-agent`.
- [x] Full structural re-check: all 10 `SKILL.md` files pass frontmatter, all 9 tools have docs-corpus entries, `pixi lock --dry-run` succeeds.
