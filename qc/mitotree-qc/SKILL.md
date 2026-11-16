---
name: mitotree-qc
description: 'Build the final mitotree-agent run report: supermatrix stats, species-tree UFBoot support, and (if present) gene-tree RF distances, gCF/sCF concordance factors, and ANI — from build/mitotree-builder + optional concordance/ani branches. Writes qc-report.md + run_summary.json/.txt. Triggers: internal — invoked by mitotree-agent orchestrator only.'
version: 1.2.0
updated: "2026-11-16"
triggers:
  - "mitotree qc"
disable-model-invocation: true
---

# mitotree-qc

> Phase 3 of `mitotree-agent`. Wraps `qc_report.py` (live-verified: gracefully reports "not found" per section when the optional `concordance/`/`ani/` branches weren't run, rather than failing the whole report).

## 0. Inputs / Outputs contract

### Inputs

| Path                                                            | Source                             | Required?                                  |
| --------------------------------------------------------------- | ---------------------------------- | ------------------------------------------ |
| `$RUN_DIR/build/species_tree.nwk` + `build/05_stats/*.txt`      | `build/mitotree-builder`           | yes                                        |
| `$RUN_DIR/concordance/rf_distances.tsv` + `concordance.cf.stat` | `concordance/mitotree-concordance` | no — section reports "not found" if absent |
| `$RUN_DIR/ani/fastani_results.tsv`                              | `ancillary/mitotree-ani`           | no — section reports "not found" if absent |

### Outputs

| Path                        | Format                                                                |
| --------------------------- | --------------------------------------------------------------------- |
| `$RUN_DIR/run_summary.json` | JSON, machine-readable                                                |
| `$RUN_DIR/run_summary.txt`  | Text, human-readable, includes a `gotree draw text` ASCII tree        |
| `$RUN_DIR/qc-report.md`     | Markdown wrapper the orchestrator's stage-detection ladder checks for |

## 0.5 Ask-User Stop Points

### SP1 — Optional branches not run

| Trigger                                         | Evidence                                 | Action                                                                                                                                                                                                            |
| ----------------------------------------------- | ---------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `concordance/` and/or `ani/` directories absent | `qc_report.py`'s per-section "not found" | **Recommend, don't block**: mention in the handoff which optional branches weren't run and how to add them, but still produce `qc-report.md` from what's available — QC should never be gated on optional stages. |

**Auto-pick when**: always — this is informational only, never a stop.

### SP2 — QC flags raised (exit code 1)

| Trigger                                                                                                                     | Evidence                             | Action                                                                                                                                                                                   |
| --------------------------------------------------------------------------------------------------------------------------- | ------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `qc_report.py` exits 1 (low UFBoot support, RF outlier genes, misleading gCF/sCF nodes, or near-identical (>99%) ANI pairs) | `run_summary.txt`'s `[WARN]` markers | Surface each `[WARN]` section to the user with the concrete node/gene/pair named — don't just say "QC flagged issues," name them (mirrors `phylo-qc`'s interpretation-layer convention). |

**Auto-pick when**: exit code 0 → report clean, no flags to surface.

## Procedure

```bash
run_dir="$RUN_DIR"
python3 qc/mitotree-qc/scripts/qc_report.py "$run_dir" --output-dir "$run_dir"
exit_code=$?
```

Then wrap `run_summary.txt` into `qc-report.md` (Markdown headers + the same content, or a thin pointer file — mirror `phylo-qc`'s existing convention for the parent skill's report format).

## Interpretation Guidelines

- **UFBoot `nodes_below_70`**: weak clades — call these out by name, don't just report a count.
- **`flag_rf_outliers`**: a gene tree whose topology disagrees strongly with the species tree — could be incomplete lineage sorting, a paralog mixed in, or a misassembled/misannotated gene for that sample.
- **`flag_misleading_concordance`**: UFBoot ≥95% but gCF<50% or sCF<33% — high bootstrap support does not always mean the underlying gene trees or alignment sites agree; this is the single most important reliability check IQ-TREE offers beyond bootstrap alone.
- **`flag_near_identical`** (ANI >99%): possible duplicate/mislabeled sample, or genuinely very recently diverged populations — don't auto-resolve, surface to the user.

## Troubleshooting — Signature library

| Signature                                                 | Likely cause                                                                                 | Fix                                                                                                                                     |
| --------------------------------------------------------- | -------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| `not found: .../build/species_tree.nwk`                   | `build/mitotree-builder` hasn't completed, or ran with a different `$RUN_DIR`                | Check `checkpoints/reroot.done.json` exists first.                                                                                      |
| `no support values found in species_tree.nwk`             | Tree was built without bootstrap (`-B` omitted per SP2 in `mitotree-builder`)                | Expected on very small taxon sets where bootstrap was intentionally skipped — not an error, `species_tree` section will just be sparse. |
| `(gotree not found on PATH...)` in the ASCII-tree section | `qc_report.py` invoked outside the pixi env                                                  | Run via `pixi run python3 qc/mitotree-qc/scripts/qc_report.py ...`, or ensure `gotree` resolves on PATH.                                |
| `KeyError` on `No_of_taxa` / `Alignment_length` columns   | AMAS output format changed between versions, or `supermatrix_stats.txt` is truncated/corrupt | Re-run `build/mitotree-builder` Step 1; check `amas --version`.                                                                         |

## Verification

- [ ] `qc-report.md` exists and names every `[WARN]` section explicitly (not just a pass/fail count).
- [ ] `run_summary.json` parses as valid JSON.
- [ ] If `concordance/`/`ani/` were run, their sections are populated (not "not found") in the report.

## What NOT to do

- Do **not** treat a missing optional branch (concordance/ANI) as a QC failure — only the flags computed from data that *is* present should affect the verdict.
- Do **not** collapse `[WARN]` sections into a bare count in the user-facing summary — name the specific node/gene/pair, matching `phylo-qc`'s existing interpretation-layer precedent (see `phylogenetics-agent/qc/phylo-qc/SKILL.md`).
- **Support-scale (v1.0.1)**: `species_tree.nwk` support values are auto-normalized — fraction-encoded trees (0–1, e.g. RAxML-NG `.supportTBE`-style Newicks) are rescaled ×100 before the 95/70 threshold checks; percent-encoded trees (IQ-TREE `.treefile`/`.contree`, this pipeline's normal output) pass through unchanged. If support means ever print as "0.5" for a clearly 50%-supported tree, the input is fraction-encoded and this normalization has been reverted or bypassed.

## Handoff

> QC complete. `qc-report.md` written. `<N>` `[WARN]` flags raised across supermatrix/species-tree/concordance/ANI sections (see report for specifics). Optional branches not run: `<list, if any>`.
