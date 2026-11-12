---
name: mitotree-debug
description: 'Failure interpreter for mitotree-agent: matches stderr/exit-code signatures from preflight, build, concordance, ani, codon-usage, or qc against a signature library and recommends a fix. Triggers: "mitotree failed", "why did the build error", "debug mitotree".'
version: 1.0.0
updated: "2026-08-26"
triggers:
  - "mitotree failed"
  - "why did the mitotree build error"
  - "debug mitotree"
disable-model-invocation: true
---

# mitotree-debug

> **v1.0.0.** Cross-phase failure interpreter for `mitotree-agent`. Each phase sub-skill has its own local signature library (see their `## Troubleshooting` sections); this skill is the aggregator for when the failing phase isn't obvious yet.

## 0. Inputs / Outputs contract

### Inputs

| Path | Required? |
| --- | --- |
| The failing command's stderr, or `$RUN_DIR`'s partial state (which checkpoints exist) | yes |

### Outputs

| Path | Format |
| --- | --- |
| `$RUN_DIR/debug-report.md` | Markdown: diagnosis + concrete fix, verdict `RECOVERABLE`/`BLOCKED` |

## Procedure

1. Identify which phase produced the failure from `$RUN_DIR`'s checkpoint state (no `params.json` → preflight; `params.json` but no `checkpoints/supermatrix.done.json` → build stage 1; etc.).
2. Match the stderr against that phase's own `## Troubleshooting` table first (`preflight/mitotree-preflight/SKILL.md`, `build/mitotree-builder/SKILL.md`, `concordance/mitotree-concordance/SKILL.md`, `ancillary/mitotree-ani/SKILL.md`, `ancillary/mitotree-codon-usage/SKILL.md`, `qc/mitotree-qc/SKILL.md`).
3. If no match there, check the cross-phase signatures below.
4. Write `debug-report.md` with the matched signature, root cause, and the exact re-run command.

## Cross-phase Signature library

| Signature | Likely cause | Fix |
| --- | --- | --- |
| `FileNotFoundError: [Errno 2] No such file or directory: 'pixi run <tool>'` | A `--<tool>-cmd`/`--<tool>` argument was passed as one space-containing string to a wrapper that treats it as a single executable path, instead of being split into argv form | Resolve the real binary path first (`pixi run which <tool>`) and pass that, or confirm the wrapper script splits multi-word commands via `.split()` before building its subprocess argv (this project's convention — see `build_supermatrix.py`, `prune_tree.py`). |
| `command not found: pixi` | Sub-skill invoked outside an activated pixi env and without `pixi run` prefixing | Run every tool invocation through `pixi run <tool>` or inside `pixi shell`, from the skill root so the manifest resolves. |
| A downstream phase's input file is silently stale (edited `params.json` after a checkpoint was written) | No hash/mtime check ties a checkpoint back to the `params.json` state that produced it | See `build/mitotree-builder` SP1 — ask the user whether to resume or rebuild rather than silently trusting a stale checkpoint. |
| Any R script (`reroot_tree.R`) exits 1 with an `ape`-specific error about tip labels | Sample label mismatch between `--fasta LABEL:FILE` and what's recorded in `params.json.sample_labels` | Labels must match exactly (case-sensitive) end to end from `preflight` → `params.json` → `build`'s `--fasta` invocation. |

## Verification

- [ ] `debug-report.md` names the specific phase and signature matched, not a generic "something failed."
- [ ] The recommended fix is a concrete command, not a vague suggestion.

## What NOT to do

- Do **not** recommend a divergence-time-dating workaround for any tree-shape/support-value complaint — that's a deliberately out-of-scope phase (see `preflight/mitotree-preflight` SP3), not a bug to route around.

## Handoff

> Debug complete. Verdict: `<RECOVERABLE|BLOCKED>`. `debug-report.md` written with the matched signature and fix.
