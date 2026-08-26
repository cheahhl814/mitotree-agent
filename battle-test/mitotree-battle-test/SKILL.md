---
name: mitotree-battle-test
description: 'Structural quality gate for mitotree-agent: frontmatter coherence across all SKILL.md files, sub-skill contract wiring, signature-library completeness, docs-corpus population for all tools, pixi.toml resolution, and (opt-in) an execution smoke test on a tiny synthetic dataset. Writes battle-test-report.md. Triggers: "battle-test mitotree-agent", "mitotree-battle-test", "is the mitotree skill ready".'
version: 1.0.0
updated: "2026-08-26"
triggers:
  - "battle-test mitotree-agent"
  - "mitotree-battle-test"
  - "is the mitotree skill ready"
---

# mitotree-battle-test

> **v1.0.0.** Quality gate for `mitotree-agent`. Mirrors `phylo-battle-test`'s structural-check pattern. Never modifies the skill — read-only checks plus one opt-in execution smoke test.

## 0. Inputs / Outputs contract

### Inputs

| Path | Required? |
| --- | --- |
| `$RUN_DIR/SKILL.md`, all sub-skill `SKILL.md` | yes |
| `$RUN_DIR/pixi.toml`, `$RUN_DIR/docs-corpus/` | yes |

### Outputs

| Path | Format |
| --- | --- |
| `$RUN_DIR/battle-test-report.md` | Markdown, verdict `PASS`/`PASS-WITH-WARNINGS`/`FAIL` |

## Procedure

### 1. Frontmatter coherence

```bash
for f in SKILL.md preflight/mitotree-preflight/SKILL.md \
         preflight/mitotree-taxon-sampling/SKILL.md \
         build/mitotree-builder/SKILL.md \
         qc/mitotree-qc/SKILL.md debug/mitotree-debug/SKILL.md \
         battle-test/mitotree-battle-test/SKILL.md \
         concordance/mitotree-concordance/SKILL.md \
         ancillary/mitotree-ani/SKILL.md ancillary/mitotree-codon-usage/SKILL.md; do
    python3 -c "
import yaml
fm = open('$f').read().split('---')[1]
d = yaml.safe_load(fm)
assert all(k in d for k in ('name','description','version','updated','triggers')), 'missing key in $f'
print('$f OK:', d['name'], d['version'])
"
done
```

### 2. Contract wiring

- `preflight/mitotree-preflight` outputs `params.json` + `preflight.md` → `build/mitotree-builder` inputs list the same two paths. ✅
- `build/mitotree-builder` outputs `build/species_tree.nwk` → `qc/mitotree-qc` inputs list it. ✅
- `qc/mitotree-qc` outputs `qc-report.md` → orchestrator's stage-detection ladder checks for it (`qc-done` stage). ✅
- `concordance/mitotree-concordance` and `ancillary/mitotree-ani` outputs (`concordance/`, `ani/`) → `qc/mitotree-qc` inputs list them as optional. ✅
- `preflight/mitotree-taxon-sampling` outputs `taxon_sampling/genes/*.genes.fasta` → consumable as `preflight/mitotree-preflight`'s `--fasta LABEL:FILE` input. ✅

### 3. Signature-library completeness

```bash
for f in preflight/mitotree-preflight/SKILL.md preflight/mitotree-taxon-sampling/SKILL.md \
         build/mitotree-builder/SKILL.md \
         qc/mitotree-qc/SKILL.md debug/mitotree-debug/SKILL.md \
         concordance/mitotree-concordance/SKILL.md \
         ancillary/mitotree-ani/SKILL.md ancillary/mitotree-codon-usage/SKILL.md; do
    N=$(grep -c '^| `' "$f")
    [ "$N" -ge 3 ] && echo "$f: $N signature rows OK" || echo "$f: only $N rows — FAIL"
done
```

### 4. Docs-corpus freshness

```bash
for tool in mafft trimal amas iqtree gotree fastani entrez-direct vsearch blast; do
    d="docs-corpus/$tool"
    [ -f "$d/README.md" ] || { echo "$tool: MISSING README.md — FAIL"; continue; }
    if [ -f "$d/cli-reference.md" ] || [ -f "$d/github-readme.md" ] || [ -f "$d/web-docs.md" ]; then
        echo "$tool: OK ($(ls "$d" | tr '\n' ' '))"
    else
        echo "$tool: index exists but zero content sources — FAIL"
    fi
done
```

### 5. Pixi parse

```bash
pixi lock --dry-run --manifest-path pixi.toml
```

### 6. Git hygiene

```bash
git -C . status --porcelain | grep -q . && echo "uncommitted changes present" || echo "clean"
```

### 7. (Opt-in, `--execution`) Tiny synthetic-dataset smoke test

Only run when explicitly requested. Uses the same 4-taxon/3-gene synthetic FASTA style already live-verified during this skill's build (see `battle-test-report.md`'s build addendum for the exact dataset and results).

```bash
# Run preflight + build/mitotree-builder against 4 synthetic per-species gene FASTAs
# (--fasta LABEL:FILE for each, --genes 12S 16S Cytb). Verify all 3 checkpoints +
# a parseable build/species_tree.nwk. IQ-TREE UFBoot on 4 taxa is expected to be
# low/uninformative (verified live) — that is not itself a smoke-test failure.
```

## Interpretation Guidelines

- **PASS**: all structural checks pass, `pixi lock --dry-run` succeeds, every tool has ≥ 1 docs-corpus source, at most 1 `WARN`.
- **PASS-WITH-WARNINGS**: 1 non-blocking issue (e.g. uncommitted `git status`, one tool with only a `web-docs.md` source).
- **FAIL**: any missing docs-corpus entry, broken contract wiring, `pixi lock` failure, or a sub-skill with < 3 signature-library rows.

## Verification

- [ ] All 10 `SKILL.md` files have valid frontmatter.
- [ ] Contract wiring verified for all phase transitions, including the optional concordance/ANI/taxon-sampling branches.
- [ ] All non-orchestrator sub-skills have ≥ 3 signature-library rows.
- [ ] All 9 tools have a populated `docs-corpus/<tool>/`.
- [ ] `pixi lock --dry-run` succeeds.
- [ ] `git status` reported (informational, not blocking unless the user asked for a clean commit).

## Handoff

> Battle-test complete. Verdict: `<verdict>`. If `PASS`/`PASS-WITH-WARNINGS`, recommend `git init` + first commit + `git push` to `github.com/cheahhl814/mitotree-agent` and deploy via `rsync -av --delete` to `~/.pi/agent/skills/mitotree-agent/`.
