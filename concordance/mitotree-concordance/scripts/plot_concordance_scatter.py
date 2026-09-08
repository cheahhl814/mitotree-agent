#!/usr/bin/env python3
"""Concordance-factor scatter plot (gCF vs sCF) from concordance.cf.stat.

Python/matplotlib port of the manuscript Fig4_concordance_scatter.R
(ggplot2) — the pixi env carries matplotlib but not ggplot2, matching the
same convention as ancillary/mitotree-ani's heatmap (Python replacement for
the original R/tidyverse script).

Figure semantics (kept faithful to the R original):
  - x = gCF (%), y = sCF (%), one point per internal branch of concordance.cf.stat
  - dashed reference lines at 33.33% ("Random" expectation: 1/3 of gene trees
    or sites) on both axes
  - points coloured by UFBoot support category from the Label column:
    < 70% (red), 70-95% (orange), >= 95% (green)
  - optional node labels: manual via --labels-tsv (ID<TAB>label), or automatic
    via --label-nodes N (labels the N most discordant nodes — lowest
    gCF+sCF) — placed with a greedy anti-overlap nudge

Inputs:  $RUN_DIR/concordance/concordance.cf.stat  (IQ-TREE --gcf/--scf output,
         '#' comment header, columns ID/Label/gCF/sCF)
Outputs: concordance/concordance_scatter.pdf + .png (+ optional labels TSV
         template when --labels-template is given)

Requires: matplotlib on PATH of the pixi env (pass --python-cmd to override).
"""
import argparse
import csv
import sys

import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt

# Same palette as the R original
PALETTE = {"< 70%": "#d73027", "70-95%": "#fc8d59", ">=95%": "#1a9641"}
RANDOM_CF = 100.0 / 3.0


def read_cf_stat(path: str) -> list[dict]:
    """Parse concordance.cf.stat: '#' comment lines, tab-separated body."""
    with open(path) as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    if not lines:
        raise SystemExit(f"ERROR: no data rows in {path}")
    return list(csv.DictReader(lines, delimiter="\t"))


def ufboot_category(label: float) -> str:
    if label < 70:
        return "< 70%"
    if label < 95:
        return "70-95%"
    return ">=95%"


def load_labels(path: str | None) -> dict[int, str]:
    if not path:
        return {}
    labels = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                print(f"WARNING: skipping malformed labels line: {line!r}", file=sys.stderr)
                continue
            labels[int(parts[0])] = parts[1]
    return labels


def auto_label_ids(rows: list[dict], n: int) -> set[int]:
    """The N most discordant nodes: lowest gCF, then lowest sCF."""
    ranked = sorted(rows, key=lambda r: (float(r["gCF"]), float(r["sCF"])))
    return {int(r["ID"]) for r in ranked[:n]}


def place_labels(points: list[dict]) -> None:
    """Greedy anti-overlap: alternate label placement left/right of each point
    in fixed-size steps, avoiding collisions with already-placed labels.
    Mutates each point dict with label_x/label_y/hjust/segment."""
    placed: list[tuple[float, float]] = []
    for i, p in enumerate(points):
        side = -1 if i % 2 == 0 else 1  # alternate left/right starting side
        lx, ly = p["_gCF"], p["_sCF"]
        nx, ny = 0.0, 0.0
        for step in range(1, 11):
            done = False
            for dy_sign in (1, -1):  # try above then below at each radius
                nx = side * 4 * step
                ny = 2.5 * step * dy_sign
                lx = min(max(p["_gCF"] + nx, 2.0), 104.0)
                ly = min(max(p["_sCF"] + ny, 2.0), 103.0)  # clamp: labels stay inside the axes
                if all(abs(lx - px) >= 16 or abs(ly - py) >= 5 for px, py in placed):
                    done = True
                    break
            if done:
                break
        p.update(label_x=lx, label_y=ly, hjust="left" if nx > 0 else "right",
                 segment=(abs(lx - p["_gCF"]) > 0.1 or abs(ly - p["_sCF"]) > 0.1))
        placed.append((lx, ly))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cf_stat", help="path to concordance.cf.stat")
    ap.add_argument("-o", "--output-prefix", required=True,
                    help="output prefix; writes <prefix>.pdf and <prefix>.png")
    ap.add_argument("--labels-tsv", help="optional TSV: node_id<TAB>clade label (manual labels)")
    ap.add_argument("--label-nodes", type=int, default=0, metavar="N",
                    help="auto-label the N most discordant nodes (lowest gCF, then sCF)")
    ap.add_argument("--labels-template-out", metavar="TSV",
                    help="write a ready-to-edit ID<TAB>label template for all nodes, then exit")
    ap.add_argument("--title", default=None, help="plot title")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    rows = read_cf_stat(args.cf_stat)
    for r in rows:
        r["_id"] = int(r["ID"])
        r["_gCF"] = float(r["gCF"])
        r["_sCF"] = float(r["sCF"])
        r["_label"] = float(r["Label"]) if r.get("Label") not in (None, "", "NA") else 0.0
        r["_cat"] = ufboot_category(r["_label"])

    if args.labels_template_out:
        with open(args.labels_template_out, "w") as f:
            f.write("# node_id\tclade_label\n")
            for r in sorted(rows, key=lambda r: r["_id"]):
                f.write(f"{r['_id']}\t\n")
        print(f"Labels template (one row per node, fill in the label column) -> {args.labels_template_out}")
        return

    labels = load_labels(args.labels_tsv)
    if args.label_nodes > 0:
        auto = auto_label_ids(rows, args.label_nodes)
        missing = auto - set(labels)
        for r in rows:
            if r["_id"] in missing:
                labels[r["_id"]] = f"node {r['_id']}"
    to_label = [r for r in rows if r["_id"] in labels]
    place_labels(to_label)

    import matplotlib.pyplot as plt  # noqa: F811 (keep pyplot import local to main)

    fig, ax = plt.subplots(figsize=(7, 6))

    # Random-expectation reference lines (1/3 of gene trees or sites)
    ax.axvline(RANDOM_CF, ls="--", lw=0.8, color="grey", alpha=0.6)
    ax.axhline(RANDOM_CF, ls="--", lw=0.8, color="grey", alpha=0.6)
    ax.annotate("Random (33%)", xy=(RANDOM_CF, 2), xytext=(2, 0), textcoords="offset points",
                fontsize=8, style="italic", color="grey", ha="left", va="bottom")

    # Points by UFBoot category (draw order: worst first so green sits on top)
    for cat in ("< 70%", "70-95%", ">=95%"):
        pts = [r for r in rows if r["_cat"] == cat]
        if pts:
            ax.scatter([p["_gCF"] for p in pts], [p["_sCF"] for p in pts],
                       s=70, alpha=0.9, color=PALETTE[cat], label=cat,
                       edgecolors="white", linewidths=0.5, zorder=3)

    # Label segments + text for labelled nodes
    for p in to_label:
        if p["segment"]:
            ax.plot([p["_gCF"], p["label_x"]], [p["_sCF"], p["label_y"]],
                    lw=0.6, color="grey", zorder=2)
        ax.text(p["label_x"], p["label_y"], labels[p["_id"]],
                fontsize=7.5, color=PALETTE[p["_cat"]], ha=p["hjust"], va="center",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="grey", lw=0.3), zorder=4)

    ax.set_xlim(0, 105)
    ax.set_ylim(0, 105)
    ax.set_xticks(range(0, 101, 20))
    ax.set_yticks(range(0, 101, 20))
    ax.set_xlabel("Gene concordance factor, gCF (%)")
    ax.set_ylabel("Site concordance factor, sCF (%)")
    if args.title:
        ax.set_title(args.title, fontsize=11)
    ax.grid(True, which="major", lw=0.3, color="#ebebeb")
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color("grey")
    leg = ax.legend(title="UFBoot support", loc="upper left", bbox_to_anchor=(0.08, 0.88),
                    fontsize=8.5, title_fontsize=9, framealpha=0.9, edgecolor="grey")
    leg.get_frame().set_linewidth(0.3)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = f"{args.output_prefix}.{ext}"
        fig.savefig(out, dpi=args.dpi if ext == "png" else None)
        print(f"Wrote {out}")
    fig.savefig  # keep reference alive

    # Console verification table (same as the R script's tail print)
    print("\nNode values (sorted by gCF):")
    print(f"{'ID':>4}  {'Label':>5}  {'gCF':>6}  {'sCF':>6}  {'gN':>4}  label")
    for r in sorted(rows, key=lambda r: r["_gCF"]):
        lbl = labels.get(r["_id"], "")
        print(f"{r['_id']:>4}  {r['_label']:>5.0f}  {r['_gCF']:>6.2f}  {r['_sCF']:>6.2f}  "
              f"{r.get('gN', 'NA'):>4}  {lbl}")


if __name__ == "__main__":
    main()