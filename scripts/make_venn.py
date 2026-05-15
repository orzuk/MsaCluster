"""Regenerate the four-way Venn over centered-diverse pair sets.

Reads `docs/fold_diversity_summary.csv` (one row per pair x method, with
`has_both_folds` and `has_both_folds_centered` columns) and builds four
sets of pair IDs, one per method, of pairs flagged as centered-diverse.
Writes `docs/figs/vend.png` (4-way Venn) and `docs/figs/vend_counts.csv`
(all 15 intersection counts) for the paper figure `figs/vend.png`.

Usage:
    python3 scripts/make_venn.py
    python3 scripts/make_venn.py --classification raw   # use absolute thresholds
"""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

METHODS = ["AF2", "ESM", "MSAT", "DDG"]


def load_method_sets(survey_csv: Path, classification: str) -> dict[str, set[str]]:
    flag_col = "has_both_folds_centered" if classification == "centered" else "has_both_folds"
    df = pd.read_csv(survey_csv)
    if flag_col not in df.columns:
        raise SystemExit(f"column {flag_col!r} missing from {survey_csv}")
    sets: dict[str, set[str]] = {}
    for m in METHODS:
        sub = df[(df["method"] == m) & (df[flag_col].astype(bool))]
        sets[m] = set(sub["pair_id"].astype(str))
    return sets


def all_intersections(sets: dict[str, set[str]]) -> pd.DataFrame:
    rows = []
    for r in range(1, len(sets) + 1):
        for combo in combinations(sets.keys(), r):
            inter = set.intersection(*[sets[m] for m in combo])
            rows.append({
                "methods": "&".join(combo),
                "n_methods": r,
                "count": len(inter),
                "pair_ids": ";".join(sorted(inter)) if inter else "",
            })
    return pd.DataFrame(rows)


def draw_venn4(sets: dict[str, set[str]], out_png: Path, title: str) -> None:
    """Four-set Venn using `venn` package if available, else UpSet-style bar plot."""
    try:
        from venn import venn  # type: ignore[import-not-found]
        fig, ax = plt.subplots(figsize=(7, 7))
        venn(sets, ax=ax, fmt="{size}")
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(out_png, dpi=200, bbox_inches="tight")
        plt.close(fig)
        return
    except ImportError:
        pass
    _draw_upset(sets, out_png, title)


def _draw_upset(sets: dict[str, set[str]], out_png: Path, title: str) -> None:
    counts = all_intersections(sets)
    counts = counts[counts["count"] > 0].sort_values("count", ascending=True)
    fig, ax = plt.subplots(figsize=(8, max(4, 0.35 * len(counts))))
    ax.barh(counts["methods"], counts["count"], color="#4C72B0")
    for i, (label, val) in enumerate(zip(counts["methods"], counts["count"])):
        ax.text(val, i, f"  {val}", va="center", fontsize=9)
    ax.set_xlabel("number of pairs")
    ax.set_title(title + "  (UpSet-style fallback; install `venn` for true 4-set Venn)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--survey-csv", type=Path, default=Path("docs/fold_diversity_summary.csv"))
    parser.add_argument("--out-dir", type=Path, default=Path("docs/figs"))
    parser.add_argument("--classification", choices=("centered", "raw"), default="centered")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    sets = load_method_sets(args.survey_csv, args.classification)

    title = f"Pairs with cluster-level fold diversity ({args.classification} classification)"
    for m in METHODS:
        print(f"  {m}: {len(sets[m])} pairs flagged")

    counts = all_intersections(sets)
    counts_csv = args.out_dir / f"vend_counts_{args.classification}.csv"
    counts.to_csv(counts_csv, index=False)
    print(f"wrote {counts_csv}")

    out_png = args.out_dir / "vend.png" if args.classification == "centered" \
              else args.out_dir / "vend_raw.png"
    draw_venn4(sets, out_png, title)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
