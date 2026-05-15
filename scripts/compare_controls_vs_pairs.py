"""Intra-cluster diversity: fold-switching pairs vs. negative controls.

Reads (after the controls aggregator pass finishes):
  - docs/foldpair_intra_diversity.csv
        per-pair mean/median/min/max pairwise TM-score across cluster prediction PDBs
  - Pipeline/Controls/<CT>/controls_diversity_summary.csv      (one row per control protein)
        per-control-protein equivalent stat for each of three control sets

Prints (and writes to docs/controls_vs_pairs_comparison.csv):
  - per-group n, mean, median, std of mean_pairwise_TM
  - Mann-Whitney U p-value: fold-switch pairs vs. each control set
  - fraction of fold-switch pairs below the 5% quantile of each control distribution

Usage (after all 249 controls finish + aggregator pass):
    source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
    for ct in PorterIDP MatchedSingleFold RandomSingleFold; do
        python3 run_controls_pipeline.py --run_mode postprocess --protein_ids ALL --control_type "$ct"
    done
    python3 scripts/compare_controls_vs_pairs.py
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

try:
    from scipy.stats import mannwhitneyu
except ImportError:
    mannwhitneyu = None  # type: ignore[assignment]

REPO = Path(".")
PAIRS_CSV = REPO / "docs" / "foldpair_intra_diversity.csv"
CONTROL_TYPES = ["PorterIDP", "MatchedSingleFold", "RandomSingleFold"]
OUT_CSV = REPO / "docs" / "controls_vs_pairs_comparison.csv"
TM_COL = "mean_pairwise_TM"


def load_pairs() -> pd.DataFrame:
    if not PAIRS_CSV.exists():
        raise SystemExit(f"missing {PAIRS_CSV}")
    df = pd.read_csv(PAIRS_CSV)
    df["group"] = "FoldSwitchPairs"
    df = df.rename(columns={"pair_id": "id"})
    return df[["id", "group", TM_COL]]


def load_controls() -> dict[str, pd.DataFrame]:
    out: dict[str, pd.DataFrame] = {}
    for ct in CONTROL_TYPES:
        csv = REPO / "Pipeline" / "Controls" / ct / "controls_diversity_summary.csv"
        if not csv.exists():
            print(f"  WARN: missing {csv}; skipping {ct}")
            continue
        df = pd.read_csv(csv)
        if TM_COL not in df.columns:
            cand = [c for c in df.columns if "TM" in c and "mean" in c.lower()]
            if not cand:
                print(f"  WARN: {ct} has no mean-pairwise-TM-like column ({df.columns.tolist()})")
                continue
            df = df.rename(columns={cand[0]: TM_COL})
        df["group"] = ct
        id_col = "protein_id" if "protein_id" in df.columns else df.columns[0]
        df = df.rename(columns={id_col: "id"})
        out[ct] = df[["id", "group", TM_COL]].dropna(subset=[TM_COL])
    return out


def summarize(label: str, vals: np.ndarray) -> dict:
    return {
        "group": label,
        "n": int(len(vals)),
        "mean": float(np.mean(vals)) if len(vals) else float("nan"),
        "median": float(np.median(vals)) if len(vals) else float("nan"),
        "std": float(np.std(vals, ddof=1)) if len(vals) > 1 else float("nan"),
        "q05": float(np.quantile(vals, 0.05)) if len(vals) else float("nan"),
        "q95": float(np.quantile(vals, 0.95)) if len(vals) else float("nan"),
    }


def main() -> None:
    pairs = load_pairs()
    controls = load_controls()

    pair_vals = pairs[TM_COL].dropna().to_numpy()
    rows = [summarize("FoldSwitchPairs", pair_vals)]
    print("\n=== Per-group summary (mean pairwise TM-score across cluster predictions) ===")
    print(f"{'group':<22} {'n':>5} {'mean':>7} {'median':>8} {'std':>7} {'q05':>7} {'q95':>7}")
    print(f"{'FoldSwitchPairs':<22} {len(pair_vals):>5d} {np.mean(pair_vals):>7.3f} "
          f"{np.median(pair_vals):>8.3f} {np.std(pair_vals, ddof=1):>7.3f} "
          f"{np.quantile(pair_vals, 0.05):>7.3f} {np.quantile(pair_vals, 0.95):>7.3f}")

    for ct, df in controls.items():
        vals = df[TM_COL].dropna().to_numpy()
        s = summarize(ct, vals)
        rows.append(s)
        print(f"{ct:<22} {len(vals):>5d} {np.mean(vals):>7.3f} "
              f"{np.median(vals):>8.3f} {np.std(vals, ddof=1):>7.3f} "
              f"{np.quantile(vals, 0.05):>7.3f} {np.quantile(vals, 0.95):>7.3f}")

    print("\n=== Pairs vs. each control: contrast tests ===")
    print(f"{'control':<22} {'mean_pairs':>11} {'mean_control':>14} "
          f"{'delta':>8} {'MW_U_p':>10} {'%pairs<q05(ctrl)':>18}")
    for ct, df in controls.items():
        cvals = df[TM_COL].dropna().to_numpy()
        if not len(cvals):
            continue
        # Mann-Whitney U, two-sided (no assumption about direction; we report mean delta)
        if mannwhitneyu is not None:
            p = mannwhitneyu(pair_vals, cvals, alternative="two-sided").pvalue
        else:
            p = float("nan")
        q05 = float(np.quantile(cvals, 0.05))
        below = float(np.mean(pair_vals < q05))
        delta = float(np.mean(pair_vals) - np.mean(cvals))
        print(f"{ct:<22} {np.mean(pair_vals):>11.3f} {np.mean(cvals):>14.3f} "
              f"{delta:>+8.3f} {p:>10.2e} {below:>18.2%}")

    out_df = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT_CSV, index=False)
    print(f"\nwrote {OUT_CSV}")


if __name__ == "__main__":
    main()
