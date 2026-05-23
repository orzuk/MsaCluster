#!/usr/bin/env python3
"""Rank pairs by strength of clade-clustering signal on the cluster tree.

Two ranking metrics, computed from docs/phylo_placement.csv:
  - best_p   = minimum (parsimony_p, nn_concordance_p, nn3_concordance_p)
               across all methods reporting that statistic for the pair.
  - best_eff = peak NN-concordance achieved by any method on the pair,
               restricted to rows that classified at least 2 F1c AND 2 F2c
               leaves (otherwise the test has no power and the score is
               vacuously 1.0 -- we exclude those).

Outputs (stdout):
  Top-20 pairs by raw-label best_p (clade structure on uncorrected labels)
  Top-20 pairs by corrected-label best_p (clade structure that survives
    the sequence-divergence residual correction).
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

DOCS = "docs"
TOP_N = 20
MIN_F1c = 2
MIN_F2c = 2


def load(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    for col in ("parsimony_p", "nn_concordance_p", "nn3_concordance_p",
                "nn_concordance"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ("n_F1c_leaves", "n_F2c_leaves"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def rank(df: pd.DataFrame) -> pd.DataFrame:
    powered = df[(df["n_F1c_leaves"] >= MIN_F1c) & (df["n_F2c_leaves"] >= MIN_F2c)].copy()
    powered["min_p"] = powered[["parsimony_p", "nn_concordance_p",
                                 "nn3_concordance_p"]].min(axis=1, skipna=True)
    pair_min_p = powered.groupby("pair_id")["min_p"].min().rename("best_p")
    pair_max_nn = powered.groupby("pair_id")["nn_concordance"].max().rename("best_NN")
    pair_methods_sig = (powered[powered["min_p"] <= 0.05]
                         .groupby("pair_id")["method"].nunique()
                         .rename("n_methods_p_le_0p05"))
    pair_F1c_max = powered.groupby("pair_id")["n_F1c_leaves"].max().rename("max_F1c")
    pair_F2c_max = powered.groupby("pair_id")["n_F2c_leaves"].max().rename("max_F2c")
    pair_method_best = (powered.sort_values("min_p")
                          .groupby("pair_id").first()["method"]
                          .rename("best_method"))
    out = pd.concat([pair_min_p, pair_max_nn, pair_methods_sig,
                     pair_F1c_max, pair_F2c_max, pair_method_best], axis=1)
    out["n_methods_p_le_0p05"] = out["n_methods_p_le_0p05"].fillna(0).astype(int)
    return out.sort_values("best_p")


def show(out: pd.DataFrame, label: str) -> None:
    print(f"\n=== {label}: Top {TOP_N} pairs by clade-clustering strength ===")
    print(f"{'pair_id':<22} {'best_p':>8} {'best_NN':>8} "
          f"{'F1c':>4} {'F2c':>4} {'#sig':>5}  best_method")
    print("-" * 80)
    for pid, row in out.head(TOP_N).iterrows():
        bp = row["best_p"]
        nn = row["best_NN"]
        f1 = int(row["max_F1c"]); f2 = int(row["max_F2c"])
        nsig = int(row["n_methods_p_le_0p05"])
        m = row["best_method"]
        bp_s = f"{bp:.4f}" if np.isfinite(bp) else "  -  "
        nn_s = f"{nn:.3f}" if np.isfinite(nn) else "  -  "
        print(f"{pid:<22} {bp_s:>8} {nn_s:>8} {f1:>4d} {f2:>4d} {nsig:>5d}  {m}")


def main() -> None:
    raw_path = os.path.join(DOCS, "phylo_placement.csv")
    cor_path = os.path.join(DOCS, "phylo_placement_corrected.csv")
    if not os.path.isfile(raw_path) or not os.path.isfile(cor_path):
        print("[error] missing phylo_placement csvs in docs/", file=sys.stderr)
        sys.exit(2)
    raw = load(raw_path)
    cor = load(cor_path)
    show(rank(raw), "RAW labels")
    show(rank(cor), "RESIDUAL-CORRECTED labels")


if __name__ == "__main__":
    main()
