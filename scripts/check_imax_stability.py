#!/usr/bin/env python3
"""
Stability of the Moran's-I normalisation (I_norm = (I - E) / (I_max - E)).

I_max is found by a capped 2-opt local search started from the spectral guess
and from the data (utils/phylo_tests.morans_I_test), so it is an *estimate* of
the attainable maximum.  For every (pair, method) cell with I_norm >= --min-inorm
this re-runs the search from many random starts with a much larger step cap
and reports the best I_max found, the re-normalised I_norm, and whether the
"sharp" call (I_norm >= 0.5) changes.

  python3 scripts/check_imax_stability.py \
      --phylo-csv docs/phylo_placement_corrected.csv \
      --survey-csv docs/fold_diversity_survey_corrected_9method.csv \
      --value-col TMdiff_residual --out docs/imax_stability.csv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.ancestral_utils import build_coarse_tree, _name_internal_nodes  # noqa: E402
from utils.phylo_tests import _patristic_matrix  # noqa: E402
from scripts.phylo_placement import load_fold_value, normalize_cluster  # noqa: E402


def two_opt_max(M, arr, max_iter):
    arr = np.array(arr, dtype=float)
    for _ in range(max_iter):
        f = M @ arr
        DZ = arr[None, :] - arr[:, None]
        DF = f[:, None] - f[None, :]
        dN = 2.0 * DZ * DF - 2.0 * (DZ ** 2) * M
        np.fill_diagonal(dN, -np.inf)
        a, b = np.unravel_index(int(np.argmax(dN)), dN.shape)
        if dN[a, b] <= 1e-12:
            break
        arr[a], arr[b] = arr[b], arr[a]
    return float(arr @ (M @ arr))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phylo-csv", required=True)
    ap.add_argument("--survey-csv", required=True)
    ap.add_argument("--value-col", default="TMdiff_residual")
    ap.add_argument("--min-inorm", type=float, default=0.3)
    ap.add_argument("--min-n", type=int, default=20)
    ap.add_argument("--n-starts", type=int, default=50)
    ap.add_argument("--max-iter", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", default="docs/imax_stability.csv")
    args = ap.parse_args()

    ph = pd.read_csv(args.phylo_csv)
    cells = ph[(ph.morans_I_norm >= args.min_inorm) & (ph.morans_n >= args.min_n)]
    print(f"{len(cells)} cells with I_norm >= {args.min_inorm} and n >= {args.min_n}")
    values = load_fold_value(args.survey_csv, args.value_col)
    rng = np.random.default_rng(args.seed)
    rows, trees = [], {}
    for _, c in cells.iterrows():
        pid, meth = c.pair_id, c.method
        if pid not in trees:
            tree, tag_to_label, _ = build_coarse_tree(pid)
            _name_internal_nodes(tree)
            trees[pid] = (tree, tag_to_label)
        tree, tag_to_label = trees[pid]
        lv = {lab: values.get((pid, meth, normalize_cluster(tag))) for tag, lab in tag_to_label.items()}
        names = [k for k, v in lv.items() if v is not None and np.isfinite(v)]
        x = np.array([lv[k] for k in names]); xc = x - x.mean(); den = float(xc @ xc)
        W = 1.0 / (_patristic_matrix(tree, names) + 1e-9); np.fill_diagonal(W, 0.0)
        rs = W.sum(axis=1, keepdims=True); rs[rs == 0] = 1.0; W = W / rs
        M = 0.5 * (W + W.T)
        I_obs = float(xc @ (W @ xc)) / den
        if abs(I_obs - float(c.morans_I)) > 1e-3:
            print(f"[warn] {pid} {meth}: recomputed I {I_obs:.4f} != stored {c.morans_I:.4f}", flush=True)
        E = -1.0 / (len(xc) - 1)
        best = max(two_opt_max(M, rng.permutation(xc), args.max_iter) for _ in range(args.n_starts)) / den
        best = max(best, two_opt_max(M, xc, args.max_iter) / den)
        stored = float(c.morans_max_perm)
        imax = max(best, stored)
        inorm_new = (I_obs - E) / (imax - E)
        rows.append({"pair_id": pid, "method": meth, "n": len(xc), "I": round(I_obs, 4),
                     "I_max_stored": stored, "I_max_multistart": round(best, 4),
                     "rel_gain": round((imax - stored) / (stored - E), 4),
                     "I_norm_stored": c.morans_I_norm, "I_norm_multistart": round(inorm_new, 4),
                     "sharp_stored": bool(c.morans_I_norm >= 0.5 and c.morans_I > 0),
                     "sharp_multistart": bool(inorm_new >= 0.5 and I_obs > 0)})
        print(f"{pid} {meth}: I_max {stored:.4f} -> {imax:.4f}; I_norm {c.morans_I_norm:.3f} -> {inorm_new:.3f}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    if len(df):
        print(f"\nmedian relative I_max gain {df.rel_gain.median():.3f} (max {df.rel_gain.max():.3f}); "
              f"sharp calls changed: {(df.sharp_stored != df.sharp_multistart).sum()} of {len(df)}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
