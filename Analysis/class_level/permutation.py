"""Multi-method aggregate permutation test for class effects.

The Fisher/Stouffer combination of per-method Kruskal p-values assumes
those tests are independent. The methods aren't strictly independent
(AF2 and AF3 share training data, rho ~ 0.68). To get a fully
distribution-free p-value that handles this correctly, we run a
permutation test:

  observed statistic T_obs = sum over methods of per-method
                              Kruskal-Wallis H on the per-pair statistic
                              stratified by trigger_class

We then permute the trigger_class labels at the *pair* level (so each
pair gets a single permuted class, and that class is broadcast to all
8 method-rows for that pair, preserving the within-pair correlation
structure across methods). The permutation null is the distribution of
T over n_perm such relabelings, and p = (#perms >= T_obs + 1) / (n_perm + 1).

The pair-level permutation correctly handles the across-method
correlation: AF2 and AF3 stay yoked to the same pair-id, so any
shared method bias is preserved in the null.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import kruskal


def _per_method_H_sum(
    stats: pd.DataFrame,
    value_col: str,
    class_col: str,
    method_col: str = "method",
    min_per_group: int = 2,
) -> float:
    """Sum of Kruskal-Wallis H statistics across methods.

    A method whose stratification is degenerate (one class with all values,
    or any class with <min_per_group observations) contributes 0 -- this is
    conservative because it shrinks H under both observed and permuted data.
    """
    total = 0.0
    for method in sorted(stats[method_col].unique()):
        sub = stats[stats[method_col] == method].dropna(subset=[value_col, class_col])
        groups = [g[value_col].values for _, g in sub.groupby(class_col)]
        if len(groups) < 2 or any(len(g) < min_per_group for g in groups):
            continue
        try:
            H, _ = kruskal(*groups)
            if np.isfinite(H):
                total += float(H)
        except Exception:
            continue
    return total


def permutation_test_class_effect(
    stats: pd.DataFrame,
    value_col: str,
    class_col: str = "trigger_class",
    method_col: str = "method",
    pair_col: str = "pair_id",
    n_perm: int = 10_000,
    seed: int = 42,
) -> dict:
    """Permute class labels at the pair level; recompute summed-H statistic.

    Returns
    -------
    dict with observed, null_mean, null_std, p_perm, n_perm
    """
    rng = np.random.default_rng(seed)
    sub = stats.dropna(subset=[value_col, class_col, pair_col]).copy()

    observed = _per_method_H_sum(sub, value_col, class_col, method_col)

    # Build the pair_id -> class mapping
    pair_class = (
        sub[[pair_col, class_col]]
        .drop_duplicates(pair_col)
        .set_index(pair_col)[class_col]
    )
    pairs = pair_class.index.to_numpy()
    # Force a real numpy object array even if pandas is using Arrow-backed
    # strings (which numpy.shuffle warns about and may not shuffle correctly).
    classes = np.asarray(pair_class.values, dtype=object)

    # For each permutation: shuffle the class assignment across pairs,
    # then broadcast to all method-rows by lookup.
    null_stats = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        perm_classes = classes.copy()
        rng.shuffle(perm_classes)
        perm_map = dict(zip(pairs, perm_classes))
        sub["_perm_class"] = sub[pair_col].map(perm_map)
        null_stats[i] = _per_method_H_sum(sub, value_col, "_perm_class", method_col)

    p_perm = (float(np.sum(null_stats >= observed)) + 1.0) / (n_perm + 1.0)
    return {
        "observed": float(observed),
        "null_mean": float(np.mean(null_stats)),
        "null_std": float(np.std(null_stats)),
        "p_perm": float(p_perm),
        "n_perm": int(n_perm),
    }
