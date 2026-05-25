"""Per-method Kruskal-Wallis omnibus tests + meta-analytic combination across
methods (Fisher and Stouffer).

Test logic. For each method m, run a Kruskal-Wallis test of the statistic
across the 5 trigger classes (omnibus: any class differs from any other).
Then combine the 8 per-method p-values using:
  - Fisher's combined test: chi^2 = -2 sum log p_m, df = 2k
  - Stouffer's Z method: Z = sum Phi^-1(1-p_m) / sqrt(k)
  - Optionally a Jonckheere-Terpstra ordered alternative test per method,
    combined the same way.

Fisher is slightly more sensitive to small-p outliers; Stouffer is more
robust. We report both. The methods are not strictly independent (AF2/AF3
share training data, rho = 0.68); Fisher is mildly anti-conservative under
weak positive correlation but the inflation is small.
"""

from __future__ import annotations

from typing import Dict

import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues, kruskal


def residualize_per_method(
    stats: pd.DataFrame,
    value_col: str,
    covariate_col: str,
    method_col: str = "method",
) -> pd.DataFrame:
    """Regress `value_col` on `covariate_col` separately for each method
    (per-method OLS), and replace value_col with the residuals. Returns
    a copy of `stats` with the residualized column.

    The per-method residualization lets us strip a covariate (e.g.
    seq_identity, MSA diversity) from each method's signal at its own
    scale before running the class comparison. Methods that have a
    near-zero covariate effect contribute essentially un-changed values;
    methods strongly confounded by the covariate get stripped.

    Returns a DataFrame with the same shape as stats; value_col is
    overwritten with residuals (NaN where the covariate or value is NaN).
    """
    sub = stats.copy()
    sub[f"_resid_{value_col}"] = np.nan
    for m, grp in sub.groupby(method_col):
        mask = grp[value_col].notna() & grp[covariate_col].notna()
        if mask.sum() < 5:
            continue
        x = grp.loc[mask, covariate_col].astype(float).values
        y = grp.loc[mask, value_col].astype(float).values
        # OLS slope and intercept
        x_mean = x.mean()
        y_mean = y.mean()
        denom = ((x - x_mean) ** 2).sum()
        if denom <= 0:
            sub.loc[grp.index[mask], f"_resid_{value_col}"] = y - y_mean
            continue
        beta = ((x - x_mean) * (y - y_mean)).sum() / denom
        alpha = y_mean - beta * x_mean
        resid = y - (alpha + beta * x)
        sub.loc[grp.index[mask], f"_resid_{value_col}"] = resid
    sub[value_col] = sub[f"_resid_{value_col}"]
    sub = sub.drop(columns=[f"_resid_{value_col}"])
    return sub


def per_method_kruskal(
    stats: pd.DataFrame,
    value_col: str,
    class_col: str = "trigger_class",
    method_col: str = "method",
    min_per_group: int = 2,
) -> pd.DataFrame:
    """Run a Kruskal-Wallis omnibus on `value_col` across trigger classes
    independently for each method. Returns one row per method with the
    test statistic, df-implied p, and the n per class.
    """
    rows = []
    for method in sorted(stats[method_col].unique()):
        sub = stats[stats[method_col] == method].dropna(subset=[value_col, class_col])
        groups = [g[value_col].values for _, g in sub.groupby(class_col)]
        ns = sub.groupby(class_col).size().to_dict()
        if len(groups) < 2 or any(len(g) < min_per_group for g in groups):
            rows.append({
                "method": method, "H": np.nan, "p": np.nan,
                "df": len(groups) - 1 if groups else 0,
                "ns": ns,
            })
            continue
        H, p = kruskal(*groups)
        rows.append({
            "method": method, "H": float(H), "p": float(p),
            "df": len(groups) - 1, "ns": ns,
        })
    return pd.DataFrame(rows)


def combine_methods(per_method: pd.DataFrame) -> Dict[str, float]:
    """Combine per-method Kruskal p-values via Fisher and Stouffer.
    Returns {fisher_stat, fisher_p, stouffer_stat, stouffer_p, k_methods}.
    """
    valid_p = per_method.dropna(subset=["p"])["p"].astype(float).values
    if len(valid_p) < 2:
        return {"fisher_stat": np.nan, "fisher_p": np.nan,
                "stouffer_stat": np.nan, "stouffer_p": np.nan,
                "k_methods": int(len(valid_p))}
    # scipy combine_pvalues: 'fisher' returns chi^2 stat; 'stouffer' returns Z
    fisher_stat, fisher_p = combine_pvalues(valid_p, method="fisher")
    stouffer_stat, stouffer_p = combine_pvalues(valid_p, method="stouffer")
    return {
        "fisher_stat": float(fisher_stat), "fisher_p": float(fisher_p),
        "stouffer_stat": float(stouffer_stat), "stouffer_p": float(stouffer_p),
        "k_methods": int(len(valid_p)),
    }


def jonckheere_terpstra_p(
    samples: list[np.ndarray],
    n_perm: int = 10_000,
    rng: np.random.Generator | None = None,
) -> tuple[float, float]:
    """Jonckheere-Terpstra one-sided ordered-alternative test (H1: ordered
    in the given sample order). Computed via permutation to avoid the
    asymptotic-normality assumption at small n.

    Returns (J statistic, one-sided permutation p-value).
    """
    if rng is None:
        rng = np.random.default_rng(0)

    def J(xs):
        s = 0
        for i in range(len(xs) - 1):
            xi = xs[i]
            for j in range(i + 1, len(xs)):
                xj = xs[j]
                # Count pairs (a, b) with a in xi, b in xj, a < b (ties = 0.5)
                a = xi[:, None]
                b = xj[None, :]
                s += int(np.sum(a < b)) + 0.5 * int(np.sum(a == b))
        return s

    obs = J(samples)

    # Permutation: pool everything, shuffle, slice back into the same group sizes
    pooled = np.concatenate(samples)
    sizes = [len(s) for s in samples]
    splits = np.cumsum(sizes)[:-1]

    count_ge = 0
    for _ in range(n_perm):
        rng.shuffle(pooled)
        new_groups = np.split(pooled, splits)
        if J(new_groups) >= obs:
            count_ge += 1

    p_one_sided = (count_ge + 1) / (n_perm + 1)
    return float(obs), float(p_one_sided)


def per_method_jonckheere(
    stats: pd.DataFrame,
    value_col: str,
    class_order: list[str],
    class_col: str = "trigger_class",
    method_col: str = "method",
    n_perm: int = 2000,
    seed: int = 42,
) -> pd.DataFrame:
    """Per-method Jonckheere-Terpstra test against the ordered alternative
    given by class_order (the order interpreted as monotone increase of the
    statistic). Returns one row per method.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for method in sorted(stats[method_col].unique()):
        sub = stats[stats[method_col] == method].dropna(subset=[value_col, class_col])
        samples = []
        present_order = []
        for cls in class_order:
            xs = sub[sub[class_col] == cls][value_col].values
            if len(xs) >= 2:
                samples.append(xs)
                present_order.append(cls)
        if len(samples) < 2:
            rows.append({"method": method, "J": np.nan, "p_jt": np.nan,
                         "class_order": present_order})
            continue
        # Use a fresh RNG slice per method to keep results reproducible per call
        method_rng = np.random.default_rng(seed + abs(hash(method)) % 10_000)
        J, p = jonckheere_terpstra_p(samples, n_perm=n_perm, rng=method_rng)
        rows.append({"method": method, "J": J, "p_jt": p,
                     "class_order": present_order})
    return pd.DataFrame(rows)
