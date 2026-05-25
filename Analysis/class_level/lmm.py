"""Linear mixed model on long-format per-pair x per-method statistics.

Model:
  response ~ C(trigger_class) + C(method) + (1 | pair_id)

This properly pools across methods while accounting for pair-level random
intercepts. The F-test on the C(trigger_class) fixed effect gives a single
class-effect p-value.

We additionally fit the interaction model
  response ~ C(trigger_class) * C(method) + (1 | pair_id)
to test whether class effects are themselves method-specific.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd


def fit_class_lmm(
    stats: pd.DataFrame,
    response_col: str,
    class_col: str = "trigger_class",
    method_col: str = "method",
    pair_col: str = "pair_id",
    include_interaction: bool = False,
    covariates: list[str] | None = None,
) -> dict:
    """Fit a mixed-effects model and return a summary dict.

    Uses statsmodels mixedlm with a random intercept for pair_id. The
    class fixed-effect significance is reported as the p-value of the
    likelihood-ratio test between the full model and the no-class
    null model.

    Returns
    -------
    dict with keys:
      response, n_obs, n_pairs, n_methods,
      llf_full, llf_null, lrt_stat, lrt_df, lrt_p,
      class_coefs (DataFrame: per-class coef, SE, z, p),
      interaction_p (only if include_interaction=True),
      converged (bool)
    """
    try:
        import statsmodels.formula.api as smf
    except ImportError as e:
        return {"error": f"statsmodels not installed: {e}"}

    drop_cols = [response_col, class_col, method_col, pair_col]
    if covariates:
        drop_cols += list(covariates)
    sub = stats.dropna(subset=drop_cols).copy()
    if sub.empty:
        return {"error": "no rows after dropping NA"}

    # Ensure a reproducible reference category (alphabetical first) so the
    # individual coefficients are interpretable; reference is conventionally
    # the "control" class -- we'll use equilibrium_or_unknown as ref.
    sub[class_col] = pd.Categorical(
        sub[class_col],
        categories=["equilibrium_or_unknown", "ligand", "mutation",
                    "oligomerization", "protein_binding"],
        ordered=False,
    )

    def _fit(formula: str):
        return smf.mixedlm(formula, sub, groups=sub[pair_col]).fit(
            method="lbfgs", reml=False
        )

    # Optional covariate terms (e.g. "seq_identity") added to both null and
    # full models -- the LRT then tests class effect *after* controlling
    # for the covariates.
    cov_term = (" + " + " + ".join(covariates)) if covariates else ""
    formula_null = f"{response_col} ~ C({method_col}){cov_term}"
    formula_full = f"{response_col} ~ C({class_col}) + C({method_col}){cov_term}"

    try:
        m_null = _fit(formula_null)
        m_full = _fit(formula_full)
    except Exception as e:
        return {"error": f"mixedlm fit failed: {e}"}

    converged = bool(getattr(m_full, "converged", True))

    lrt_stat = 2.0 * (m_full.llf - m_null.llf)
    # df = number of class levels - 1
    n_class = sub[class_col].nunique()
    lrt_df = n_class - 1
    from scipy.stats import chi2
    lrt_p = float(1.0 - chi2.cdf(lrt_stat, df=lrt_df)) if lrt_df > 0 else np.nan

    # Per-class coefficient table from the full model
    class_coefs_rows = []
    for name in m_full.params.index:
        if name.startswith(f"C({class_col})"):
            class_coefs_rows.append({
                "term": name,
                "coef": float(m_full.params[name]),
                "se": float(m_full.bse[name]),
                "z": float(m_full.tvalues[name]),
                "p": float(m_full.pvalues[name]),
            })

    out = {
        "response": response_col,
        "n_obs": int(len(sub)),
        "n_pairs": int(sub[pair_col].nunique()),
        "n_methods": int(sub[method_col].nunique()),
        "llf_full": float(m_full.llf),
        "llf_null": float(m_null.llf),
        "lrt_stat": float(lrt_stat),
        "lrt_df": int(lrt_df),
        "lrt_p": lrt_p,
        "class_coefs": class_coefs_rows,
        "converged": converged,
        "reference_class": "equilibrium_or_unknown",
    }

    if include_interaction:
        formula_int = f"{response_col} ~ C({class_col}) * C({method_col})"
        try:
            m_int = _fit(formula_int)
            int_lrt = 2.0 * (m_int.llf - m_full.llf)
            n_methods = sub[method_col].nunique()
            int_df = (n_class - 1) * (n_methods - 1)
            int_p = float(1.0 - chi2.cdf(int_lrt, df=int_df)) if int_df > 0 else np.nan
            out["interaction_lrt"] = float(int_lrt)
            out["interaction_df"] = int(int_df)
            out["interaction_p"] = int_p
        except Exception as e:
            out["interaction_error"] = str(e)

    return out
