"""Build the per-pair x per-method statistics table that downstream tests
consume. One row per (pair_id, method); columns are the candidate per-pair
statistics that could differ across trigger classes.

Inputs (CSV paths default to repo docs/):
  fold_diversity_survey.csv             - per-cluster TMdiff_centered per method
  fold_diversity_survey_corrected.csv   - same plus residual + regression_beta
  triggers_from_pdb.csv                 - per-pair trigger_class + seq_identity

Output columns (per (pair_id, method)):
  pair_id, method, trigger_class, seq_identity, n_clusters,
  stddev_centered     - within-pair stddev of TMdiff_centered
  stddev_residual     - within-pair stddev of TMdiff_residual (post-correction)
  max_abs_centered    - max |TMdiff_centered| within the pair
  max_abs_residual    - max |TMdiff_residual| within the pair
  abs_beta            - |regression slope| (PID confound magnitude)
"""

from __future__ import annotations

import os
from typing import Optional

import numpy as np
import pandas as pd


def build_per_pair_per_method_stats(
    survey_csv: str,
    corrected_csv: Optional[str] = None,
    triggers_csv: Optional[str] = None,
    min_clusters: int = 5,
) -> pd.DataFrame:
    """Aggregate per-cluster preference values to per-(pair, method) statistics.

    A pair-method cell is kept only if it has >= min_clusters non-missing
    TMdiff_centered values.
    """
    survey = pd.read_csv(survey_csv)
    # TMdiff_centered may be present for some methods (TM-based) but encoded
    # differently for DDG / S4PRED via TMdiff_centered. Trust the survey.
    survey = survey.dropna(subset=["TMdiff_centered"]).copy()
    survey["abs_centered"] = survey["TMdiff_centered"].abs()

    agg_centered = (
        survey.groupby(["pair_id", "method"])
        .agg(
            n_clusters=("TMdiff_centered", "size"),
            stddev_centered=("TMdiff_centered", "std"),
            max_abs_centered=("abs_centered", "max"),
        )
        .reset_index()
    )

    out = agg_centered

    if corrected_csv and os.path.isfile(corrected_csv):
        corr = pd.read_csv(corrected_csv)
        if "TMdiff_residual" in corr.columns:
            corr = corr.dropna(subset=["TMdiff_residual"]).copy()
            corr["abs_residual"] = corr["TMdiff_residual"].abs()
            has_beta = "regression_beta" in corr.columns
            agg_dict = dict(
                stddev_residual=("TMdiff_residual", "std"),
                max_abs_residual=("abs_residual", "max"),
            )
            if has_beta:
                # regression_beta is the same for all clusters within a
                # pair-method, so first() pulls one row's value.
                agg_dict["regression_beta"] = ("regression_beta", "first")
            agg_resid = (
                corr.groupby(["pair_id", "method"])
                .agg(**agg_dict)
                .reset_index()
            )
            if has_beta:
                agg_resid["abs_beta"] = agg_resid["regression_beta"].abs()
                agg_resid = agg_resid.drop(columns=["regression_beta"])
            else:
                # No regression_beta column -> abs_beta is unavailable.
                # Downstream tests that need abs_beta will be skipped.
                print("[build_stats] NOTE: corrected CSV has TMdiff_residual "
                      "but no regression_beta column; abs_beta tests skipped.")
            out = out.merge(agg_resid, on=["pair_id", "method"], how="left")

    if triggers_csv and os.path.isfile(triggers_csv):
        trig = pd.read_csv(triggers_csv)
        keep = ["pair_id", "trigger_class"]
        if "seq_identity" in trig.columns:
            keep.append("seq_identity")
        trig = trig[keep].drop_duplicates("pair_id")
        out = out.merge(trig, on="pair_id", how="left")

    # Filter: drop pair-method cells with too few clusters
    out = out[out["n_clusters"] >= min_clusters].reset_index(drop=True)

    # Drop pairs we couldn't classify (unresolved_no_pdb etc.)
    if "trigger_class" in out.columns:
        out = out[
            out["trigger_class"].isin(
                ["mutation", "ligand", "oligomerization",
                 "protein_binding", "equilibrium_or_unknown"]
            )
        ].reset_index(drop=True)

    return out


def per_class_summary(stats: pd.DataFrame, value_col: str) -> pd.DataFrame:
    """Median, IQR, n per (method, trigger_class) for a given statistic column."""
    g = stats.groupby(["method", "trigger_class"])[value_col]
    summary = pd.DataFrame({
        "n_pairs": g.size(),
        "median": g.median(),
        "q1": g.quantile(0.25),
        "q3": g.quantile(0.75),
    }).reset_index()
    return summary
