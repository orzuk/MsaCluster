"""Orchestrator: run Fisher + LMM + permutation tests on 3 candidate
statistics (stddev_centered, stddev_residual, abs_beta) and save outputs
to docs/class_level/.

Usage:
    python3 -m Analysis.class_level.run_all
    python3 -m Analysis.class_level.run_all --n-perm 5000
    python3 -m Analysis.class_level.run_all --include-interaction
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import pandas as pd

from .build_stats import build_per_pair_per_method_stats, per_class_summary
from .combined_tests import (
    combine_methods,
    per_method_jonckheere,
    per_method_kruskal,
    residualize_per_method,
)
from .lmm import fit_class_lmm
from .permutation import permutation_test_class_effect


STATS_TO_TEST = [
    "stddev_centered",
    "stddev_residual",
    "abs_beta",
    "max_abs_centered",
    "max_abs_residual",
]

# Hypothesised increasing-signal order (Chakravarty 2023 + class-specific
# biology, see paper Methods Sec. 2.11 / SI Sec. supp-trigger-classes).
CLASS_ORDER_ASCENDING = [
    "equilibrium_or_unknown",
    "oligomerization",
    "mutation",
    "ligand",
    "protein_binding",
]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--survey", default="docs/fold_diversity_survey.csv")
    ap.add_argument("--corrected", default="docs/fold_diversity_survey_corrected.csv")
    ap.add_argument("--triggers", default="docs/triggers_from_pdb.csv")
    ap.add_argument("--outdir", default="docs/class_level")
    ap.add_argument("--n-perm", type=int, default=10_000,
                    help="Permutations for class-effect permutation test")
    ap.add_argument("--n-perm-jt", type=int, default=2000,
                    help="Permutations for per-method Jonckheere-Terpstra test")
    ap.add_argument("--min-clusters", type=int, default=5,
                    help="Drop pair-method cells with fewer clusters than this")
    ap.add_argument("--include-interaction", action="store_true",
                    help="Also fit class * method interaction in the LMM")
    ap.add_argument("--high-identity-threshold", type=float, default=0.95,
                    help="When --include-identity-controls is set, the "
                         "restricted-to-high-identity analysis uses this "
                         "minimum seq_identity cutoff (default 0.95).")
    ap.add_argument("--skip-identity-controls", action="store_true",
                    help="Skip the seq_identity-controlled re-analyses "
                         "(residualization + high-identity filter + LMM "
                         "with seq_identity covariate). Use to keep the "
                         "old-style class-only output.")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[class_level] loading inputs")
    stats = build_per_pair_per_method_stats(
        survey_csv=args.survey,
        corrected_csv=args.corrected,
        triggers_csv=args.triggers,
        min_clusters=args.min_clusters,
    )
    print(f"[class_level] per-pair x per-method stats: {stats.shape}")
    print(f"[class_level] classes: "
          f"{stats['trigger_class'].value_counts().to_dict()}")
    stats.to_csv(outdir / "per_pair_per_method_stats.csv", index=False)

    # Per-class summaries
    summary_frames = []
    for s in STATS_TO_TEST:
        if s not in stats.columns:
            continue
        s_summary = per_class_summary(stats, s)
        s_summary["statistic"] = s
        summary_frames.append(s_summary)
    if summary_frames:
        per_class_summary_df = pd.concat(summary_frames, ignore_index=True)
        per_class_summary_df.to_csv(outdir / "per_class_summary.csv", index=False)

    def _run_tests_on(stats_df: pd.DataFrame, mode_tag: str,
                      include_lmm_seq_identity_covariate: bool = False) -> dict:
        """Run the full Kruskal + JT + LMM + permutation battery on a
        given stats DataFrame (which may be the raw stats, residualised
        stats, or filtered stats). Returns dict {statistic: results}.
        """
        out: dict = {}
        for s in STATS_TO_TEST:
            if s not in stats_df.columns:
                continue
            if stats_df[s].notna().sum() < 20:
                # Not enough non-NaN observations to run anything meaningful
                continue
            per_method_kw = per_method_kruskal(stats_df, s)
            per_method_kw["statistic"] = s
            per_method_kw["test"] = "kruskal_wallis"
            per_method_kw["mode"] = mode_tag
            per_method_jt = per_method_jonckheere(
                stats_df, s, CLASS_ORDER_ASCENDING,
                n_perm=args.n_perm_jt, seed=42,
            )
            per_method_jt["statistic"] = s
            per_method_jt["test"] = "jonckheere_terpstra"
            per_method_jt["mode"] = mode_tag
            per_method_rows.append(per_method_kw)
            per_method_rows.append(per_method_jt)

            combined_kw = combine_methods(per_method_kw)
            jt_for_combine = per_method_jt.rename(columns={"p_jt": "p"})[
                ["method", "p"]
            ]
            combined_jt = combine_methods(jt_for_combine)

            covariates = (
                ["seq_identity"]
                if include_lmm_seq_identity_covariate
                and "seq_identity" in stats_df.columns
                else None
            )
            lmm_result = fit_class_lmm(
                stats_df, response_col=s,
                include_interaction=args.include_interaction,
                covariates=covariates,
            )
            perm_result = permutation_test_class_effect(
                stats_df, value_col=s, n_perm=args.n_perm, seed=42,
            )
            out[s] = {
                "mode": mode_tag,
                "n_pairs": int(stats_df["pair_id"].nunique()),
                "n_obs": int(stats_df[s].notna().sum()),
                "combined_kruskal_fisher_p": combined_kw["fisher_p"],
                "combined_kruskal_stouffer_p": combined_kw["stouffer_p"],
                "combined_jonckheere_fisher_p": combined_jt["fisher_p"],
                "combined_jonckheere_stouffer_p": combined_jt["stouffer_p"],
                "lmm_class_lrt_p": lmm_result.get("lrt_p"),
                "lmm_class_lrt_df": lmm_result.get("lrt_df"),
                "lmm_covariates": ",".join(covariates) if covariates else "",
                "permutation_p": perm_result["p_perm"],
                "permutation_observed_sumH": perm_result["observed"],
                "permutation_null_mean": perm_result["null_mean"],
                "lmm_class_coefs": lmm_result.get("class_coefs", []),
                "lmm_converged": lmm_result.get("converged"),
                "lmm_interaction_p": lmm_result.get("interaction_p"),
                "lmm_error": lmm_result.get("error"),
            }
            print(f"[class_level] {mode_tag}/{s}: "
                  f"Fisher_p={combined_kw['fisher_p']:.3g}  "
                  f"perm_p={perm_result['p_perm']:.3g}  "
                  f"lmm_lrt_p={lmm_result.get('lrt_p')}")
        return out

    # === Mode 1: raw, no seq_identity controls ===
    print(f"\n[class_level] ====================================")
    print(f"[class_level] MODE 1: raw (no seq_identity controls)")
    print(f"[class_level] ====================================")
    per_method_rows: list = []
    results_raw = _run_tests_on(stats, mode_tag="raw")

    # === Mode 2: per-method residualization on seq_identity ===
    results_resid: dict = {}
    results_filtered: dict = {}
    results_lmm_covariate: dict = {}
    if not args.skip_identity_controls and "seq_identity" in stats.columns:
        print(f"\n[class_level] ====================================")
        print(f"[class_level] MODE 2: residualised on seq_identity")
        print(f"[class_level] ====================================")
        stats_resid = stats.copy()
        for s in STATS_TO_TEST:
            if s in stats_resid.columns:
                stats_resid = residualize_per_method(
                    stats_resid, s, covariate_col="seq_identity",
                )
        results_resid = _run_tests_on(stats_resid, mode_tag="resid_on_seq_identity")

        # === Mode 3: filtered to high seq_identity ===
        print(f"\n[class_level] ====================================")
        thr = args.high_identity_threshold
        print(f"[class_level] MODE 3: filtered to seq_identity >= {thr}")
        print(f"[class_level] ====================================")
        stats_high = stats[stats["seq_identity"] >= thr].copy()
        print(f"[class_level] kept {stats_high['pair_id'].nunique()} pairs "
              f"out of {stats['pair_id'].nunique()} "
              f"({len(stats_high)} pair-method cells)")
        results_filtered = _run_tests_on(
            stats_high, mode_tag=f"identity_ge_{thr}"
        )

        # === Mode 4: LMM with seq_identity as covariate (raw data) ===
        # This re-uses raw stats but adds seq_identity to the LMM null and
        # full formulas. We re-run only the LMM (Kruskal/perm are unchanged
        # by adding an LMM covariate, so we skip re-running those).
        print(f"\n[class_level] ====================================")
        print(f"[class_level] MODE 4: LMM with seq_identity covariate")
        print(f"[class_level] ====================================")
        for s in STATS_TO_TEST:
            if s not in stats.columns or stats[s].notna().sum() < 20:
                continue
            lmm_r = fit_class_lmm(
                stats, response_col=s,
                include_interaction=args.include_interaction,
                covariates=["seq_identity"],
            )
            results_lmm_covariate[s] = {
                "mode": "lmm_seq_identity_covariate",
                "lmm_class_lrt_p": lmm_r.get("lrt_p"),
                "lmm_class_lrt_df": lmm_r.get("lrt_df"),
                "lmm_class_coefs": lmm_r.get("class_coefs", []),
                "lmm_converged": lmm_r.get("converged"),
                "lmm_interaction_p": lmm_r.get("interaction_p"),
                "lmm_error": lmm_r.get("error"),
            }
            print(f"[class_level] lmm_cov/{s}: "
                  f"lmm_lrt_p={lmm_r.get('lrt_p')}  "
                  f"converged={lmm_r.get('converged')}")

    results = {
        **{f"raw__{k}": v for k, v in results_raw.items()},
        **{f"resid__{k}": v for k, v in results_resid.items()},
        **{f"highid__{k}": v for k, v in results_filtered.items()},
        **{f"lmmcov__{k}": v for k, v in results_lmm_covariate.items()},
    }

    # Save outputs
    if per_method_rows:
        per_method_df = pd.concat(per_method_rows, ignore_index=True, sort=False)
        per_method_df.to_csv(outdir / "per_method_tests.csv", index=False)

    summary_rows = []
    for key, r in results.items():
        # key is "<mode_prefix>__<statistic>", split it back out
        if "__" in key:
            mode_prefix, statistic = key.split("__", 1)
        else:
            mode_prefix, statistic = "", key
        row = {"mode_prefix": mode_prefix, "statistic": statistic}
        for k, v in r.items():
            if isinstance(v, list):
                continue
            row[k] = v
        summary_rows.append(row)
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(outdir / "summary.csv", index=False)

    with open(outdir / "summary.json", "w") as f:
        json.dump({k: v for k, v in results.items()}, f, indent=2, default=str)

    print(f"\n[class_level] outputs written to {outdir}/")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
