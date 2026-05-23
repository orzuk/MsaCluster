#!/usr/bin/env python3
"""Permutation-null test for per-pair cross-method concordance and
centered-diversity statistics produced by ``scripts/fold_diversity_survey.py``.

Why
---
The current concordance results report per-pair statistics like
``mean_concordance``, ``has_both_folds_centered``, and pairwise sign-agreement
percentages between methods (AF2~ESM, AF2~DDG, etc.). The manuscript currently
attaches a normal-approximation binomial p-value (``_binom_p_two_sided`` in
``fold_diversity_survey.py``) for the per-method-pair sign-agreement, but no
proper null for the *aggregate* statistics (mean concordance across method-pairs
and the ``has_both_folds_centered`` diversity call).

This script provides an empirical permutation-null:
   if you shuffle which cluster gets which TMdiff value within a pair (per
   method, independently), how often do you observe a concordance / diversity
   call as extreme as the real data?

It is purely a re-statistical computation on the existing per-cluster CSV
written by ``fold_diversity_survey.py`` — no AF/ESM/MSAT/DDG re-run.

Algorithm (per pair p)
----------------------
1. For each method m present for p, gather the per-cluster TMdiff_max vector
   over ShallowMsa clusters (DeepMsa is excluded, just as in the original
   centered-classification logic).
2. Apply ``_apply_centered_classification`` to compute observed
   ``TMdiff_centered`` and ``pref_centered``. This reproduces the centering
   used in ``fold_diversity_survey.py``.
3. Compute observed test statistics:
     * has_both_folds_centered (any F1 AND any F2 among ShallowMsa rows)
     * pairwise sign-agreement on TMdiff_centered for every method-pair
       (AF2~ESM, AF2~MSAT, AF2~CCMpred, AF2~DDG, ESM~MSAT, …)
     * mean_concordance = mean over present method-pairs of agree / (agree+disagree)
4. Repeat N_PERM times: independently shuffle each method's TMdiff_max vector
   over its clusters, re-derive ``TMdiff_centered`` (same per-pair median
   subtraction — the median is invariant under shuffling so centered values
   are equal to a permutation of the observed centered vector) and
   ``pref_centered``, then recompute the test statistics.
5. Empirical one-sided p-value per statistic =
       (#perms with statistic >= observed) / N_PERM.
   * For the boolean ``has_both_folds_centered`` we report
     P[shuffle == True] when observed is True, and 1.0 when observed is False
     (i.e. the permutation cannot make a False call "more extreme").

Output
------
``docs/fold_diversity_concordance_perm.csv`` with one row per pair:

    pair_id, n_clusters, n_perm,
    observed_mean_concordance, perm_mean_concordance_p,
    observed_<A>_<B>, perm_<A>_<B>_p, observed_<A>_<B>_n  [for each method-pair],
    observed_has_both_folds_centered, perm_has_both_folds_centered_p

Sanity check
------------
After running, the script prints the top-15 pairs by ascending
``perm_mean_concordance_p``. The fold-switching pairs ``2namA_1uxmK`` and
``2n54B_2hdmA`` are expected to rank near the top with p < 0.01.

Usage
-----
    python3 scripts/permutation_null_concordance.py
    python3 scripts/permutation_null_concordance.py --n-perm 1000 \
            --input  docs/fold_diversity_survey.csv \
            --output docs/fold_diversity_concordance_perm.csv \
            --seed 12345

Runs in seconds-to-minutes on 93 pairs * 1000 permutations.
"""

from __future__ import annotations

import argparse
import os
import sys
from itertools import combinations

import numpy as np
import pandas as pd

# Reuse the centering logic from fold_diversity_survey.py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fold_diversity_survey import _apply_centered_classification  # noqa: E402


# Method-pair ordering for column naming. We include every method that may
# appear in the survey CSV; pairs not actually present in a given pair_id
# are simply skipped.
ALL_METHODS = ["AF2", "ESM", "MSAT", "CCMpred", "DDG"]


# ---------------------------------------------------------------------------
# Statistic computation
# ---------------------------------------------------------------------------

def _classify_centered(tmdiff_centered: np.ndarray, delta: float) -> np.ndarray:
    """Vectorized analog of the ``pref_centered`` block in
    ``_apply_centered_classification``: F1 if > delta, F2 if < -delta, Amb
    otherwise. NaNs map to 'Amb'.
    """
    out = np.full(tmdiff_centered.shape, "Amb", dtype="<U3")
    finite = np.isfinite(tmdiff_centered)
    out[finite & (tmdiff_centered > delta)] = "F1"
    out[finite & (tmdiff_centered < -delta)] = "F2"
    return out


def _has_both_folds_centered(centered_by_method: dict) -> bool:
    """Pair-level diversity call: any method has at least one F1 AND any
    method (possibly the same) has at least one F2 — matching the
    summarize_pair logic which is per-method, but here we aggregate across
    methods only for a pair-level call.

    NOTE: ``summarize_pair`` actually computes ``has_both_folds_centered``
    PER METHOD. The pair-level concordance CSV does not currently store an
    aggregated flag. We expose two variants:

      * ``has_both_folds_centered_any``: True if ANY method shows both folds
      * ``has_both_folds_centered_all``: True if EVERY present method shows
        both folds

    The most natural interpretation of the manuscript "has_both_folds_centered"
    pair-level summary is the per-method maximum (any), which is the strictest
    possible "the pair shows fold diversity in at least one independent
    method". That's the variant we return as the primary statistic.
    """
    if not centered_by_method:
        return False
    for prefs in centered_by_method.values():
        has_f1 = bool(np.any(prefs == "F1"))
        has_f2 = bool(np.any(prefs == "F2"))
        if has_f1 and has_f2:
            return True
    return False


def _pairwise_concordance(
    centered_values_by_method: dict,
    cluster_ids_by_method: dict,
):
    """For every ordered method-pair (A, B) with both methods present,
    compute (agree_count, disagree_count, agree_pct) on per-cluster sign of
    TMdiff_centered. Skips (cluster) rows where either side is exactly 0 or
    NaN. Matches ``compute_pair_concordance`` in fold_diversity_survey.py.

    Returns dict {(A, B): (agree, disagree, n, pct)} for every method-pair
    with at least one shared cluster.
    """
    result = {}
    methods = [m for m in ALL_METHODS if m in centered_values_by_method]
    for A, B in combinations(methods, 2):
        clA = cluster_ids_by_method[A]
        clB = cluster_ids_by_method[B]
        valA = centered_values_by_method[A]
        valB = centered_values_by_method[B]
        # Map cluster -> value for each method
        mapA = dict(zip(clA, valA))
        mapB = dict(zip(clB, valB))
        shared = set(mapA) & set(mapB)
        agree = disagree = 0
        for cl in shared:
            a = mapA[cl]
            b = mapB[cl]
            if not np.isfinite(a) or not np.isfinite(b):
                continue
            if a == 0 or b == 0:
                continue
            if (a > 0) == (b > 0):
                agree += 1
            else:
                disagree += 1
        n = agree + disagree
        if n > 0:
            result[(A, B)] = (agree, disagree, n, agree / n)
    return result


def _aggregate_stats(
    centered_values_by_method: dict,
    cluster_ids_by_method: dict,
    delta: float,
):
    """Compute pair-level test statistics from a (possibly permuted) set of
    per-method centered TMdiff vectors aligned to per-method cluster lists.

    Returns dict with:
      * mean_concordance (float or NaN)
      * pairwise[(A, B)] = (agree, disagree, n, pct) for each present pair
      * has_both_folds_centered (bool, ANY-method variant)
    """
    pairwise = _pairwise_concordance(centered_values_by_method, cluster_ids_by_method)
    if pairwise:
        mean_pct = float(np.mean([v[3] for v in pairwise.values()]))
    else:
        mean_pct = float("nan")

    # has_both_folds_centered (per-method, any)
    classified = {
        m: _classify_centered(np.asarray(v, dtype=float), delta)
        for m, v in centered_values_by_method.items()
    }
    diversity = _has_both_folds_centered(classified)

    return {
        "mean_concordance": mean_pct,
        "pairwise": pairwise,
        "has_both_folds_centered": diversity,
    }


# ---------------------------------------------------------------------------
# Per-pair permutation driver
# ---------------------------------------------------------------------------

def _build_pair_inputs(detail_rows: pd.DataFrame, delta: float):
    """For one pair's per-cluster rows (DeepMsa included), apply the
    centered-classification logic to each method to derive TMdiff_centered
    over ShallowMsa clusters, then return per-method aligned lists usable
    for permutation:
        cluster_ids_by_method: {method: [cluster_id, ...]} ShallowMsa only
        centered_values_by_method: {method: np.array TMdiff_centered}
        tmdiff_max_by_method: {method: np.array TMdiff_max} ShallowMsa only

    The pair-median used for centering is always recomputed in the same way
    as in the original (median over ShallowMsa rows with finite TMdiff_max).
    """
    cluster_ids = {}
    centered_vals = {}
    raw_vals = {}
    medians = {}

    for method, mdf in detail_rows.groupby("method"):
        # Convert to dict-rows so the existing helper works unmodified.
        rows = []
        for _, r in mdf.iterrows():
            tdv = r["TMdiff_max"]
            if pd.isna(tdv):
                tdv = float("nan")
            rows.append({
                "pair_id": r["pair_id"],
                "cluster": r["cluster"],
                "method": method,
                "TMdiff_max": float(tdv) if not pd.isna(r["TMdiff_max"]) else float("nan"),
            })
        _apply_centered_classification(rows, delta)

        # Keep only ShallowMsa rows with finite centered values for permutation.
        cl, cv, rv = [], [], []
        for rr in rows:
            if rr["cluster"] == "DeepMsa":
                continue
            tc = rr["TMdiff_centered"]
            tm = rr["TMdiff_max"]
            if tc is None or (isinstance(tc, float) and not np.isfinite(tc)):
                continue
            cl.append(rr["cluster"])
            cv.append(float(tc))
            rv.append(float(tm))
        if cl:
            cluster_ids[method] = cl
            centered_vals[method] = np.asarray(cv, dtype=float)
            raw_vals[method] = np.asarray(rv, dtype=float)
            medians[method] = rows[0].get("pair_median_used", 0.0)

    return cluster_ids, centered_vals, raw_vals, medians


def permute_pair(
    pair_id: str,
    detail_rows: pd.DataFrame,
    delta: float,
    n_perm: int,
    seed: int,
):
    """Run permutation null for one pair.

    Returns a dict suitable for one row of the output CSV.
    """
    cluster_ids, centered_obs, _raw, _medians = _build_pair_inputs(detail_rows, delta)
    if len(centered_obs) < 2:
        # Not enough methods for any concordance; still report observed.
        observed = _aggregate_stats(centered_obs, cluster_ids, delta)
        return _format_row(
            pair_id, cluster_ids, observed,
            perm_mean_p=float("nan"),
            perm_pairwise_p={},
            perm_diversity_p=float("nan"),
            n_perm=n_perm,
        )

    observed = _aggregate_stats(centered_obs, cluster_ids, delta)

    # Permutation sampling — independent shuffle per method on its centered
    # vector. The cluster_ids list stays fixed; the values are permuted.
    rng = np.random.default_rng(seed)
    method_keys = list(centered_obs.keys())

    geq_mean = 0
    geq_diversity = 0
    geq_pairwise = {k: 0 for k in observed["pairwise"]}

    obs_mean = observed["mean_concordance"]
    obs_diversity = bool(observed["has_both_folds_centered"])
    obs_pairwise_pct = {k: v[3] for k, v in observed["pairwise"].items()}

    # Pre-extract centered vectors for fast shuffling.
    base_vectors = {m: centered_obs[m].copy() for m in method_keys}

    for _ in range(n_perm):
        perm_centered = {}
        for m in method_keys:
            v = base_vectors[m].copy()
            rng.shuffle(v)
            perm_centered[m] = v
        stats = _aggregate_stats(perm_centered, cluster_ids, delta)

        # mean_concordance — count permutations >= observed
        if not np.isnan(stats["mean_concordance"]) and not np.isnan(obs_mean):
            if stats["mean_concordance"] >= obs_mean:
                geq_mean += 1

        # has_both_folds_centered — count permutations matching observed=True
        if obs_diversity and stats["has_both_folds_centered"]:
            geq_diversity += 1

        # pairwise — count permutations with pct >= observed
        for k, obs_pct in obs_pairwise_pct.items():
            v = stats["pairwise"].get(k)
            if v is not None and v[3] >= obs_pct:
                geq_pairwise[k] += 1

    perm_mean_p = geq_mean / n_perm if not np.isnan(obs_mean) else float("nan")
    if obs_diversity:
        perm_diversity_p = geq_diversity / n_perm
    else:
        # Permutations cannot make a False call "more extreme"; report 1.0
        # to indicate non-significant.
        perm_diversity_p = 1.0
    perm_pairwise_p = {k: geq_pairwise[k] / n_perm for k in obs_pairwise_pct}

    return _format_row(
        pair_id, cluster_ids, observed,
        perm_mean_p=perm_mean_p,
        perm_pairwise_p=perm_pairwise_p,
        perm_diversity_p=perm_diversity_p,
        n_perm=n_perm,
    )


def _format_row(
    pair_id: str,
    cluster_ids_by_method: dict,
    observed: dict,
    perm_mean_p: float,
    perm_pairwise_p: dict,
    perm_diversity_p: float,
    n_perm: int,
):
    n_clusters = max((len(v) for v in cluster_ids_by_method.values()), default=0)
    row = {
        "pair_id": pair_id,
        "n_clusters": n_clusters,
        "n_methods_present": len(cluster_ids_by_method),
        "n_perm": n_perm,
        "observed_mean_concordance": (
            None if np.isnan(observed["mean_concordance"])
            else round(observed["mean_concordance"], 4)
        ),
        "perm_mean_concordance_p": (
            None if (isinstance(perm_mean_p, float) and np.isnan(perm_mean_p))
            else round(perm_mean_p, 4)
        ),
        "observed_has_both_folds_centered": bool(observed["has_both_folds_centered"]),
        "perm_has_both_folds_centered_p": (
            None if (isinstance(perm_diversity_p, float) and np.isnan(perm_diversity_p))
            else round(perm_diversity_p, 4)
        ),
    }
    # Pairwise stats — fixed column ordering for stable schema across pairs.
    for A, B in combinations(ALL_METHODS, 2):
        key = (A, B)
        col_pct = f"observed_{A}_{B}"
        col_n = f"observed_{A}_{B}_n"
        col_p = f"perm_{A}_{B}_p"
        if key in observed["pairwise"]:
            agree, disagree, n, pct = observed["pairwise"][key]
            row[col_pct] = round(pct, 4)
            row[col_n] = n
            row[col_p] = round(perm_pairwise_p.get(key, float("nan")), 4)
        else:
            row[col_pct] = None
            row[col_n] = 0
            row[col_p] = None
    return row


# ---------------------------------------------------------------------------
# Corrected-mode (TMdiff_residual) logic
# Merged from former scripts/permutation_null_corrected.py (2026-05-23).
# Used when --mode corrected: reads TMdiff_residual column from
# fold_diversity_survey_corrected.csv (produced by seq_divergence_correction.py)
# and computes empirical permutation null + BH FDR on the per-pair mean
# pairwise sign-concordance of residuals.
# ---------------------------------------------------------------------------

def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN-safe, monotone)."""
    p = np.asarray(pvals, dtype=float)
    mask = np.isfinite(p)
    out = np.full_like(p, np.nan, dtype=float)
    if not mask.any():
        return out
    pv = p[mask]
    n = len(pv)
    order = np.argsort(pv)
    ranks = np.arange(1, n + 1)
    adj_sorted = pv[order] * n / ranks
    adj_sorted = np.minimum.accumulate(adj_sorted[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(adj_sorted, 0.0, 1.0)
    out[mask] = adj
    return out


def _mean_concordance_from_residuals(per_method: dict) -> tuple:
    """Mean pairwise sign-agreement across method-pairs on the residual scale.

    Returns (mean_A, mean_n_overlap_clusters).
    """
    present = [m for m in ALL_METHODS if m in per_method and len(per_method[m]) > 0]
    if len(present) < 2:
        return float("nan"), 0
    agreements: list = []
    ns: list = []
    for a, b in combinations(present, 2):
        va, vb = per_method[a], per_method[b]
        n = min(len(va), len(vb))
        if n < 1:
            continue
        sa = np.sign(va[:n])
        sb = np.sign(vb[:n])
        mask = (sa != 0) & (sb != 0)
        if mask.sum() < 1:
            continue
        agreements.append((sa[mask] == sb[mask]).mean())
        ns.append(int(mask.sum()))
    if not agreements:
        return float("nan"), 0
    return float(np.mean(agreements)), int(np.mean(ns))


def _run_corrected_mode(args) -> None:
    """Corrected-mode path: read TMdiff_residual, permute within method,
    write per-pair p-values + BH-adjusted variant. Mirrors former
    scripts/permutation_null_corrected.py behavior.
    """
    if not os.path.isfile(args.input):
        print(f"[error] corrected survey CSV not found: {args.input}", file=sys.stderr)
        print("        Run scripts/seq_divergence_correction.py first.", file=sys.stderr)
        sys.exit(2)

    df = pd.read_csv(args.input)
    needed = {"pair_id", "method", "TMdiff_residual"}
    if not needed.issubset(df.columns):
        raise SystemExit(
            f"missing columns in {args.input}: need {needed - set(df.columns)}"
        )
    cluster_col = "cluster_norm" if "cluster_norm" in df.columns else "cluster"
    df = df[df[cluster_col] != "DeepMsa"].copy()

    rng = np.random.default_rng(args.seed)
    rows = []
    for pid, sub in df.groupby("pair_id"):
        per_method_obs: dict = {}
        for m in ALL_METHODS:
            ms = sub[sub["method"] == m]
            vals = ms["TMdiff_residual"].dropna().to_numpy(dtype=float)
            if len(vals):
                per_method_obs[m] = vals
        if len(per_method_obs) < 2:
            continue
        obs_mean, obs_n = _mean_concordance_from_residuals(per_method_obs)
        if np.isnan(obs_mean):
            continue
        n_clusters = max(len(v) for v in per_method_obs.values())

        n_ge = 0
        for _ in range(args.n_perm):
            permuted = {m: rng.permutation(v) for m, v in per_method_obs.items()}
            perm_mean, _ = _mean_concordance_from_residuals(permuted)
            if not np.isnan(perm_mean) and perm_mean >= obs_mean:
                n_ge += 1
        p_emp = (n_ge + 1) / (args.n_perm + 1)

        rows.append({
            "pair_id": pid,
            "n_clusters": int(n_clusters),
            "n_methods_present": len(per_method_obs),
            "observed_mean_concordance_corrected": obs_mean,
            "observed_mean_n_clusters": obs_n,
            "n_perm": args.n_perm,
            "perm_mean_concordance_p": p_emp,
        })

    out = pd.DataFrame(rows).sort_values("perm_mean_concordance_p")
    out["p_bh"] = bh_adjust(out["perm_mean_concordance_p"].to_numpy())
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    out.to_csv(args.output, index=False)
    # also write the BH variant (matches old permutation_null_corrected.py)
    bh_out = args.output.replace(".csv", "_bh.csv")
    if bh_out == args.output:
        bh_out = args.output + ".bh.csv"
    out.to_csv(bh_out, index=False)
    print(f"wrote {args.output} and {bh_out}: {len(out)} pairs, "
          f"N_PERM={args.n_perm}")

    passing = out[out["p_bh"] <= args.alpha]
    print(f"\n=== Pairs passing post-correction BH FDR <= {args.alpha} ===")
    for _, r in passing.iterrows():
        print(f"  {r['pair_id']:<18}  "
              f"obs={r['observed_mean_concordance_corrected']:.3f}  "
              f"raw_p={r['perm_mean_concordance_p']:.4f}  "
              f"BH_p={r['p_bh']:.4f}  n_clu={int(r['n_clusters'])}")
    if not len(passing):
        print("  (none)")


def _submit_self_sbatch(args, script_path: str) -> None:
    """Re-invoke this script under sbatch with --run-job-mode local.

    Mirrors the run_foldswitch_pipeline.py pattern: --run-job-mode sbatch
    submits a single sbatch job that runs the full analysis inline.
    """
    import shlex
    import subprocess
    jobs_dir = "Pipeline/FoldPairs/jobs"
    os.makedirs(jobs_dir, exist_ok=True)
    log = os.path.join(jobs_dir, f"permutation_null_{args.mode}.out")
    inner_args = []
    for k, v in vars(args).items():
        if k == "run_job_mode":
            continue
        if v is None or v is False:
            continue
        flag = "--" + k.replace("_", "-")
        if v is True:
            inner_args.append(flag)
        else:
            inner_args.extend([flag, str(v)])
    inner = " ".join([shlex.quote(sys.executable), script_path,
                      "--run-job-mode", "local"] + inner_args)
    sbatch_cmd = [
        "sbatch", "--parsable",
        f"--job-name=perm_null_{args.mode}",
        f"--cpus-per-task={args.sbatch_cpus}",
        f"--mem={args.sbatch_mem}",
        f"--time={args.sbatch_time}",
        f"--output={log}",
        "--wrap", shlex.quote(inner),
    ]
    if args.sbatch_partition:
        sbatch_cmd.insert(1, f"--partition={args.sbatch_partition}")
    cmd_str = " ".join(sbatch_cmd)
    print(f"[sbatch] {cmd_str}")
    jid = subprocess.check_output(cmd_str, shell=True, text=True).strip()
    print(f"[sbatch] submitted permutation_null ({args.mode}) job {jid}")
    print(f"[sbatch] log: {log}")


# ---------------------------------------------------------------------------
# Main CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Permutation-null test for fold_diversity_survey concordance. "
                    "--mode raw uses TMdiff_centered (pre-correction); "
                    "--mode corrected uses TMdiff_residual (post-correction)."
    )
    parser.add_argument(
        "--mode", choices=["raw", "corrected"], default="raw",
        help="raw = pre-correction TMdiff_centered + centered classification (default); "
             "corrected = post-correction TMdiff_residual from "
             "fold_diversity_survey_corrected.csv (+ BH FDR adjustment). "
             "Replaces the former scripts/permutation_null_corrected.py."
    )
    parser.add_argument(
        "--input", type=str, default=None,
        help="Per-cluster detail CSV. Default depends on --mode: "
             "raw → docs/fold_diversity_survey.csv; "
             "corrected → docs/fold_diversity_survey_corrected.csv"
    )
    parser.add_argument(
        "--output", type=str, default=None,
        help="Output CSV path. Default depends on --mode: "
             "raw → docs/fold_diversity_concordance_perm.csv; "
             "corrected → docs/fold_diversity_concordance_corrected_perm.csv"
    )
    parser.add_argument(
        "--n-perm", type=int, default=1000,
        help="Number of permutations per pair (default: 1000)"
    )
    parser.add_argument(
        "--delta", type=float, default=0.05,
        help="TM threshold for F1/F2/Amb classification (raw mode only; "
             "default: 0.05; matches fold_diversity_survey.py)"
    )
    parser.add_argument(
        "--seed", type=int, default=12345,
        help="Base random seed (default: 12345)"
    )
    parser.add_argument(
        "--top", type=int, default=15,
        help="How many pairs to print in top-by-significance table (default: 15)"
    )
    parser.add_argument(
        "--alpha", type=float, default=0.05,
        help="BH FDR threshold for highlighting (corrected mode; default 0.05)"
    )
    # sbatch self-submission (matches run_foldswitch_pipeline.py pattern)
    parser.add_argument(
        "--run-job-mode", choices=["local", "sbatch"], default="local",
        help="local = run inline (default); sbatch = self-submit one sbatch job"
    )
    parser.add_argument("--sbatch-partition", default=None)
    parser.add_argument("--sbatch-cpus", default=2)
    parser.add_argument("--sbatch-mem", default="8G")
    parser.add_argument("--sbatch-time", default="04:00:00")
    args = parser.parse_args()

    # sbatch self-submission
    if args.run_job_mode == "sbatch":
        _submit_self_sbatch(args, os.path.abspath(__file__))
        return

    # Resolve defaults from --mode
    if args.input is None:
        args.input = (os.path.join("docs", "fold_diversity_survey_corrected.csv")
                      if args.mode == "corrected"
                      else os.path.join("docs", "fold_diversity_survey.csv"))
    if args.output is None:
        args.output = (os.path.join("docs", "fold_diversity_concordance_corrected_perm.csv")
                       if args.mode == "corrected"
                       else os.path.join("docs", "fold_diversity_concordance_perm.csv"))

    # Dispatch by mode
    if args.mode == "corrected":
        _run_corrected_mode(args)
        return

    if not os.path.isfile(args.input):
        print(f"[error] input CSV not found: {args.input}", file=sys.stderr)
        sys.exit(2)

    print(f"Reading per-cluster detail: {args.input}")
    df = pd.read_csv(args.input)
    print(f"  {len(df)} rows; "
          f"{df['pair_id'].nunique()} pairs; "
          f"{df['method'].nunique()} methods "
          f"({sorted(df['method'].unique())})")

    pair_ids = sorted(df["pair_id"].unique())
    print(f"Running {args.n_perm} permutations per pair on {len(pair_ids)} pairs "
          f"(seed={args.seed}, delta={args.delta})")

    out_rows = []
    for i, pid in enumerate(pair_ids, 1):
        sub = df[df["pair_id"] == pid]
        seed = args.seed + abs(hash(pid)) % 10000
        row = permute_pair(pid, sub, args.delta, args.n_perm, seed)
        out_rows.append(row)
        if i % 10 == 0 or i == len(pair_ids):
            print(f"  [{i:>3d}/{len(pair_ids)}] done")

    out_df = pd.DataFrame(out_rows)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    out_df.to_csv(args.output, index=False)
    print(f"\nWrote {args.output} ({len(out_df)} rows)")

    # --- Sanity-check report ---
    print(f"\n{'='*78}")
    print(f"Top-{args.top} pairs by ascending perm_mean_concordance_p (most significant)")
    print(f"{'='*78}")
    sig = out_df.dropna(subset=["perm_mean_concordance_p"]).copy()
    sig = sig.sort_values(["perm_mean_concordance_p", "observed_mean_concordance"],
                          ascending=[True, False])
    if sig.empty:
        print("  (no pairs with computable mean_concordance)")
    else:
        print(f"  {'pair':<22s} {'n_cl':>4s} {'mean':>5s} {'p_mean':>7s} "
              f"{'div?':>4s} {'p_div':>7s}")
        for _, r in sig.head(args.top).iterrows():
            div_obs = "YES" if r["observed_has_both_folds_centered"] else "no"
            p_div = r.get("perm_has_both_folds_centered_p")
            p_div_s = "    --" if p_div is None or pd.isna(p_div) else f"{p_div:>7.4f}"
            mean_s = "  --" if pd.isna(r["observed_mean_concordance"]) else f"{r['observed_mean_concordance']:>5.2f}"
            p_mean_s = "    --" if pd.isna(r["perm_mean_concordance_p"]) else f"{r['perm_mean_concordance_p']:>7.4f}"
            print(f"  {r['pair_id']:<22s} {int(r['n_clusters']):>4d} {mean_s} "
                  f"{p_mean_s} {div_obs:>4s} {p_div_s}")

    # Highlight expected fold-switching pairs.
    EXPECTED = ["2namA_1uxmK", "2n54B_2hdmA"]
    found = sig[sig["pair_id"].isin(EXPECTED)]
    if not found.empty:
        print(f"\nExpected fold-switching pairs (sanity check):")
        for _, r in found.iterrows():
            p_mean = r["perm_mean_concordance_p"]
            p_mean_s = "--" if pd.isna(p_mean) else f"{p_mean:.4f}"
            mean_s = "--" if pd.isna(r["observed_mean_concordance"]) else f"{r['observed_mean_concordance']:.2f}"
            print(f"  {r['pair_id']:<22s} mean={mean_s}  p_mean={p_mean_s}")
    print()


if __name__ == "__main__":
    main()
