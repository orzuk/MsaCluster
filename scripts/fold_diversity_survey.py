#!/usr/bin/env python3
"""Survey fold-preference diversity across clusters for all pairs.

For each pair, reads existing df_af.csv, df_esm.csv, df_cmap.csv and computes
per-cluster fold-preference metrics from three independent methods:
  - AF2:  TM-score toward fold1 vs fold2
  - ESM:  TM-score toward fold1 vs fold2 (10 samples per cluster)
  - MSAT: contact map recall for fold1 vs fold2 full contact maps

Outputs:
  docs/fold_diversity_survey.csv   — one row per (pair, cluster, method)
  docs/fold_diversity_summary.csv  — one row per pair, aggregated stats
  stdout summary of which method shows the most fold differentiation

Usage:
  python scripts/fold_diversity_survey.py
  python scripts/fold_diversity_survey.py --pairs 1dzlA_5keqF,2qqjA_4qdsA
  python scripts/fold_diversity_survey.py --delta 0.05
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config import DATA_DIR, PAIR_DIR_RE, TABLES_RES


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_MSA_TAG_RE = re.compile(r"(?:Shallow)?Msa_(\d+)", re.IGNORECASE)


def _normalize_cluster(s):
    """Normalize cluster tag -> ShallowMsa_NNN or DeepMsa."""
    s = str(s).strip()
    if s.lower() in ("deep", "deepmsa", "msa_deep", "deep_msa", "query"):
        return "DeepMsa"
    m = _MSA_TAG_RE.fullmatch(s)
    if m:
        return f"ShallowMsa_{int(m.group(1)):03d}"
    if s.isdigit():
        return f"ShallowMsa_{int(s):03d}"
    return s


def _assign_pref(tm1, tm2, delta):
    """Assign fold preference from TM-score pair."""
    if np.isnan(tm1) or np.isnan(tm2):
        return "Amb"
    if tm1 - tm2 > delta:
        return "F1"
    if tm2 - tm1 > delta:
        return "F2"
    return "Amb"


# ---------------------------------------------------------------------------
# AF2 loader
# ---------------------------------------------------------------------------

def load_af2_diversity(pair_id, delta):
    """Load df_af.csv and compute per-cluster fold preference.

    Handles both old format (score_pdb1/score_pdb2) and new format
    (TMscore_fold1/TMscore_fold2).

    Returns list of dicts, one per cluster.
    """
    csv_path = os.path.join(DATA_DIR, pair_id, "Analysis", "df_af.csv")
    if not os.path.isfile(csv_path):
        return []

    df = pd.read_csv(csv_path)
    if df.empty:
        return []

    # Normalize column names
    if "TMscore_fold1" not in df.columns and "score_pdb1" in df.columns:
        df = df.rename(columns={"score_pdb1": "TMscore_fold1",
                                "score_pdb2": "TMscore_fold2"})

    if "TMscore_fold1" not in df.columns:
        return []

    # Filter AF2 only (new format has 'model' column)
    if "model" in df.columns:
        df = df[df["model"].astype(str).str.upper() == "AF2"]

    # Parse cluster tags
    col = "cluster_num" if "cluster_num" in df.columns else "cluster"
    if col not in df.columns:
        return []
    df["_tag"] = df[col].apply(_normalize_cluster)

    df["TM1"] = pd.to_numeric(df["TMscore_fold1"], errors="coerce")
    df["TM2"] = pd.to_numeric(df["TMscore_fold2"], errors="coerce")

    rows = []
    for tag, grp in df.groupby("_tag"):
        tm1_max = grp["TM1"].max()
        tm2_max = grp["TM2"].max()
        tm1_mean = grp["TM1"].mean()
        tm2_mean = grp["TM2"].mean()
        n_models = len(grp)
        # Per-model fold assignment (argmax)
        n_f1 = (grp["TM1"] > grp["TM2"]).sum()
        n_f2 = (grp["TM2"] > grp["TM1"]).sum()
        pref = _assign_pref(tm1_max, tm2_max, delta)
        rows.append({
            "pair_id": pair_id,
            "cluster": tag,
            "method": "AF2",
            "TM1_max": round(tm1_max, 4),
            "TM2_max": round(tm2_max, 4),
            "TM1_mean": round(tm1_mean, 4),
            "TM2_mean": round(tm2_mean, 4),
            "TMdiff_max": round(tm1_max - tm2_max, 4),
            "TMdiff_mean": round(tm1_mean - tm2_mean, 4),
            "n_models": n_models,
            "n_toward_f1": int(n_f1),
            "n_toward_f2": int(n_f2),
            "vote_frac_f1": round(n_f1 / n_models, 3) if n_models > 0 else 0,
            "pref": pref,
        })
    return rows


# ---------------------------------------------------------------------------
# ESM loader
# ---------------------------------------------------------------------------

def load_esm_diversity(pair_id, delta):
    """Load df_esm.csv and compute per-cluster fold preference.

    ESM has ~10 samples per cluster, so we use VOTE-BASED preference:
    the fold preferred by the majority of individual samples determines
    the cluster's preference. This avoids the max-aggregation problem
    where a single outlier sample dominates.
    """
    csv_path = os.path.join(DATA_DIR, pair_id, "Analysis", "df_esm.csv")
    if not os.path.isfile(csv_path):
        return []

    df = pd.read_csv(csv_path)
    if df.empty:
        return []

    # Normalize column names
    if "TMscore_fold1" not in df.columns and "score_pdb1" in df.columns:
        df = df.rename(columns={"score_pdb1": "TMscore_fold1",
                                "score_pdb2": "TMscore_fold2"})

    if "TMscore_fold1" not in df.columns:
        return []

    col = "cluster_num" if "cluster_num" in df.columns else "cluster"
    if col not in df.columns:
        return []
    df["_tag"] = df[col].apply(_normalize_cluster)

    df["TM1"] = pd.to_numeric(df["TMscore_fold1"], errors="coerce")
    df["TM2"] = pd.to_numeric(df["TMscore_fold2"], errors="coerce")

    rows = []
    for tag, grp in df.groupby("_tag"):
        tm1_max = grp["TM1"].max()
        tm2_max = grp["TM2"].max()
        tm1_mean = grp["TM1"].mean()
        tm2_mean = grp["TM2"].mean()
        n_models = len(grp)
        n_f1 = int((grp["TM1"] > grp["TM2"]).sum())
        n_f2 = int((grp["TM2"] > grp["TM1"]).sum())
        vote_frac = n_f1 / n_models if n_models > 0 else 0.5

        # Vote-based preference: use mean TMdiff, not max
        mean_diff = tm1_mean - tm2_mean
        if mean_diff > delta:
            pref = "F1"
        elif mean_diff < -delta:
            pref = "F2"
        else:
            pref = "Amb"

        rows.append({
            "pair_id": pair_id,
            "cluster": tag,
            "method": "ESM",
            "TM1_max": round(tm1_max, 4),
            "TM2_max": round(tm2_max, 4),
            "TM1_mean": round(tm1_mean, 4),
            "TM2_mean": round(tm2_mean, 4),
            "TMdiff_max": round(tm1_max - tm2_max, 4),
            "TMdiff_mean": round(mean_diff, 4),
            "n_models": n_models,
            "n_toward_f1": n_f1,
            "n_toward_f2": n_f2,
            "vote_frac_f1": round(vote_frac, 3),
            "pref": pref,
        })
    return rows


# ---------------------------------------------------------------------------
# MSAT (contact map) loader
# ---------------------------------------------------------------------------

def load_msat_diversity(pair_id, delta, msat_metric="t_recall"):
    """Load df_cmap.csv and compute per-cluster fold preference from contacts.

    msat_metric controls which contact recall columns to use:
      "t_recall"    — full-fold recall (t1_recall vs t2_recall). Uses all
                      contacts including shared ones. More robust but less
                      discriminating (shared contacts compress the range).
      "uniq_recall" — fold-unique recall (uniq1_recall vs uniq2_recall).
                      Only contacts unique to each fold. More discriminating
                      but sparser (may be NaN for some clusters).
                      Requires the cmap_analysis.py Layer 2 fix (2026-02-24).
    """
    csv_path = os.path.join(DATA_DIR, pair_id, "Analysis", "df_cmap.csv")
    if not os.path.isfile(csv_path):
        return []

    df = pd.read_csv(csv_path)
    if df.empty:
        return []

    if msat_metric == "uniq_recall":
        col1, col2 = "uniq1_recall", "uniq2_recall"
        metric_label = "uniq_recall"
    else:
        col1, col2 = "t1_recall", "t2_recall"
        metric_label = "t_recall"

    if col1 not in df.columns or col2 not in df.columns:
        return []

    # Parse cluster tags
    if "cluster" in df.columns:
        cluster_col = "cluster"
    elif "file" in df.columns:
        # Extract cluster from filename: msa_t__ShallowMsa_000.npy
        df["cluster"] = df["file"].apply(
            lambda f: re.search(r"((?:Shallow|Deep)Msa[_\d]*)", str(f)).group(1)
            if re.search(r"((?:Shallow|Deep)Msa[_\d]*)", str(f)) else str(f)
        )
        cluster_col = "cluster"
    else:
        return []

    df["_tag"] = df[cluster_col].apply(_normalize_cluster)
    df["R1"] = pd.to_numeric(df[col1], errors="coerce")
    df["R2"] = pd.to_numeric(df[col2], errors="coerce")

    rows = []
    for tag, grp in df.groupby("_tag"):
        # Drop rows where recall is NaN (missing data, not zero signal)
        valid = grp.dropna(subset=["R1", "R2"])
        if valid.empty:
            continue

        r1 = valid["R1"].max()
        r2 = valid["R2"].max()
        r1_mean = valid["R1"].mean()
        r2_mean = valid["R2"].mean()

        # Fold preference from contact recall difference
        diff = r1 - r2
        if abs(diff) > delta:
            pref = "F1" if diff > 0 else "F2"
        else:
            pref = "Amb"

        rows.append({
            "pair_id": pair_id,
            "cluster": tag,
            "method": "MSAT",
            "TM1_max": round(r1, 4),  # reusing TM column names for uniformity
            "TM2_max": round(r2, 4),
            "TM1_mean": round(r1_mean, 4),
            "TM2_mean": round(r2_mean, 4),
            "TMdiff_max": round(r1 - r2, 4),
            "TMdiff_mean": round(r1_mean - r2_mean, 4),
            "n_models": len(valid),
            "n_toward_f1": int((valid["R1"] > valid["R2"]).sum()),
            "n_toward_f2": int((valid["R2"] > valid["R1"]).sum()),
            "vote_frac_f1": round(
                (valid["R1"] > valid["R2"]).sum() / len(valid), 3
            ) if len(valid) > 0 else 0,
            "pref": pref,
            "t1_recall": round(r1, 4) if not np.isnan(r1) else np.nan,
            "t2_recall": round(r2, 4) if not np.isnan(r2) else np.nan,
            "metric_used": metric_label,
        })
    return rows


# ---------------------------------------------------------------------------
# CCMpred (contact map) loader
# ---------------------------------------------------------------------------

def load_ccmpred_diversity(pair_id, delta, msat_metric="t_recall"):
    """Load df_cmap_ccmpred.csv and compute per-cluster fold preference.

    Uses the same metric logic as load_msat_diversity but reads from
    df_cmap_ccmpred.csv (CCMpred contact predictions) instead of
    df_cmap.csv (MSA-Transformer predictions).
    """
    csv_path = os.path.join(DATA_DIR, pair_id, "Analysis", "df_cmap_ccmpred.csv")
    if not os.path.isfile(csv_path):
        return []

    df = pd.read_csv(csv_path)
    if df.empty:
        return []

    if msat_metric == "uniq_recall":
        col1, col2 = "uniq1_recall", "uniq2_recall"
        metric_label = "uniq_recall"
    else:
        col1, col2 = "t1_recall", "t2_recall"
        metric_label = "t_recall"

    if col1 not in df.columns or col2 not in df.columns:
        return []

    # Parse cluster tags
    if "cluster" in df.columns:
        cluster_col = "cluster"
    elif "file" in df.columns:
        df["cluster"] = df["file"].apply(
            lambda f: re.search(r"((?:Shallow|Deep)Msa[_\d]*)", str(f)).group(1)
            if re.search(r"((?:Shallow|Deep)Msa[_\d]*)", str(f)) else str(f)
        )
        cluster_col = "cluster"
    else:
        return []

    df["_tag"] = df[cluster_col].apply(_normalize_cluster)
    df["R1"] = pd.to_numeric(df[col1], errors="coerce")
    df["R2"] = pd.to_numeric(df[col2], errors="coerce")

    rows = []
    for tag, grp in df.groupby("_tag"):
        valid = grp.dropna(subset=["R1", "R2"])
        if valid.empty:
            continue

        r1 = valid["R1"].max()
        r2 = valid["R2"].max()
        r1_mean = valid["R1"].mean()
        r2_mean = valid["R2"].mean()

        diff = r1 - r2
        if abs(diff) > delta:
            pref = "F1" if diff > 0 else "F2"
        else:
            pref = "Amb"

        rows.append({
            "pair_id": pair_id,
            "cluster": tag,
            "method": "CCMpred",
            "TM1_max": round(r1, 4),
            "TM2_max": round(r2, 4),
            "TM1_mean": round(r1_mean, 4),
            "TM2_mean": round(r2_mean, 4),
            "TMdiff_max": round(r1 - r2, 4),
            "TMdiff_mean": round(r1_mean - r2_mean, 4),
            "n_models": len(valid),
            "n_toward_f1": int((valid["R1"] > valid["R2"]).sum()),
            "n_toward_f2": int((valid["R2"] > valid["R1"]).sum()),
            "vote_frac_f1": round(
                (valid["R1"] > valid["R2"]).sum() / len(valid), 3
            ) if len(valid) > 0 else 0,
            "pref": pref,
            "t1_recall": round(r1, 4) if not np.isnan(r1) else np.nan,
            "t2_recall": round(r2, 4) if not np.isnan(r2) else np.nan,
            "metric_used": metric_label,
        })
    return rows


# ---------------------------------------------------------------------------
# Per-pair summary
# ---------------------------------------------------------------------------

def summarize_pair(pair_id, cluster_rows):
    """Compute pair-level diversity metrics from cluster-level rows."""
    if not cluster_rows:
        return None

    df = pd.DataFrame(cluster_rows)

    # Exclude DeepMsa from diversity computation (it's the full MSA)
    shallow = df[df["cluster"] != "DeepMsa"]
    if shallow.empty:
        shallow = df  # fallback

    method = cluster_rows[0]["method"]
    n_clusters = len(shallow)
    prefs = shallow["pref"].value_counts().to_dict()

    # Spread: std of TMdiff across clusters
    tm_diffs = shallow["TMdiff_max"].values
    spread = float(np.std(tm_diffs)) if len(tm_diffs) > 1 else 0.0

    # Range of TMdiff
    diff_range = float(np.max(tm_diffs) - np.min(tm_diffs)) if len(tm_diffs) > 1 else 0.0

    # Does ANY cluster prefer a different fold from the majority?
    has_f1 = prefs.get("F1", 0) > 0
    has_f2 = prefs.get("F2", 0) > 0
    has_diversity = has_f1 and has_f2

    # Best TM1 and TM2 across clusters
    best_tm1 = shallow["TM1_max"].max()
    best_tm2 = shallow["TM2_max"].max()

    return {
        "pair_id": pair_id,
        "method": method,
        "n_clusters": n_clusters,
        "n_f1": prefs.get("F1", 0),
        "n_f2": prefs.get("F2", 0),
        "n_amb": prefs.get("Amb", 0),
        "has_both_folds": has_diversity,
        "TMdiff_spread": round(spread, 4),
        "TMdiff_range": round(diff_range, 4),
        "TMdiff_mean": round(float(np.mean(tm_diffs)), 4),
        "best_score_f1": round(best_tm1, 4),
        "best_score_f2": round(best_tm2, 4),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def discover_pairs():
    """Find all pair directories."""
    pairs = []
    if not os.path.isdir(DATA_DIR):
        return pairs
    for name in sorted(os.listdir(DATA_DIR)):
        if PAIR_DIR_RE.match(name):
            pairs.append(name)
    return pairs


def main():
    parser = argparse.ArgumentParser(description="Fold-diversity survey across methods")
    parser.add_argument("--pairs", type=str, default="ALL",
                        help="Comma-separated pair IDs, or ALL")
    parser.add_argument("--delta", type=float, default=0.05,
                        help="Threshold for F1/F2/Amb assignment")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Output directory (default: docs/)")
    parser.add_argument("--msat-metric", type=str, default="t_recall",
                        choices=["t_recall", "uniq_recall"],
                        help="MSAT metric: t_recall (full-fold, default) or "
                             "uniq_recall (fold-unique contacts)")
    args = parser.parse_args()

    if args.pairs.upper() == "ALL":
        pairs = discover_pairs()
    else:
        pairs = [p.strip() for p in args.pairs.split(",") if p.strip()]

    out_dir = args.output_dir or TABLES_RES
    os.makedirs(out_dir, exist_ok=True)

    print(f"Scanning {len(pairs)} pairs for fold diversity...")
    print(f"Delta threshold: {args.delta}")
    print(f"MSAT metric:     {args.msat_metric}")
    print()

    all_cluster_rows = []
    all_summaries = []

    n_af2 = n_esm = n_msat = n_ccmpred = 0

    for pair_id in pairs:
        # AF2
        af_rows = load_af2_diversity(pair_id, args.delta)
        if af_rows:
            n_af2 += 1
            all_cluster_rows.extend(af_rows)
            s = summarize_pair(pair_id, af_rows)
            if s:
                all_summaries.append(s)

        # ESM
        esm_rows = load_esm_diversity(pair_id, args.delta)
        if esm_rows:
            n_esm += 1
            all_cluster_rows.extend(esm_rows)
            s = summarize_pair(pair_id, esm_rows)
            if s:
                all_summaries.append(s)

        # MSAT
        msat_rows = load_msat_diversity(pair_id, args.delta,
                                        msat_metric=args.msat_metric)
        if msat_rows:
            n_msat += 1
            all_cluster_rows.extend(msat_rows)
            s = summarize_pair(pair_id, msat_rows)
            if s:
                all_summaries.append(s)

        # CCMpred
        cc_rows = load_ccmpred_diversity(pair_id, args.delta,
                                         msat_metric=args.msat_metric)
        if cc_rows:
            n_ccmpred += 1
            all_cluster_rows.extend(cc_rows)
            s = summarize_pair(pair_id, cc_rows)
            if s:
                all_summaries.append(s)

    # --- Save detailed CSV ---
    df_detail = pd.DataFrame(all_cluster_rows)
    detail_path = os.path.join(out_dir, "fold_diversity_survey.csv")
    if not df_detail.empty:
        df_detail.to_csv(detail_path, index=False)
        print(f"Saved per-cluster detail: {detail_path} ({len(df_detail)} rows)")
    else:
        print("[warn] No data found.")
        return

    # --- Save summary CSV ---
    df_summary = pd.DataFrame(all_summaries)
    summary_path = os.path.join(out_dir, "fold_diversity_summary.csv")
    df_summary.to_csv(summary_path, index=False)
    print(f"Saved pair-level summary: {summary_path} ({len(df_summary)} rows)")

    # --- Print analysis ---
    print(f"\n{'='*70}")
    print(f"FOLD DIVERSITY SURVEY RESULTS")
    print(f"{'='*70}")
    print(f"Pairs with AF2 data:     {n_af2}")
    print(f"Pairs with ESM data:     {n_esm}")
    print(f"Pairs with MSAT data:    {n_msat}")
    print(f"Pairs with CCMpred data: {n_ccmpred}")

    for method in ["AF2", "ESM", "MSAT", "CCMpred"]:
        mdf = df_summary[df_summary["method"] == method]
        if mdf.empty:
            continue

        print(f"\n--- {method} ({len(mdf)} pairs) ---")

        # How many pairs show diversity?
        n_diverse = mdf["has_both_folds"].sum()
        print(f"  Pairs with BOTH F1 and F2 clusters: {n_diverse}/{len(mdf)} "
              f"({100*n_diverse/len(mdf):.0f}%)")

        # Preference distribution
        total_f1 = mdf["n_f1"].sum()
        total_f2 = mdf["n_f2"].sum()
        total_amb = mdf["n_amb"].sum()
        total = total_f1 + total_f2 + total_amb
        print(f"  Total clusters: F1={total_f1} F2={total_f2} Amb={total_amb} "
              f"({total} total)")
        if total > 0:
            print(f"  Fractions:      F1={100*total_f1/total:.1f}% "
                  f"F2={100*total_f2/total:.1f}% "
                  f"Amb={100*total_amb/total:.1f}%")

        # Spread statistics
        print(f"  TMdiff spread:  mean={mdf['TMdiff_spread'].mean():.4f}  "
              f"median={mdf['TMdiff_spread'].median():.4f}  "
              f"max={mdf['TMdiff_spread'].max():.4f}")
        print(f"  TMdiff range:   mean={mdf['TMdiff_range'].mean():.4f}  "
              f"median={mdf['TMdiff_range'].median():.4f}  "
              f"max={mdf['TMdiff_range'].max():.4f}")

        # Per-pair breakdown
        print(f"\n  Per-pair breakdown:")
        print(f"  {'pair':<20s} {'clusters':>8s} {'F1':>4s} {'F2':>4s} "
              f"{'Amb':>4s} {'diverse':>8s} {'spread':>8s} {'range':>8s}")
        for _, row in mdf.iterrows():
            div_mark = " ***" if row["has_both_folds"] else ""
            print(f"  {row['pair_id']:<20s} {row['n_clusters']:>8d} "
                  f"{row['n_f1']:>4d} {row['n_f2']:>4d} {row['n_amb']:>4d} "
                  f"{'YES' if row['has_both_folds'] else 'no':>8s} "
                  f"{row['TMdiff_spread']:>8.4f} {row['TMdiff_range']:>8.4f}"
                  f"{div_mark}")

    # --- Cross-method comparison ---
    if len(df_summary["method"].unique()) > 1:
        print(f"\n{'='*70}")
        print("CROSS-METHOD COMPARISON")
        print(f"{'='*70}")

        # Build a clean pair x method table
        pairs_seen = df_summary["pair_id"].unique()
        print(f"\n  {'pair':<20s}", end="")
        for m in ["AF2", "ESM", "MSAT", "CCMpred"]:
            print(f"  {m:>12s}", end="")
        print()

        for pid in sorted(pairs_seen):
            print(f"  {pid:<20s}", end="")
            for m in ["AF2", "ESM", "MSAT", "CCMpred"]:
                row = df_summary[(df_summary["pair_id"] == pid) &
                                 (df_summary["method"] == m)]
                if row.empty:
                    print(f"  {'--':>12s}", end="")
                else:
                    r = row.iloc[0]
                    label = f"{r['n_f1']}F1/{r['n_f2']}F2/{r['n_amb']}A"
                    print(f"  {label:>12s}", end="")
            print()

    print(f"\n{'='*70}")
    print("CONCLUSION")
    print(f"{'='*70}")
    print("If most pairs show 'uniform' (no diversity) across ALL methods,")
    print("then fold preference is genuinely not separable by MSA clustering.")
    print("If MSAT shows more diversity than AF2, it suggests co-evolutionary")
    print("signal exists but AF2 ignores it (consistent with Porter et al.).")
    print("If ESM shows more diversity than AF2, it suggests single-sequence")
    print("composition differs between clusters even when AF2 doesn't respond.")
    print()


if __name__ == "__main__":
    main()
