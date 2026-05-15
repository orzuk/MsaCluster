"""Sequence-divergence regression correction for per-cluster fold preference.

Addresses the training-set-proximity confound (Section 7A of RESULTS_STATUS,
\\zuk L466 of paper): per-cluster centered TM-diff may be a linear function of
cluster sequence-similarity-to-query (clusters closer to AF2 training data get
higher-confidence "default-fold" predictions). For each pair x method we fit

    TMdiff_centered_c = beta * pid_to_query_c + alpha + r_c

and re-classify clusters using the residual r_c. Pairs whose signal survives
correction are robust against training proximity; pairs whose signal washes
out had a confounded apparent signal.

Inputs (read from cwd):
  - docs/fold_diversity_survey.csv       per-cluster, has TMdiff_centered
  - Pipeline/FoldPairs/<pair>/output_msa_cluster/ShallowMsa_*.a3m
  - Pipeline/FoldPairs/<pair>/output_get_msa/DeepMsa.a3m (query = first seq)

Outputs:
  - docs/per_cluster_pid_to_query.csv         cached PID per (pair, cluster)
  - docs/fold_diversity_survey_corrected.csv  per-cluster + residual + pref_corrected
  - docs/fold_diversity_concordance_corrected.csv  per-pair, residual-sign concordance
  - prints a before/after comparison for high-confidence candidate pairs

Usage:
    python3 scripts/seq_divergence_correction.py
    python3 scripts/seq_divergence_correction.py --rebuild-pid-cache
    python3 scripts/seq_divergence_correction.py --delta 0.05
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd

REPO = Path(".")
PIPELINE_DIR = REPO / "Pipeline" / "FoldPairs"
SURVEY_CSV = REPO / "docs" / "fold_diversity_survey.csv"
PID_CACHE = REPO / "docs" / "per_cluster_pid_to_query.csv"
OUT_SURVEY = REPO / "docs" / "fold_diversity_survey_corrected.csv"
OUT_CONCORD = REPO / "docs" / "fold_diversity_concordance_corrected.csv"

CANDIDATE_PAIRS = ["2n54B_2hdmA", "2namA_1uxmK", "1jfkA_2nxqB", "3kdsG_2ce7C"]
METHODS = ["AF2", "ESM", "MSAT", "DDG"]
_MSA_TAG_RE = re.compile(r"(\d{1,4})")


# ---------------------------------------------------------------------------
# Per-cluster PID-to-query (cluster sequence-similarity to deep-MSA seed)
# ---------------------------------------------------------------------------

def normalize_cluster(s: str) -> str:
    s = str(s).strip()
    s = re.sub(r"^(msa_t_)?(clusters_|cmap_)?", "", s)
    if s.startswith("DeepMsa") or s.lower() in ("deep", "deepmsa", "query"):
        return "DeepMsa"
    m = re.search(r"ShallowMsa_(\d+)", s)
    if m:
        return f"ShallowMsa_{int(m.group(1)):03d}"
    if s.isdigit():
        return f"ShallowMsa_{int(s):03d}"
    return s


def load_a3m_seqs(path: Path) -> list[str]:
    """Return raw a3m sequences in order. First entry is the query."""
    out, cur = [], []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    out.append("".join(cur))
                    cur = []
            else:
                cur.append(line)
    if cur:
        out.append("".join(cur))
    return out


def gapaware_pid(query: str, seq: str) -> float:
    """Percent identity over match columns (upper-case, non-gap) of query."""
    nm = nt = 0
    for q, s in zip(query, seq):
        if q == "-" or q.islower():
            continue
        if s == "-" or s.islower():
            continue
        nt += 1
        if q.upper() == s.upper():
            nm += 1
    return nm / nt if nt > 0 else float("nan")


def compute_pair_pids(pair_id: str) -> list[dict]:
    pair_dir = PIPELINE_DIR / pair_id
    deep = pair_dir / "output_get_msa" / "DeepMsa.a3m"
    if not deep.exists():
        return []
    deep_seqs = load_a3m_seqs(deep)
    if not deep_seqs:
        return []
    query = deep_seqs[0]
    msa_dir = pair_dir / "output_msa_cluster"
    if not msa_dir.exists():
        return []
    rows = []
    for f in sorted(msa_dir.glob("ShallowMsa_*.a3m")):
        cid_num = re.search(r"ShallowMsa_(\d+)", f.stem)
        if not cid_num:
            continue
        cid = f"ShallowMsa_{int(cid_num.group(1)):03d}"
        seqs = load_a3m_seqs(f)
        if len(seqs) <= 1:
            continue
        pids = [gapaware_pid(query, s) for s in seqs[1:]]
        pids = [p for p in pids if not np.isnan(p)]
        if not pids:
            continue
        rows.append({
            "pair_id": pair_id,
            "cluster": cid,
            "n_seqs": len(pids),
            "mean_pid_to_query": float(np.mean(pids)),
            "median_pid_to_query": float(np.median(pids)),
        })
    return rows


def build_pid_cache() -> pd.DataFrame:
    if not PIPELINE_DIR.exists():
        raise SystemExit(f"missing {PIPELINE_DIR} — run from repo root")
    pairs = sorted(d.name for d in PIPELINE_DIR.iterdir()
                   if d.is_dir() and d.name != "jobs")
    rows = []
    for i, p in enumerate(pairs):
        r = compute_pair_pids(p)
        rows.extend(r)
        if (i + 1) % 10 == 0 or i == len(pairs) - 1:
            print(f"  [{i+1}/{len(pairs)}] {p}: {len(r)} clusters")
    df = pd.DataFrame(rows)
    PID_CACHE.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(PID_CACHE, index=False)
    print(f"wrote {PID_CACHE}: {len(df)} rows")
    return df


# ---------------------------------------------------------------------------
# Regression correction
# ---------------------------------------------------------------------------

def apply_regression(survey: pd.DataFrame, pid_df: pd.DataFrame, delta: float) -> pd.DataFrame:
    """Per (pair, method) regress TMdiff_centered on mean_pid_to_query, classify residual."""
    survey = survey.copy()
    survey["cluster_norm"] = survey["cluster"].map(normalize_cluster)
    pid_df = pid_df.copy()
    pid_df["cluster_norm"] = pid_df["cluster"].map(normalize_cluster)
    merged = survey.merge(
        pid_df[["pair_id", "cluster_norm", "mean_pid_to_query"]],
        on=["pair_id", "cluster_norm"], how="left",
    )
    merged["TMdiff_residual"] = np.nan
    merged["pref_corrected"] = "Amb"
    merged["regression_beta"] = np.nan
    merged["regression_alpha"] = np.nan

    for (pid, method), sub in merged.groupby(["pair_id", "method"]):
        sub_valid = sub.dropna(subset=["TMdiff_centered", "mean_pid_to_query"])
        sub_valid = sub_valid[sub_valid["cluster_norm"] != "DeepMsa"]
        if len(sub_valid) < 3:
            merged.loc[sub.index, "TMdiff_residual"] = sub["TMdiff_centered"]
            continue
        x = sub_valid["mean_pid_to_query"].to_numpy(dtype=float)
        y = sub_valid["TMdiff_centered"].to_numpy(dtype=float)
        if np.std(x) < 1e-9:
            merged.loc[sub.index, "TMdiff_residual"] = sub["TMdiff_centered"]
            continue
        beta, alpha = np.polyfit(x, y, 1)
        resid = sub["TMdiff_centered"].to_numpy(dtype=float) - (
            beta * sub["mean_pid_to_query"].to_numpy(dtype=float) + alpha
        )
        merged.loc[sub.index, "TMdiff_residual"] = resid
        merged.loc[sub.index, "regression_beta"] = beta
        merged.loc[sub.index, "regression_alpha"] = alpha

    def _classify(v):
        if pd.isna(v):
            return "Amb"
        if v > delta:
            return "F1"
        if v < -delta:
            return "F2"
        return "Amb"

    merged["pref_corrected"] = merged["TMdiff_residual"].map(_classify)
    return merged


# ---------------------------------------------------------------------------
# Concordance on residual signs
# ---------------------------------------------------------------------------

def compute_concordance(corrected: pd.DataFrame) -> pd.DataFrame:
    """Per-pair: mean pairwise sign-agreement of TMdiff_residual across methods."""
    rows = []
    for pid, sub in corrected.groupby("pair_id"):
        sub = sub[sub["cluster_norm"] != "DeepMsa"]
        pivot = sub.pivot_table(
            index="cluster_norm", columns="method", values="TMdiff_residual", aggfunc="first",
        )
        present_methods = [m for m in METHODS if m in pivot.columns]
        if len(present_methods) < 2:
            continue
        pair_pcts = {}
        for a, b in combinations(present_methods, 2):
            both = pivot[[a, b]].dropna()
            if len(both) < 1:
                continue
            sa = np.sign(both[a].to_numpy())
            sb = np.sign(both[b].to_numpy())
            mask = (sa != 0) & (sb != 0)
            if mask.sum() < 1:
                continue
            agree = (sa[mask] == sb[mask]).mean()
            pair_pcts[f"{a}~{b}"] = (float(agree), int(mask.sum()))
        if not pair_pcts:
            continue
        vals = [v[0] for v in pair_pcts.values()]
        ns = [v[1] for v in pair_pcts.values()]
        row = {
            "pair_id": pid,
            "mean_concordance_corrected": float(np.mean(vals)),
            "mean_n_clusters_corrected": float(np.mean(ns)),
            "n_method_pairs": len(pair_pcts),
        }
        for k, (a, n) in pair_pcts.items():
            row[f"{k}_pct_corrected"] = a
            row[f"{k}_n_corrected"] = n
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def global_summary(corrected: pd.DataFrame, before_concord: pd.DataFrame,
                   after_concord: pd.DataFrame, delta: float) -> None:
    """All-pairs summaries: per-method diversity, biggest movers, top-10 corrected."""
    # (1) Per-method centered-diverse counts: before vs after
    print("\n=== Per-method centered-diverse counts (all pairs): before / after ===")
    print(f"{'method':<8} {'coverage':>10} {'before':>10} {'after':>10} {'delta':>8}")
    for m in METHODS:
        sub = corrected[corrected["method"] == m]
        coverage = sub["pair_id"].nunique()
        per_pair = sub.groupby("pair_id").agg(
            before_F1=("pref_centered", lambda s: int((s == "F1").sum())),
            before_F2=("pref_centered", lambda s: int((s == "F2").sum())),
            after_F1=("pref_corrected", lambda s: int((s == "F1").sum())),
            after_F2=("pref_corrected", lambda s: int((s == "F2").sum())),
        )
        before_div = int(((per_pair["before_F1"] > 0) & (per_pair["before_F2"] > 0)).sum())
        after_div = int(((per_pair["after_F1"] > 0) & (per_pair["after_F2"] > 0)).sum())
        print(f"{m:<8} {coverage:>10d} {before_div:>10d} {after_div:>10d} {after_div - before_div:>+8d}")

    # (2) Biggest concordance movers
    bi = before_concord.set_index("pair_id")["mean_concordance"]
    ai = after_concord.set_index("pair_id")["mean_concordance_corrected"]
    joined = pd.DataFrame({"before": bi, "after": ai}).dropna()
    joined["delta"] = joined["after"] - joined["before"]

    print("\n=== Pairs with concordance DROP > 0.10 (signal eroded by correction) ===")
    dropped = joined[joined["delta"] < -0.10].sort_values("delta")
    for pid, row in dropped.iterrows():
        print(f"  {pid:<18}  {row['before']:.3f} -> {row['after']:.3f}  ({row['delta']:+.3f})")
    if len(dropped) == 0:
        print("  (none)")

    print("\n=== Pairs with concordance GAIN > 0.10 (signal sharpened by correction) ===")
    gained = joined[joined["delta"] > 0.10].sort_values("delta", ascending=False)
    for pid, row in gained.iterrows():
        print(f"  {pid:<18}  {row['before']:.3f} -> {row['after']:.3f}  ({row['delta']:+.3f})")
    if len(gained) == 0:
        print("  (none)")

    # (3) Top-10 by corrected concordance
    print("\n=== Top-10 pairs by CORRECTED mean concordance ===")
    top = joined.sort_values("after", ascending=False).head(10)
    for pid, row in top.iterrows():
        flag = "  <-- HIGH CONFIDENCE" if row["after"] >= 0.65 else ""
        print(f"  {pid:<18}  before={row['before']:.3f}  after={row['after']:.3f}{flag}")


def before_after_table(corrected: pd.DataFrame, before_concord: pd.DataFrame,
                       after_concord: pd.DataFrame, pairs: list[str]) -> None:
    print("\n=== Before / after regression-correction (candidate pairs) ===")
    print(f"{'pair_id':<18} {'method':<6} "
          f"{'n_F1c→':>8} {'n_F2c→':>8} {'n_F1c*':>8} {'n_F2c*':>8} {'beta':>8}")
    for p in pairs:
        sub = corrected[corrected["pair_id"] == p]
        for m in METHODS:
            ms = sub[sub["method"] == m]
            if not len(ms):
                continue
            n_f1_before = int((ms["pref_centered"] == "F1").sum())
            n_f2_before = int((ms["pref_centered"] == "F2").sum())
            n_f1_after = int((ms["pref_corrected"] == "F1").sum())
            n_f2_after = int((ms["pref_corrected"] == "F2").sum())
            beta = ms["regression_beta"].dropna()
            beta_s = f"{beta.iloc[0]:+.3f}" if len(beta) else "  --  "
            print(f"{p:<18} {m:<6} "
                  f"{n_f1_before:>8d} {n_f2_before:>8d} {n_f1_after:>8d} {n_f2_after:>8d} {beta_s:>8}")
    print("\n=== Mean cross-method concordance: before vs after ===")
    print(f"{'pair_id':<18} {'before':>8} {'after':>8} {'delta':>8}")
    before_idx = before_concord.set_index("pair_id")
    after_idx = after_concord.set_index("pair_id")
    for p in pairs:
        b = before_idx.loc[p, "mean_concordance"] if p in before_idx.index else np.nan
        a = after_idx.loc[p, "mean_concordance_corrected"] if p in after_idx.index else np.nan
        d = (a - b) if (not np.isnan(a) and not np.isnan(b)) else np.nan
        print(f"{p:<18} {b:>8.3f} {a:>8.3f} {d:>+8.3f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rebuild-pid-cache", action="store_true",
                    help="rebuild docs/per_cluster_pid_to_query.csv from a3m files")
    ap.add_argument("--delta", type=float, default=0.05,
                    help="classification threshold on TM-score residual (default 0.05)")
    ap.add_argument("--candidate-pairs", default=",".join(CANDIDATE_PAIRS),
                    help="comma-separated pair IDs to print before/after for")
    args = ap.parse_args()

    if args.rebuild_pid_cache or not PID_CACHE.exists():
        print("Building per-cluster PID-to-query cache...")
        pid_df = build_pid_cache()
    else:
        pid_df = pd.read_csv(PID_CACHE)
        print(f"Loaded {PID_CACHE}: {len(pid_df)} rows")

    survey = pd.read_csv(SURVEY_CSV)
    print(f"Loaded {SURVEY_CSV}: {len(survey)} rows")

    if "TMdiff_centered" not in survey.columns:
        raise SystemExit("survey CSV missing TMdiff_centered column — re-run fold_diversity_survey.py first")

    corrected = apply_regression(survey, pid_df, args.delta)
    corrected.to_csv(OUT_SURVEY, index=False)
    print(f"wrote {OUT_SURVEY}: {len(corrected)} rows")

    after_concord = compute_concordance(corrected)
    after_concord.to_csv(OUT_CONCORD, index=False)
    print(f"wrote {OUT_CONCORD}: {len(after_concord)} rows")

    before_concord_path = REPO / "docs" / "fold_diversity_concordance.csv"
    if before_concord_path.exists():
        before_concord = pd.read_csv(before_concord_path)
    else:
        before_concord = pd.DataFrame(columns=["pair_id", "mean_concordance"])

    pairs = [p.strip() for p in args.candidate_pairs.split(",") if p.strip()]
    global_summary(corrected, before_concord, after_concord, args.delta)
    before_after_table(corrected, before_concord, after_concord, pairs)


if __name__ == "__main__":
    main()
