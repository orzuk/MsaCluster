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
PID_CACHE = REPO / "docs" / "per_cluster_pid_symmetric.csv"   # F1+F2 symmetric
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


# --- F2-chain extraction (symmetric correction) ---------------------------

def _get_chain_seq_from_pdb(pdb_path: Path, chain_id: str) -> str:
    """Extract amino-acid sequence of a given chain from a PDB file.
    Returns "" on failure. Used for getting the F2-truth chain sequence
    so the regression can use both PID-to-F1 and PID-to-F2 as predictors
    (symmetric correction)."""
    try:
        from Bio.PDB import PDBParser, PPBuilder
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("s", str(pdb_path))
        ppb = PPBuilder()
        for chain in structure[0]:
            if chain.id == chain_id:
                seq = ""
                for pp in ppb.build_peptides(chain):
                    seq += str(pp.get_sequence())
                return seq
    except Exception:
        pass
    return ""


def _get_F2_chain_seq(pair_id: str) -> str:
    """For pair "<pdb1><chain1>_<pdb2><chain2>", return F2 (second) chain seq."""
    if "_" not in pair_id:
        return ""
    foldA, foldB = pair_id.split("_", 1)
    # Allow chain IDs of 1 character (most common)
    if len(foldB) < 5:
        return ""
    pdb_id, chain_id = foldB[:4].lower(), foldB[4:]
    for cand in (PIPELINE_DIR / pair_id / f"{pdb_id}.pdb",
                 PIPELINE_DIR / pair_id / f"{pdb_id.upper()}.pdb"):
        if cand.is_file():
            return _get_chain_seq_from_pdb(cand, chain_id)
    return ""


def _build_F1_F2_column_map(seq_F1_ungap: str, seq_F2_ungap: str
                            ) -> list[int | None]:
    """Pairwise-align F1 to F2 ungapped sequences. Returns a list whose
    i-th entry is the F2 position aligned to F1 position i (or None if F1
    position i aligns to a gap in F2)."""
    if not seq_F1_ungap or not seq_F2_ungap:
        return []
    try:
        from Bio import Align
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -1
        alignment = aligner.align(seq_F1_ungap, seq_F2_ungap)[0]
        s1, s2 = str(alignment).strip().split("\n")[0], \
                 str(alignment).strip().split("\n")[2]
        mapping: list[int | None] = []
        f1_pos = f2_pos = -1
        for c1, c2 in zip(s1, s2):
            if c1 != "-":
                f1_pos += 1
            if c2 != "-":
                f2_pos += 1
            if c1 != "-":
                mapping.append(f2_pos if c2 != "-" else None)
        return mapping
    except Exception:
        return []


def _pid_to_F2_via_map(cluster_seq: str, query_F1: str,
                       F2_seq: str, f1_to_f2_map: list) -> float:
    """PID of cluster_seq to F2 chain via the F1<->F2 column mapping.
    cluster_seq and query_F1 are column-aligned a3m strings; F2_seq is
    ungapped."""
    nm = nt = 0
    f1_pos = -1
    for q, s in zip(query_F1, cluster_seq):
        if q == "-" or q.islower():
            continue
        f1_pos += 1
        if s == "-" or s.islower():
            continue
        if f1_pos >= len(f1_to_f2_map):
            break
        f2_pos = f1_to_f2_map[f1_pos]
        if f2_pos is None or f2_pos >= len(F2_seq):
            continue
        nt += 1
        if s.upper() == F2_seq[f2_pos].upper():
            nm += 1
    return nm / nt if nt > 0 else float("nan")


def compute_pair_pids(pair_id: str) -> list[dict]:
    """Compute per-cluster mean PID to BOTH F1 chain (= MSA seed) and F2
    chain (= second truth structure's chain). Symmetric correction needs
    both as predictors in the regression."""
    pair_dir = PIPELINE_DIR / pair_id
    deep = pair_dir / "output_get_msa" / "DeepMsa.a3m"
    if not deep.exists():
        return []
    deep_seqs = load_a3m_seqs(deep)
    if not deep_seqs:
        return []
    query = deep_seqs[0]                     # F1 chain (a3m, with gaps)
    # F1 ungapped (lowercase = insertions, kept as upper-only for alignment)
    query_ungap = "".join(c for c in query if c not in "-" and not c.islower())
    seq_F2 = _get_F2_chain_seq(pair_id)
    f1_to_f2 = _build_F1_F2_column_map(query_ungap, seq_F2) if seq_F2 else []
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
        # PID to F1 (original behavior): column-wise to query
        pids_F1 = [gapaware_pid(query, s) for s in seqs[1:]]
        pids_F1 = [p for p in pids_F1 if not np.isnan(p)]
        # PID to F2 (new): via F1->F2 column mapping
        if seq_F2 and f1_to_f2:
            pids_F2 = [_pid_to_F2_via_map(s, query, seq_F2, f1_to_f2)
                       for s in seqs[1:]]
            pids_F2 = [p for p in pids_F2 if not np.isnan(p)]
        else:
            pids_F2 = []
        if not pids_F1:
            continue
        rows.append({
            "pair_id": pair_id,
            "cluster": cid,
            "n_seqs": len(pids_F1),
            "mean_pid_to_F1": float(np.mean(pids_F1)),
            "median_pid_to_F1": float(np.median(pids_F1)),
            "mean_pid_to_F2": float(np.mean(pids_F2)) if pids_F2 else float("nan"),
            "median_pid_to_F2": float(np.median(pids_F2)) if pids_F2 else float("nan"),
            # Legacy alias for compatibility with downstream consumers
            "mean_pid_to_query": float(np.mean(pids_F1)),
            "median_pid_to_query": float(np.median(pids_F1)),
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

def apply_regression(survey: pd.DataFrame, pid_df: pd.DataFrame,
                     delta_tm: float, delta_ddg: float) -> pd.DataFrame:
    """Per (pair, method) regress TMdiff_centered on BOTH mean_pid_to_F1 AND
    mean_pid_to_F2 (symmetric correction). Falls back to single-predictor
    PID-to-F1 only when F2 chain sequence can't be extracted.

    Model (symmetric):
        TMdiff_centered_c = alpha + beta1 * pid_to_F1_c + beta2 * pid_to_F2_c + r_c

    Note on scales: AF2/ESM/MSAT TMdiff_centered is on TM-score scale [-1, 1];
    DDG TMdiff_centered is actually centered Delta-Delta-G in kcal/mol, typical
    range much wider. We therefore apply delta_tm to AF2/ESM/MSAT residuals and
    delta_ddg to DDG residuals, matching the original survey thresholds.
    """
    survey = survey.copy()
    survey["cluster_norm"] = survey["cluster"].map(normalize_cluster)
    pid_df = pid_df.copy()
    pid_df["cluster_norm"] = pid_df["cluster"].map(normalize_cluster)
    # Carry both F1 and F2 PID columns (back-compat: fall back to legacy
    # "mean_pid_to_query" if new columns aren't present in cache)
    pid_cols = ["pair_id", "cluster_norm"]
    for c in ("mean_pid_to_F1", "mean_pid_to_F2", "mean_pid_to_query"):
        if c in pid_df.columns:
            pid_cols.append(c)
    merged = survey.merge(pid_df[pid_cols],
                          on=["pair_id", "cluster_norm"], how="left")
    # Unify "PID to F1" column: prefer new column, fall back to legacy
    if "mean_pid_to_F1" not in merged.columns:
        merged["mean_pid_to_F1"] = merged.get("mean_pid_to_query", np.nan)
    if "mean_pid_to_F2" not in merged.columns:
        merged["mean_pid_to_F2"] = np.nan
    merged["TMdiff_residual"] = np.nan
    merged["pref_corrected"] = "Amb"
    merged["regression_beta_F1"] = np.nan
    merged["regression_beta_F2"] = np.nan
    merged["regression_alpha"] = np.nan
    merged["regression_n_predictors"] = 0

    for (_, method), sub in merged.groupby(["pair_id", "method"]):
        sub_valid = sub.dropna(subset=["TMdiff_centered", "mean_pid_to_F1"])
        sub_valid = sub_valid[sub_valid["cluster_norm"] != "DeepMsa"]
        if len(sub_valid) < 3:
            merged.loc[sub.index, "TMdiff_residual"] = sub["TMdiff_centered"]
            continue

        y = sub_valid["TMdiff_centered"].to_numpy(dtype=float)
        x1 = sub_valid["mean_pid_to_F1"].to_numpy(dtype=float)
        x2 = sub_valid["mean_pid_to_F2"].to_numpy(dtype=float)
        has_F2 = not np.isnan(x2).all()

        if has_F2 and np.std(x2) >= 1e-9 and np.std(x1) >= 1e-9:
            # Symmetric 2-predictor regression via least squares
            X = np.column_stack([np.ones(len(y)), x1, x2])
            try:
                coefs, *_ = np.linalg.lstsq(X, y, rcond=None)
                alpha, beta1, beta2 = float(coefs[0]), float(coefs[1]), float(coefs[2])
            except Exception:
                merged.loc[sub.index, "TMdiff_residual"] = sub["TMdiff_centered"]
                continue
            full_x1 = sub["mean_pid_to_F1"].to_numpy(dtype=float)
            full_x2 = sub["mean_pid_to_F2"].to_numpy(dtype=float)
            pred = alpha + beta1 * full_x1 + beta2 * full_x2
            resid = sub["TMdiff_centered"].to_numpy(dtype=float) - pred
            merged.loc[sub.index, "TMdiff_residual"] = resid
            merged.loc[sub.index, "regression_alpha"] = alpha
            merged.loc[sub.index, "regression_beta_F1"] = beta1
            merged.loc[sub.index, "regression_beta_F2"] = beta2
            merged.loc[sub.index, "regression_n_predictors"] = 2
        else:
            # Fallback: single-predictor F1 only (legacy behavior)
            if np.std(x1) < 1e-9:
                merged.loc[sub.index, "TMdiff_residual"] = sub["TMdiff_centered"]
                continue
            beta1, alpha = np.polyfit(x1, y, 1)
            full_x1 = sub["mean_pid_to_F1"].to_numpy(dtype=float)
            resid = sub["TMdiff_centered"].to_numpy(dtype=float) - (
                beta1 * full_x1 + alpha
            )
            merged.loc[sub.index, "TMdiff_residual"] = resid
            merged.loc[sub.index, "regression_alpha"] = float(alpha)
            merged.loc[sub.index, "regression_beta_F1"] = float(beta1)
            merged.loc[sub.index, "regression_n_predictors"] = 1

    def _classify(v, d):
        if pd.isna(v):
            return "Amb"
        if v > d:
            return "F1"
        if v < -d:
            return "F2"
        return "Amb"

    def _row_classify(row):
        d = delta_ddg if row["method"] == "DDG" else delta_tm
        return _classify(row["TMdiff_residual"], d)

    merged["pref_corrected"] = merged.apply(_row_classify, axis=1)
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

def threshold_sensitivity(survey: pd.DataFrame, pid_df: pd.DataFrame,
                          deltas: list[float], delta_ddg: float) -> None:
    """Run the regression correction at multiple TM-score thresholds.

    Reports per-method centered-diverse pair counts (raw vs corrected)
    and Venn 4-way intersection size at each delta. Concordance itself
    is sign-based and delta-invariant, so the high-confidence pair
    count by mean_concordance >= 0.65 is unchanged across deltas; the
    diversity counts that feed Table 1 / Venn are not.
    """
    from itertools import combinations
    print("\n=== Threshold sensitivity (per-method centered-diverse counts) ===")
    print(f"{'delta':<7} {'method':<6} {'raw':>5} {'corrected':>10}   {'corr-raw':>9}")
    venn_sizes = {}
    for d in deltas:
        corrected = apply_regression(survey, pid_df, d, delta_ddg)
        sets = {}
        for m in METHODS:
            sub = corrected[corrected["method"] == m]
            per_pair = sub.groupby("pair_id").agg(
                raw_F1=("pref_centered", lambda s: int((s == "F1").sum())),
                raw_F2=("pref_centered", lambda s: int((s == "F2").sum())),
                cor_F1=("pref_corrected", lambda s: int((s == "F1").sum())),
                cor_F2=("pref_corrected", lambda s: int((s == "F2").sum())),
            )
            raw_div = int(((per_pair["raw_F1"] > 0) & (per_pair["raw_F2"] > 0)).sum())
            cor_div = int(((per_pair["cor_F1"] > 0) & (per_pair["cor_F2"] > 0)).sum())
            sets[m] = set(per_pair[(per_pair["cor_F1"] > 0) & (per_pair["cor_F2"] > 0)].index)
            print(f"{d:<7.2f} {m:<6} {raw_div:>5d} {cor_div:>10d}   {cor_div - raw_div:>+9d}")
        venn = set.intersection(*sets.values()) if sets else set()
        venn_sizes[d] = (len(venn), sorted(venn))
        print(f"  4-way intersection at delta={d:.2f}: {len(venn)} pairs  {sorted(venn)}")
    print("\n=== Venn 4-way intersection size by delta ===")
    for d, (n, pairs) in venn_sizes.items():
        print(f"  delta={d:.2f}: {n} pairs  {pairs}")


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
    ap.add_argument("--delta-ddg", type=float, default=2.0,
                    help="classification threshold on DDG residual in kcal/mol (default 2.0)")
    ap.add_argument("--delta-sweep", action="store_true",
                    help="run threshold-sensitivity sweep at delta in {0.03, 0.05, 0.07}")
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

    corrected = apply_regression(survey, pid_df, args.delta, args.delta_ddg)
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
    if args.delta_sweep:
        threshold_sensitivity(survey, pid_df, [0.03, 0.05, 0.07], args.delta_ddg)


if __name__ == "__main__":
    main()
