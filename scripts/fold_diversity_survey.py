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

try:
    from utils.trigger_utils import bound_fold, condition_spec
except Exception:                    # trigger CSV optional at import time
    bound_fold = lambda pid: None    # noqa: E731
    condition_spec = lambda pid: {}  # noqa: E731


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_MSA_TAG_RE = re.compile(r"(?:Shallow)?Msa_(\d+)", re.IGNORECASE)


def _normalize_cluster(s):
    """Normalize cluster tag -> ShallowMsa_NNN or DeepMsa.

    Handles several observed naming conventions across methods/files:
      - "000", "0", 0                       -> "ShallowMsa_000"
      - "ShallowMsa_000", "ShallowMsa_0"    -> "ShallowMsa_000"
      - "DeepMsa", "DeepMsa_DeepMsa"        -> "DeepMsa"
      - "clusters_ShallowMsa_007"           -> "ShallowMsa_007"
        (MSAT cmap CSVs use this prefix from msa_t_clusters_ShallowMsa_*.npz)
      - "msa_t_clusters_ShallowMsa_007"     -> "ShallowMsa_007"
      - "msa_t_DeepMsa_DeepMsa"             -> "DeepMsa"
    """
    s = str(s).strip()
    # Strip common file-derived prefixes
    s = re.sub(r"^(msa_t_)?(clusters_|cmap_)?", "", s)
    # Handle DeepMsa-flavored variants up front
    if s.startswith("DeepMsa") or s.lower() in ("deep", "deepmsa", "msa_deep", "deep_msa", "query"):
        return "DeepMsa"
    # Standard "<digits>" or "<dd>" forms via existing _MSA_TAG_RE
    m = _MSA_TAG_RE.fullmatch(s)
    if m:
        return f"ShallowMsa_{int(m.group(1)):03d}"
    # Fallback: extract ShallowMsa_NNN from anywhere in the string
    m2 = re.search(r"ShallowMsa_(\d+)", s)
    if m2:
        return f"ShallowMsa_{int(m2.group(1)):03d}"
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


def _apply_centered_classification(rows, delta):
    """Add per-pair-CENTERED classification fields to per-cluster rows.

    Motivation
    ----------
    The raw `pref` (F1/F2/Amb based on absolute TMdiff_max) can be confounded
    by systematic method bias — e.g., AF2 trained more on one fold's family,
    ESM defaulting to a single-sequence prior, MSAT favoring whichever fold's
    contacts survive alignment intersection, ThermoMPNN's intrinsic stability
    differential between two backbones. In all those cases EVERY cluster
    looks like it prefers the same fold in absolute terms, washing out the
    "which clusters lean MORE F1 than others" signal we actually want.

    Fix: per pair (and per method), compute the median TMdiff_max across
    ShallowMsa clusters, subtract it, and reclassify with the same delta.

    A cluster that's "F1_centered" leans MORE toward F1 than the typical
    cluster in this family, even if all clusters absolutely prefer F2.

    The spread/range stats are centering-invariant and unchanged.

    Adds fields per row:
      TMdiff_centered (float)        — TMdiff_max minus per-pair shallow median
      pref_centered   ("F1"|"F2"|"Amb")
      pair_median_used (float)       — the centering value, for transparency
    """
    if not rows:
        return rows

    shallow_diffs = [r["TMdiff_max"] for r in rows
                     if r.get("cluster") != "DeepMsa"
                     and not (r.get("TMdiff_max") is None or
                              (isinstance(r["TMdiff_max"], float) and np.isnan(r["TMdiff_max"])))]
    if len(shallow_diffs) >= 2:
        pair_median = float(np.median(shallow_diffs))
    else:
        pair_median = 0.0   # too few clusters to center meaningfully

    for r in rows:
        td = r.get("TMdiff_max")
        if td is None or (isinstance(td, float) and np.isnan(td)):
            r["TMdiff_centered"] = float("nan")
            r["pref_centered"] = "Amb"
        else:
            # DeepMsa stays uncentered (it's the family-wide reference).
            centered = float(td) - (pair_median if r.get("cluster") != "DeepMsa" else 0.0)
            r["TMdiff_centered"] = round(centered, 4)
            if centered > delta:
                r["pref_centered"] = "F1"
            elif centered < -delta:
                r["pref_centered"] = "F2"
            else:
                r["pref_centered"] = "Amb"
        r["pair_median_used"] = round(pair_median, 4)
    return rows


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
# AF3 loader (uses same df_af.csv; filtered by model=="AF3")
# ---------------------------------------------------------------------------

def load_af3_diversity(pair_id, delta):
    """Load df_af.csv and compute per-cluster AF3 fold preference.

    AF2 and AF3 predictions are aggregated into the same df_af.csv by
    Analysis/AF_analysis.py via the recursive glob
    output_AF/AF*/*.pdb. They are distinguished by the `model` column
    (values "AF2" and "AF3"). This loader is the AF3 counterpart of
    load_af2_diversity and shares the same downstream centered-
    classification / concordance pipeline.

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

    # Filter AF3 rows only (df_af.csv may contain both AF2 and AF3)
    if "model" not in df.columns:
        return []  # legacy CSV with no model column => no AF3 rows
    df = df[df["model"].astype(str).str.upper() == "AF3"]
    if df.empty:
        return []

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
        n_f1 = (grp["TM1"] > grp["TM2"]).sum()
        n_f2 = (grp["TM2"] > grp["TM1"]).sum()
        pref = _assign_pref(tm1_max, tm2_max, delta)
        rows.append({
            "pair_id": pair_id,
            "cluster": tag,
            "method": "AF3",
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

def load_boltz2_diversity(pair_id, delta):
    """Load df_boltz2.csv (apo Boltz-2) and compute per-cluster fold preference.

    Schema (one row per ShallowMsa cluster, produced by run_Boltz2.py or
    scripts/recover_boltz2_csv.py):
        pair_id, cluster, n_models, rank1_tm_fold1, rank1_tm_fold2,
        rank1_TMdiff, mean_tm_fold1, mean_tm_fold2,
        partner_yaml, mode, ligand_ccds

    We treat the rank-1 TM-scores as the per-cluster prediction strength
    toward each fold, parallel to load_af2_diversity. With n_models=1 per
    cluster the rank1 and mean values coincide for tm_fold1/tm_fold2.

    Holo (df_boltz2_holo.csv) is NOT commingled into the apo concordance:
    the condition-aware contrast is computed separately via
    load_boltz2_holo_diversity + compute_condition_switch under the
    --condition-switch flag (writes fold_diversity_condition_switch.csv).
    """
    csv_path = os.path.join(DATA_DIR, pair_id, "Analysis", "df_boltz2.csv")
    if not os.path.isfile(csv_path):
        return []
    try:
        df = pd.read_csv(csv_path)
    except (pd.errors.EmptyDataError, pd.errors.ParserError):
        return []
    if df.empty:
        return []

    if "rank1_tm_fold1" not in df.columns or "rank1_tm_fold2" not in df.columns:
        return []

    col = "cluster_num" if "cluster_num" in df.columns else "cluster"
    if col not in df.columns:
        return []
    df["_tag"] = df[col].apply(_normalize_cluster)

    df["TM1"] = pd.to_numeric(df["rank1_tm_fold1"], errors="coerce")
    df["TM2"] = pd.to_numeric(df["rank1_tm_fold2"], errors="coerce")

    rows = []
    for tag, grp in df.groupby("_tag"):
        tm1_max = grp["TM1"].max()
        tm2_max = grp["TM2"].max()
        tm1_mean = grp["TM1"].mean()
        tm2_mean = grp["TM2"].mean()
        n_models = int(grp["n_models"].iloc[0]) if "n_models" in grp.columns else len(grp)
        # Per-model fold assignment (argmax). Each cluster has one model in
        # the rank-1 row, so the vote is binary on that single prediction.
        n_f1 = (grp["TM1"] > grp["TM2"]).sum()
        n_f2 = (grp["TM2"] > grp["TM1"]).sum()
        pref = _assign_pref(tm1_max, tm2_max, delta)
        rows.append({
            "pair_id": pair_id,
            "cluster": tag,
            "method": "Boltz2",
            "TM1_max": round(tm1_max, 4) if pd.notna(tm1_max) else float("nan"),
            "TM2_max": round(tm2_max, 4) if pd.notna(tm2_max) else float("nan"),
            "TM1_mean": round(tm1_mean, 4) if pd.notna(tm1_mean) else float("nan"),
            "TM2_mean": round(tm2_mean, 4) if pd.notna(tm2_mean) else float("nan"),
            "TMdiff_max": round(tm1_max - tm2_max, 4)
                          if (pd.notna(tm1_max) and pd.notna(tm2_max)) else float("nan"),
            "TMdiff_mean": round(tm1_mean - tm2_mean, 4)
                           if (pd.notna(tm1_mean) and pd.notna(tm2_mean)) else float("nan"),
            "n_models": n_models,
            "n_toward_f1": int(n_f1),
            "n_toward_f2": int(n_f2),
            "vote_frac_f1": round(n_f1 / max(n_models, 1), 3) if n_models > 0 else 0,
            "pref": pref,
        })
    return rows


# ---------------------------------------------------------------------------
# Condition-aware (holo) loaders + Delta_switch scoring
# ---------------------------------------------------------------------------
#
# The methods that can encode a folding condition (AF3, Boltz-2) are run a
# second time WITH the physiological trigger (ligand / homo-oligomer copies /
# partner), producing a "holo" prediction alongside the default "apo" one.
# Fold preference is then scored not as a single static TMdiff but as the
# CHANGE in preference between conditions -- the difference-in-differences
#
#     Delta_switch = b(holo) - b(apo),   b(c) = TM(P_c, F_bound) - TM(P_c, F_free)
#
# oriented toward the fold that the trigger is known to stabilize (F_bound,
# from utils.trigger_utils.bound_fold). Delta_switch > 0 means the trigger
# pushed the prediction toward the bound fold, i.e. the method reproduced the
# switch. See the appendix "Folding-condition dependence" table in the paper.


def _load_boltz2_like(csv_path, pair_id, delta, method):
    """Shared Boltz-2 rank1 reader used for both apo and holo CSVs."""
    try:
        df = pd.read_csv(csv_path)
    except (pd.errors.EmptyDataError, pd.errors.ParserError):
        return []
    if df.empty or "rank1_tm_fold1" not in df.columns or "rank1_tm_fold2" not in df.columns:
        return []
    col = "cluster_num" if "cluster_num" in df.columns else "cluster"
    if col not in df.columns:
        return []
    df["_tag"] = df[col].apply(_normalize_cluster)
    df["TM1"] = pd.to_numeric(df["rank1_tm_fold1"], errors="coerce")
    df["TM2"] = pd.to_numeric(df["rank1_tm_fold2"], errors="coerce")
    rows = []
    for tag, grp in df.groupby("_tag"):
        tm1 = grp["TM1"].max()
        tm2 = grp["TM2"].max()
        if pd.isna(tm1) or pd.isna(tm2):
            continue
        rows.append({
            "pair_id": pair_id, "cluster": tag, "method": method,
            "TM1_max": round(float(tm1), 4), "TM2_max": round(float(tm2), 4),
            "TMdiff_max": round(float(tm1 - tm2), 4),
            "pref": _assign_pref(tm1, tm2, delta),
        })
    return rows


def load_boltz2_holo_diversity(pair_id, delta):
    """Per-cluster holo Boltz-2 rows (df_boltz2_holo.csv), method 'Boltz2_holo'."""
    csv_path = os.path.join(DATA_DIR, pair_id, "Analysis", "df_boltz2_holo.csv")
    return _load_boltz2_like(csv_path, pair_id, delta, method="Boltz2_holo")


def load_af3_holo_diversity(pair_id, delta):
    """Per-cluster holo AF3 rows.

    AF3 holo predictions are scored the same way as apo, but stored either as
    a dedicated ``df_af_holo.csv`` or as rows of ``df_af.csv`` carrying a
    ``condition == 'holo'`` column. We try both so this works regardless of
    how the (cluster-side) aggregation lands them. Method label 'AF3_holo'.
    """
    ana = os.path.join(DATA_DIR, pair_id, "Analysis")
    # (a) dedicated holo CSV
    holo_csv = os.path.join(ana, "df_af_holo.csv")
    src = None
    if os.path.isfile(holo_csv):
        src = pd.read_csv(holo_csv)
    else:
        af_csv = os.path.join(ana, "df_af.csv")
        if os.path.isfile(af_csv):
            df = pd.read_csv(af_csv)
            if "condition" in df.columns:
                src = df[(df.get("model", "").astype(str).str.upper() == "AF3")
                         & (df["condition"].astype(str).str.lower() == "holo")]
    if src is None or src.empty:
        return []
    if "TMscore_fold1" not in src.columns and "score_pdb1" in src.columns:
        src = src.rename(columns={"score_pdb1": "TMscore_fold1",
                                  "score_pdb2": "TMscore_fold2"})
    if "TMscore_fold1" not in src.columns:
        return []
    col = "cluster_num" if "cluster_num" in src.columns else "cluster"
    if col not in src.columns:
        return []
    src = src.copy()
    src["_tag"] = src[col].apply(_normalize_cluster)
    src["TM1"] = pd.to_numeric(src["TMscore_fold1"], errors="coerce")
    src["TM2"] = pd.to_numeric(src["TMscore_fold2"], errors="coerce")
    rows = []
    for tag, grp in src.groupby("_tag"):
        tm1 = grp["TM1"].max()
        tm2 = grp["TM2"].max()
        if pd.isna(tm1) or pd.isna(tm2):
            continue
        rows.append({
            "pair_id": pair_id, "cluster": tag, "method": "AF3_holo",
            "TM1_max": round(float(tm1), 4), "TM2_max": round(float(tm2), 4),
            "TMdiff_max": round(float(tm1 - tm2), 4),
            "pref": _assign_pref(tm1, tm2, delta),
        })
    return rows


def compute_condition_switch(pair_id, method, apo_rows, holo_rows, delta):
    """Delta_switch fold-preference rows for one (pair, method) with both
    apo and holo predictions.

    For every cluster present in BOTH conditions:
      TMdiff_apo/holo  = TM(F1) - TM(F2)                    (signed toward F1)
      dswitch_f1       = TMdiff_holo - TMdiff_apo           (more F1 under holo)
      dswitch_bound    = s * dswitch_f1,  s = +1 if bound fold is F1 else -1
                         (>0  => trigger pushed prediction toward the bound
                          fold, i.e. the method reproduces the switch)
      pref_switch      = toward_bound | away_from_bound | Amb    (|.| vs delta)
    Plus condition-matched TMs (holo-to-bound, apo-to-free).

    ``bound`` is utils.trigger_utils.bound_fold(pair_id); when it is unknown
    (None -- e.g. oligomerization pairs whose stoichiometry the CSV lacks) the
    bound-oriented columns are left as NaN but dswitch_f1 is still emitted.
    """
    bf = bound_fold(pair_id)
    apo = {r["cluster"]: r for r in apo_rows}
    holo = {r["cluster"]: r for r in holo_rows}
    s = None if bf not in (1, 2) else (1.0 if bf == 1 else -1.0)
    out = []
    for cl in sorted(set(apo) & set(holo)):
        a, h = apo[cl], holo[cl]
        da, dh = a.get("TMdiff_max"), h.get("TMdiff_max")
        if da is None or dh is None or np.isnan(da) or np.isnan(dh):
            continue
        dsw_f1 = round(float(dh - da), 4)
        dsw_bound = round(s * dsw_f1, 4) if s is not None else float("nan")
        if s is None or np.isnan(dsw_bound):
            pref = "NA"
        elif dsw_bound > delta:
            pref = "toward_bound"
        elif dsw_bound < -delta:
            pref = "away_from_bound"
        else:
            pref = "Amb"
        # condition-matched TMs
        tm_bound = tm_free = float("nan")
        if bf in (1, 2):
            tm_bound = h.get("TM1_max") if bf == 1 else h.get("TM2_max")
            tm_free = a.get("TM2_max") if bf == 1 else a.get("TM1_max")
        out.append({
            "pair_id": pair_id, "cluster": cl, "method": f"{method}_switch",
            "bound_fold": bf,
            "TMdiff_apo": da, "TMdiff_holo": dh,
            "dswitch_f1": dsw_f1, "dswitch_bound": dsw_bound,
            "pref_switch": pref,
            "holo_tm_bound": (round(float(tm_bound), 4)
                              if tm_bound == tm_bound else float("nan")),
            "apo_tm_free": (round(float(tm_free), 4)
                            if tm_free == tm_free else float("nan")),
        })
    return out


def load_bioemu_diversity(pair_id, delta):
    """Load df_bioemu.csv (BioEmu equilibrium ensemble) -> per-cluster preference.

    Schema (one row per ShallowMsa cluster, from run_bioemu.py):
        pair_id, cluster, n_samples, mean_tm_fold1, mean_tm_fold2,
        mean_TMdiff, best_tm_fold1, best_tm_fold2, frac_toward_f1

    Unlike the point-estimate predictors (AF2/AF3/Boltz-2), BioEmu samples the
    equilibrium ENSEMBLE, so the natural per-cluster signal is the ENSEMBLE MEAN
    preference (the population balance between the folds). We therefore map the
    ensemble means onto the standard TM1/TM2 fields and use mean_TMdiff as the
    per-cluster signal (TMdiff_max), so BioEmu flows through the centered
    classification + concordance exactly like the other methods. frac_toward_f1
    (fraction of the ensemble closer to fold 1) is carried as the vote fraction.
    """
    csv_path = os.path.join(DATA_DIR, pair_id, "Analysis", "df_bioemu.csv")
    if not os.path.isfile(csv_path):
        return []
    try:
        df = pd.read_csv(csv_path)
    except (pd.errors.EmptyDataError, pd.errors.ParserError):
        return []
    if df.empty or "mean_tm_fold1" not in df.columns or "mean_tm_fold2" not in df.columns:
        return []

    col = "cluster_num" if "cluster_num" in df.columns else "cluster"
    if col not in df.columns:
        return []
    df["_tag"] = df[col].apply(_normalize_cluster)
    df["TM1"] = pd.to_numeric(df["mean_tm_fold1"], errors="coerce")
    df["TM2"] = pd.to_numeric(df["mean_tm_fold2"], errors="coerce")

    rows = []
    for tag, grp in df.groupby("_tag"):
        tm1 = grp["TM1"].mean()
        tm2 = grp["TM2"].mean()
        n_samples = int(grp["n_samples"].iloc[0]) if "n_samples" in grp.columns else 1
        frac_f1 = float(grp["frac_toward_f1"].iloc[0]) if "frac_toward_f1" in grp.columns else float("nan")
        best1 = pd.to_numeric(grp.get("best_tm_fold1"), errors="coerce").max() if "best_tm_fold1" in grp.columns else tm1
        best2 = pd.to_numeric(grp.get("best_tm_fold2"), errors="coerce").max() if "best_tm_fold2" in grp.columns else tm2
        pref = _assign_pref(tm1, tm2, delta)
        rows.append({
            "pair_id": pair_id,
            "cluster": tag,
            "method": "BioEmu",
            # ensemble MEAN is the per-cluster signal -> use it as the "max"
            # (centering) field so BioEmu behaves like the other TM methods.
            "TM1_max": round(tm1, 4) if pd.notna(tm1) else float("nan"),
            "TM2_max": round(tm2, 4) if pd.notna(tm2) else float("nan"),
            "TM1_mean": round(tm1, 4) if pd.notna(tm1) else float("nan"),
            "TM2_mean": round(tm2, 4) if pd.notna(tm2) else float("nan"),
            "TMdiff_max": round(tm1 - tm2, 4)
                          if (pd.notna(tm1) and pd.notna(tm2)) else float("nan"),
            "TMdiff_mean": round(tm1 - tm2, 4)
                           if (pd.notna(tm1) and pd.notna(tm2)) else float("nan"),
            "best_tm_fold1": round(float(best1), 4) if pd.notna(best1) else float("nan"),
            "best_tm_fold2": round(float(best2), 4) if pd.notna(best2) else float("nan"),
            "n_models": n_samples,
            "vote_frac_f1": round(frac_f1, 3) if pd.notna(frac_f1) else float("nan"),
            "pref": pref,
        })
    return rows


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


def load_esmfold2_diversity(pair_id, delta):
    """Load ESMFold2 per-cluster TM-scores (the v2 single-sequence baseline).

    Reads docs/esmfold2_tmscores_<pair>.csv (written by
    scripts/score_esmfold2.py): ONE prediction per cluster (the medoid),
    scored vs both truth folds (columns TM_F1, TM_F2). This REPLACES the old
    ESMFold-v1 loader (`load_esm_diversity`, which read df_esm.csv). The
    method label is kept as ``"ESM"`` so the ~8 downstream scripts that
    hard-code that name keep working, but the data is now ESMFold2.

    Pairs whose ESMFold2 scoring failed (e.g. all-NaN TM) yield no rows -
    we do NOT silently fall back to stale v1 ESMFold.
    """
    csv_path = os.path.join(TABLES_RES, f"esmfold2_tmscores_{pair_id}.csv")
    if not os.path.isfile(csv_path) or os.path.getsize(csv_path) == 0:
        return []
    df = pd.read_csv(csv_path)
    if df.empty or "TM_F1" not in df.columns or "TM_F2" not in df.columns \
            or "cluster" not in df.columns:
        return []
    df["_tag"] = df["cluster"].apply(_normalize_cluster)
    df["TM1"] = pd.to_numeric(df["TM_F1"], errors="coerce")
    df["TM2"] = pd.to_numeric(df["TM_F2"], errors="coerce")

    rows = []
    for tag, grp in df.groupby("_tag"):
        g = grp.dropna(subset=["TM1", "TM2"])
        if g.empty:
            continue                       # scoring failed for this cluster
        # One prediction per cluster; if >1 row (e.g. multi-chain truth),
        # take the best score toward each fold.
        tm1 = float(g["TM1"].max())
        tm2 = float(g["TM2"].max())
        diff = tm1 - tm2
        n_models = len(g)
        n_f1 = int((g["TM1"] > g["TM2"]).sum())
        n_f2 = int((g["TM2"] > g["TM1"]).sum())
        vote_frac = n_f1 / n_models if n_models else 0.5
        if diff > delta:
            pref = "F1"
        elif diff < -delta:
            pref = "F2"
        else:
            pref = "Amb"
        rows.append({
            "pair_id": pair_id,
            "cluster": tag,
            "method": "ESM",
            "TM1_max": round(tm1, 4),
            "TM2_max": round(tm2, 4),
            "TM1_mean": round(tm1, 4),
            "TM2_mean": round(tm2, 4),
            "TMdiff_max": round(diff, 4),
            "TMdiff_mean": round(diff, 4),
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
    if not os.path.isfile(csv_path) or os.path.getsize(csv_path) == 0:
        return []

    try:
        df = pd.read_csv(csv_path)
    except (pd.errors.EmptyDataError, pd.errors.ParserError):
        return []
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
    if not os.path.isfile(csv_path) or os.path.getsize(csv_path) == 0:
        return []

    try:
        df = pd.read_csv(csv_path)
    except (pd.errors.EmptyDataError, pd.errors.ParserError):
        return []
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
# DDG (ThermoMPNN conformation-bias) loader
# ---------------------------------------------------------------------------

def load_ddg_diversity(pair_id, delta_ddg):
    """Load df_ddg.csv and compute per-cluster fold preference from DDG bias.

    df_ddg.csv schema (one row per ShallowMsa cluster):
      pair_id, cluster, n_seqs, bias_mean, bias_std, frac_F1_preferring,
      sum_F1_mean, sum_F2_mean

    bias_mean = mean over the sampled cluster sequences of
                ΣDDG(F2 backbone) − ΣDDG(F1 backbone), in kcal/mol.
      bias_mean > 0  → cluster prefers F1 (F2 destabilized more by mutations)
      bias_mean < 0  → cluster prefers F2

    *** Centering for relative-classification ***
    The absolute bias_mean is biased by the intrinsic stability difference
    between the two backbones — if F1 is generally more rigid, EVERY cluster
    will look "F1-preferring" in absolute terms. The fold-switching signal
    of interest is whether DIFFERENT CLUSTERS lean differently, not whether
    they all lean the same way. So we classify each ShallowMsa cluster's
    bias relative to the per-pair MEDIAN (over ShallowMsa clusters):

        centered_bias = bias_mean - median(bias_mean over ShallowMsa clusters)
        F1   if centered_bias >  delta_ddg
        F2   if centered_bias < -delta_ddg
        Amb  otherwise

    This is centering-invariant for spread/range (the existing diversity
    statistics) and removes the systematic-stability confound.

    `delta_ddg` is the kcal/mol threshold (default 1.0; DDG is on a
    different scale than TM-score's 0.05).

    Row-format mapping for compatibility with TM-score loaders:
      TM1_max = -sum_F1_mean   (higher = more stable on F1)
      TM2_max = -sum_F2_mean
      TMdiff_max = bias_mean   ( = TM1_max − TM2_max ; NOT centered )
    """
    csv_path = os.path.join(DATA_DIR, pair_id, "Analysis", "df_ddg.csv")
    if not os.path.isfile(csv_path) or os.path.getsize(csv_path) == 0:
        return []

    try:
        df = pd.read_csv(csv_path)
    except (pd.errors.EmptyDataError, pd.errors.ParserError):
        return []
    if df.empty:
        return []

    needed = {"cluster", "bias_mean", "n_seqs", "sum_F1_mean", "sum_F2_mean"}
    if not needed.issubset(df.columns):
        return []

    rows = []
    for _, r in df.iterrows():
        tag = _normalize_cluster(r["cluster"])
        bias = float(r["bias_mean"])
        n_seqs = int(r.get("n_seqs", 0) or 0)
        frac_f1 = float(r.get("frac_F1_preferring", 0) or 0)
        sum_f1 = float(r["sum_F1_mean"])
        sum_f2 = float(r["sum_F2_mean"])
        bias_std = float(r.get("bias_std", 0) or 0)

        # Raw classification (absolute kcal/mol)
        if np.isnan(bias):
            pref_raw = "Amb"
        elif bias > delta_ddg:
            pref_raw = "F1"
        elif bias < -delta_ddg:
            pref_raw = "F2"
        else:
            pref_raw = "Amb"

        rows.append({
            "pair_id": pair_id,
            "cluster": tag,
            "method": "DDG",
            "TM1_max": round(-sum_f1, 4),
            "TM2_max": round(-sum_f2, 4),
            "TM1_mean": round(-sum_f1, 4),
            "TM2_mean": round(-sum_f2, 4),
            "TMdiff_max": round(bias, 4),         # raw bias; centered version
            "TMdiff_mean": round(bias, 4),         # added by helper afterwards
            "n_models": n_seqs,
            "n_toward_f1": int(round(frac_f1 * n_seqs)),
            "n_toward_f2": int(round((1 - frac_f1) * n_seqs)),
            "vote_frac_f1": round(frac_f1, 3),
            "pref": pref_raw,
            # DDG-specific extras
            "ddg_bias_mean": round(bias, 4),
            "ddg_bias_std": round(bias_std, 4),
            "ddg_sum_F1_mean": round(sum_f1, 4),
            "ddg_sum_F2_mean": round(sum_f2, 4),
            "metric_used": f"ddg_kcal_per_mol(delta={delta_ddg})",
        })
    return rows


# ---------------------------------------------------------------------------
# Per-pair summary
# ---------------------------------------------------------------------------

def summarize_pair(pair_id, cluster_rows):
    """Compute pair-level diversity metrics from cluster-level rows.

    Reports BOTH classifications side by side:
      - n_f1 / n_f2 / n_amb / has_both_folds          (RAW: absolute threshold)
      - n_f1_centered / n_f2_centered / n_amb_centered / has_both_folds_centered
        (CENTERED: per-pair median subtracted; removes systematic method bias)

    Spread/range stats are centering-invariant so they apply to both.
    """
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

    # Centered prefs (added by _apply_centered_classification)
    if "pref_centered" in shallow.columns:
        prefs_c = shallow["pref_centered"].value_counts().to_dict()
    else:
        prefs_c = {}

    # Spread: std of TMdiff across clusters (centering-invariant)
    tm_diffs = shallow["TMdiff_max"].values
    spread = float(np.std(tm_diffs)) if len(tm_diffs) > 1 else 0.0
    diff_range = float(np.max(tm_diffs) - np.min(tm_diffs)) if len(tm_diffs) > 1 else 0.0

    has_f1 = prefs.get("F1", 0) > 0
    has_f2 = prefs.get("F2", 0) > 0
    has_diversity = has_f1 and has_f2

    has_f1_c = prefs_c.get("F1", 0) > 0
    has_f2_c = prefs_c.get("F2", 0) > 0
    has_diversity_c = has_f1_c and has_f2_c

    best_tm1 = shallow["TM1_max"].max()
    best_tm2 = shallow["TM2_max"].max()

    pair_median = (cluster_rows[0].get("pair_median_used")
                   if "pair_median_used" in cluster_rows[0] else None)

    return {
        "pair_id": pair_id,
        "method": method,
        "n_clusters": n_clusters,
        # raw classification
        "n_f1": prefs.get("F1", 0),
        "n_f2": prefs.get("F2", 0),
        "n_amb": prefs.get("Amb", 0),
        "has_both_folds": has_diversity,
        # centered (relative-to-pair-median) classification
        "n_f1_centered": prefs_c.get("F1", 0),
        "n_f2_centered": prefs_c.get("F2", 0),
        "n_amb_centered": prefs_c.get("Amb", 0),
        "has_both_folds_centered": has_diversity_c,
        "pair_median_used": pair_median,
        # spread/range — centering-invariant, apply to both
        "TMdiff_spread": round(spread, 4),
        "TMdiff_range": round(diff_range, 4),
        "TMdiff_mean": round(float(np.mean(tm_diffs)), 4),
        "best_score_f1": round(best_tm1, 4),
        "best_score_f2": round(best_tm2, 4),
    }


# ---------------------------------------------------------------------------
# Cross-method concordance (per-cluster sign agreement on TMdiff_centered)
# ---------------------------------------------------------------------------

def _binom_p_two_sided(k, n):
    """One-sided P[X>=k | n trials, p=0.5] via normal approximation; cheap.
    Returns 1.0 for n<5 (too small to call). Used for ranking, not formal test.
    """
    if n < 5:
        return 1.0
    from math import erfc, sqrt
    # P(X >= k) under Binom(n, 0.5) — normal approx with continuity correction
    z = (k - 0.5 - n/2) / (sqrt(n)/2)
    return 0.5 * erfc(z / sqrt(2))


def compute_pair_concordance(pair_id, cluster_rows_by_method):
    """For one pair, compute per-cluster sign agreement of TMdiff_centered
    across every pair of methods. Returns one dict (one row per pair).

    cluster_rows_by_method: dict {"AF2": [rows...], "ESM": [...], ...}

    Concordance for method-pair (A, B):
      - For each cluster present in both A's and B's rows where both
        TMdiff_centered values are nonzero, count "agree" if signs match
      - agree_pct = agree / (agree + disagree)
      - n = clusters compared
      - p = approximate one-sided binomial p-value vs null p=0.5

    Schema (per pair): pair_id, n_methods,
                       AF2_ESM_pct, AF2_ESM_n, AF2_MSAT_pct, ..., DDG__pct,
                       mean_concordance, mean_concordance_n,
                       max_concordance_methods (e.g. "AF2~DDG"),
                       max_concordance_pct
    """
    methods = ["AF2", "AF3", "ESM", "MSAT", "CCMpred", "DDG", "Boltz2", "BioEmu"]
    # Build per-method dict cluster -> centered value
    by_m = {}
    for m in methods:
        by_m[m] = {}
        for r in cluster_rows_by_method.get(m, []):
            if r.get("cluster") == "DeepMsa":
                continue
            v = r.get("TMdiff_centered")
            if v is None or (isinstance(v, float) and np.isnan(v)):
                continue
            by_m[m][r["cluster"]] = float(v)

    out = {"pair_id": pair_id, "n_methods": sum(1 for m in methods if by_m[m])}

    pair_pcts = []
    method_pairs_present = []
    for i, A in enumerate(methods):
        for B in methods[i+1:]:
            if not by_m[A] or not by_m[B]:
                continue
            shared = set(by_m[A]) & set(by_m[B])
            agree = disagree = 0
            for cl in shared:
                a, b = by_m[A][cl], by_m[B][cl]
                if a == 0 or b == 0:
                    continue
                if (a > 0) == (b > 0):
                    agree += 1
                else:
                    disagree += 1
            n = agree + disagree
            pct = (agree / n) if n > 0 else float("nan")
            key = f"{A}_{B}"
            out[f"{key}_pct"] = round(pct, 3) if n > 0 else None
            out[f"{key}_n"] = n
            out[f"{key}_p"] = round(_binom_p_two_sided(agree, n), 4) if n > 0 else None
            if n > 0:
                pair_pcts.append((pct, n, f"{A}~{B}"))
                method_pairs_present.append(f"{A}~{B}")

    if pair_pcts:
        # Mean concordance across method-pairs
        mean_pct = sum(p for p, _, _ in pair_pcts) / len(pair_pcts)
        out["mean_concordance"] = round(mean_pct, 3)
        out["mean_n_clusters"] = round(sum(n for _, n, _ in pair_pcts) / len(pair_pcts), 1)
        # Best (highest concordance) method-pair
        best = max(pair_pcts, key=lambda x: x[0])
        out["max_concordance_methods"] = best[2]
        out["max_concordance_pct"] = round(best[0], 3)
    else:
        out["mean_concordance"] = None
        out["mean_n_clusters"] = None
        out["max_concordance_methods"] = None
        out["max_concordance_pct"] = None

    return out


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
    parser.add_argument("--delta-ddg", type=float, default=1.0,
                        help="kcal/mol threshold for DDG F1/F2/Amb assignment "
                             "(default 1.0; DDG is on a different scale than TM-score)")
    parser.add_argument("--condition-switch", action="store_true",
                        help="Also score condition-aware fold preference "
                             "(Delta_switch) for methods with both apo and "
                             "holo predictions (AF3, Boltz-2). Reads "
                             "df_boltz2_holo.csv / df_af holo rows and writes "
                             "docs/fold_diversity_condition_switch.csv plus a "
                             "per-pair summary. Does not alter the apo survey.")
    parser.add_argument("--min-clusters-for-concordance", type=int, default=10,
                        help="Min mean_n_clusters per pair for inclusion in the "
                             "printed concordance ranking (default 10). Pairs with "
                             "fewer clusters trivially hit 100%% — exclude them. "
                             "Affects only the printed top-25; the CSV always has all pairs.")
    args = parser.parse_args()

    if args.pairs.upper() == "ALL":
        pairs = discover_pairs()
    else:
        pairs = [p.strip() for p in args.pairs.split(",") if p.strip()]

    out_dir = args.output_dir or TABLES_RES
    os.makedirs(out_dir, exist_ok=True)

    print(f"Scanning {len(pairs)} pairs for fold diversity...")
    print(f"Delta threshold (TM): {args.delta}")
    print(f"Delta threshold (DDG, kcal/mol): {args.delta_ddg}")
    print(f"MSAT metric:     {args.msat_metric}")
    print()

    all_cluster_rows = []
    all_summaries = []
    all_concordance = []   # one row per pair (cross-method concordance)
    all_switch_rows = []       # per (pair, cluster, method) Delta_switch rows
    all_switch_summary = []    # per (pair, method) Delta_switch summary

    n_af2 = n_af3 = n_esm = n_msat = n_ccmpred = n_ddg = n_boltz2 = 0

    for pair_id in pairs:
        per_method_rows = {}   # used for cross-method concordance below

        # AF2
        af_rows = load_af2_diversity(pair_id, args.delta)
        if af_rows:
            _apply_centered_classification(af_rows, args.delta)
            n_af2 += 1
            all_cluster_rows.extend(af_rows)
            per_method_rows["AF2"] = af_rows
            s = summarize_pair(pair_id, af_rows)
            if s:
                all_summaries.append(s)

        # AF3 (uses same df_af.csv as AF2, filtered by model=="AF3")
        af3_rows = load_af3_diversity(pair_id, args.delta)
        if af3_rows:
            _apply_centered_classification(af3_rows, args.delta)
            n_af3 += 1
            all_cluster_rows.extend(af3_rows)
            per_method_rows["AF3"] = af3_rows
            s = summarize_pair(pair_id, af3_rows)
            if s:
                all_summaries.append(s)

        # ESM slot = ESMFold2 (v2 single-sequence baseline; replaces ESMFold v1)
        esm_rows = load_esmfold2_diversity(pair_id, args.delta)
        if esm_rows:
            _apply_centered_classification(esm_rows, args.delta)
            n_esm += 1
            all_cluster_rows.extend(esm_rows)
            per_method_rows["ESM"] = esm_rows
            s = summarize_pair(pair_id, esm_rows)
            if s:
                all_summaries.append(s)

        # MSAT
        msat_rows = load_msat_diversity(pair_id, args.delta,
                                        msat_metric=args.msat_metric)
        if msat_rows:
            _apply_centered_classification(msat_rows, args.delta)
            n_msat += 1
            all_cluster_rows.extend(msat_rows)
            per_method_rows["MSAT"] = msat_rows
            s = summarize_pair(pair_id, msat_rows)
            if s:
                all_summaries.append(s)

        # CCMpred
        cc_rows = load_ccmpred_diversity(pair_id, args.delta,
                                         msat_metric=args.msat_metric)
        if cc_rows:
            _apply_centered_classification(cc_rows, args.delta)
            n_ccmpred += 1
            all_cluster_rows.extend(cc_rows)
            per_method_rows["CCMpred"] = cc_rows
            s = summarize_pair(pair_id, cc_rows)
            if s:
                all_summaries.append(s)

        # DDG (ThermoMPNN conformation-bias)
        ddg_rows = load_ddg_diversity(pair_id, args.delta_ddg)
        if ddg_rows:
            _apply_centered_classification(ddg_rows, args.delta_ddg)
            n_ddg += 1
            all_cluster_rows.extend(ddg_rows)
            per_method_rows["DDG"] = ddg_rows
            s = summarize_pair(pair_id, ddg_rows)
            if s:
                all_summaries.append(s)

        # Boltz-2 (apo per-cluster predictions; reads df_boltz2.csv)
        boltz2_rows = load_boltz2_diversity(pair_id, args.delta)
        if boltz2_rows:
            _apply_centered_classification(boltz2_rows, args.delta)
            n_boltz2 += 1
            all_cluster_rows.extend(boltz2_rows)
            per_method_rows["Boltz2"] = boltz2_rows
            s = summarize_pair(pair_id, boltz2_rows)
            if s:
                all_summaries.append(s)

        # BioEmu (per-cluster equilibrium ensemble; reads df_bioemu.csv)
        bioemu_rows = load_bioemu_diversity(pair_id, args.delta)
        if bioemu_rows:
            _apply_centered_classification(bioemu_rows, args.delta)
            all_cluster_rows.extend(bioemu_rows)
            per_method_rows["BioEmu"] = bioemu_rows
            s = summarize_pair(pair_id, bioemu_rows)
            if s:
                all_summaries.append(s)

        # Condition-aware Delta_switch (apo vs holo) for AF3 and Boltz-2
        if args.condition_switch:
            for method, apo_rows, holo_loader in (
                ("AF3", af3_rows, load_af3_holo_diversity),
                ("Boltz2", boltz2_rows, load_boltz2_holo_diversity),
            ):
                if not apo_rows:
                    continue
                holo_rows = holo_loader(pair_id, args.delta)
                if not holo_rows:
                    continue
                sw = compute_condition_switch(pair_id, method, apo_rows,
                                              holo_rows, args.delta)
                if not sw:
                    continue
                all_switch_rows.extend(sw)
                bound = [r for r in sw if r["pref_switch"] == "toward_bound"]
                away = [r for r in sw if r["pref_switch"] == "away_from_bound"]
                vals = [r["dswitch_bound"] for r in sw
                        if r["dswitch_bound"] == r["dswitch_bound"]]
                all_switch_summary.append({
                    "pair_id": pair_id, "method": method,
                    "trigger_class": condition_spec(pair_id).get("trigger_class", ""),
                    "bound_fold": sw[0]["bound_fold"],
                    "n_clusters": len(sw),
                    "n_toward_bound": len(bound),
                    "n_away_from_bound": len(away),
                    "mean_dswitch_bound": round(float(np.mean(vals)), 4) if vals else float("nan"),
                    "frac_toward_bound": round(len(bound) / len(sw), 3) if sw else 0.0,
                })

        # Cross-method concordance for THIS pair
        if len(per_method_rows) >= 2:
            all_concordance.append(compute_pair_concordance(pair_id, per_method_rows))

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

    # --- Save cross-method concordance CSV ---
    df_concord = pd.DataFrame(all_concordance)
    concord_path = os.path.join(out_dir, "fold_diversity_concordance.csv")
    if not df_concord.empty:
        df_concord.to_csv(concord_path, index=False)
        print(f"Saved cross-method concordance: {concord_path} "
              f"({len(df_concord)} rows)")

    # --- Save condition-aware Delta_switch CSVs (opt-in) ---
    if args.condition_switch:
        if all_switch_rows:
            sw_path = os.path.join(out_dir, "fold_diversity_condition_switch.csv")
            pd.DataFrame(all_switch_rows).to_csv(sw_path, index=False)
            sws_path = os.path.join(out_dir,
                                    "fold_diversity_condition_switch_summary.csv")
            pd.DataFrame(all_switch_summary).to_csv(sws_path, index=False)
            n_pairs_sw = len({r["pair_id"] for r in all_switch_rows})
            print(f"Saved condition-switch (Delta_switch): {sw_path} "
                  f"({len(all_switch_rows)} cluster rows, {n_pairs_sw} pairs) "
                  f"+ summary {sws_path}")
        else:
            print("[condition-switch] no pair had both apo and holo predictions "
                  "(need df_boltz2_holo.csv or AF3 holo rows); nothing written.")

    # --- Print analysis ---
    print(f"\n{'='*70}")
    print(f"FOLD DIVERSITY SURVEY RESULTS")
    print(f"{'='*70}")
    print(f"Pairs with AF2 data:     {n_af2}")
    print(f"Pairs with AF3 data:     {n_af3}")
    print(f"Pairs with ESM data:     {n_esm}")
    print(f"Pairs with MSAT data:    {n_msat}")
    print(f"Pairs with CCMpred data: {n_ccmpred}")
    print(f"Pairs with DDG data:     {n_ddg}")
    print(f"Pairs with Boltz2 data:  {n_boltz2}")

    for method in ["AF2", "ESM", "MSAT", "CCMpred", "DDG"]:
        mdf = df_summary[df_summary["method"] == method]
        if mdf.empty:
            continue

        print(f"\n--- {method} ({len(mdf)} pairs) ---")

        # How many pairs show diversity? (raw absolute classification)
        n_diverse = mdf["has_both_folds"].sum()
        print(f"  RAW (absolute):       diverse pairs = {n_diverse}/{len(mdf)} "
              f"({100*n_diverse/len(mdf):.0f}%)")

        # Centered (relative-to-pair-median) diversity
        if "has_both_folds_centered" in mdf.columns:
            n_diverse_c = mdf["has_both_folds_centered"].sum()
            print(f"  CENTERED (relative):  diverse pairs = {n_diverse_c}/{len(mdf)} "
                  f"({100*n_diverse_c/len(mdf):.0f}%)")

        # Preference distribution (raw)
        total_f1 = mdf["n_f1"].sum()
        total_f2 = mdf["n_f2"].sum()
        total_amb = mdf["n_amb"].sum()
        total = total_f1 + total_f2 + total_amb
        print(f"  Raw clusters:      F1={total_f1} F2={total_f2} Amb={total_amb} "
              f"({total} total)")
        if total > 0:
            print(f"  Raw fractions:     F1={100*total_f1/total:.1f}% "
                  f"F2={100*total_f2/total:.1f}% "
                  f"Amb={100*total_amb/total:.1f}%")
        if "n_f1_centered" in mdf.columns:
            tf1c = mdf["n_f1_centered"].sum()
            tf2c = mdf["n_f2_centered"].sum()
            tac = mdf["n_amb_centered"].sum()
            tt = tf1c + tf2c + tac
            if tt > 0:
                print(f"  Centered fractions: F1={100*tf1c/tt:.1f}% "
                      f"F2={100*tf2c/tt:.1f}% "
                      f"Amb={100*tac/tt:.1f}%")

        # Spread statistics
        print(f"  TMdiff spread:  mean={mdf['TMdiff_spread'].mean():.4f}  "
              f"median={mdf['TMdiff_spread'].median():.4f}  "
              f"max={mdf['TMdiff_spread'].max():.4f}")
        print(f"  TMdiff range:   mean={mdf['TMdiff_range'].mean():.4f}  "
              f"median={mdf['TMdiff_range'].median():.4f}  "
              f"max={mdf['TMdiff_range'].max():.4f}")

        # Per-pair breakdown — show raw + centered side by side
        print(f"\n  Per-pair breakdown (RAW | CENTERED):")
        print(f"  {'pair':<20s} {'clusters':>8s} | "
              f"{'F1':>3s} {'F2':>3s} {'Am':>3s} {'div':>4s} | "
              f"{'F1c':>3s} {'F2c':>3s} {'Amc':>3s} {'divc':>4s} | "
              f"{'spread':>8s} {'range':>8s}")
        for _, row in mdf.iterrows():
            d_raw = "YES" if row.get("has_both_folds") else "no"
            d_cen = "YES" if row.get("has_both_folds_centered") else "no"
            star = ""
            if row.get("has_both_folds_centered") and not row.get("has_both_folds"):
                star = " *new*"   # centered only — would have been missed by raw
            elif row.get("has_both_folds"):
                star = " ***"
            print(f"  {row['pair_id']:<20s} {row['n_clusters']:>8d} | "
                  f"{row.get('n_f1', 0):>3d} {row.get('n_f2', 0):>3d} "
                  f"{row.get('n_amb', 0):>3d} {d_raw:>4s} | "
                  f"{row.get('n_f1_centered', 0):>3d} {row.get('n_f2_centered', 0):>3d} "
                  f"{row.get('n_amb_centered', 0):>3d} {d_cen:>4s} | "
                  f"{row['TMdiff_spread']:>8.4f} {row['TMdiff_range']:>8.4f}"
                  f"{star}")

    # --- Cross-method comparison ---
    if len(df_summary["method"].unique()) > 1:
        print(f"\n{'='*70}")
        print("CROSS-METHOD COMPARISON")
        print(f"{'='*70}")

        # Build a clean pair x method table
        pairs_seen = df_summary["pair_id"].unique()
        print(f"\n  {'pair':<20s}", end="")
        for m in ["AF2", "ESM", "MSAT", "CCMpred", "DDG"]:
            print(f"  {m:>12s}", end="")
        print()

        for pid in sorted(pairs_seen):
            print(f"  {pid:<20s}", end="")
            for m in ["AF2", "ESM", "MSAT", "CCMpred", "DDG"]:
                row = df_summary[(df_summary["pair_id"] == pid) &
                                 (df_summary["method"] == m)]
                if row.empty:
                    print(f"  {'--':>12s}", end="")
                else:
                    r = row.iloc[0]
                    label = f"{r['n_f1']}F1/{r['n_f2']}F2/{r['n_amb']}A"
                    print(f"  {label:>12s}", end="")
            print()

    # --- Cross-method concordance (per-cluster sign agreement) ---
    if not df_concord.empty:
        print(f"\n{'='*70}")
        print("CROSS-METHOD CONCORDANCE — per-cluster TMdiff_centered sign agreement")
        print(f"{'='*70}")
        print("Independent of diversity-call thresholds: even if no single method")
        print("crosses its own delta, sign-agreement across orthogonal methods")
        print("(different training sets / inductive biases) is hard to dismiss")
        print("as method-specific bias. ~50% = chance; >70% = significant.\n")

        MIN_CLUSTERS = args.min_clusters_for_concordance

        # Filter pairs by min cluster count for reliable concordance estimates.
        # Tiny pairs (1-4 clusters) hit 100% concordance trivially — exclude them.
        cd_all = df_concord.dropna(subset=["mean_concordance"])
        cd = cd_all[cd_all["mean_n_clusters"] >= MIN_CLUSTERS].copy()
        cd = cd.sort_values("mean_concordance", ascending=False)

        n_filtered_out = len(cd_all) - len(cd)
        if n_filtered_out > 0:
            print(f"  ({n_filtered_out} pairs excluded with mean_n_clusters < "
                  f"{MIN_CLUSTERS} — too few clusters for reliable concordance)")

        # Top 25
        print(f"  Top 25 pairs by mean cross-method concordance "
              f"(mean_n_clusters >= {MIN_CLUSTERS}):")
        print(f"  {'pair':<18s} {'mean':>5s} {'A~E':>5s} {'A~M':>5s} {'A~D':>5s} "
              f"{'E~M':>5s} {'E~D':>5s} {'M~D':>5s}  best_pair (pct, n)")
        for _, r in cd.head(25).iterrows():
            def _fmt(p):
                v = r.get(f"{p}_pct")
                return f"{v:.2f}" if v is not None and not pd.isna(v) else " -- "
            best_lab = r.get("max_concordance_methods") or "-"
            best_pct = r.get("max_concordance_pct")
            best_pct_s = f"{best_pct:.2f}" if best_pct is not None and not pd.isna(best_pct) else "--"
            print(f"  {r['pair_id']:<18s} {r['mean_concordance']:>5.2f} "
                  f"{_fmt('AF2_ESM'):>5s} {_fmt('AF2_MSAT'):>5s} {_fmt('AF2_DDG'):>5s} "
                  f"{_fmt('ESM_MSAT'):>5s} {_fmt('ESM_DDG'):>5s} {_fmt('MSAT_DDG'):>5s}  "
                  f"{best_lab} ({best_pct_s}, n={int(r.get('mean_n_clusters') or 0)})")

        n_above = int((cd["mean_concordance"] > 0.65).sum())
        n_total = len(cd)
        print(f"\n  Pairs with mean concordance > 0.65: {n_above} / {n_total} "
              f"({100*n_above/max(n_total,1):.0f}%) — strong cross-method signal.")
        print("  See docs/fold_diversity_concordance.csv for the full table.")

    print(f"\n{'='*70}")
    print("CONCLUSION")
    print(f"{'='*70}")
    print("Three complementary signals per pair:")
    print("  1. Diversity counts (raw / centered): does any method cross its threshold?")
    print("  2. Spread/range of TMdiff: how much do clusters vary in fold preference?")
    print("  3. Cross-method concordance: do orthogonal methods agree per-cluster on")
    print("     direction, even if no single method crosses its threshold?")
    print()
    print("Strongest fold-switching candidates satisfy 2+ of these signals.")
    print("Cross-method concordance is least sensitive to method-specific")
    print("training bias and threshold choice.")
    print()


if __name__ == "__main__":
    main()
