#!/usr/bin/env python3
"""False-positive calibration of the divergence-corrected concordance test.

The concordance test asks whether seven predictors agree on which fold each
cluster of a family leans toward more often than chance. Its null shuffles
each method's values among the clusters it covers, so it asks whether the
methods share *any* per-cluster signal. Agreement driven by a shared sequence
property rather than by fold would also break that null, and the reported
33 of 87 families would then be an overcount.

This replaces every method's per-cluster values with a fold-free pseudo-method
and runs the identical chain: divergence regression on identity to both
references, residual signs, the same permutation null, the same BH
correction, the same n_clusters >= 3 filter. Each pseudo-method is a random
linear projection of the cluster's amino-acid composition, so the seven
pseudo-methods read the same alignment through different lenses and share
exactly the kind of non-structural signal the test must not mistake for fold
preference. Nothing in them can know either fold. The number of families they
flag is the false-positive rate.

Usage:
  concordance_fp_calibration.py COMP_CSV OUT_CSV [--reps 20] [--n-perm 1000]
"""
from __future__ import annotations

import argparse
import os
import sys
from itertools import combinations
from multiprocessing import Pool

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "scripts"))
from permutation_null_concordance import (          # noqa: E402
    _mean_concordance_from_residuals, bh_adjust)

SURVEY = os.path.join(REPO, "docs", "fold_diversity_survey_corrected_9method.csv")
PID = os.path.join(REPO, "docs", "per_cluster_pid_symmetric.csv")
NO_CORRECT = False
FAM_DIR = ""
METHODS = ["AF2", "AF3", "ESM", "Boltz2", "BioEmu", "DDG", "S4PRED"]
MIN_CLUSTERS = 3
AA20 = "ACDEFGHIKLMNPQRSTVWY"
FEATS = (["gap_frac", "log_size", "net_charge", "hydrophobicity"]
         + [f"f_{a}" for a in AA20])


def regress_out(y: np.ndarray, x1: np.ndarray, x2: np.ndarray) -> np.ndarray:
    """The divergence correction of seq_divergence_correction.py, verbatim in
    behaviour: two-predictor least squares on identity to both references,
    single-predictor fallback, and no correction when it cannot be fit."""
    ok = np.isfinite(y) & np.isfinite(x1)
    if ok.sum() < 3:
        return y
    yv, x1v, x2v = y[ok], x1[ok], x2[ok]
    # np.std of a vector holding any NaN is NaN, and NaN >= 1e-9 is False, so
    # a partially missing second reference falls through to the single
    # predictor branch. seq_divergence_correction.py behaves the same way.
    has_f2 = not np.isnan(x2v).all()
    if has_f2 and np.std(x2v) >= 1e-9 and np.std(x1v) >= 1e-9:
        X = np.column_stack([np.ones(len(yv)), x1v, x2v])
        try:
            a, b1, b2 = np.linalg.lstsq(X, yv, rcond=None)[0]
        except Exception:
            return y
        return y - (a + b1 * x1 + b2 * x2)
    if np.std(x1v) < 1e-9:
        return y
    b1, a = np.polyfit(x1v, yv, 1)
    return y - (b1 * x1 + a)


def build_pair_tables():
    """Per pair: cluster list, the (cluster x method) presence mask of the real
    data, the standardized composition matrix, and the two identity vectors."""
    sv = pd.read_csv(SURVEY, usecols=["pair_id", "cluster_norm", "method",
                                      "TMdiff_residual"])
    sv = sv[(sv.cluster_norm != "DeepMsa") & sv.method.isin(METHODS)]
    sv = sv.dropna(subset=["TMdiff_residual"])
    comp = pd.read_csv(COMP_CSV)
    comp["cluster_norm"] = comp["cluster"]
    pid = pd.read_csv(PID)
    pid["cluster_norm"] = pid["cluster"]

    tables = []
    for pair, sub in sv.groupby("pair_id"):
        pres = (sub.assign(v=1)
                   .pivot_table(index="cluster_norm", columns="method",
                                values="v", aggfunc="first"))
        present = [m for m in METHODS if m in pres.columns]
        if len(present) < 2:
            continue
        mask = np.isfinite(pres[present].to_numpy(dtype=float))
        clusters = list(pres.index)

        c = comp[comp.pair_id == pair].set_index("cluster_norm")
        c = c.reindex(clusters)
        Z = c[FEATS].to_numpy(dtype=float)
        if not np.isfinite(Z).all():
            col = np.nanmean(Z, axis=0)
            idx = np.where(~np.isfinite(Z))
            if len(idx[0]):
                Z[idx] = np.take(col, idx[1])
        if not np.isfinite(Z).all():
            continue
        sd = Z.std(axis=0)
        keep = sd > 1e-12
        if keep.sum() < 2:
            continue
        Z = (Z[:, keep] - Z[:, keep].mean(axis=0)) / sd[keep]

        p = pid[pid.pair_id == pair].set_index("cluster_norm").reindex(clusters)
        x1 = p["mean_pid_to_F1"].to_numpy(dtype=float)
        x2 = p["mean_pid_to_F2"].to_numpy(dtype=float)
        tables.append(dict(pair=pair, present=present, mask=mask, Z=Z,
                           x1=x1, x2=x2))
    return tables


def one_replicate(rep: int) -> dict:
    rng_global = np.random.default_rng([SEED, rep])
    rows = []
    for t in TABLES:
        Z, mask, present = t["Z"], t["mask"], t["present"]
        n_clu, n_m = mask.shape
        vals = np.full((n_clu, n_m), np.nan)
        for j in range(n_m):
            rows_j = np.flatnonzero(mask[:, j])
            if len(rows_j) == 0:
                continue
            w = rng_global.standard_normal(Z.shape[1])
            v = Z[rows_j] @ w
            # Centering: the survey subtracts the per-pair median before the
            # regression, so the pseudo-method gets the same treatment.
            v = v - np.median(v)
            vals[rows_j, j] = (v if NO_CORRECT else
                               regress_out(v, t["x1"][rows_j], t["x2"][rows_j]))

        signs = np.nan_to_num(np.sign(vals), nan=0.0)
        obs, obs_n = _mean_concordance_from_residuals(signs)
        if np.isnan(obs):
            continue
        n_clusters = int(mask.sum(axis=0).max())

        rng = np.random.default_rng([SEED, rep, zlib_hash(t["pair"])])
        col_rows = [np.flatnonzero(mask[:, j]) for j in range(n_m)]
        n_ge = 0
        perm = signs.copy()
        for _ in range(N_PERM):
            for j, rj in enumerate(col_rows):
                perm[rj, j] = rng.permutation(signs[rj, j])
            pm, _ = _mean_concordance_from_residuals(perm)
            if not np.isnan(pm) and pm >= obs:
                n_ge += 1
        rows.append(dict(pair_id=t["pair"], n_clusters=n_clusters,
                         observed=obs,
                         p=(n_ge + 1) / (N_PERM + 1)))
    out = pd.DataFrame(rows)
    out["p_bh"] = bh_adjust(out["p"].to_numpy())
    out["rep"] = rep
    ev = out[out.n_clusters >= MIN_CLUSTERS]
    if FAM_DIR:
        ev.to_csv(os.path.join(FAM_DIR, f"fam_rep{rep:03d}.csv"), index=False)
    return dict(rep=rep, n_eval=len(ev), hits=int((ev.p_bh < 0.05).sum()),
                hits_raw=int((ev.p < 0.05).sum()),
                median_obs=float(ev.observed.median()))


def zlib_hash(s: str) -> int:
    import zlib
    return zlib.crc32(s.encode()) % (2 ** 31)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("comp_csv")
    ap.add_argument("out_csv")
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261024)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--fam-dir", default="",
                    help="directory for the per-family rows of every replicate")
    ap.add_argument("--no-correct", action="store_true",
                    help="skip the divergence regression, to show how "
                         "much of the false-positive rate it removes")
    a = ap.parse_args()

    COMP_CSV, N_PERM, SEED = a.comp_csv, a.n_perm, a.seed
    NO_CORRECT = a.no_correct
    FAM_DIR = a.fam_dir
    if FAM_DIR:
        os.makedirs(FAM_DIR, exist_ok=True)
    TABLES = build_pair_tables()
    print(f"pairs with >=2 methods and composition: {len(TABLES)}", flush=True)

    with Pool(a.jobs) as pool:
        res = pool.map(one_replicate, range(a.reps))
    df = pd.DataFrame(res).sort_values("rep")
    df.to_csv(a.out_csv, index=False)
    print(df.to_string(index=False))
    print(f"\nfold-free hits out of ~{int(df.n_eval.median())} families: "
          f"mean {df.hits.mean():.1f}  median {df.hits.median():.0f}  "
          f"range {df.hits.min()}-{df.hits.max()}")
    print(f"observed in real data: 33")
