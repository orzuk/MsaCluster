#!/usr/bin/env python3
"""Negative-trait baseline for phylogenetic fold-preference signal (referee #6).

Question: is the Moran's I of per-cluster fold preference larger than what
OTHER per-cluster, sequence-derived quantities show on the same coarse
cluster trees?  Sequence-derived traits are expected to carry phylogenetic
signal (the tree IS built from the sequences), so they give an empirical
"how much autocorrelation does an arbitrary sequence property have here"
reference for the fold-preference numbers in docs/phylo_placement*.csv.

Per cluster (Pipeline/FoldPairs/<pair>/output_msa_cluster/<tag>.a3m; row 0 is
the query/medoid row and is skipped; members = rows 1..):
  pid_to_query   mean gap-aware identity of members to the pair's query
                 (DeepMsa.a3m row 0 as in scripts/seq_divergence_correction.py;
                 falls back to the cluster a3m row 0), over query match columns
                 (``gapaware_pid``).  This is the training-proximity confound
                 that the paper's correction regresses out -> strong
                 autocorrelation expected (sanity check).
  hydrophobicity mean over members of the mean Kyte-Doolittle score over
                 aligned (upper-case) residues
  gap_frac       mean over members of the fraction of match columns that are
                 '-' (a3m: '-' = deletion at a match column; lower case =
                 insertion, not counted)
  log_size       log(number of member sequences)
  net_charge     mean over members of (#K + #R - #D - #E) over aligned residues
  random_proj_k  mean over members of the mean of a FIXED random per-residue
                 score (20 amino acids ~ N(0,1), seeded; k = 1..--n-random)
                 -- a generic sequence-composition control

Each trait -> Moran's I on the pair's coarse tree via
``utils.phylo_tests.run_test('morans', ...)`` with the paper's settings
(n_perm = 1000, rng = default_rng(seed + 7), seed 12345), so I / I_norm / p
are directly comparable to docs/phylo_placement*.csv.

Outputs
-------
--out            per (pair, trait): n, I, I_norm, z, p           (long CSV)
--out-values     per (pair, cluster) trait values                (for re-use)
--out-summary    per group (trait or fold-pref method x label set): #families,
                 fraction p < 0.05, median I, median I_norm
--out-families   per family: best fold-pref I (max over methods, corrected
                 and centered labels), median fold-pref I, each trait's I,
                 and whether best fold-pref I > hydrophobicity / random I
--fig            boxplots of I and I_norm by trait / method
                 (docs/figures/global/negative_traits_morans.png)

Usage (repo root)
-----------------
python3 scripts/negative_traits_phylo.py --n-jobs 8
python3 scripts/negative_traits_phylo.py --max-pairs 3          # quick test
python3 scripts/negative_traits_phylo.py --summary-only         # re-summarise --out
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache-" + os.environ.get("USER", "x"))
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import argparse
import multiprocessing as mp
import signal
import sys
import time
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "scripts"))

from utils.phylo_tests import run_test  # noqa: E402
from seq_divergence_correction import gapaware_pid, load_a3m_seqs  # noqa: E402

AA20 = "ACDEFGHIKLMNPQRSTVWY"
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
CHARGE = {"K": 1.0, "R": 1.0, "D": -1.0, "E": -1.0}
BASE_TRAITS = ("pid_to_query", "hydrophobicity", "gap_frac", "log_size", "net_charge")
FP_METHODS = ("AF2", "AF3", "Boltz2", "BioEmu", "ESM", "MSAT", "CCMpred", "DDG", "S4PRED")


# ---------------------------------------------------------------------------
# Per-cluster trait computation
# ---------------------------------------------------------------------------

def random_scores(k: int, seed: int) -> Dict[str, float]:
    """Fixed random per-residue score table #k (N(0,1) over the 20 aa)."""
    rs = np.random.default_rng([seed, k])
    return dict(zip(AA20, rs.standard_normal(len(AA20))))


def _mean_score(aligned: str, table: Dict[str, float]) -> float:
    v = [table[c] for c in aligned if c in table]
    return float(np.mean(v)) if v else np.nan


def _match_cols(s: str) -> str:
    """Drop a3m insertions (lower case) so rows are column-aligned to the
    query's match columns before the position-wise ``gapaware_pid``.
    No-op for a3m files already stripped of insertions (AF-Cluster output)."""
    return "".join(c for c in s if not c.islower())


def cluster_traits(seqs: Sequence[str], query: Optional[str],
                   rp_tables: Sequence[Dict[str, float]]) -> Optional[dict]:
    """Traits of one cluster from its a3m rows (row 0 skipped).

    ``query`` is the reference for identity (None -> row 0). Returns None
    when the cluster has no member rows.
    """
    if len(seqs) <= 1:
        return None
    members = seqs[1:]
    ref = _match_cols(query if query else seqs[0])
    pid, hyd, gap, chg = [], [], [], []
    rps = [[] for _ in rp_tables]
    for s in members:
        aligned = "".join(c for c in s if c.isupper())          # match-col residues
        n_match = len(aligned) + s.count("-")
        p = gapaware_pid(ref, _match_cols(s))
        if np.isfinite(p):
            pid.append(p)
        hyd.append(_mean_score(aligned, KYTE_DOOLITTLE))
        gap.append(s.count("-") / n_match if n_match else np.nan)
        chg.append(sum(CHARGE.get(c, 0.0) for c in aligned))
        for k, tab in enumerate(rp_tables):
            rps[k].append(_mean_score(aligned, tab))
    nm = lambda v: float(np.nanmean(v)) if np.isfinite(v).any() else np.nan  # noqa: E731
    out = {
        "n_members": len(members),
        "pid_to_query": nm(np.asarray(pid, float)) if pid else np.nan,
        "hydrophobicity": nm(np.asarray(hyd, float)),
        "gap_frac": nm(np.asarray(gap, float)),
        "log_size": float(np.log(len(members))),
        "net_charge": nm(np.asarray(chg, float)),
    }
    for k, v in enumerate(rps, 1):
        out[f"random_proj_{k}"] = nm(np.asarray(v, float))
    return out


def trait_names(n_random: int) -> List[str]:
    return list(BASE_TRAITS) + [f"random_proj_{k}" for k in range(1, n_random + 1)]


# ---------------------------------------------------------------------------
# Per-pair analysis
# ---------------------------------------------------------------------------

def morans_rows(pair_id: str, tree, label_values: Dict[str, Dict[str, float]],
                traits: Sequence[str], n_perm: int, seed: int) -> List[dict]:
    """Moran's I (paper settings) for each trait. ``label_values`` is
    {leaf_label: {trait: value}}."""
    rows = []
    for t in traits:
        vals = {lab: d.get(t) for lab, d in label_values.items()}
        r = run_test("morans", tree, vals, n_perm=min(n_perm, 1000),
                     rng=np.random.default_rng(seed + 7))
        rows.append({"pair_id": pair_id, "trait": t, "n": r["morans_n"],
                     "morans_I": r["morans_I"], "morans_I_norm": r["morans_I_norm"],
                     "morans_z": r["morans_z"], "morans_p": r["morans_p"],
                     "morans_note": r["morans_note"]})
    return rows


def _alarm(signum, frame):
    raise TimeoutError("timed out")


def analyze_pair(job) -> dict:
    """Worker: coarse tree + traits + Moran's for one pair."""
    pair_id, n_random, rp_seed, n_perm, seed, timeout = job
    from config import DATA_DIR
    try:
        from utils.ancestral_utils import build_coarse_tree, _name_internal_nodes
        signal.signal(signal.SIGALRM, _alarm)
        signal.alarm(timeout)
        try:
            tree, tag_to_label, _ = build_coarse_tree(pair_id)
            _name_internal_nodes(tree)
        finally:
            signal.alarm(0)
    except Exception as e:  # noqa: BLE001
        return {"pair_id": pair_id, "error": f"{type(e).__name__}: {e}"}
    pair_dir = os.path.join(DATA_DIR, pair_id)
    deep = os.path.join(pair_dir, "output_get_msa", "DeepMsa.a3m")
    query = None
    if os.path.exists(deep):
        ds = load_a3m_seqs(deep)
        query = ds[0] if ds else None
    rp_tables = [random_scores(k, rp_seed) for k in range(1, n_random + 1)]
    values, label_values = [], {}
    for tag, lab in tag_to_label.items():
        f = os.path.join(pair_dir, "output_msa_cluster", f"{tag}.a3m")
        if not os.path.exists(f):
            continue
        tr = cluster_traits(load_a3m_seqs(f), query, rp_tables)
        if tr is None:
            continue
        label_values[lab] = tr
        values.append({"pair_id": pair_id, "cluster": tag, "label": lab, **tr})
    if len(label_values) < 3:
        return {"pair_id": pair_id, "error": f"only {len(label_values)} clusters with traits"}
    rows = morans_rows(pair_id, tree, label_values, trait_names(n_random), n_perm, seed)
    n_leaves = sum(1 for _ in tree.iter_leaves())
    for r in rows:
        r["n_leaves"] = n_leaves
    return {"pair_id": pair_id, "rows": rows, "values": values,
            "query_source": "DeepMsa" if query else "cluster_row0"}


# ---------------------------------------------------------------------------
# Comparison with fold-preference Moran's I
# ---------------------------------------------------------------------------

def load_fold_pref(paths: Dict[str, str], pairs: Optional[set]) -> pd.DataFrame:
    parts = []
    for label_set, p in paths.items():
        if not os.path.exists(p):
            print(f"[warn] {p} missing; skipping {label_set} fold-pref rows")
            continue
        d = pd.read_csv(p, usecols=lambda c: c in (
            "pair_id", "method", "morans_I", "morans_I_norm", "morans_p"))
        d["label_set"] = label_set
        parts.append(d)
    if not parts:
        return pd.DataFrame(columns=["pair_id", "method", "morans_I",
                                     "morans_I_norm", "morans_p", "label_set"])
    fp = pd.concat(parts, ignore_index=True).dropna(subset=["morans_I"])
    if pairs is not None:
        fp = fp[fp.pair_id.isin(pairs)]
    return fp


def summarise(tr: pd.DataFrame, fp: pd.DataFrame, alpha: float = 0.05):
    """(group summary, per-family comparison)."""
    tr = tr.dropna(subset=["morans_I"])
    long = pd.concat([
        tr.assign(group=tr.trait, kind="trait")[
            ["pair_id", "group", "kind", "morans_I", "morans_I_norm", "morans_p"]],
        fp.assign(group="FP_" + fp.label_set + ":" + fp.method, kind="fold_pref")[
            ["pair_id", "group", "kind", "morans_I", "morans_I_norm", "morans_p"]],
    ], ignore_index=True)
    summ = long.groupby(["kind", "group"]).agg(
        n_families=("pair_id", "nunique"),
        frac_sig=("morans_p", lambda p: float((p < alpha).mean())),
        median_I=("morans_I", "median"), mean_I=("morans_I", "mean"),
        median_I_norm=("morans_I_norm", "median"),
    ).reset_index().sort_values(["kind", "median_I"], ascending=[False, False])

    wide = tr.pivot_table(index="pair_id", columns="trait", values="morans_I")
    wide.columns = [f"I_{c}" for c in wide.columns]
    fam = wide
    for ls, sub in fp.groupby("label_set"):
        g = sub.groupby("pair_id")
        best = g.morans_I.max().rename(f"best_fp_I_{ls}")
        arg = sub.loc[g.morans_I.idxmax(), ["pair_id", "method"]].set_index(
            "pair_id").method.rename(f"best_fp_method_{ls}")
        med = g.morans_I.median().rename(f"median_fp_I_{ls}")
        nsig = g.morans_p.apply(lambda p: int((p < alpha).sum())).rename(f"n_fp_sig_{ls}")
        fam = fam.join([best, arg, med, nsig], how="left")
        for t in [c for c in wide.columns if c in ("I_hydrophobicity", "I_random_proj_1",
                                                   "I_pid_to_query", "I_net_charge")]:
            b = fam[f"best_fp_I_{ls}"]
            fam[f"best_fp_{ls}_gt_{t[2:]}"] = np.where(
                b.notna() & fam[t].notna(), (b > fam[t]).astype(float), np.nan)
    return summ, fam.reset_index()


def print_report(summ: pd.DataFrame, fam: pd.DataFrame) -> None:
    with pd.option_context("display.width", 180, "display.max_rows", 100):
        print("\n=== Moran's I by trait / fold-preference method ===")
        print(summ.round(4).to_string(index=False))
        cmp_cols = [c for c in fam.columns if "_gt_" in c]
        if cmp_cols:
            print("\n=== fraction of families where the BEST fold-preference I "
                  "exceeds the trait I (best = max over methods; biased upward) ===")
            for c in cmp_cols:
                v = fam[c].dropna()
                print(f"  {c:45s} {v.mean():.3f}  (n={len(v)})")


def plot(tr: pd.DataFrame, fp: pd.DataFrame, fig_path: str, label_set: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tr = tr.dropna(subset=["morans_I"])
    fps = fp[fp.label_set == label_set] if label_set in set(fp.label_set) else fp
    traits = list(tr.groupby("trait").morans_I.median().sort_values(ascending=False).index)
    meths = [m for m in FP_METHODS if m in set(fps.method)]
    groups = [("trait", t) for t in traits] + [("fp", m) for m in meths]
    c_trait, c_fp = "#52514e", "#2a78d6"
    fig, axes = plt.subplots(2, 1, figsize=(max(9, 0.7 * len(groups) + 2), 8.5), sharex=True)
    for ax, col, ylab in ((axes[0], "morans_I", "Moran's I"),
                          (axes[1], "morans_I_norm", "Moran's I_norm")):
        data, colors = [], []
        for kind, g in groups:
            if kind == "trait":
                v = tr.loc[tr.trait == g, col].dropna().values
                colors.append(c_trait)
            else:
                v = fps.loc[fps.method == g, col].dropna().values
                colors.append(c_fp)
            data.append(v)
        bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=True,
                        flierprops=dict(marker="o", ms=3, alpha=0.5, mec="none"),
                        medianprops=dict(color="white", lw=2))
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_edgecolor(c)
        for k, c in enumerate(colors):
            for part in ("whiskers", "caps"):
                for ln in bp[part][2 * k: 2 * k + 2]:
                    ln.set_color(c)
            bp["fliers"][k].set_markerfacecolor(c)
        ax.axhline(0, color="#8a8984", lw=0.8, ls="--")
        ax.set_ylabel(ylab)
        ax.grid(axis="y", color="#e4e3df", lw=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        # fraction significant above each box
        if col == "morans_I":
            top = ax.get_ylim()[1]
            for k, (kind, g) in enumerate(groups, 1):
                sub = tr[tr.trait == g] if kind == "trait" else fps[fps.method == g]
                ax.text(k, top, f"{(sub.morans_p < 0.05).mean():.0%}", ha="center",
                        va="bottom", fontsize=7, color="#52514e")
    axes[1].set_xticks(range(1, len(groups) + 1))
    axes[1].set_xticklabels([g for _, g in groups], rotation=45, ha="right")
    from matplotlib.patches import Patch
    axes[0].legend(handles=[Patch(color=c_trait, label="sequence-derived trait"),
                            Patch(color=c_fp, label=f"fold preference ({label_set})")],
                   frameon=False, loc="upper right", fontsize=8)
    axes[0].set_title("Phylogenetic autocorrelation (Moran's I) of per-cluster traits on "
                      f"the coarse cluster trees ({tr.pair_id.nunique()} families; "
                      "% = families with p < 0.05)", fontsize=10, pad=14)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(fig_path)), exist_ok=True)
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)
    print(f"[fig] {fig_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def discover_pairs() -> List[str]:
    from config import DATA_DIR, PAIR_DIR_RE
    if not os.path.isdir(DATA_DIR):
        return []
    return sorted(d for d in os.listdir(DATA_DIR)
                  if PAIR_DIR_RE.match(d) and os.path.isdir(os.path.join(DATA_DIR, d)))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pairs", default="ALL")
    ap.add_argument("--max-pairs", type=int, default=0)
    ap.add_argument("--dry-run", action="store_true",
                    help="list pairs and exit (use --max-pairs for a small real run)")
    ap.add_argument("--n-random", type=int, default=3, help="# random-projection traits")
    ap.add_argument("--rp-seed", type=int, default=777)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--n-jobs", type=int, default=1)
    ap.add_argument("--tree-timeout", type=int, default=300)
    ap.add_argument("--fp-corrected", default=os.path.join(REPO, "docs", "phylo_placement_corrected.csv"))
    ap.add_argument("--fp-centered", default=os.path.join(REPO, "docs", "phylo_placement.csv"))
    ap.add_argument("--plot-label-set", choices=["corrected", "centered"], default="corrected")
    ap.add_argument("--out", default=os.path.join(REPO, "docs", "negative_traits_morans.csv"))
    ap.add_argument("--out-values", default=os.path.join(REPO, "docs", "negative_traits_values.csv"))
    ap.add_argument("--out-summary", default=os.path.join(REPO, "docs", "negative_traits_summary.csv"))
    ap.add_argument("--out-families", default=os.path.join(REPO, "docs", "negative_traits_families.csv"))
    ap.add_argument("--fig", default=os.path.join(REPO, "docs", "figures", "global",
                                                   "negative_traits_morans.png"))
    ap.add_argument("--summary-only", action="store_true",
                    help="skip computation; summarise/plot an existing --out")
    args = ap.parse_args(argv)
    t0 = time.time()

    if args.summary_only:
        tr = pd.read_csv(args.out)
    else:
        pairs = (discover_pairs() if args.pairs.upper() == "ALL"
                 else [p.strip() for p in args.pairs.split(",") if p.strip()])
        if args.max_pairs:
            pairs = pairs[: args.max_pairs]
        print(f"[pairs] {len(pairs)} pairs; traits: {trait_names(args.n_random)}")
        if args.dry_run:
            print("\n".join(pairs))
            return
        jobs = [(p, args.n_random, args.rp_seed, args.n_perm, args.seed, args.tree_timeout)
                for p in pairs]
        if args.n_jobs > 1:
            with mp.Pool(args.n_jobs) as pool:
                res = list(pool.imap_unordered(analyze_pair, jobs, chunksize=1))
        else:
            res = [analyze_pair(j) for j in jobs]
        rows, vals = [], []
        for r in sorted(res, key=lambda d: d["pair_id"]):
            if "error" in r:
                print(f"  [skip] {r['pair_id']}: {r['error']}")
                continue
            rows += r["rows"]
            vals += r["values"]
            print(f"  {r['pair_id']}: {len(r['values'])} clusters (query={r['query_source']})")
        if not rows:
            sys.exit("[error] no pairs analysed")
        tr = pd.DataFrame(rows)
        for p in (args.out, args.out_values):
            os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
        tr.to_csv(args.out, index=False)
        pd.DataFrame(vals).to_csv(args.out_values, index=False)
        print(f"[out] {args.out}\n[out] {args.out_values}")

    fp = load_fold_pref({"corrected": args.fp_corrected, "centered": args.fp_centered},
                        set(tr.pair_id))
    summ, fam = summarise(tr, fp)
    summ.to_csv(args.out_summary, index=False)
    fam.to_csv(args.out_families, index=False)
    print(f"[out] {args.out_summary}\n[out] {args.out_families}")
    print_report(summ, fam)
    plot(tr, fp, args.fig, args.plot_label_set)
    print(f"[done] {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
