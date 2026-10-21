#!/usr/bin/env python3
"""Ancestral gain/loss test of the 4-state fold capability on the coarse trees.

Per pair (default mode):
  tree    ``utils.ancestral_utils.build_coarse_tree`` (default representative
          strategy ``largest_piece_medoid``) + ``_name_internal_nodes`` --
          the same coarse tree as scripts/phylo_placement.py.
  labels  4-state F1 / F2 / both / none per cluster, exactly as the per-pair
          tree figures (utils/protein_plot_utils.py gain/loss markers):
          ``assign_fold_preference_4state(load_consensus_cluster_tm(pair),
          tau=0.5)`` -- consensus max-TM to each fold over AF2 / AF3 /
          ESMFold2; both = TM1 >= tau and TM2 >= tau, none = neither. Leaves
          without TM data are "none".
  ASR     Fitch parsimony (``utils.ancestral_tests.fitch_reconstruct``; same
          algorithm and tie-break as the Fitch path of
          ``ancestral_utils.run_asr``). Transitions = edges whose
          reconstructed parent/child states differ ("none" is a state).
          Events: gain = into "both", loss = out of "both", swap = F1<->F2,
          other = into/out of "none".
  test    testable = >= 2 distinct informative states among F1/F2/both.
          Label shuffle over leaves (--n-perm 1000, seed 42):
          p = (1 + #{perm transitions <= observed}) / (n_perm + 1);
          BH q over testable pairs.
  trait   Sequence-trait control: K_d = --n-trait (50) random sequence traits
  control (random per-residue score projections; tables k = 1..K_d of
          scripts/negative_traits_phylo.random_scores, rp_seed 777 -- the
          first K_d traits of scripts/trait_null_phylo.py) are each
          discretised by rank into the SAME states with EXACTLY the same
          state counts as the observed fold labels (top n_1 values -> most
          frequent state, ...; ``ancestral_tests.discretise_like``) and
          Fitch-counted. p_trait_anc = (1 + #{trait transitions <= observed
          fold transitions}) / (K_d + 1): small = fold state changes more
          rarely along the tree than generic sequence features do.
  timeout --timeout seconds per pair (default $ANCESTRAL_PAIR_TIMEOUT or 120).

Outputs
-------
docs/ancestral_summary.csv (one row per pair; overwritten). Columns:
  pair_id, state_mode (abs4), tau, rep_strategy, recon_method (fitch),
  n_leaves, n_F1, n_F2, n_both, n_none, n_other_state, informative_states,
  testable, root_state, pattern (ancestral_utils.classify_pattern),
  n_transitions, n_gain, n_loss, n_swap, n_other,
  p_value, q_value (BH over testable pairs), n_perm, perm_null_median,
  K_trait, n_trait_leaves, p_trait_anc, trait_trans_median, trait_trans_min,
  frac_trait_le_obs, secs, error,
  legacy aliases: n_clusters (= n_leaves), n_f1, n_f2, n_amb (= 0),
  n_fold_transitions (= n_swap), method (= recon_method).
docs/ancestral_events.csv (one row per changed edge): pair_id, parent,
  child, from_state, to_state, event, child_is_leaf, n_leaves_below.

--legacy runs the previous per-pair ``ancestral_utils.analyze_pair``
(pastml if available, figures, old summary columns).

Usage:
  python scripts/run_ancestral_analysis.py --pairs ALL --n-jobs 8
  python scripts/run_ancestral_analysis.py --pairs 2qqjA_4qdsA,1fzpD_2frhA
  python scripts/run_ancestral_analysis.py --pairs ALL --summary-only   # re-summarise CSV
  python scripts/run_ancestral_analysis.py --pairs ALL --legacy --delta 0.05
"""

import os
# Headless rendering env vars MUST be set BEFORE importing matplotlib / ete3.
# Without these, on a headless compute node ete3.TreeStyle hangs waiting on Qt
# and matplotlib fails to find a writable config dir.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache-orzuk")
os.environ.setdefault("PYTHONUNBUFFERED", "1")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import argparse
import multiprocessing as mp
import signal
import sys
import time
import zlib
from contextlib import contextmanager

import numpy as np
import pandas as pd


@contextmanager
def _per_pair_timeout(seconds: int):
    """SIGALRM-based per-pair timeout. Skips a single hung pair without
    killing the whole run. Linux-only; uses SIGALRM."""
    def _handler(signum, frame):
        raise TimeoutError(f"per-pair timeout after {seconds}s")
    old = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)

# Allow running from repo root
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "scripts"))

print(f"[ancestral] importing modules...", flush=True)
_t0 = time.time()
from config import DATA_DIR, PAIR_DIR_RE, TABLES_RES
from utils.ancestral_tests import ancestral_core, bh_adjust, EVENT_TYPES
print(f"[ancestral] imports done in {time.time()-_t0:.1f}s", flush=True)


def discover_eligible_pairs() -> list:
    """Find pairs that have both a tree file and df_af.csv."""
    eligible = []
    if not os.path.isdir(DATA_DIR):
        print(f"[error] DATA_DIR not found: {DATA_DIR}")
        return eligible

    for name in sorted(os.listdir(DATA_DIR)):
        if not PAIR_DIR_RE.match(name):
            continue
        pair_dir = os.path.join(DATA_DIR, name)
        tree_file = os.path.join(pair_dir, "output_phytree", "DeepMsa_tree.nwk")
        af_csv = os.path.join(pair_dir, "Analysis", "df_af.csv")
        cluster_dir = os.path.join(pair_dir, "output_msa_cluster")
        if (os.path.isfile(tree_file) and os.path.isfile(af_csv)
                and os.path.isdir(cluster_dir)):
            eligible.append(name)
    return eligible


# ---------------------------------------------------------------------------
# 4-state gain/loss analysis (default)
# ---------------------------------------------------------------------------

def analyze_pair_4state(job):
    """Worker: one pair -> (summary row, event rows). Never raises."""
    (pair_id, tau, n_perm, n_trait, rp_seed, timeout, render) = job
    t0 = time.time()
    row = {"pair_id": pair_id, "state_mode": "abs4", "tau": tau,
           "recon_method": "fitch", "error": ""}
    try:
        with _per_pair_timeout(timeout):
            from utils.ancestral_utils import (
                build_coarse_tree, _name_internal_nodes, load_consensus_cluster_tm,
                assign_fold_preference_4state, classify_pattern, DEFAULT_REP_STRATEGY)
            row["rep_strategy"] = DEFAULT_REP_STRATEGY
            tree, tag_to_label, _ = build_coarse_tree(pair_id)
            _name_internal_nodes(tree)
            prefs = assign_fold_preference_4state(load_consensus_cluster_tm(pair_id), tau)
            leaf_states = {tag_to_label[t]: s for t, s in prefs.items() if t in tag_to_label}
            traits = None
            if n_trait > 0:
                from negative_traits_phylo import random_score_matrix
                from trait_null_phylo import cluster_random_traits
                R = random_score_matrix(n_trait, rp_seed)
                rt = cluster_random_traits(os.path.join(DATA_DIR, pair_id),
                                           list(tag_to_label), R)
                traits = {tag_to_label[t]: v for t, v in rt.items()}
            summ, events, states = ancestral_core(
                tree, leaf_states, traits, n_perm=n_perm, seed=42,
                trait_seed=zlib.crc32(pair_id.encode()), default="none")
            node_states = {nd.name: s for nd, s in zip(tree.traverse("preorder"), states)}
            row["pattern"] = classify_pattern(tree, node_states)
            row.update(summ)
            if render:
                from utils.ancestral_utils import render_ancestral_tree
                fig = os.path.join(DATA_DIR, pair_id, "output_figs",
                                   f"{pair_id}_ancestral_tree.png")
                try:
                    render_ancestral_tree(tree, node_states, pair_id, fig, None)
                except Exception as e:  # noqa: BLE001
                    print(f"  [warn] render failed for {pair_id}: {e}", flush=True)
        ev_rows = [{"pair_id": pair_id, **e} for e in events]
    except Exception as e:  # noqa: BLE001  (incl. TimeoutError)
        row["error"] = f"{type(e).__name__}: {e}"
        ev_rows = []
    row["secs"] = round(time.time() - t0, 1)
    return row, ev_rows


SUMMARY_COLS = [
    "pair_id", "state_mode", "tau", "rep_strategy", "recon_method",
    "n_leaves", "n_F1", "n_F2", "n_both", "n_none", "n_other_state",
    "informative_states", "testable", "root_state", "pattern",
    "n_transitions"] + [f"n_{e}" for e in EVENT_TYPES] + [
    "p_value", "q_value", "n_perm", "perm_null_median",
    "K_trait", "n_trait_leaves", "p_trait_anc", "trait_trans_median",
    "trait_trans_min", "frac_trait_le_obs", "secs", "error",
    "n_clusters", "n_f1", "n_f2", "n_amb", "n_fold_transitions", "method"]


def finalize_summary(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for c in SUMMARY_COLS:
        if c not in df.columns:
            df[c] = np.nan
    ok = df["error"].fillna("") == ""
    test = ok & (df["testable"] == True)  # noqa: E712
    df["q_value"] = np.nan
    df.loc[test, "q_value"] = bh_adjust(df.loc[test, "p_value"].values)
    df["n_clusters"] = df["n_leaves"]
    df["n_f1"], df["n_f2"] = df["n_F1"], df["n_F2"]
    df.loc[ok, "n_amb"] = 0
    df["n_fold_transitions"] = df["n_swap"]
    df["method"] = df["recon_method"]
    return df[SUMMARY_COLS]


def print_summary(df: pd.DataFrame, alpha: float = 0.05) -> None:
    ok = df[df["error"].fillna("") == ""]
    test = ok[ok["testable"] == True]  # noqa: E712
    print(f"\n{'='*60}\nAncestral 4-state gain/loss summary")
    print(f"  pairs attempted         {len(df)}")
    print(f"  reconstructed           {len(ok)}   (errors/timeouts: {len(df) - len(ok)})")
    print(f"  testable (>=2 of F1/F2/both) {len(test)}")
    print(f"  label-shuffle p <= {alpha}  {int((test.p_value <= alpha).sum())}")
    print(f"  label-shuffle q <= {alpha}  {int((test.q_value <= alpha).sum())}")
    pt_ = test["p_trait_anc"]
    print(f"  p_trait_anc <= {alpha}     {int((pt_ <= alpha).sum())}"
          f"   (of {int(np.isfinite(pt_.astype(float)).sum())} with traits)")
    both = test[(test.p_value <= alpha) & (test.p_trait_anc <= alpha)]
    print(f"  both p <= {alpha}          {len(both)}"
          + (": " + ", ".join(both.pair_id) if len(both) else ""))
    tot = {e: int(ok[f"n_{e}"].fillna(0).sum()) for e in EVENT_TYPES}
    print(f"  total events (reconstructed pairs): gains {tot['gain']}, losses "
          f"{tot['loss']}, swaps {tot['swap']}, other {tot['other']} "
          f"(transitions {int(ok['n_transitions'].fillna(0).sum())})")
    if len(ok):
        print("\nPattern distribution:")
        for pat, c in ok["pattern"].value_counts().items():
            print(f"  {pat:25s} {c:3d}")
        print("Root state distribution:")
        for st, c in ok["root_state"].value_counts().items():
            print(f"  {st:5s} {c:3d}")


def run_4state(args, pairs):
    timeout = args.timeout
    jobs = [(p, args.tau, args.n_perm, args.n_trait, args.rp_seed, timeout,
             args.render_figures) for p in pairs]
    rows, events = [], []
    t0 = time.time()

    def _log(k, r):
        msg = (r["error"] if r["error"] else
               f"{r['n_leaves']} leaves F1={r['n_F1']} F2={r['n_F2']} "
               f"both={r['n_both']} none={r['n_none']} | trans={r['n_transitions']} "
               f"(g{r['n_gain']}/l{r['n_loss']}/s{r['n_swap']}/o{r['n_other']})"
               + (f" | p={r['p_value']:.3f} p_trait={r['p_trait_anc']:.3f}"
                  if r.get("testable") else " | untestable"))
        print(f"[{k}/{len(jobs)}] {r['pair_id']}: {msg} ({r['secs']}s)", flush=True)

    if args.n_jobs > 1:
        with mp.Pool(args.n_jobs) as pool:
            for k, (r, ev) in enumerate(pool.imap_unordered(analyze_pair_4state, jobs,
                                                            chunksize=1), 1):
                rows.append(r)
                events += ev
                _log(k, r)
    else:
        for k, j in enumerate(jobs, 1):
            r, ev = analyze_pair_4state(j)
            rows.append(r)
            events += ev
            _log(k, r)
    print(f"\nProcessed {len(pairs)} pairs in {time.time()-t0:.1f}s")
    df = finalize_summary(pd.DataFrame(rows).sort_values("pair_id"))
    ev = pd.DataFrame(events, columns=["pair_id", "parent", "child", "from_state",
                                       "to_state", "event", "child_is_leaf",
                                       "n_leaves_below"])
    return df, ev.sort_values(["pair_id", "event", "child"]) if len(ev) else ev


# ---------------------------------------------------------------------------
# Legacy path (previous behaviour)
# ---------------------------------------------------------------------------

def run_legacy(args, pairs):
    from utils.ancestral_utils import analyze_pair
    results = []
    for i, pair_id in enumerate(pairs):
        t_pair = time.time()
        print(f"\n[{i+1}/{len(pairs)}] {pair_id} ... starting "
              f"(timeout {args.timeout}s)", flush=True)
        try:
            with _per_pair_timeout(args.timeout):
                r = analyze_pair(pair_id, delta=args.delta, n_perm=args.n_perm)
            results.append(r)
            print(f"  {r['n_clusters']} clusters | pattern={r['pattern']} | "
                  f"root={r['root_state']} | transitions={r['n_fold_transitions']}"
                  + (f" | p={r['p_value']:.3f}" if r['p_value'] is not None else "")
                  + f"  ({time.time()-t_pair:.1f}s)", flush=True)
        except TimeoutError:
            print(f"  TIMEOUT after {args.timeout}s; skipped", flush=True)
            results.append({"pair_id": pair_id, "error": f"timeout after {args.timeout}s"})
        except Exception as e:
            print(f"  [FAILED] {e}  ({time.time()-t_pair:.1f}s)", flush=True)
            results.append({"pair_id": pair_id, "error": str(e)})
    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(
        description="Ancestral 4-state fold-capability gain/loss test")
    parser.add_argument("--pairs", type=str, default="ALL",
                        help="Comma-separated pair IDs, or ALL")
    parser.add_argument("--delta", type=float, default=0.05,
                        help="(--legacy only) F1/F2/Amb threshold on the AF2 TM diff")
    parser.add_argument("--tau", type=float, default=0.5,
                        help="absolute TM threshold for the 4-state labels")
    parser.add_argument("--n-perm", type=int, default=1000,
                        help="Number of label permutations")
    parser.add_argument("--n-trait", type=int, default=50,
                        help="K_d random sequence traits for the trait control (0 = off)")
    parser.add_argument("--rp-seed", type=int, default=777)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--timeout", type=int,
                        default=int(os.environ.get("ANCESTRAL_PAIR_TIMEOUT", "120")),
                        help="per-pair timeout in seconds")
    parser.add_argument("--render-figures", action="store_true",
                        help="also render <pair>_ancestral_tree.png (Fitch states)")
    parser.add_argument("--summary-only", action="store_true",
                        help="skip computation; re-summarise an existing --output-csv")
    parser.add_argument("--legacy", action="store_true",
                        help="previous analyze_pair pipeline (pastml/figures/old columns)")
    parser.add_argument("--output-csv", type=str, default=None,
                        help="Path for summary CSV (default: docs/ancestral_summary.csv)")
    parser.add_argument("--events-csv", type=str, default=None,
                        help="Per-branch events CSV (default: docs/ancestral_events.csv)")
    args = parser.parse_args()

    csv_path = args.output_csv or os.path.join(TABLES_RES, "ancestral_summary.csv")
    ev_path = args.events_csv or os.path.join(TABLES_RES, "ancestral_events.csv")
    if args.summary_only:
        df = pd.read_csv(csv_path)
        if "testable" in df.columns:
            print_summary(df)
        else:
            print(df.describe(include="all").T.to_string())
        return

    if args.pairs.upper() == "ALL":
        pairs = discover_eligible_pairs()
        print(f"Found {len(pairs)} eligible pairs (have tree + df_af.csv)")
    else:
        pairs = [p.strip() for p in args.pairs.split(",") if p.strip()]
    if not pairs:
        print("[error] No pairs to process.")
        return

    os.makedirs(os.path.dirname(os.path.abspath(csv_path)), exist_ok=True)
    if args.legacy:
        df = run_legacy(args, pairs)
        df.to_csv(csv_path, index=False)
        print(f"\nSaved summary: {csv_path}")
        return

    df, ev = run_4state(args, pairs)
    df.to_csv(csv_path, index=False)
    ev.to_csv(ev_path, index=False)
    print(f"\nSaved summary: {csv_path}\nSaved events:  {ev_path}")
    print_summary(df)


if __name__ == "__main__":
    main()
