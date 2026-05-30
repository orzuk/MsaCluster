#!/usr/bin/env python3
"""porter_93_distribution.py -- show the joint TM x SS-switch x seq-id
distribution for the existing 93 Porter pairs.

Reads docs/summary_final_res_all_pairs_df.csv (per config.py) which
already has TM-score for each Porter pair. Merges in seq_id from
docs/triggers_from_pdb.csv if available. Reports distributions and
identifies the right thresholds for re-mining.

Usage:
  python3 scripts/porter_93_distribution.py
  python3 scripts/porter_93_distribution.py --summary <other.csv>
"""

import argparse
import sys
from pathlib import Path
import pandas as pd


def find_col(df, candidates):
    """Find the first column whose name (case-insensitive substring) matches
    any candidate. Returns column name or None."""
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        cand_low = cand.lower()
        for low, orig in cols_lower.items():
            if cand_low in low:
                return orig
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default="docs/summary_final_res_all_pairs_df.csv")
    ap.add_argument("--triggers", default="docs/triggers_from_pdb.csv")
    args = ap.parse_args()

    if not Path(args.summary).is_file():
        print(f"ERROR: not found: {args.summary}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(args.summary)
    print(f"Loaded {len(df)} rows from {args.summary}")
    print(f"Columns ({len(df.columns)}): {list(df.columns)[:25]}...")
    print()

    # Try to auto-detect the relevant columns
    pair_col = find_col(df, ["pair_id", "pair", "name"])
    # PAIR_TM (the actual column name in summary_final_res_all_pairs_df.csv)
    # comes first so we don't accidentally pick AF2Clust_TM1 etc.
    tm_col = find_col(df, ["pair_tm", "tm_score_truth_to_truth", "tm_truth",
                            "tm_f1_f2", "tm_pair", "tm_score_pair", "pdb_tm"])
    if tm_col is None:
        # Fallback: any column with "tm" that looks like the pair score
        tm_candidates = [c for c in df.columns if "tm" in c.lower() and "score" in c.lower()]
        if tm_candidates:
            tm_col = tm_candidates[0]
    ss_col = find_col(df, ["ss_switch", "ss_diff", "ss_change", "ss_identity"])
    # If SS column not in summary, try the detailed CSV
    if ss_col is None:
        det_path = Path(args.summary).parent / "detailed_final_res_all_pairs_df.csv"
        if det_path.is_file():
            det = pd.read_csv(det_path)
            det_ss = find_col(det, ["ss_switch", "ss_diff", "ss_change", "ss_identity"])
            det_pair = find_col(det, ["pair_id", "pair", "name"])
            if det_ss and det_pair:
                # If multiple rows per pair, take mean
                agg = det.groupby(det_pair)[det_ss].mean().reset_index()
                df = df.merge(agg, left_on=pair_col, right_on=det_pair, how="left")
                ss_col = det_ss
                print(f"Merged in SS-switch column '{det_ss}' from {det_path.name}")
                print()
    # Also try the existing triggers CSV
    if ss_col is None:
        trig_path = Path(args.triggers)
        if trig_path.is_file():
            trig = pd.read_csv(trig_path)
            trig_ss = find_col(trig, ["ss_switch", "ss_diff", "ss_change", "ss_identity"])
            if trig_ss:
                trig_pair = find_col(trig, ["pair_id", "pair", "name"])
                if trig_pair:
                    df = df.merge(trig[[trig_pair, trig_ss]], left_on=pair_col, right_on=trig_pair, how="left")
                    ss_col = trig_ss
                    print(f"Merged in SS-switch column '{trig_ss}' from triggers CSV")

    print(f"Auto-detected columns:")
    print(f"  pair_id col:    {pair_col}")
    print(f"  TM-score col:   {tm_col}")
    print(f"  SS-switch col:  {ss_col}")
    print()

    # Merge in seq-id from triggers_from_pdb.csv if available
    seq_id_col = None
    if Path(args.triggers).is_file():
        trig = pd.read_csv(args.triggers)
        sid_col_in_trig = find_col(trig, ["seq_id", "sequence_identity"])
        pair_in_trig = find_col(trig, ["pair_id", "pair", "name"])
        if sid_col_in_trig and pair_in_trig:
            df = df.merge(
                trig[[pair_in_trig, sid_col_in_trig]],
                left_on=pair_col, right_on=pair_in_trig, how="left"
            )
            seq_id_col = sid_col_in_trig
            print(f"Merged in seq-id column '{sid_col_in_trig}' from {args.triggers}")
            print()

    print("=" * 65)
    print("PORTER-93 STRUCTURAL DISTRIBUTION")
    print("=" * 65)

    if tm_col:
        tm = pd.to_numeric(df[tm_col], errors="coerce").dropna()
        # Normalize to 0-1 if it's in 0-100
        if tm.max() > 1.5:
            tm = tm / 100.0
            print(f"  (TM was 0-100 scale; rescaled to 0-1)")
        print(f"\nTM-score (n={len(tm)}):")
        print(f"  min={tm.min():.3f}, 5%={tm.quantile(0.05):.3f}, "
              f"25%={tm.quantile(0.25):.3f}, median={tm.median():.3f}, "
              f"75%={tm.quantile(0.75):.3f}, 95%={tm.quantile(0.95):.3f}, "
              f"max={tm.max():.3f}")
        for thr in [0.3, 0.5, 0.7, 0.85, 0.95]:
            covered = (tm <= thr).sum()
            print(f"  TM <= {thr}: {covered} / {len(tm)} pairs ({covered/len(tm)*100:.0f}%)")

    if ss_col:
        ss = pd.to_numeric(df[ss_col], errors="coerce").dropna()
        if ss.max() > 1.5:
            ss = ss / 100.0
            print(f"  (SS was 0-100 scale; rescaled to 0-1)")
        print(f"\nSS-switch fraction (n={len(ss)}):")
        print(f"  min={ss.min():.3f}, 5%={ss.quantile(0.05):.3f}, "
              f"25%={ss.quantile(0.25):.3f}, median={ss.median():.3f}, "
              f"75%={ss.quantile(0.75):.3f}, 95%={ss.quantile(0.95):.3f}, "
              f"max={ss.max():.3f}")
        for thr in [0.05, 0.10, 0.15, 0.20, 0.30]:
            covered = (ss >= thr).sum()
            print(f"  ss_switch >= {thr}: {covered} / {len(ss)} pairs ({covered/len(ss)*100:.0f}%)")

    if seq_id_col:
        sid = pd.to_numeric(df[seq_id_col], errors="coerce").dropna()
        if sid.max() > 1.5:
            sid = sid / 100.0
            print(f"  (seq-id was 0-100 scale; rescaled to 0-1)")
        print(f"\nseq-id (n={len(sid)}):")
        print(f"  min={sid.min():.3f}, 5%={sid.quantile(0.05):.3f}, "
              f"25%={sid.quantile(0.25):.3f}, median={sid.median():.3f}, "
              f"75%={sid.quantile(0.75):.3f}, 95%={sid.quantile(0.95):.3f}, "
              f"max={sid.max():.3f}")
        for thr in [0.50, 0.70, 0.85, 0.95]:
            covered = (sid >= thr).sum()
            print(f"  seq_id >= {thr}: {covered} / {len(sid)} pairs ({covered/len(sid)*100:.0f}%)")

    # The key joint question: what filter recovers >= 90% of Porter?
    if tm_col and ss_col:
        print("\n" + "=" * 65)
        print("JOINT FILTER COVERAGE (would mining find this Porter pair?)")
        print("=" * 65)
        tmv = pd.to_numeric(df[tm_col], errors="coerce")
        ssv = pd.to_numeric(df[ss_col], errors="coerce")
        if tmv.max() > 1.5:  tmv = tmv / 100.0
        if ssv.max() > 1.5:  ssv = ssv / 100.0
        joint = pd.DataFrame({"tm": tmv, "ss": ssv}).dropna()
        scenarios = [
            ("current mining (tm in [0.3, 0.5] AND ss>=0.10)",
                (joint.tm >= 0.3) & (joint.tm <= 0.5) & (joint.ss >= 0.10)),
            ("drop TM upper bound:  tm>=0.3 AND ss>=0.10",
                (joint.tm >= 0.3) & (joint.ss >= 0.10)),
            ("no TM filter:         ss>=0.10",
                (joint.ss >= 0.10)),
            ("no TM filter:         ss>=0.15",
                (joint.ss >= 0.15)),
            ("no TM filter:         ss>=0.20",
                (joint.ss >= 0.20)),
            ("tighten ss only:      ss>=0.25",
                (joint.ss >= 0.25)),
            ("exclude duplicates:   tm<0.99 AND ss>=0.10",
                (joint.tm < 0.99) & (joint.ss >= 0.10)),
        ]
        for name, mask in scenarios:
            n = int(mask.sum())
            print(f"  {name:55s} : {n} / {len(joint)} ({n/len(joint)*100:.0f}%)")

    # Show the extreme cases
    if tm_col and ss_col and pair_col:
        print("\nHigh-TM Porter pairs (regional fold-switchers):")
        tmv = pd.to_numeric(df[tm_col], errors="coerce")
        if tmv.max() > 1.5:  tmv = tmv / 100.0
        ssv = pd.to_numeric(df[ss_col], errors="coerce")
        if ssv.max() > 1.5:  ssv = ssv / 100.0
        sub = df.assign(_tm=tmv, _ss=ssv).dropna(subset=["_tm", "_ss"])
        top = sub.nlargest(10, "_tm")[[pair_col, "_tm", "_ss"]]
        top.columns = [pair_col, "TM", "SS_switch"]
        print(top.to_string(index=False))
        print("\nLow-TM Porter pairs (KaiB-style fold-switchers):")
        bot = sub.nsmallest(10, "_tm")[[pair_col, "_tm", "_ss"]]
        bot.columns = [pair_col, "TM", "SS_switch"]
        print(bot.to_string(index=False))


if __name__ == "__main__":
    main()
