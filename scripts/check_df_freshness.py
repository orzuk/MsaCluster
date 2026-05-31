#!/usr/bin/env python3
"""Flag pairs whose per-pair Analysis CSVs are STALE vs the fine clustering.

The survey reads Pipeline/FoldPairs/<pair>/Analysis/df_*.csv. If a df has
fewer distinct clusters than the fine clustering (output_msa_cluster), it
wasn't regenerated after re-clustering (the KaiB bug). Lists the stale pairs
so you can re-run postprocess only on those:

    python3 Analysis/postprocess_unified.py --pairs <stale list> --force_rerun

Usage:
    python3 scripts/check_df_freshness.py            # checks df_af (+df_ddg,df_cmap)
"""
import glob
import os
import re

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAIRS = os.path.join(ROOT, "Pipeline", "FoldPairs")
# per-pair CSVs that carry per-cluster rows (fine-resolution methods)
CHECK = {"df_af.csv": "AF2/AF3", "df_ddg.csv": "DDG"}


def n_clusters(csv):
    if not os.path.isfile(csv):
        return None
    try:
        d = pd.read_csv(csv)
    except Exception:
        return -1
    ccol = next((c for c in ("cluster", "cluster_num") if c in d.columns), None)
    if ccol is None:
        return -2
    cl = d[ccol].astype(str)
    cl = cl[~cl.str.contains("Deep", case=False, na=False)]
    # normalize to the ShallowMsa index so '000' and 'ShallowMsa_000' collapse
    idx = cl.str.extract(r"(\d+)", expand=False)
    return idx.dropna().nunique()


def main():
    stale = []
    print(f"{'pair':16}{'fine':>6}{'df_af':>7}{'df_ddg':>8}  flag")
    for d in sorted(glob.glob(os.path.join(PAIRS, "*", ""))):
        pair = os.path.basename(os.path.normpath(d))
        if pair in ("jobs", "_s4pred_work"):
            continue
        fine = len(glob.glob(os.path.join(d, "output_msa_cluster", "ShallowMsa_*.a3m")))
        if fine == 0:
            continue
        naf = n_clusters(os.path.join(d, "Analysis", "df_af.csv"))
        nddg = n_clusters(os.path.join(d, "Analysis", "df_ddg.csv"))
        flag = ""
        # stale if df_af present but has notably fewer clusters than fine
        if naf is not None and naf >= 0 and naf < 0.9 * fine:
            flag = "STALE df_af"
            stale.append(pair)
        print(f"{pair:16}{fine:>6}{str(naf):>7}{str(nddg):>8}  {flag}")
    print(f"\n{len(stale)} pair(s) with STALE df_af:")
    print(" ".join(stale) if stale else "(none)")
    if stale:
        print("\nFix:\n  python3 Analysis/postprocess_unified.py --pairs "
              + " ".join(stale) + " --force_rerun")


if __name__ == "__main__":
    main()
