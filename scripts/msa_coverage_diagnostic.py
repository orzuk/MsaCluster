#!/usr/bin/env python3
"""MSA coverage diagnostic for the four candidate pairs.

For XCL1, SOD1, RelB, Selecase report:
  - total Deep MSA seq count
  - n_clusters
  - per-cluster (n_seqs, mean PID-to-query, median PID-to-query)

Interpretation: PID-to-query tells us how diverged each cluster's seqs are
from the query. For XCL1 the literature shift (Dishman 2021) sits at the
C-chemokine / CC-CXC chemokine split, i.e. roughly PID ~25-35% to query.
If all clusters are PID>40% (close C-chemokine orthologs), our MSA does
not actually cross Dishman's branch and we shouldn't expect to see the
clade-level shift in the per-cluster tree.

Reads docs/per_cluster_pid_to_query.csv (pre-computed).
"""

from __future__ import annotations

import os
import sys

import pandas as pd

DOCS = "docs"
PAIRS = {
    "2n54B_2hdmA": "XCL1",
    "5jytA_2qkeE": "KaiB",
    "6c6sD_2ougC": "RfaH",
    "4qhhA_4qhfA": "Selecase",
    "2namA_1uxmK": "SOD1",
    "1zk9A_3jv6A": "RelB",
    "1nrjB_2gedB": "SR-beta",
}


def main() -> None:
    csv = os.path.join(DOCS, "per_cluster_pid_to_query.csv")
    if not os.path.isfile(csv):
        print(f"[error] missing {csv}", file=sys.stderr)
        sys.exit(2)
    df = pd.read_csv(csv)
    print(f"Loaded {len(df)} rows from {csv}")
    print()

    for pid, name in PAIRS.items():
        sub = df[df["pair_id"] == pid].copy()
        if sub.empty:
            print(f"=== {name} ({pid}): no rows ===")
            continue
        sub = sub.sort_values("cluster")
        total_seqs = int(sub["n_seqs"].sum())
        n_clusters = len(sub)
        mean_pid_all = float(sub["mean_pid_to_query"].mean())
        min_pid = float(sub["mean_pid_to_query"].min())
        max_pid = float(sub["mean_pid_to_query"].max())
        n_close = int((sub["mean_pid_to_query"] >= 0.40).sum())
        n_mid = int(((sub["mean_pid_to_query"] >= 0.25)
                     & (sub["mean_pid_to_query"] < 0.40)).sum())
        n_far = int((sub["mean_pid_to_query"] < 0.25).sum())

        print(f"=== {name} ({pid}) ===")
        print(f"  total seqs across clusters: {total_seqs}")
        print(f"  n_clusters:                 {n_clusters}")
        print(f"  cluster mean-PID range:     [{min_pid:.3f}, {max_pid:.3f}]")
        print(f"  mean cluster mean-PID:      {mean_pid_all:.3f}")
        print(f"  clusters by PID band:")
        print(f"     close (>=0.40, likely orthologs/close paralogs): {n_close}")
        print(f"     mid   (0.25 - 0.40, distant homologs/sister clade): {n_mid}")
        print(f"     far   (<0.25, very distant / questionable): {n_far}")
        print()
        # detailed per-cluster lines
        for _, r in sub.iterrows():
            cl = r["cluster"]
            ns = int(r["n_seqs"])
            mp = float(r["mean_pid_to_query"])
            md = float(r["median_pid_to_query"])
            print(f"    {cl:<22} n={ns:>5d}  mean_PID={mp:.3f}  median_PID={md:.3f}")
        print()


if __name__ == "__main__":
    main()
