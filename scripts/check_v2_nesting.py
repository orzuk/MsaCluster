#!/usr/bin/env python3
"""Check whether coarse clusters are unions of fine clusters across all
viable pairs.

If every fine cluster is fully contained in exactly one coarse cluster
(frac_overlap = 1.0 for all rows), then the independent tree-cut
approach happens to produce nested clusterings -- equivalent in spirit
to bottom-up agglomeration from fine clusters. We can then stick with
the simpler tree-cut implementation.

If some fine clusters span multiple coarse clusters (frac_overlap < 1.0),
the merging step has broken nesting and we should switch to building
coarse FROM fine via UPGMA agglomeration.

For each viable pair this script:
  1. Reads output_msa_cluster/ShallowMsa_*.a3m (fine memberships)
  2. Reads output_msa_cluster_coarse/CoarseMsa_*.a3m (coarse memberships)
  3. Computes per-fine-cluster frac_overlap with its best-matching coarse
  4. Reports global nesting quality

Usage:
    python3 scripts/check_v2_nesting.py
"""
from __future__ import annotations

import sys
from pathlib import Path


def read_headers(a3m_path: Path) -> set[str]:
    ids = set()
    with a3m_path.open() as f:
        for line in f:
            if line.startswith(">"):
                tok = line[1:].split()[0]
                if tok:
                    ids.add(tok)
    return ids


def check_pair(pair_dir: Path) -> dict:
    """Return per-pair nesting summary."""
    fine_dir = pair_dir / "output_msa_cluster"
    coarse_dir = pair_dir / "output_msa_cluster_coarse"
    if not fine_dir.is_dir() or not coarse_dir.is_dir():
        return {"pair": pair_dir.name, "status": "missing_dir"}

    fine_files = sorted(fine_dir.glob("ShallowMsa_*.a3m"))
    coarse_files = sorted(coarse_dir.glob("CoarseMsa_*.a3m"))
    if not fine_files or not coarse_files:
        return {"pair": pair_dir.name, "status": "missing_files",
                "n_fine": len(fine_files), "n_coarse": len(coarse_files)}

    coarse_members = {c.stem: read_headers(c) for c in coarse_files}

    fracs = []
    impure_fines = 0
    for f in fine_files:
        f_ids = read_headers(f)
        if not f_ids:
            continue
        best_overlap = max(len(f_ids & c_ids) for c_ids in coarse_members.values())
        frac = best_overlap / len(f_ids)
        fracs.append(frac)
        if frac < 0.999:
            impure_fines += 1

    return {
        "pair": pair_dir.name,
        "status": "ok",
        "n_fine": len(fine_files),
        "n_coarse": len(coarse_files),
        "n_impure_fines": impure_fines,
        "frac_clean":  sum(1 for f in fracs if f >= 0.999) / len(fracs),
        "min_frac":    min(fracs),
        "median_frac": sorted(fracs)[len(fracs) // 2],
    }


def main():
    viable_path = Path("docs/v2_viable_pairs.txt")
    if viable_path.is_file():
        pairs = [p.strip() for p in viable_path.read_text().splitlines() if p.strip()]
        print(f"Checking {len(pairs)} pairs from docs/v2_viable_pairs.txt")
    else:
        pairs = sorted(p.name for p in Path("Pipeline/FoldPairs").glob("*/")
                       if p.name not in ("_s4pred_work", "jobs"))
        print(f"docs/v2_viable_pairs.txt not found; checking all "
              f"{len(pairs)} pair dirs")

    results = []
    perfect = 0
    broken = []
    for pid in pairs:
        r = check_pair(Path("Pipeline/FoldPairs") / pid)
        results.append(r)
        if r["status"] != "ok":
            continue
        if r["n_impure_fines"] == 0:
            perfect += 1
        else:
            broken.append((pid, r["n_impure_fines"], r["n_fine"],
                           r["min_frac"]))

    ok_count = sum(1 for r in results if r["status"] == "ok")
    print()
    print(f"=== Nesting check across {ok_count} pairs ===")
    print(f"Perfectly nested (all fine clusters fully in one coarse): "
          f"{perfect}/{ok_count}")
    print(f"Has at least one impure fine cluster: {len(broken)}/{ok_count}")
    if broken:
        broken.sort(key=lambda x: x[3])  # by worst min_frac
        print()
        print("Pairs with broken nesting (lowest frac first):")
        print(f"  {'pair':<22s}  impure/total  min_frac")
        for pid, imp, tot, mf in broken[:20]:
            print(f"  {pid:<22s}  {imp:3d}/{tot:<3d}      {mf:.3f}")
        if len(broken) > 20:
            print(f"  ... and {len(broken)-20} more")

    # Verdict
    print()
    if perfect == ok_count:
        print("VERDICT: nesting is perfect across all pairs -- the "
              "independent tree-cut approach is equivalent to bottom-up "
              "agglomeration. Keep as-is.")
    elif perfect / ok_count >= 0.95:
        print(f"VERDICT: {100*perfect/ok_count:.0f}% pairs are perfectly "
              f"nested; the rest have minor merging artifacts. Acceptable "
              f"in practice -- document the imperfect cases in the SI.")
    else:
        print(f"VERDICT: only {100*perfect/ok_count:.0f}% pairs are "
              f"perfectly nested -- nesting is materially broken. Switch "
              f"to bottom-up agglomeration.")


if __name__ == "__main__":
    main()
