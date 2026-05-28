#!/usr/bin/env python3
"""Build coarse clusters by walking up the phylogenetic tree from fine clusters.

This GUARANTEES coarse clusters are unions of fine clusters: we simply
re-cut the same tree at a higher threshold. The merge order is determined
by the tree's branch lengths (Newick branch lengths from DeepMsa_tree.nwk),
not by any re-computed sequence distance.

Algorithm (per pair):
  1. Load output_phytree/DeepMsa_tree.nwk (Bio.Phylo).
  2. For each fine cluster (output_msa_cluster/ShallowMsa_NNN.a3m),
     find the MRCA of its leaves in the tree.
  3. For each pair of fine clusters, the height at which they MERGE
     in the tree = root-distance of LCA(MRCA_i, MRCA_j).
  4. Sort all K_fine choose 2 merge events by height (ascending).
  5. Process events in order via union-find: union two clusters
     unless already merged. Stop when active group count = K_target.
  6. K_target = clip(N_DeepMsa / 200, K_min=2, K_max=30).
  7. Each resulting group = a coarse cluster = union of member
     sequences from its constituent fine clusters.
  8. Write CoarseMsa_NNN.a3m + fine_to_coarse.csv mapping.

Speed: K_fine <= 100, so K_fine*(K_fine-1)/2 <= 4950 merge events per
pair. Each LCA lookup uses cached root-paths -> O(tree_depth) per call.
Total: a few seconds per pair, ~5 minutes for 82 pairs.

Usage:
    python3 scripts/build_coarse_from_fine.py --foldpair_ids ALL
    python3 scripts/build_coarse_from_fine.py --pair 5jytA_2qkeE
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

try:
    from Bio import Phylo
except ImportError:
    sys.exit("Biopython is required: pip install biopython")


def read_a3m_lines(a3m_path: Path) -> list[tuple[str, str]]:
    """Return [(header_first_token, raw_alignment), ...] preserving inserts."""
    out = []
    name, chunks = None, []
    with a3m_path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            out.append((name, "".join(chunks)))
    return out


def count_deep_msa_seqs(deep_path: Path) -> int:
    if not deep_path.is_file():
        return 0
    n = 0
    with deep_path.open() as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def build_coarse_for_pair(pair_dir: Path,
                           k_min: int = 2,
                           k_max: int = 30,
                           seqs_per_cluster: int = 200) -> dict:
    fine_dir = pair_dir / "output_msa_cluster"
    coarse_dir = pair_dir / "output_msa_cluster_coarse"
    tree_path = pair_dir / "output_phytree" / "DeepMsa_tree.nwk"
    deep_path = pair_dir / "output_get_msa" / "DeepMsa.a3m"

    fine_files = sorted(fine_dir.glob("ShallowMsa_*.a3m"))
    if not fine_files:
        return {"status": "no_fine_clusters", "n_fine": 0}
    if not tree_path.is_file():
        return {"status": "no_tree"}

    # ---- target K_coarse ----
    n_deep = count_deep_msa_seqs(deep_path)
    k_target = max(k_min, min(k_max, n_deep // seqs_per_cluster))

    n_fine = len(fine_files)
    fine_id_to_idx = {f.stem: i for i, f in enumerate(fine_files)}
    # Cache fine cluster members (set of leaf names)
    fine_members = []
    for f in fine_files:
        members = set()
        with f.open() as fh:
            for line in fh:
                if line.startswith(">"):
                    tok = line[1:].split()[0]
                    if tok:
                        members.add(tok)
        fine_members.append(members)

    # If we already have <= K_target fine clusters, no merging needed
    if n_fine <= k_target:
        return _emit_1to1(pair_dir, fine_files, fine_members, coarse_dir,
                          fine_dir, k_target, n_deep)

    # ---- parse tree, find MRCA of each fine cluster ----
    tree = Phylo.read(str(tree_path), "newick")
    leaf_by_name = {l.name: l for l in tree.get_terminals() if l.name}

    # Precompute depth (root-distance) of every clade for fast LCA-height
    # via cached root paths.
    parent_of: dict[int, object] = {}
    def _walk_with_parent(clade, par=None):
        parent_of[id(clade)] = par
        for c in clade.clades:
            _walk_with_parent(c, clade)
    _walk_with_parent(tree.root)

    def root_path(clade):
        """Return list of clade nodes from root to clade, inclusive."""
        path = []
        cur = clade
        while cur is not None:
            path.append(cur)
            cur = parent_of.get(id(cur))
        path.reverse()
        return path

    # Depth of each clade = sum of branch lengths from root
    depth_of: dict[int, float] = {}
    def _walk_depth(clade, depth=0.0):
        depth_of[id(clade)] = depth
        for c in clade.clades:
            bl = c.branch_length if c.branch_length is not None else 0.0
            _walk_depth(c, depth + bl)
    _walk_depth(tree.root)

    # Find MRCA per fine cluster
    mrcas = []
    skipped = 0
    for members in fine_members:
        leaves_in_tree = [leaf_by_name[n] for n in members if n in leaf_by_name]
        if not leaves_in_tree:
            mrcas.append(None)
            skipped += 1
            continue
        if len(leaves_in_tree) == 1:
            mrcas.append(leaves_in_tree[0])
        else:
            mrcas.append(tree.common_ancestor(leaves_in_tree))
    if skipped == n_fine:
        return {"status": "no_mrca_found"}

    # Precompute root paths for each MRCA (None entries kept as None)
    paths = [root_path(m) if m is not None else None for m in mrcas]

    def lca_depth(i: int, j: int) -> float:
        """Return root-depth of LCA(mrca_i, mrca_j) using cached paths."""
        pi, pj = paths[i], paths[j]
        if pi is None or pj is None:
            return float("inf")
        # find divergence point
        L = min(len(pi), len(pj))
        last_common = pi[0]
        for k in range(L):
            if pi[k] is pj[k]:
                last_common = pi[k]
            else:
                break
        return depth_of[id(last_common)]

    # ---- pairwise merge heights ----
    events = []
    for i in range(n_fine):
        if paths[i] is None:
            continue
        for j in range(i + 1, n_fine):
            if paths[j] is None:
                continue
            events.append((lca_depth(i, j), i, j))
    events.sort()

    # ---- union-find ----
    parent = list(range(n_fine))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
            return True
        return False

    # Count active groups (excluding fine clusters with no MRCA)
    valid_idx = [i for i in range(n_fine) if paths[i] is not None]
    active = len(valid_idx)
    for h, i, j in events:
        if active <= k_target:
            break
        if union(i, j):
            active -= 1
        if h == float("inf"):
            break

    # ---- collect groups ----
    groups: dict[int, list[int]] = {}
    for i in valid_idx:
        r = find(i)
        groups.setdefault(r, []).append(i)
    # Stable ordering: sort coarse cluster ids by min fine-cluster index
    sorted_groups = sorted(groups.values(), key=lambda lst: min(lst))

    return _emit_groups(pair_dir, fine_files, fine_members, coarse_dir,
                         fine_dir, sorted_groups, k_target, n_deep)


def _emit_1to1(pair_dir, fine_files, fine_members, coarse_dir, fine_dir,
               k_target, n_deep):
    """Edge case: K_fine <= K_target, no merging. 1:1 fine->coarse."""
    coarse_dir.mkdir(parents=True, exist_ok=True)
    for old in coarse_dir.glob("CoarseMsa_*.a3m"):
        old.unlink()
    rows = []
    for i, fp in enumerate(fine_files):
        coarse_path = coarse_dir / f"CoarseMsa_{i:03d}.a3m"
        # Read fine, write same content as coarse (1:1)
        coarse_path.write_text(fp.read_text())
        rows.append({
            "fine_id":   fp.stem,
            "coarse_id": coarse_path.stem,
            "n_fine":    len(fine_members[i]),
            "n_coarse":  len(fine_members[i]),
        })
    map_path = fine_dir / "fine_to_coarse.csv"
    with map_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["fine_id", "coarse_id",
                                            "n_fine", "n_coarse"])
        w.writeheader()
        w.writerows(rows)
    return {"status": "ok_1to1",
            "n_fine":   len(fine_files),
            "n_coarse": len(fine_files),
            "k_target": k_target,
            "n_deep":   n_deep}


def _emit_groups(pair_dir, fine_files, fine_members, coarse_dir, fine_dir,
                  groups, k_target, n_deep):
    coarse_dir.mkdir(parents=True, exist_ok=True)
    for old in coarse_dir.glob("CoarseMsa_*.a3m"):
        old.unlink()
    rows = []
    for c_idx, fine_indices in enumerate(groups):
        coarse_path = coarse_dir / f"CoarseMsa_{c_idx:03d}.a3m"
        # union sequences across all fine clusters in this group
        seen = set()
        with coarse_path.open("w") as fh:
            for fi in fine_indices:
                lines_pairs = read_a3m_lines(fine_files[fi])
                for hdr, aln in lines_pairs:
                    if hdr in seen:
                        continue
                    seen.add(hdr)
                    fh.write(f">{hdr}\n{aln}\n")
        n_coarse = len(seen)
        for fi in fine_indices:
            rows.append({
                "fine_id":   fine_files[fi].stem,
                "coarse_id": coarse_path.stem,
                "n_fine":    len(fine_members[fi]),
                "n_coarse":  n_coarse,
            })
    map_path = fine_dir / "fine_to_coarse.csv"
    with map_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["fine_id", "coarse_id",
                                            "n_fine", "n_coarse"])
        w.writeheader()
        w.writerows(rows)
    return {"status":   "ok",
            "n_fine":   len(fine_files),
            "n_coarse": len(groups),
            "k_target": k_target,
            "n_deep":   n_deep}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair", help="Single pair id")
    ap.add_argument("--foldpair_ids", default=None,
                    help="Comma-separated pair ids, or 'ALL'")
    ap.add_argument("--k-min", type=int, default=2)
    ap.add_argument("--k-max", type=int, default=30)
    ap.add_argument("--seqs-per-cluster", type=int, default=200)
    args = ap.parse_args()

    if args.pair:
        pair_list = [args.pair]
    elif args.foldpair_ids:
        if args.foldpair_ids == "ALL":
            pair_list = sorted(p.name for p in Path("Pipeline/FoldPairs").glob("*/")
                               if p.name not in ("_s4pred_work", "jobs"))
        else:
            pair_list = [s.strip() for s in args.foldpair_ids.split(",") if s.strip()]
    else:
        sys.exit("Must give --pair or --foldpair_ids")

    print(f"Tree-walk coarse-from-fine on {len(pair_list)} pair(s): "
          f"K_min={args.k_min}, K_max={args.k_max}, "
          f"seqs/cluster={args.seqs_per_cluster}")

    ok = skip = 0
    for pid in pair_list:
        pair_dir = Path("Pipeline/FoldPairs") / pid
        if not pair_dir.is_dir():
            print(f"[{pid}] SKIP: no pair dir")
            skip += 1
            continue
        r = build_coarse_for_pair(pair_dir, args.k_min, args.k_max,
                                   args.seqs_per_cluster)
        if r["status"].startswith("ok"):
            ok += 1
            print(f"[{pid}] {r['status']}: n_fine={r['n_fine']} -> "
                  f"n_coarse={r['n_coarse']} (K_target={r['k_target']}, "
                  f"N_deep={r['n_deep']})")
        else:
            skip += 1
            print(f"[{pid}] SKIP: {r['status']}")

    print(f"\nDone. {ok} pairs built, {skip} skipped.")


if __name__ == "__main__":
    main()
