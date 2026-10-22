"""Choose one representative tree leaf per MSA cluster (for the coarse tree).

``utils.ancestral_utils.build_coarse_tree`` prunes the per-pair sequence tree
to ONE leaf per MSA cluster. Clusters are usually not monophyletic on that
tree (they split into several "pieces"), so which leaf is kept decides where
the cluster sits on the coarse tree. Strategies:

  first                  the first tree leaf (in leaf order) of the cluster
                         -- the historical behaviour, arbitrary piece.
  largest_piece_medoid   (default) split the cluster's on-tree members into
                         maximal single-cluster clades ("pieces"), take the
                         piece with the most members (ties: first in tree
                         preorder), and keep its medoid leaf -- the member
                         minimising the summed patristic distance to the
                         other members of the piece (ties: leaf order).

Leaves that belong to no cluster are ignored, which is equivalent to pruning
them before labelling (as in ``scripts/cluster_monophyly_census.py``).

Cost: one postorder + one preorder pass over the tree for the piece labels,
plus a linear two-pass sum-of-distances inside each chosen piece. No
all-pairs distances over the whole tree.

ete3 is NOT imported: any tree with ``name``, ``dist``, ``children``, ``up``,
``is_leaf`` and ``traverse('postorder'|'preorder')`` works.
"""

from __future__ import annotations

from typing import Dict, Hashable, List, Tuple

REP_STRATEGIES = ("first", "largest_piece_medoid")
DEFAULT_REP_STRATEGY = "largest_piece_medoid"

_MIXED = object()   # node whose cluster-assigned leaves span >1 cluster


def _dist(n) -> float:
    d = getattr(n, "dist", None)
    return float(d) if d is not None else 0.0


def cluster_pieces(tree, leaf_to_tag: Dict[str, Hashable]):
    """Maximal single-cluster clades of every cluster.

    Returns ``{tag: [piece_root_node, ...]}`` in tree preorder. A piece root is
    a node all of whose cluster-assigned leaves belong to ``tag`` while its
    parent's do not (unassigned leaves are ignored). Also returns the
    per-node count of assigned leaves underneath (``n_assigned``).
    """
    label, n_assigned = {}, {}
    for n in tree.traverse("postorder"):
        if n.is_leaf():
            tag = leaf_to_tag.get(n.name)
            label[n] = tag            # None = unassigned leaf
            n_assigned[n] = 1 if tag is not None else 0
        else:
            labs = set()
            cnt = 0
            for c in n.children:
                cnt += n_assigned[c]
                if label[c] is not None:
                    labs.add(label[c])
            n_assigned[n] = cnt
            if not labs:
                label[n] = None
            elif len(labs) == 1 and _MIXED not in labs:
                label[n] = labs.pop()
            else:
                label[n] = _MIXED
    pieces: Dict[Hashable, list] = {}
    for n in tree.traverse("preorder"):
        lab = label[n]
        if lab is None or lab is _MIXED:
            continue
        if n.up is None or label.get(n.up) != lab:
            pieces.setdefault(lab, []).append(n)
    return pieces, n_assigned


def piece_medoid(piece_root, member_names: set) -> str:
    """Member leaf under ``piece_root`` minimising the summed patristic
    distance to the other members (linear time; ties -> leaf order)."""
    nodes = list(piece_root.traverse("postorder"))
    cnt, down = {}, {}
    for n in nodes:
        if n.is_leaf():
            cnt[n] = 1 if n.name in member_names else 0
            down[n] = 0.0
        else:
            cnt[n] = sum(cnt[c] for c in n.children)
            down[n] = sum(down[c] + cnt[c] * _dist(c) for c in n.children)
    k = cnt[piece_root]
    total = {piece_root: down[piece_root]}
    best, best_val = None, None
    for n in reversed(nodes):          # reverse postorder: parents first
        if n is not piece_root:
            # moving from parent to n: cnt[n] members get closer by dist(n),
            # the other k - cnt[n] get farther by dist(n)
            total[n] = total[n.up] + _dist(n) * (k - 2 * cnt[n])
    for n in piece_root.traverse("preorder"):   # leaf order for ties
        if n.is_leaf() and n.name in member_names:
            v = total[n]
            if best_val is None or v < best_val - 1e-12:
                best, best_val = n.name, v
    return best


def choose_cluster_reps(tree, cluster_to_leaves: Dict[Hashable, List[str]],
                        strategy: str = DEFAULT_REP_STRATEGY
                        ) -> Tuple[Dict[Hashable, str], Dict[Hashable, dict]]:
    """Pick one representative leaf per cluster.

    Parameters
    ----------
    tree : ete3-like tree (not modified)
    cluster_to_leaves : {tag: [leaf_name, ...]} on-tree members, leaf order
    strategy : one of ``REP_STRATEGIES``

    Returns
    -------
    reps : {tag: leaf_name}
    info : {tag: {"n_on_tree", "n_pieces", "rep_piece_size",
                  "rep_piece_frac"}} -- the piece containing the chosen leaf
           and the fraction of the cluster's on-tree members in it.
    """
    if strategy not in REP_STRATEGIES:
        raise ValueError(f"unknown rep strategy {strategy!r}; "
                         f"choose from {REP_STRATEGIES}")
    leaf_to_tag = {}
    for tag, leaves in cluster_to_leaves.items():
        for lf in leaves:
            leaf_to_tag[lf] = tag
    pieces, n_assigned = cluster_pieces(tree, leaf_to_tag)

    reps, info = {}, {}
    for tag, leaves in cluster_to_leaves.items():
        if not leaves:
            continue
        tag_pieces = pieces.get(tag, [])
        if strategy == "first":
            rep = leaves[0]
            piece = next((p for p in tag_pieces
                          if any(lf.name == rep for lf in _leaves(p))), None)
        else:
            # largest piece; max() keeps the first (preorder) on ties
            piece = max(tag_pieces, key=lambda p: n_assigned[p])
            rep = piece_medoid(piece, set(leaves))
        size = n_assigned[piece] if piece is not None else 1
        reps[tag] = rep
        info[tag] = {"n_on_tree": len(leaves), "n_pieces": len(tag_pieces),
                     "rep_piece_size": size,
                     "rep_piece_frac": size / len(leaves)}
    return reps, info


def _leaves(node):
    for n in node.traverse("preorder"):
        if n.is_leaf():
            yield n
