"""Fast Fitch reconstruction, gain/loss event typing and null tests for the
4-state (F1 / F2 / both / none) fold-capability character on the COARSE
cluster tree. ete3-free (works on any ete3-like tree: ``name``,
``children``, ``up``, ``traverse``, ``is_leaf``, ``is_root``,
``iter_leaves``), so it is unit-testable without ete3.

The reconstruction is the same Fitch parsimony as
``utils.ancestral_utils._fitch_bottom_up`` / ``_fitch_top_down`` (the Fitch
path of ``run_asr``): bottom-up intersection-else-union sets, top-down
"keep the parent's state if allowed, else the first allowed state in
``FITCH_PREF`` order". The number of transitions is the number of edges
whose reconstructed parent and child states differ (what
``run_permutation_test`` counts). Nodes are indexed once (``index_tree``) and
states are bitmasks, so 1000 permutations take well under a second.

Event types on an edge parent -> child (4-state character):
  gain   child == "both"  (gain of fold-switching capability)
  loss   parent == "both" (loss of fold-switching capability)
  swap   {parent, child} == {F1, F2}
  other  any other change (into / out of "none", or legacy "Amb")
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# Tie-break order for the Fitch top-down pass. utils.ancestral_utils imports
# this, so the figures and the tests resolve ties identically.
FITCH_PREF = ("both", "F1", "F2", "Amb", "none")
INFORMATIVE = ("F1", "F2", "both")
EVENT_TYPES = ("gain", "loss", "swap", "other")


# ---------------------------------------------------------------------------
# Tree indexing
# ---------------------------------------------------------------------------

def index_tree(tree) -> dict:
    """Flat preorder index of ``tree``.

    Returns {names, parent (int array; -1 at root), children (list of
    lists), postorder (list of node indices), leaves (node indices of the
    leaves in preorder), leaf_names, n_desc (#leaves below each node)}.
    Unnamed nodes get "N<preorder index>" (as ``run_asr`` keys them).
    """
    nodes = list(tree.traverse("preorder"))
    pos = {id(nd): i for i, nd in enumerate(nodes)}
    names = [nd.name if nd.name else f"N{i}" for i, nd in enumerate(nodes)]
    parent = np.full(len(nodes), -1, dtype=int)
    children: List[List[int]] = [[] for _ in nodes]
    for i, nd in enumerate(nodes):
        for c in nd.children:
            j = pos[id(c)]
            parent[j] = i
            children[i].append(j)
    post: List[int] = []
    stack = [(0, False)]
    while stack:                                   # iterative postorder
        i, done = stack.pop()
        if done:
            post.append(i)
            continue
        stack.append((i, True))
        for j in reversed(children[i]):
            stack.append((j, False))
    leaves = [i for i in range(len(nodes)) if not children[i]]
    n_desc = np.zeros(len(nodes), dtype=int)
    for i in post:
        n_desc[i] = 1 if not children[i] else sum(n_desc[j] for j in children[i])
    return {"names": names, "parent": parent, "children": children,
            "postorder": post, "leaves": leaves,
            "leaf_names": [names[i] for i in leaves], "n_desc": n_desc}


def _alphabet(states: Sequence[str]) -> List[str]:
    extra = sorted({s for s in states if s not in FITCH_PREF})
    return list(FITCH_PREF) + extra


# ---------------------------------------------------------------------------
# Fitch reconstruction
# ---------------------------------------------------------------------------

def fitch_reconstruct(ti: dict, leaf_states: Dict[str, str],
                      default: str = "Amb") -> List[str]:
    """Fitch states for every node of ``ti`` (preorder index order).

    ``leaf_states`` maps leaf name -> state; leaves absent get ``default``
    (as ``run_asr``). Ties: parent's state if allowed, else the first
    allowed state in FITCH_PREF order (then alphabetical for extra states).
    """
    leaf_sts = [leaf_states.get(ti["names"][i], default) for i in ti["leaves"]]
    alpha = _alphabet(leaf_sts)
    code = {s: 1 << k for k, s in enumerate(alpha)}
    n = len(ti["names"])
    sets = [0] * n
    for i, s in zip(ti["leaves"], leaf_sts):
        sets[i] = code[s]
    children = ti["children"]
    for i in ti["postorder"]:
        ch = children[i]
        if not ch:
            continue
        inter, uni = -1, 0
        for j in ch:
            inter &= sets[j]
            uni |= sets[j]
        sets[i] = inter if inter else uni
    order = [code[s] for s in alpha]
    st = [0] * n
    parent = ti["parent"]
    for i in range(n):                              # preorder: parent first
        if not children[i]:
            st[i] = sets[i]
            continue
        p = parent[i]
        if p >= 0 and (st[p] & sets[i]):
            st[i] = st[p]
        else:
            st[i] = next(b for b in order if b & sets[i])
    inv = {b: s for s, b in code.items()}
    return [inv[b] for b in st]


def count_changes(ti: dict, states: Sequence[str]) -> int:
    par = ti["parent"]
    return int(sum(1 for i in range(1, len(states))
                   if par[i] >= 0 and states[i] != states[par[i]]))


def event_type(parent_state: str, child_state: str) -> Optional[str]:
    """gain / loss / swap / other for a changed edge; None if unchanged."""
    if parent_state == child_state:
        return None
    if child_state == "both":
        return "gain"
    if parent_state == "both":
        return "loss"
    if {parent_state, child_state} == {"F1", "F2"}:
        return "swap"
    return "other"


def fitch_events(ti: dict, states: Sequence[str]) -> List[dict]:
    """One dict per changed edge: parent, child, from_state, to_state,
    event, child_is_leaf, n_leaves_below (child's clade size)."""
    out = []
    par = ti["parent"]
    for i in range(len(states)):
        p = par[i]
        if p < 0:
            continue
        ev = event_type(states[p], states[i])
        if ev is None:
            continue
        out.append({"parent": ti["names"][p], "child": ti["names"][i],
                    "from_state": states[p], "to_state": states[i], "event": ev,
                    "child_is_leaf": not ti["children"][i],
                    "n_leaves_below": int(ti["n_desc"][i])})
    return out


def leaf_transitions(ti: dict, leaf_names: Sequence[str],
                     leaf_state_list: Sequence[str], default: str = "none") -> int:
    """Fitch transitions for states given as a list aligned to ``leaf_names``."""
    return count_changes(ti, fitch_reconstruct(
        ti, dict(zip(leaf_names, leaf_state_list)), default=default))


# ---------------------------------------------------------------------------
# Null tests
# ---------------------------------------------------------------------------

def fitch_permutation_test(ti: dict, leaf_states: Dict[str, str],
                           n_perm: int = 1000, seed: int = 42,
                           default: str = "none"
                           ) -> Tuple[int, float, np.ndarray]:
    """Label-shuffle test of clade structure of a discrete character.

    Shuffles the observed leaf states over the leaves of the tree (state
    frequencies fixed) and recounts Fitch transitions.
    p = (1 + #{perm transitions <= observed}) / (n_perm + 1): low p = fewer
    changes than random placement = the states are clustered on the tree.
    Returns (observed transitions, p, null transition counts).
    """
    names = ti["leaf_names"]
    obs_states = [leaf_states.get(nm, default) for nm in names]
    obs = leaf_transitions(ti, names, obs_states, default)
    rng = np.random.default_rng(seed)
    arr = np.array(obs_states, dtype=object)
    null = np.empty(n_perm, dtype=int)
    for k in range(n_perm):
        null[k] = leaf_transitions(ti, names, list(rng.permutation(arr)), default)
    p = (1 + int((null <= obs).sum())) / (n_perm + 1)
    return obs, p, null


def state_order_by_frequency(labels: Sequence[str]) -> List[Tuple[str, int]]:
    """[(state, count)] most frequent first; ties broken by FITCH_PREF order
    (then alphabetically) so the assignment is deterministic."""
    alpha = _alphabet(labels)
    rank = {s: k for k, s in enumerate(alpha)}
    cnt: Dict[str, int] = {}
    for s in labels:
        cnt[s] = cnt.get(s, 0) + 1
    return sorted(cnt.items(), key=lambda kv: (-kv[1], rank[kv[0]]))


def discretise_like(values: Sequence[float], labels: Sequence[str],
                    rng: Optional[np.random.Generator] = None) -> List[str]:
    """Discretise a continuous trait into the observed label distribution.

    The leaves are ranked by ``values`` (descending) and the top n_1 get the
    most frequent observed state, the next n_2 the second most frequent, etc.
    (n_s = observed count of state s), so the returned labels have EXACTLY
    the observed label-frequency vector. Value ties are broken at random;
    leaves with a missing (NaN) value are inserted at uniformly random rank
    positions (no systematic state for them). ``rng`` defaults to seed 0.
    """
    v = np.asarray(values, dtype=float)
    n = len(v)
    if n != len(labels):
        raise ValueError("values and labels differ in length")
    rng = np.random.default_rng(0) if rng is None else rng
    fin = np.flatnonzero(np.isfinite(v))
    jitter = rng.random(len(fin))
    order = list(fin[np.lexsort((jitter, -v[fin]))])   # by -v, ties random
    for j in np.flatnonzero(~np.isfinite(v)):
        order.insert(int(rng.integers(0, len(order) + 1)), int(j))
    out = [""] * n
    pos = 0
    for s, c in state_order_by_frequency(labels):
        for j in order[pos:pos + c]:
            out[j] = s
        pos += c
    return out


def trait_transition_null(ti: dict, leaf_states: Dict[str, str],
                          trait_values: np.ndarray, obs_transitions: int,
                          seed: int = 0, default: str = "none"
                          ) -> Tuple[float, np.ndarray]:
    """Sequence-trait control for the fold-state transition count.

    ``trait_values`` is (K, n_leaves) aligned to ``ti['leaf_names']``. Each
    trait is discretised to the observed state frequencies
    (``discretise_like``) and its Fitch transitions counted.
    p_trait_anc = (1 + #{trait transitions <= observed fold transitions})
    / (K + 1): small = fold state changes more rarely on the tree than a
    generic sequence feature with the same state frequencies.
    Returns (p_trait_anc, per-trait transition counts).
    """
    names = ti["leaf_names"]
    labels = [leaf_states.get(nm, default) for nm in names]
    tv = np.atleast_2d(np.asarray(trait_values, dtype=float))
    rng = np.random.default_rng(seed)
    tr = np.array([leaf_transitions(ti, names, discretise_like(row, labels, rng), default)
                   for row in tv], dtype=int)
    p = (1 + int((tr <= obs_transitions).sum())) / (len(tr) + 1)
    return p, tr


def bh_adjust(pvals) -> np.ndarray:
    """Benjamini-Hochberg q-values (NaN-safe, monotone)."""
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    mask = np.isfinite(p)
    if not mask.any():
        return out
    pv = p[mask]
    m = len(pv)
    order = np.argsort(pv)
    adj = pv[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    res = np.empty(m)
    res[order] = np.clip(adj, 0.0, 1.0)
    out[mask] = res
    return out


# ---------------------------------------------------------------------------
# Per-pair core (tree + labels + optional sequence traits)
# ---------------------------------------------------------------------------

def ancestral_core(tree, leaf_states: Dict[str, str],
                   trait_by_leaf: Optional[Dict[str, np.ndarray]] = None,
                   n_perm: int = 1000, seed: int = 42, trait_seed: int = 0,
                   default: str = "none") -> Tuple[dict, List[dict], List[str]]:
    """Fitch reconstruction + label-shuffle test + sequence-trait control.

    ``leaf_states`` {leaf: F1/F2/both/none}; leaves of ``tree`` absent get
    ``default``. ``trait_by_leaf`` {leaf: (K,) random sequence-trait values};
    leaves absent are NaN (random rank in ``discretise_like``).
    Testable = >= 2 distinct informative states (F1 / F2 / both) among the
    leaves; the permutation test and the trait control run only then.
    Returns (summary dict, per-branch event list, node states in
    ``index_tree`` preorder).
    """
    ti = index_tree(tree)
    names = ti["leaf_names"]
    labs = [leaf_states.get(nm, default) for nm in names]
    lab_map = dict(zip(names, labs))
    states = fitch_reconstruct(ti, lab_map, default=default)
    events = fitch_events(ti, states)
    cnt = {s: labs.count(s) for s in ("F1", "F2", "both", "none")}
    inf = sorted({s for s in labs if s in INFORMATIVE})
    out = {"n_leaves": len(names), "n_F1": cnt["F1"], "n_F2": cnt["F2"],
           "n_both": cnt["both"], "n_none": cnt["none"],
           "n_other_state": len(labs) - sum(cnt.values()),
           "informative_states": ";".join(inf), "testable": len(inf) >= 2,
           "root_state": states[0], "n_transitions": len(events)}
    for ev in EVENT_TYPES:
        out[f"n_{ev}"] = sum(1 for e in events if e["event"] == ev)
    out.update({"p_value": np.nan, "perm_null_median": np.nan, "n_perm": 0,
                "K_trait": 0, "n_trait_leaves": 0, "p_trait_anc": np.nan,
                "trait_trans_median": np.nan, "trait_trans_min": np.nan,
                "frac_trait_le_obs": np.nan})
    if not out["testable"]:
        return out, events, states
    obs, p, null = fitch_permutation_test(ti, lab_map, n_perm, seed, default)
    out.update({"p_value": p, "perm_null_median": float(np.median(null)),
                "n_perm": int(n_perm)})
    if trait_by_leaf:
        K = len(next(iter(trait_by_leaf.values())))
        X = np.full((K, len(names)), np.nan)
        for j, nm in enumerate(names):
            if nm in trait_by_leaf:
                X[:, j] = trait_by_leaf[nm]
        p_t, tr = trait_transition_null(ti, lab_map, X, obs,
                                        seed=trait_seed, default=default)
        out.update({"K_trait": K,
                    "n_trait_leaves": int(sum(nm in trait_by_leaf for nm in names)),
                    "p_trait_anc": p_t,
                    "trait_trans_median": float(np.median(tr)),
                    "trait_trans_min": int(tr.min()),
                    "frac_trait_le_obs": float((tr <= obs).mean())})
    return out, events, states
