"""Check best_split_batch reproduces utils.phylo_tests' scan (argmax split, sign)."""
import sys, importlib.util
sys.path.insert(0, "/mnt/c/Code/Github/MsaCluster")
import numpy as np
from utils.phylo_tests import _scan_core, average_ranks

spec = importlib.util.spec_from_file_location(
    "tan", "/tmp/claude-1000/-mnt-c-Code-Github-MsaCluster/fea976ea-5c76-4e3b-a35d-95918917c7f5/scratchpad/sharedclade/tree_aware_null.py")
tan = importlib.util.module_from_spec(spec); spec.loader.exec_module(tan)

rng = np.random.default_rng(0)
for trial in range(50):
    n, S = 30, 25
    M = rng.random((S, n)) < 0.4
    keep = (M.sum(1) >= 2) & ((n - M.sum(1)) >= 2)
    M = M[keep]
    x = rng.normal(size=n)
    if trial % 3 == 0:
        x[rng.choice(n, 4, replace=False)] = np.nan
    fin = np.isfinite(x)
    # reference: _scan_core on the non-missing leaves only (as run_clade_shift)
    names = np.flatnonzero(fin)
    Msub = M[:, names]
    r, tt = average_ranks(x[names])
    zc, _ = _scan_core(Msub, r[None, :], np.ones((1, len(names)), bool),
                       np.array([tt]), 2, None)
    zc = zc[0]
    a = np.abs(zc); a[~np.isfinite(a)] = -np.inf
    ref_best = int(np.argmax(a)); ref_sgn = np.sign(zc[ref_best])
    b, s, ok = tan.best_split_batch(M, x, np.arange(n)[None, :], 2)
    assert ok[0], trial
    assert b[0] == ref_best, (trial, b[0], ref_best, a[b[0]], a[ref_best])
    assert s[0] == ref_sgn, (trial, s[0], ref_sgn)
print("best_split_batch matches utils.phylo_tests scan on 50 random cases")

# agreement_stats sanity: identical clades -> clique = K
n, S, K = 20, 10, 3
M = np.zeros((S, n), bool)
for i in range(S):
    M[i, :i + 3] = True
best = np.zeros((K, 4), int); sgn = np.ones((K, 4)); ok = np.ones((K, 4), bool)
masks = np.ones((K, n), bool)
kd, knd, e = tan.agreement_stats(M, best, sgn, ok, masks)
assert (kd == K).all() and (knd == K).all() and (e == 3).all(), (kd, knd, e)
sgn[2] = -1
kd, knd, e = tan.agreement_stats(M, best, sgn, ok, masks)
assert (kd == 2).all() and (knd == K).all(), (kd, knd)
best[2] = S - 1            # a different, much larger clade
kd, knd, e = tan.agreement_stats(M, best, sgn, ok, masks)
assert (kd == 2).all() and (knd == 2).all(), (kd, knd)
print("agreement_stats OK (direction-aware and direction-agnostic cliques)")
