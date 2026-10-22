"""Global (across-family) test: is the number of families with k>=2,3,4
agreeing methods, and the total number of agreeing method pairs, higher
than the random-tree null expects?"""
import sys; sys.path.insert(0, "/tmp/claude-1000/-mnt-c-Code-Github-MsaCluster/fea976ea-5c76-4e3b-a35d-95918917c7f5/scratchpad/sharedclade")
import itertools, collections
import numpy as np
import analyse_shared_clade as A
IGNORE_DIR = __import__("os").environ.get("NODIR") == "1"

B = 2000
rng = np.random.default_rng(11)

for label, cfg in A.DATASETS.items():
    placement, ncl = A.load_placement(cfg["placement"])
    uni = A.load_universe(cfg["survey"], cfg["value_col"])
    obs_k, null_k, obs_e, null_e = [], [], 0, np.zeros(B)
    for pair, recs in placement.items():
        valid = {m: recs[m] for m in A.METHODS7
                 if m in recs and recs[m]["size"] > 0 and recs[m]["T"] is not None}
        ms = list(valid); k = len(ms)
        if k < 2:
            obs_k.append(0); null_k.append(np.zeros(B, int)); continue
        U = {m: uni.get((pair, m), set()) for m in ms}
        allu = sorted(set().union(*U.values())); idx = {t: i for i, t in enumerate(allu)}
        N = len(allu)
        adj = {}
        ne = 0
        for a, b in itertools.combinations(ms, 2):
            j = A.jaccard(valid[a]["members"], valid[b]["members"])
            sd = True if IGNORE_DIR else (valid[a]["direction"] == valid[b]["direction"] != 0)
            adj[(a, b)] = adj[(b, a)] = bool(sd and j >= 0.5)
            ne += int(adj[(a, b)])
        obs_k.append(len(A.max_clique(ms, adj))); obs_e += ne
        if N < 4:
            null_k.append(np.zeros(B, int)); continue
        masks = [np.zeros(N, bool) for _ in ms]
        for mi, m in enumerate(ms):
            for t in U[m]:
                masks[mi][idx[t]] = True
        sizes = [valid[m]["size"] for m in ms]
        dirs = [valid[m]["direction"] for m in ms]
        pl = list(itertools.combinations(range(k), 2)); pidx = {p: i for i, p in enumerate(pl)}
        sts, nes = [], []
        ntree, per = 50, B // 50
        for _ in range(ntree):
            M = A._random_tree_splits(N, rng)
            if M.shape[0] == 0:
                sts.append(np.ones(per, int)); nes.append(np.zeros(per, int)); continue
            ssz = M.sum(1)
            pick = np.zeros((k, per, N), bool)
            for mi, sz in enumerate(sizes):
                tol = max(1, int(round(0.25 * sz)))
                el = np.flatnonzero(np.abs(ssz - sz) <= tol)
                if el.size == 0:
                    d = np.abs(ssz - sz); el = np.flatnonzero(d == d.min())
                pick[mi] = M[rng.choice(el, size=per)] & masks[mi][None, :]
            edges = np.zeros((per, len(pl)), bool)
            for pi, (i, j) in enumerate(pl):
                if (not IGNORE_DIR) and (dirs[i] != dirs[j] or dirs[i] == 0):
                    continue
                inter = (pick[i] & pick[j]).sum(1).astype(float)
                union = (pick[i] | pick[j]).sum(1).astype(float)
                edges[:, pi] = np.where(union > 0, inter / union, 0.0) >= 0.5
            st = np.ones(per, int)
            for size in range(k, 1, -1):
                for sub in itertools.combinations(range(k), size):
                    need = [pidx[p] for p in itertools.combinations(sub, 2)]
                    ok = edges[:, need].all(1)
                    st = np.where(ok & (st < size), size, st)
            sts.append(st); nes.append(edges.sum(1))
        st = np.concatenate(sts); nn = np.concatenate(nes)
        null_k.append(st[:B]); null_e += nn[:B]
    NK = np.vstack(null_k)                 # (families, B)
    obs_k = np.array(obs_k)
    print(f"=== {label}")
    print(f"  agreeing method PAIRS: obs {obs_e}  null mean {null_e.mean():.1f} "
          f"sd {null_e.std():.1f}  p={(1+(null_e>=obs_e).sum())/(1+B):.4f}")
    for t in (2, 3, 4, 5):
        o = int((obs_k >= t).sum()); nd = (NK >= t).sum(0)
        print(f"  families with k>={t}: obs {o}  null mean {nd.mean():.1f} "
              f"sd {nd.std():.1f}  p={(1+(nd>=o).sum())/(1+B):.4f}")
