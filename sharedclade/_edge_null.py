"""Expected pairwise agreement rate under the random-tree null."""
import itertools, csv, collections
import sys; sys.path.insert(0, "/tmp/claude-1000/-mnt-c-Code-Github-MsaCluster/fea976ea-5c76-4e3b-a35d-95918917c7f5/scratchpad/sharedclade")
import numpy as np
import analyse_shared_clade as A

rng = np.random.default_rng(7)
out = []
for label, cfg in A.DATASETS.items():
    placement, ncl = A.load_placement(cfg["placement"])
    uni = A.load_universe(cfg["survey"], cfg["value_col"])
    n_edge = n_tot = 0
    exp_edge = 0.0
    for pair, recs in placement.items():
        valid = {m: recs[m] for m in A.METHODS7
                 if m in recs and recs[m]["size"] > 0 and recs[m]["T"] is not None}
        ms = list(valid)
        if len(ms) < 2:
            continue
        U = {m: uni.get((pair, m), set()) for m in ms}
        allu = sorted(set().union(*U.values()))
        idx = {t: i for i, t in enumerate(allu)}
        N = len(allu)
        if N < 4:
            continue
        masks = {m: np.zeros(N, bool) for m in ms}
        for m in ms:
            for t in U[m]:
                masks[m][idx[t]] = True
        # observed
        for a, b in itertools.combinations(ms, 2):
            j = A.jaccard(valid[a]["members"], valid[b]["members"])
            sd = valid[a]["direction"] == valid[b]["direction"] != 0
            n_tot += 1
            n_edge += int(sd and j >= 0.5)
        # null: same-direction indicator kept, clade placement randomised
        ntree, per = 60, 60
        acc = collections.Counter(); den = collections.Counter()
        for _ in range(ntree):
            M = A._random_tree_splits(N, rng)
            if M.shape[0] == 0:
                continue
            ssz = M.sum(axis=1)
            pick = {}
            for m in ms:
                sz = valid[m]["size"]
                tol = max(1, int(round(0.25 * sz)))
                el = np.flatnonzero(np.abs(ssz - sz) <= tol)
                if el.size == 0:
                    d = np.abs(ssz - sz); el = np.flatnonzero(d == d.min())
                pick[m] = M[rng.choice(el, size=per)] & masks[m][None, :]
            for a, b in itertools.combinations(ms, 2):
                sd = valid[a]["direction"] == valid[b]["direction"] != 0
                if not sd:
                    den[(a, b)] += per; continue
                inter = (pick[a] & pick[b]).sum(1).astype(float)
                union = (pick[a] | pick[b]).sum(1).astype(float)
                jac = np.where(union > 0, inter / union, 0.0)
                acc[(a, b)] += int((jac >= 0.5).sum()); den[(a, b)] += per
        exp_edge += sum(acc[k] / den[k] for k in den)
    print(f"{label}: observed pairwise edge rate {n_edge/n_tot:.3f} "
          f"({n_edge}/{n_tot});  random-tree expectation {exp_edge/n_tot:.3f}")
