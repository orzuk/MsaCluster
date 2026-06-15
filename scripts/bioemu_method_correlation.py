#!/usr/bin/env python3
"""Method x method agreement of per-cluster fold preference, pooled over all
(pair, cluster). Answers: is BioEmu redundant with AF2, or additive?

Uses TMdiff_centered from docs/fold_diversity_survey.csv. Reports, for each
method pair: sign-concordance (fraction agreeing on fold direction; ~0.50 =
chance, >0.70 = strong) and Spearman rho (if scipy present). BioEmu's row is
printed first; the full matrix follows so you can see where BioEmu lands
relative to AF2~AF3.
"""
import csv
import itertools
from collections import defaultdict

try:
    from scipy.stats import spearmanr
    def rho(a, b):
        return spearmanr(a, b).correlation
except Exception:
    def rho(a, b):
        return float("nan")

CSV = "docs/fold_diversity_survey.csv"

val = defaultdict(dict)        # (pair,cluster) -> {method: TMdiff_centered}
methods = set()
with open(CSV) as f:
    for r in csv.DictReader(f):
        try:
            x = float(r["TMdiff_centered"])
        except (ValueError, KeyError, TypeError):
            continue
        val[(r["pair_id"], r["cluster"])][r["method"]] = x
        methods.add(r["method"])
methods = sorted(methods)


def vecs(ma, mb):
    a, b = [], []
    for d in val.values():
        if ma in d and mb in d:
            a.append(d[ma]); b.append(d[mb])
    return a, b


def signconc(a, b):
    n = agree = 0
    for x, y in zip(a, b):
        if x == 0 or y == 0:
            continue
        n += 1
        agree += (x > 0) == (y > 0)
    return (agree / n if n else float("nan")), n


print(f"methods present: {methods}\n")

if "BioEmu" in methods:
    print("=== BioEmu vs each method (per-cluster fold-pref agreement) ===")
    res = []
    for m in methods:
        if m == "BioEmu":
            continue
        a, b = vecs("BioEmu", m)
        if len(a) < 10:
            continue
        sc, n = signconc(a, b)
        res.append((sc, rho(a, b), n, m))
    for sc, r, n, m in sorted(res, reverse=True):
        print(f"  BioEmu ~ {m:9s}  sign-conc={sc:.3f}  rho={r:+.3f}  n={n}")
    print()

print("=== full method-pair table (n >= 20) ===")
rows = []
for ma, mb in itertools.combinations(methods, 2):
    a, b = vecs(ma, mb)
    if len(a) < 20:
        continue
    sc, n = signconc(a, b)
    rows.append((sc, rho(a, b), n, ma, mb))
for sc, r, n, ma, mb in sorted(rows, reverse=True):
    star = " <== BioEmu" if "BioEmu" in (ma, mb) else ""
    print(f"  {ma:9s} ~ {mb:9s}  sign-conc={sc:.3f}  rho={r:+.3f}  n={n}{star}")
