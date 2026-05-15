#!/bin/bash
# diag_phylo_bh.sh — BH-correct phylogenetic-placement p-values per method
# Reads:  docs/phylo_placement.csv
# Writes: docs/phylo_placement_bh.csv
# Prints: method counts, BH-significant rows (q<0.10), nominal p<0.05 counts, top 15.

set -e
cd "$(dirname "$0")/.."

python3 <<'PYEOF'
import pandas as pd
import numpy as np

d = pd.read_csv('docs/phylo_placement.csv')
print('Methods:', sorted(d['method'].dropna().unique()))
print('Rows per method:')
print(d.groupby('method').size())

def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out

for col in ['parsimony_p', 'nn_concordance_p', 'nn3_concordance_p']:
    d[col + '_bh'] = np.nan
    for m, g in d.groupby('method'):
        mask = g[col].notna()
        if mask.sum() > 0:
            d.loc[g.index[mask], col + '_bh'] = bh(g.loc[mask, col].values)

print('\nBH-significant rows (q<0.10) on ANY p-column:')
flt = (d[['parsimony_p_bh', 'nn_concordance_p_bh', 'nn3_concordance_p_bh']] < 0.10).any(axis=1)
print(d.loc[flt, ['pair_id', 'method', 'n_leaves',
                  'parsimony_p', 'parsimony_p_bh',
                  'nn_concordance', 'nn_concordance_p', 'nn_concordance_p_bh']].to_string(index=False))

print('\nNominal p<0.05 counts per method (no correction):')
for col in ['parsimony_p', 'nn_concordance_p', 'nn3_concordance_p']:
    print(f'  {col}:')
    for m, g in d.groupby('method'):
        n_sig = (g[col] < 0.05).sum()
        n_total = g[col].notna().sum()
        print(f'    {m}: {n_sig}/{n_total}')

print('\nTop 15 by nn_concordance_p:')
print(d.sort_values('nn_concordance_p').head(15)[
    ['pair_id', 'method', 'n_leaves', 'nn_concordance', 'nn_concordance_p', 'parsimony_p']
].to_string(index=False))

d.to_csv('docs/phylo_placement_bh.csv', index=False)
print('\nWrote docs/phylo_placement_bh.csv')
PYEOF
