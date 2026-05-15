#!/bin/bash
# diag_ancestral_summary.sh — Summarize ancestral-reconstruction results
# Reads:  docs/ancestral_summary.csv
# Prints: column layout, pattern × root state, per-method counts,
#         nominal p<0.05 hits, all non-uniform rows.

set -e
cd "$(dirname "$0")/.."

echo "== phylo header =="
head -1 docs/phylo_placement.csv
echo
echo "== ancestral header =="
head -1 docs/ancestral_summary.csv
echo

python3 <<'PYEOF'
import pandas as pd

d = pd.read_csv('docs/ancestral_summary.csv')
print('Rows total:', len(d))
print('Columns:', list(d.columns))

print('\nPattern x root_state:')
print(d.groupby(['pattern', 'root_state']).size().unstack(fill_value=0))

print('\nMethod breakdown:')
print(d.groupby('method').size())

print('\nNominal p<0.05:')
sig_cols = ['pair_id', 'method', 'pattern', 'root_state', 'p_value',
            'n_clusters', 'n_transitions', 'n_fold_transitions']
sig_cols = [c for c in sig_cols if c in d.columns]
print(d[d['p_value'] < 0.05][sig_cols].to_string(index=False))

print('\nNon-uniform patterns (regardless of p):')
nu_cols = ['pair_id', 'method', 'pattern', 'root_state', 'p_value',
           'n_clusters', 'n_transitions']
nu_cols = [c for c in nu_cols if c in d.columns]
print(d[d['pattern'] != 'uniform'][nu_cols].head(40).to_string(index=False))
PYEOF
