# Class-level analysis

Tests whether trigger classes (mutation / ligand / oligomerization /
protein_binding / equilibrium_or_unknown) differ in per-cluster
fold-preference variability, aggregated across 8 prediction methods.

Approach: aggregation happens on test statistics, not on the data.
For each (statistic, method) we run a per-class Kruskal-Wallis omnibus;
the 8 per-method p-values are combined via Fisher and Stouffer.
Complemented by a pair-level permutation test (sum of per-method H)
and a linear mixed model with pair as random intercept.

## Four analysis modes

1. **Raw** -- no controls
2. **Residualised on seq_identity** -- per-method OLS on seq_identity,
   then Kruskal/permutation on the residuals
3. **Identity-filtered** -- restrict to pairs with seq_identity >= 0.95
   (drops the low-identity outliers in the mutation class)
4. **LMM with seq_identity covariate** -- statsmodels mixedlm with
   seq_identity in both null and full formulas; LRT tests class effect
   after controlling for identity

## Key finding

The raw class effect (Fisher p ~ 1e-6 for stddev/max statistics) is
driven by the mutation class containing low-identity dataset members,
not by trigger-class biology. After per-method residualisation on
seq_identity, the variance/range class effects drop to non-significance
(Fisher p ~ 0.1). After restricting to identity >= 0.95 (which removes
the 9-65% identity outliers in the mutation class), the same pattern
holds.

The only statistic where class effect strengthens after the identity
control is |regression_beta| (the PID-to-query regression slope
magnitude per pair-method). This is a methodological finding rather
than biological: the training-proximity confound itself is
class-specific.

Bottom line: Chakravarty 2023's prediction that triggered switchers
should show more clade-level signal than equilibrium switchers is
NOT supported at our family-level MSA scale.

## Files

- `per_pair_per_method_stats.csv` -- one row per (pair, method) with
  the candidate statistics
- `per_class_summary.csv` -- per-class median/IQR/n per (statistic, method)
- `per_method_tests.csv` -- per-method Kruskal-Wallis and Jonckheere-Terpstra
  results, one row per (statistic, method, test, mode)
- `summary.csv` -- per-(mode, statistic) Fisher/Stouffer combined p,
  permutation p, LMM LRT p
- `summary.json` -- same as summary.csv but with full per-class
  coefficient tables from the LMM
- `scatter_seqid_<statistic>.png` -- per-method scatter of seq_identity
  vs statistic, coloured by trigger_class
- `summary_pvalues.png` -- -log10 p-value summary bar plot across modes

## How to reproduce

```bash
python3 -m Analysis.class_level.run_all          # default 10K perms
python3 -m Analysis.class_level.visualize        # produce PNGs
```

Requires: pandas, scipy, matplotlib. statsmodels optional (only used
for the LMM mode; Kruskal+Fisher+permutation work without it).
