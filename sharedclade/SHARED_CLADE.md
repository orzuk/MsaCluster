# Do independent fold-preference methods localise the SAME branch?

Analysis of **existing** clade-shift scan results (no new predictions; CPU-only, local).

## 0. Data and definitions

**Input used (primary).** `round3/clade_shift_results/v07/phylo_placement_v07.csv`
(raw/centered) and `..._corrected_v07.csv` (corrected). I verified that for all
623 per-method rows (89 pairs x 7 methods) the columns `clade_T`, `clade_p`,
`clade_best_members`, `clade_best_direction`, `clade_best_size` are **byte-identical**
to `round2/.../medoidrep/phylo_placement_medoidrep.csv`: v07 is the same
medoid-representative scan. The two differ only in the `COMBINED` row, which in
v07 is computed over exactly the 7 methods
(`clade_methods = AF2,AF3,Boltz2,BioEmu,ESM,DDG,S4PRED` for 84/89 pairs, 6 methods
for 3, empty for 2), whereas round2/medoidrep mixes 7-, 8- and 9-method COMBINED
rows. **I therefore used v07.** MSAT and CCMpred were excluded throughout
(they were run on mismatched clusters: e.g. XCL1 has `clade_n_leaves = 21` for
MSAT/CCMpred vs 100 for the other seven).

**Leaf universes.** Per (pair, method) the set of clusters with a finite value was
taken from `final_docs/docs/fold_diversity_survey.csv` (`TMdiff_centered`) and
`fold_diversity_survey_corrected_9method.csv` (`TMdiff_residual`). These match
`clade_n_leaves` exactly for 622/623 rows, so the member sets and the universes are
consistent. Only BioEmu is ever short (18 pairs; e.g. 4/100 for `1qlnA_1h38D`).

**Agreement.** Two methods "agree" in a family iff
Jaccard(`clade_best_members`) >= 0.5 **and** `clade_best_direction` is equal.
The family statistic *k* is the size of the **largest clique** of mutually agreeing
methods (exhaustive over the <= 2^7 subsets). Sensitivity at Jaccard 0.3 and 0.8 is
reported in the CSV (`n_agree_j03`, `n_agree_j08`). Direction +1 means the clade
leans to **F1** (the first PDB chain in the pair id), per `run_clade_shift`'s
sign convention.

Where several cliques tie at the maximum size, the one reported is the first in
enumeration order over `AF2, AF3, Boltz2, BioEmu, ESM, DDG, S4PRED`; the full
pairwise matrix is in `shared_clade_pairwise.csv` if a different tie-break matters.

A useful property of the scan: `clade_split_matrix` canonicalises bipartitions and
all methods of a pair share the same split set, so two methods that pick the *same*
node report *identical* member sets (Jaccard exactly 1). Jaccard between 0.5 and 1
means nested or near-nested nodes.

**Outputs.** `shared_clade_per_family.csv` (178 rows = 89 pairs x {raw, corrected}),
`shared_clade_pairwise.csv` (3618 rows = 21 method pairs x 89 x 2),
`tree_aware_null.py` (follow-up script for the cluster), `analyse_shared_clade.py`
(this analysis), `_global_null.py`, `_edge_null.py` (null calibration),
`_test_tree_aware_null.py` (unit test of the follow-up script's scan maths).

---

## 1. How much do the methods agree?

| statistic | raw | corrected |
|---|---|---|
| families with k >= 2 | 78 / 89 | 70 / 89 |
| families with k >= 3 | 34 | 29 |
| families with k >= 4 | 12 | 12 |
| families with k >= 5 | 4 | 3 |
| families with k >= 2 at Jaccard 0.3 / 0.8 | 82 / 63 | 79 / 61 |
| families with k >= 3 at Jaccard 0.3 / 0.8 | 42 / 23 | 38 / 17 |
| overall pairwise agreement rate | 0.143 (259/1809) | 0.123 (223/1809) |
| mean pairwise Jaccard | 0.247 | 0.210 |
| pairs choosing an **identical** clade | 205 (11.3%) | 173 (9.6%) |
| **same-direction rate over all pairs** | **0.514** | **0.505** |

The raw agreement is modest and the *direction* component carries no information
at all on its own (51% = chance).

### Pairwise agreement rate by method (raw; J >= 0.5 and same direction)

|        | AF2 | AF3 | Boltz2 | BioEmu | ESM | DDG | S4PRED |
|---|---|---|---|---|---|---|---|
| AF2    |  -  |0.31|0.22|0.27|0.09|0.08|0.08|
| AF3    |0.31|  - |0.22|0.24|0.20|0.05|0.10|
| Boltz2 |0.22|0.22|  - |0.13|0.21|0.08|0.11|
| BioEmu |0.27|0.24|0.13|  - |0.12|0.13|0.08|
| ESM    |0.09|0.20|0.21|0.12|  - |0.06|0.13|
| DDG    |0.08|0.05|0.08|0.13|0.06|  - |0.10|
| S4PRED |0.08|0.10|0.11|0.08|0.13|0.10|  - |

Agreement concentrates in the structure-predictor block (AF2-AF3 highest, 0.31).
DDG and S4PRED agree with everything at roughly the 0.05-0.13 level. This matters:
a 4-method clique made of AF2+AF3+BioEmu+Boltz2 is close to **one** piece of
evidence, not four, because those predictors see the same MSA and are strongly
correlated.

### The branch is shared more often than the direction

| | raw | corrected |
|---|---|---|
| same-direction rate, all method pairs | 0.514 | 0.505 |
| ... conditional on Jaccard >= 0.5 | 0.617 | 0.592 |
| ... conditional on an **identical** clade | 0.688 | 0.618 |
| identical clade, **structure-structure** pair -> same direction | 0.79 | 0.79 |
| identical clade, **structure-sequence** pair -> same direction | 0.53 | 0.43 |
| identical clade, **DDG-S4PRED** pair -> same direction | 0.25 | 0.20 |

This is one of the clearest results here: **when a structure predictor and a
sequence-based method land on the same branch, which fold they say that branch
leans to is a coin flip.** Examples (raw): alpha1-antitrypsin `1kctA_3t1pA`, where
AF3, ESM and DDG all pick exactly the same 19-cluster clade but DDG's sign is
opposite; SARS-CoV spike `5wrgC_1wyyB`, where ESM, DDG and S4PRED all pick exactly
the same 19-cluster clade and DDG is again opposite. Direction-agnostically
("same branch, any direction") the counts are much higher: raw 83 / 52 / 18
families with k >= 2 / 3 / 4, corrected 80 / 42 / 23.

### Top families (ranked by k, then by the largest p in the clique)

#### raw

| rank | pair | family | k | methods | clade size | frac | leans | max p | p label-free | p rand-tree | strong | repl. |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2k42A_1ceeB | - | 6 | AF2,AF3,Boltz2,BioEmu,ESM,S4PRED | 49 | 0.49 | F2 | 1e-4 | 1e-4 | 0.707 | yes | no |
| 2 | 1eboE_5fhcJ | - | 6 | AF2,AF3,Boltz2,BioEmu,DDG,S4PRED | 22 | 0.50 | F1 | 3e-4 | 1e-4 | 0.484 | yes | **yes** |
| 3 | 5jytA_2qkeE | KaiB | 5 | AF2,AF3,Boltz2,BioEmu,DDG | 45 | 0.45 | F2 | 4e-4 | 1e-4 | 0.416 | yes | no |
| 4 | 1svfC_4wsgC | - | 5 | AF2,AF3,Boltz2,ESM,S4PRED | 2 | 0.50 | F2 | 0.329 | 7e-4 | 0.692 | yes | yes |
| 5 | 3lowA_3m1bF | - | 4 | AF2,AF3,Boltz2,ESM | 45 | 0.45 | F2 | 1e-4 | 1e-4 | 0.442 | no | no |
| 6 | 4qhhA_4qhfA | - | 4 | AF2,AF3,BioEmu,S4PRED | 52 | 0.61 | F2 | 1e-4 | 2e-4 | 0.881 | yes | **yes** |
| 7 | 2namA_1uxmK | **SOD1** | 4 | Boltz2,ESM,DDG,S4PRED | 21 | 0.21 | F2 | 3e-4 | 1e-4 | **0.014** | yes | no |
| 8 | 2mwfA_2nntA | - | 4 | AF2,AF3,Boltz2,ESM | 5 | 0.05 | F1 | 0.010 | 1e-4 | **0.003** | no | yes |
| 9 | 3hdfA_3hdeA | - | 4 | AF2,AF3,Boltz2,ESM | 9 | 0.09 | F1 | 0.443 | 1e-4 | **0.005** | no | no |
| 10 | 2pbkB_3njqA | - | 4 | AF2,AF3,Boltz2,BioEmu | 3 | 0.19 | F1 | 0.482 | 1e-4 | 0.101 | no | yes |
| 11 | 5i2sA_5i2mA | - | 4 | AF2,Boltz2,BioEmu,ESM | 2 | 0.33 | F2 | 0.497 | 3e-3 | 0.379 | no | yes |
| 12 | 4yhdG_7ahlE | - | 4 | AF2,AF3,Boltz2,BioEmu | 3 | 0.20 | F1 | 0.859 | 1e-4 | 0.047 | no | yes |
| 13 | 1nqjB_1nqdA | - | 3 | AF2,AF3,BioEmu | 41 | 0.41 | F2 | 1e-4 | 1e-4 | 0.745 | no | yes |
| 14 | 4zrbH_4zrbC | - | 3 | AF3,BioEmu,ESM | 33 | 0.33 | F2 | 1e-4 | 1e-4 | 0.182 | no | yes |
| 15 | 4h2aA_3j9cA | - | 3 | AF3,Boltz2,ESM | 34 | 0.57 | F2 | 2e-4 | 2e-4 | 0.783 | no | no |

#### corrected

| rank | pair | family | k | methods | clade size | frac | leans | max p | p label-free | p rand-tree | strong | repl. |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1eboE_5fhcJ | - | 6 | AF3,Boltz2,BioEmu,ESM,DDG,S4PRED | 22 | 0.50 | F1 | 0.341 | 1e-4 | 0.355 | yes | **yes** |
| 2 | 5jytA_2qkeE | **KaiB** | 5 | AF2,AF3,Boltz2,BioEmu,DDG | 9 | 0.09 | F2 | 0.009 | 1e-4 | **0.003** | yes | no |
| 3 | 1svfC_4wsgC | - | 5 | AF2,AF3,Boltz2,ESM,S4PRED | 2 | 0.50 | F2 | 0.329 | 1.4e-3 | 0.710 | yes | yes |
| 4 | 3m7pA_3ejhA | - | 4 | AF2,AF3,Boltz2,BioEmu | 10 | 0.10 | F1 | 1.5e-3 | 1e-4 | **0.005** | no | yes |
| 5 | 4qhhA_4qhfA | - | 4 | AF2,AF3,BioEmu,S4PRED | 52 | 0.61 | F2 | 4e-3 | 2e-4 | 0.926 | yes | **yes** |
| 6 | 5ly6B_5aoeB | **pneumolysin** | 4 | AF2,AF3,Boltz2,BioEmu | 24 | 0.42 | F1 | 0.095 | 1e-4 | 0.248 | no | **yes** |
| 7 | 2lv1A_6lniA | - | 4 | AF2,AF3,BioEmu,DDG | 3 | 0.12 | F2 | 0.109 | 1e-4 | **0.004** | yes | **yes** |
| 8 | 2kkwA_2n0aD | - | 4 | Boltz2,ESM,DDG,S4PRED | 14 | 0.44 | F1 | 0.389 | 1e-4 | 0.167 | yes | **yes** |
| 9 | 2lelA_2k0qA | - | 4 | AF2,BioEmu,ESM,S4PRED | 3 | 0.15 | F2 | 0.443 | 1e-4 | 0.096 | yes | no |
| 10 | 2kb8A_6vw2A | - | 4 | AF2,Boltz2,BioEmu,ESM | 7 | 0.35 | F2 | 0.605 | 1e-4 | 0.276 | no | no |
| 11 | 4yhdG_7ahlE | - | 4 | AF2,AF3,Boltz2,BioEmu | 3 | 0.20 | F1 | 0.849 | 3e-4 | 0.067 | no | yes |
| 12 | 3zxgB_5ec5P | - | 4 | AF2,AF3,Boltz2,ESM | 2 | 0.40 | F2 | 0.933 | 1.2e-3 | 0.173 | no | no |
| 13 | 4zrbH_4zrbC | - | 3 | AF3,BioEmu,ESM | 33 | 0.33 | F2 | 4e-4 | 1e-4 | 0.181 | no | yes |
| 14 | 1mnmD_1mnmC | - | 3 | AF3,Boltz2,ESM | 40 | 0.40 | F2 | 4.8e-3 | 1e-4 | 0.660 | no | yes |
| 15 | 1miqB_1qs8B | - | 3 | Boltz2,ESM,DDG | 20 | 0.20 | F2 | 0.022 | 1e-4 | **0.086** | yes | **yes** |

"repl." = the consensus clade of this family is the same (Jaccard >= 0.5, same
direction, k >= 2 in both) under raw and corrected labels. 30/89 families replicate.

---

## 2. Nulls

### 2a. The requested label-free null (10,000 draws per family)

For each family, each method's best clade was replaced by a uniformly random subset
of the **same size** drawn from that method's own leaf universe, with the observed
directions held fixed, and the agreement statistic recomputed.

Result: **median p = 1e-4 (raw) and 2e-4 (corrected); 87% of families with k >= 2
reach p < 0.05 under both labels.**

> **This null is wrong and it overstates significance badly.** Random subsets are
> not clades. A tree on *n* leaves offers only ~*n*-2 clades, they are nested, and
> clades of similar size overlap heavily, so two methods scanning the same tree are
> *far* more likely to land on overlapping sets than two uniform random subsets are.
> Nothing in section 1 should be defended with these p-values. They are reported in
> `shared_clade_per_family.csv` as `p_null_labelfree` only for completeness.

### 2b. A better stop-gap: the random-tree null

Still tree-free in the sense that the real topology is unavailable locally, but it
respects the *structure* of the problem: for each family I generate a random
(coalescent-shape) tree on the same leaf set, enumerate its clades with the same
canonicalisation as `clade_split_matrix`, and give each method a uniformly random
clade whose size is within 25% of its observed size (nearest size if none),
intersected with that method's own leaf universe. Directions are held fixed.
10,000 draws per family (100 trees x 100 draws). Column `p_null_randtree`.

| | raw | corrected |
|---|---|---|
| median p over families with k >= 2 | 0.59 | 0.56 |
| fraction p < 0.05 | 0.103 (8/78) | 0.114 (8/70) |

So of the 78 / 70 families with "agreement", only **8 and 8** survive a null that
merely knows how big the clades are and that they come from *some* tree, against
~4 and ~3.5 expected at the 5% level. Per family, the evidence is weak.

The expected *pairwise* agreement rate under this null is **0.100** (raw) and
**0.077** (corrected), against observed 0.143 and 0.123 — a real but modest 1.4x
and 1.6x excess.

### 2c. The aggregate signal is real

Because the per-family statistic is a coarse integer, the across-family counts are
much more powerful than any one family. Summing the same random-tree null over all
89 families (2000 joint draws):

| statistic | raw obs | raw null mean (sd) | raw p | corr obs | corr null mean (sd) | corr p |
|---|---|---|---|---|---|---|
| agreeing method pairs | 259 | 180.2 (12.4) | **5e-4** | 223 | 141.1 (10.3) | **5e-4** |
| families k >= 2 | 78 | 64.0 (3.3) | **5e-4** | 70 | 58.7 (3.3) | **5e-4** |
| families k >= 3 | 34 | 22.7 (2.7) | **5e-4** | 29 | 16.8 (2.6) | **5e-4** |
| families k >= 4 | 12 | 7.3 (1.7) | **0.012** | 12 | 4.9 (1.5) | **5e-4** |
| families k >= 5 | 4 | 2.7 (1.1) | 0.244 | 3 | 1.6 (0.9) | 0.131 |

Direction-agnostically (same branch, either direction) the excess is larger in
absolute terms: 420 vs 337.9 agreeing pairs (raw, p = 5e-4) and 377 vs 264.2
(corrected, p = 5e-4); families with k >= 3 are 52 vs 41.5 (p = 1.5e-3) and
42 vs 31.3 (p = 5e-4).

**Reading:** methods *do* co-localise more than chance, consistently, across the
whole panel. No individual family is carrying that signal.

### 2d. The proper null, for the cluster

`tree_aware_null.py` is the follow-up the main session should run on HUJI Moriah
(the cluster trees are not available locally). It rebuilds each pair's coarse tree
with `utils.ancestral_utils.build_coarse_tree(rep_strategy="largest_piece_medoid")`
— the same trees the medoid-representative scan used — permutes the per-cluster
values across the leaves, rescans all splits, re-derives each method's best clade
and sign, and recomputes the same clique statistic. Two nulls:

* **joint** — one leaf permutation applied to all methods. Keeps the cross-method
  correlation within a cluster (the reason AF2 and AF3 agree so often) and destroys
  only the association with the tree. **This is the number to quote.**
* **independent** — one permutation per method (what the task sketched); more
  liberal, report as a sensitivity check.

```
# on the cluster, from the repository root
python3 sharedclade/tree_aware_null.py --pairs ALL --labels centered \
    --survey-csv docs/fold_diversity_survey.csv \
    --n-perm 10000 --out clade_agreement_null_centered.csv
python3 sharedclade/tree_aware_null.py --pairs ALL --labels corrected \
    --survey-csv docs/fold_diversity_survey_corrected_9method.csv \
    --n-perm 10000 --out clade_agreement_null_corrected.csv
```

Its scan maths (`best_split_batch`) is unit-tested against
`utils.phylo_tests._scan_core` on 50 random cases (`_test_tree_aware_null.py`);
only the tree-loading step needs the cluster.

---

## 3. Robustness

For each family the clique was recomputed inside four restricted method sets, and
counted as "independent support" only when that restricted clique has >= 2 members
*and* all of them still match the family's consensus clade (Jaccard >= 0.5, same
direction). Columns `support_drop_DDG`, `support_drop_S4PRED`, `support_struct5`,
`support_DDG_S4PRED` in `shared_clade_per_family.csv`.

| restriction | raw families supported | corrected |
|---|---|---|
| drop DDG | 68 | 60 |
| drop S4PRED | 64 | 60 |
| the five structure predictors only | 59 | 53 |
| {DDG, S4PRED} only | 5 | 4 |

The asymmetry is the point. The agreement signal lives almost entirely in the
structure-predictor block; DDG and S4PRED agree with each other on the same clade
in only 5 (raw: `1eboE_5fhcJ`, **SOD1** `2namA_1uxmK`, `3lqcA_1xntA`, `3zxgB_5ec5P`,
`4ydqB_4twaA`) and 4 (corrected: `1eboE_5fhcJ`, `2kkwA_2n0aD`, `4ydqB_4twaA`,
`4dxrA_4dxtA`) families out of 89.

### Strong class (>= 1 structure predictor **and** >= 1 of DDG/S4PRED in the clique)

**28 families, raw; 28 corrected.** At k >= 3 this drops to 14 and 14; at k >= 4 to
6 and 7. Those with k >= 3 that also **replicate across raw and corrected**:

| pair | family | k (raw) | clique (raw) | clade | leans |
|---|---|---|---|---|---|
| 1eboE_5fhcJ | - | 6 | AF2,AF3,Boltz2,BioEmu,DDG,S4PRED | 22/44 | F1 |
| 1svfC_4wsgC | - | 5 | AF2,AF3,Boltz2,ESM,S4PRED | 2/4 | F2 |
| 4qhhA_4qhfA | - | 4 | AF2,AF3,BioEmu,S4PRED | 52/85 | F2 |
| 5ly6B_5aoeB | **pneumolysin** | 3 | AF3,BioEmu,DDG | 24/57 | F1 |
| 3lqcA_1xntA | - | 3 | Boltz2,DDG,S4PRED | 22/87 | F2 |
| 3g0hA_3ewsB | - | 3 | AF2,AF3,S4PRED | 33/100 | F1 |
| 2kkwA_2n0aD | - | 3 | Boltz2,BioEmu,DDG | 14/32 | F1 |

Caveats on these seven: `1svfC_4wsgC` has 4 clusters (a 2-of-4 split is nearly
forced; rand-tree p = 0.69), `4qhhA_4qhfA`'s clade is 61% of the family
(rand-tree p = 0.88), and none of the seven has rand-tree p < 0.05. `3lqcA_1xntA`
(p = 0.27) and `2kkwA_2n0aD` (p = 0.17) are the least compromised.

## 4. Direction, size and identity of the consensus clades

* Consensus clades are **large**: among families with k >= 3 the median consensus
  clade is 38% (raw) / 33% (corrected) of the family's clusters; only 9 / 11 of them
  are <= 20%. Near-root bipartitions are both easy to detect and easy to agree on
  by accident, which is exactly what the random-tree null punishes (see 2b: the
  high-k families with 40-60% clades have rand-tree p of 0.4-0.93, while the
  small-clade families — `2mwfA_2nntA` 5/100, `3m7pA_3ejhA` 10/100, KaiB-corrected
  9/100, `2lv1A_6lniA` 3/26 — are the ones with p < 0.01).
* Direction is balanced: among k >= 2 families, 44 lean F1 / 34 lean F2 (raw),
  36 / 34 (corrected).
* The `COMBINED` row's best clade matches the consensus clade (Jaccard >= 0.5) in
  34/78 (raw) and 38/70 (corrected) families — the combined scan is dominated by
  whichever methods have the largest |z|, not by the consensus.
* 30/89 families have the same consensus clade under raw and corrected labels
  (`same_clade_both_labels`). The sequence-divergence correction frequently moves
  the best clade (XCL1, SOD1, KaiB, CWC2 and alpha1-antitrypsin all move).

## 5. Families of interest

k = direction-aware clique; k(nd) = clique ignoring direction (same branch, either
fold). "p rt" = random-tree null. "strong" = structure + DDG/S4PRED in the clique.

| family | pair | labels | k | clique | k(nd) | clique (no dir) | consensus clade | leans | max p | p rt | strong | repl. |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| XCL1 | 2n54B_2hdmA | raw | 3 | AF2,AF3,BioEmu | 3 | AF2,AF3,BioEmu | 6/100 | F1 | 5e-4 | **0.018** | no | no |
| XCL1 | 2n54B_2hdmA | corrected | 2 | AF2,DDG | 2 | AF2,DDG | 33/100 | F2 | 5e-4 | 0.745 | yes | no |
| KaiB | 5jytA_2qkeE | raw | 5 | AF2,AF3,Boltz2,BioEmu,DDG | 6 | +S4PRED | 45/100 | F2 | 4e-4 | 0.416 | yes | no |
| KaiB | 5jytA_2qkeE | corrected | 5 | AF2,AF3,Boltz2,BioEmu,DDG | 6 | +S4PRED | 9/100 | F2 | 0.009 | **0.003** | yes | no |
| SOD1 | 2namA_1uxmK | raw | 4 | Boltz2,ESM,DDG,S4PRED | 4 | same | 21/100 | F2 | 3e-4 | **0.014** | yes | no |
| SOD1 | 2namA_1uxmK | corrected | 3 | AF2,AF3,BioEmu | 3 | same | 5/100 | F1 | 0.027 | **0.028** | no | no |
| CWC2 | 5lj3M_3tp2A | raw | 2 | AF2,Boltz2 | 2 | AF2,Boltz2 | 29/100 | F2 | 0.066 | 0.593 | no | no |
| CWC2 | 5lj3M_3tp2A | corrected | 2 | AF2,DDG | 2 | AF2,DDG | 65/100 | F1 | 1.7e-3 | 0.980 | yes | no |
| pneumolysin | 5ly6B_5aoeB | raw | 3 | AF3,BioEmu,DDG | 3 | same | 24/57 | F1 | 8e-4 | 0.581 | yes | **yes** |
| pneumolysin | 5ly6B_5aoeB | corrected | 4 | AF2,AF3,Boltz2,BioEmu | 4 | same | 24/57 | F1 | 0.095 | 0.248 | no | **yes** |
| SARS-CoV spike | 5wrgC_1wyyB | raw | 2 | ESM,S4PRED | 3 | ESM,DDG,S4PRED | 19/35 | F2 | 0.287 | 0.840 | yes | **yes** |
| SARS-CoV spike | 5wrgC_1wyyB | corrected | 3 | BioEmu,ESM,S4PRED | 4 | +DDG | 19/35 | F2 | 0.248 | 0.067 | yes | **yes** |
| alpha1-antitrypsin | 1kctA_3t1pA | raw | 2 | AF3,ESM | 3 | AF3,ESM,DDG | 19/100 | F1 | 1e-4 | 0.614 | no | no |
| alpha1-antitrypsin | 1kctA_3t1pA | corrected | 2 | AF2,DDG | 2 | AF2,DDG | 78/100 | F1 | 0.579 | 0.774 | yes | no |

Notes per family (raw labels unless stated):

* **XCL1 (2n54B_2hdmA).** AF2, AF3 and BioEmu pick the *identical* 6-cluster clade
  (ShallowMsa_005/033/037/038/051/075) leaning F1, with p = 5e-4, 2e-4, 4e-4 —
  the cleanest small-clade agreement in the set and rand-tree p = 0.018. But
  Boltz2 (8 clusters, F2), ESM (13, F1), DDG (33, F2) and S4PRED (81, F2) all go
  elsewhere, and the COMBINED row follows DDG/S4PRED to a 37-cluster F2 clade.
  Under corrected labels the AF2/AF3/BioEmu clade evaporates entirely (AF2 moves to
  a 22-cluster F2 clade; BioEmu and ESM move to a 4-cluster clade). **Not
  reproducible across label definitions, and not strong class.**
* **KaiB (5jytA_2qkeE).** The single most unanimous family: raw, six of seven
  methods choose the *same* 45-cluster bipartition (AF2/AF3/Boltz2/BioEmu/DDG all
  signed F2, S4PRED same branch but F1); ESM alone picks a 6-cluster clade. But
  45/100 is a near-root split, and the rand-tree null gives p = 0.42: with six
  methods and a half-of-the-tree clade, unanimity is cheap. The *corrected* run is
  more interesting — the same five methods agree on a **9-cluster** clade
  (ShallowMsa_000/032/033/035/042/052/061/...) leaning F2, with rand-tree p = 0.003,
  and again S4PRED picks that identical clade with the opposite sign. KaiB is the
  best candidate in the whole panel: structure predictors + DDG, a small clade,
  and it survives the conservative null — though the raw and corrected clades are
  *different* clades, which has to be acknowledged.
* **SOD1 (2namA_1uxmK).** Raw: ESM and DDG pick the identical 21-cluster F2 clade
  (p = 1e-4 both), S4PRED a nested 16-cluster F2 clade, Boltz2 a nested 12-cluster
  F2 clade — k = 4, **strong class**, rand-tree p = 0.014. Meanwhile AF2, AF3 and
  BioEmu all pick the identical *5*-cluster F1 clade (015/023/042/063/072). So SOD1
  has two competing, oppositely-signed consensus clades, one per method class.
  Under corrected labels the AF2/AF3/BioEmu 5-cluster F1 clade becomes the winner
  (k = 3, rand-tree p = 0.028, COMBINED agrees, p = 4e-4) and the ESM/DDG clade
  dissolves. Both signals are individually interesting; they are not the same
  branch, and the sequence-divergence correction flips which one wins.
* **CWC2 (5lj3M_3tp2A).** No real agreement. Raw: AF3 and BioEmu share an identical 17-cluster
  F1 clade and AF2/Boltz2 share an F2 clade (both cliques have size 2; the table
  reports AF2/Boltz2 only because of the enumeration order); DDG's 65-cluster F1 clade drives COMBINED. k = 2 at best, rand-tree
  p = 0.59 / 0.98.
* **Pneumolysin (5ly6B_5aoeB).** The one named family that **replicates**: a
  24-cluster F1 clade (ShallowMsa_023/024/026/027/028/032/033/...) is picked by
  AF3, BioEmu and DDG (raw; all p <= 8e-4) and by AF2, AF3, Boltz2 and BioEmu
  (corrected), with the same identity and direction in both. Rand-tree p is only
  0.58 / 0.25, and note that ESM picks a *different* 24-cluster clade
  (ShallowMsa_000/001/004/...) with the same F1 sign in raw, and DDG flips sign between raw and corrected.
* **SARS-CoV spike (5wrgC_1wyyB).** ESM, DDG and S4PRED all pick the *identical*
  19-cluster clade (of 35) in both label sets, and so does COMBINED; BioEmu joins it
  in corrected. This is a genuinely shared branch — but DDG's sign is F1 while
  ESM/S4PRED say F2, so the direction-aware clique drops to 2-3. Direction-agnostic
  k(nd) = 3 (raw) and 4 (corrected), rand-tree p = 0.067 corrected. The clearest
  "same branch, opposite fold" case.
* **alpha1-antitrypsin (1kctA_3t1pA).** Raw: AF3, ESM and DDG pick the identical
  19-cluster clade (006/007/011/018/025/027/029/...), AF3 and ESM at p = 1e-4 and
  DDG at p = 1e-4 — but DDG's sign is opposite, so k = 2 and k(nd) = 3. Under
  corrected labels the picture breaks apart (ESM and Boltz2 move to a different
  19-cluster clade, AF2 and DDG to a 78-cluster clade). Another "same branch,
  opposite direction" family, not replicating.

---

## 6. Honest assessment

**Is this a real signal?**

*At the panel level, yes, modestly.* Against the random-tree null — which already
grants the methods matched clade sizes, matched directions and a realistic nested
clade set — the observed 259 agreeing method pairs beat the expected 180 (p = 5e-4),
and the counts of families with k >= 3 and k >= 4 are significantly inflated under
both label sets. That excess is not explained by clade sizes or by tree shape. It
is, however, partly explained by something that is not independent evidence: the
agreement is concentrated in the AF2/AF3/Boltz2/BioEmu/ESM block (pairwise rates
0.13-0.31) whose members share an MSA and are mutually correlated, versus
0.05-0.13 for anything involving DDG or S4PRED. The `joint` null in
`tree_aware_null.py` is designed precisely to remove that confound, and it will
almost certainly shrink the excess.

*At the single-family level, no, not with these data.* Only 8/78 (raw) and 8/70
(corrected) families reach p < 0.05 under the random-tree null, barely above the
~4 expected; none of the 7 "strong class + replicating" families does. The
requested label-free null says 87% of families are significant, and that number
should not be used for anything.

**What is actually worth reporting as biology.**

1. **Direction is the weak link, not location.** Structure predictors and
   sequence-based methods land on the *same* branch at a rate well above chance
   (identical clade in ~10% of all method pairs, and conditional direction
   concordance 0.79 within the structure block), but when they do, they disagree
   about *which fold* that branch leans to roughly half the time (0.53 raw,
   0.43 corrected for structure-vs-sequence; 0.25/0.20 for DDG-vs-S4PRED). SARS-CoV
   spike and alpha1-antitrypsin are clean illustrations. Any claim of the form
   "method X and method Y agree that clade C switched toward fold F" needs to show
   the signs, not just the overlap.
2. **KaiB (corrected) is the single best candidate**: five methods including DDG on
   a 9-of-100-cluster clade, rand-tree p = 0.003, S4PRED on the same branch with
   the opposite sign.
3. **SOD1 is a two-signal family**, not a one-signal family: a large F2 clade found
   by Boltz2/ESM/DDG/S4PRED and a small F1 clade found by AF2/AF3/BioEmu, with the
   sequence-divergence correction deciding which one dominates.
4. **Pneumolysin** is the only named family whose consensus clade survives the
   raw -> corrected change intact.
5. **Big clades are not evidence.** Several of the highest-*k* families
   (`2k42A_1ceeB` 49/100, `4qhhA_4qhfA` 52/85, KaiB-raw 45/100, `1eboE_5fhcJ` 22/44)
   agree on a near-root bipartition, where agreement is nearly automatic
   (rand-tree p 0.42-0.93). The small-clade agreements are the informative ones.

**Recommended next step.** Run `tree_aware_null.py` on the cluster with both label
sets (10,000 permutations, `joint` and `independent`), then re-rank families by the
`joint` p-value. Everything above should be treated as provisional until that
number exists; my expectation is that the panel-level excess survives at reduced
magnitude and that KaiB-corrected, SOD1-raw and the small-clade families
(`2mwfA_2nntA`, `3m7pA_3ejhA`, `2lv1A_6lniA`) are the ones left standing.
