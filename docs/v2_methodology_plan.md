# v2 Paper Methodology Plan

**Status:** Planning phase — implementation pending user approval of details.
**Author:** Or Zuk + Claude collaboration, 2026-05-28.

This document defines the methodology for the v2 (revised) paper on
evolutionary signatures of fold-switching, replacing the v1 default-K
clustering with a hierarchical multi-resolution scheme validated on a
6-pair pilot (KaiB, XCL1, RelB, SOD1, RfaH, Selecase).

---

## 1. Motivation

The v1 paper used a single clustering at the default tree-cut depth
(typically K=10–22 clusters per pair, with `min_output_size=200` and
`min_neff=50`). This was tuned for methods that ingest the full
cluster MSA (e.g., CCMpred, MSA Transformer). It missed
clade-correlated fold-switch signal that AF2/AF3 can detect when fed
shallower MSAs — the Wayment-Steele 2024 "AF-Cluster" regime.

The v2 pilot (K=100 + medoid + top-10 nearest-by-Hamming sampling on
6 pairs) showed:

- **5 of 6 pilot pairs** improve std(ΔTM_centered) by 1.7–8.3× over v1.
- **4 of 5 improvers** pass phylogenetic-placement p_Fitch < 0.05.
- **RfaH** (the lone non-improver) is the expected biological null:
  its switch is context-gated by NusG/RNAP binding rather than
  clade-correlated, so MSA-based methods cannot detect it.

The challenge: AF2, AF3, ESMFold v1, ESMFold2, S4PRED, and DDG all
work well with shallow MSAs (or single sequences), but CCMpred, MSA
Transformer, and Boltz-2 require deeper MSAs (50+ effective
sequences, in some cases ≥128). A single fine-grained clustering
breaks the deep-MSA methods.

**Solution:** Two-resolution hierarchical tree-cut clustering. Build
one tree per pair; cut it at two depths; assign each method to its
appropriate resolution.

---

## 2. Clustering Method

### 2.1 Common preprocessing (unchanged from v1)

1. Build deep MSA via HHblits → `Pipeline/FoldPairs/<pair>/output_get_msa/DeepMsa.a3m`
2. Build phylogenetic tree (UPGMA on Hamming distances) →
   `Pipeline/FoldPairs/<pair>/output_phytree/DeepMsa_tree.nwk`

### 2.2 Fine-resolution clustering

Used by methods that work with shallow MSAs or single sequences.

- **K_target_fine:** `clip(N_seqs / 30, 20, 100)` per pair
- **min_output_size:** 10 (matches top-10 AF input)
- **min_neff:** 5
- **Per-cluster representative:** medoid (sequence minimizing
  sum-of-Hamming-distances to other cluster members)
- **Output directory:** `Pipeline/FoldPairs/<pair>/output_msa_cluster/`
- **File naming:** `ShallowMsa_NNN.a3m` (N = 3-digit cluster index)
- **Per-cluster size range:** typically 10–50 sequences

### 2.3 Coarse-resolution clustering

Used by methods that need deep MSA input.

- **K_target_coarse:** `clip(N_seqs / 200, 10, 30)` per pair
- **min_output_size:** 100 (CCMpred / MSAT depth requirement)
- **min_neff:** 30
- **Per-cluster representative:** medoid
- **Output directory:** `Pipeline/FoldPairs/<pair>/output_msa_cluster_coarse/`
- **File naming:** `CoarseMsa_NNN.a3m`
- **Per-cluster size range:** typically 100–500 sequences

### 2.4 Fine ↔ coarse mapping

Each fine cluster is fully contained in exactly one coarse cluster
(both are tree cuts; the coarse cut is at a higher node). A mapping
file `output_msa_cluster/fine_to_coarse.csv` is generated after both
clusterings:

| fine_cluster_id | coarse_cluster_id | n_seqs_fine | n_seqs_coarse | medoid_fine | medoid_coarse |
|---|---|---|---|---|---|
| ShallowMsa_000 | CoarseMsa_002 | 23 | 187 | UniRef_xxx | UniRef_yyy |
| ... | ... | ... | ... | ... | ... |

This mapping enables cross-resolution concordance analysis (§5.2).

---

## 3. Per-Method Resolution Assignment

Initial assignment (revisitable based on per-method results):

| Method | Resolution | Justification |
|---|---|---|
| **AF2** (per-cluster, top-10 nearest as MSA) | Fine | Wayment-Steele shallow-MSA regime; pilot-validated |
| **AF3** | Fine | Same |
| **ESMFold2** (Biohub, May 2026) | Fine | Single-sequence; uses cluster medoid only. **Replaces ESMFold v1**: the latest single-sequence predictor is enough; running both versions of the same architecture family adds nothing |
| **S4PRED** (secondary structure) | Fine | Single-sequence; orthogonal to fold-switching axis |
| **DDG** (ThermoMPNN-based) | **Fine** | Code already does its own top-10 sampling per cluster (`cluster_sample_n=10`); structure-based scoring, no deep-MSA need |
| **Boltz-2** | **Fine** | Reads cluster a3m directly; at fine resolution gets ~10–50 seq MSA which IS the shallow-MSA regime |
| **CCMpred** (DCA contacts) | Coarse | Requires Neff ≥ ~50 for coevolution |
| **MSA Transformer** | Coarse | Trained on MSAs ≥ 128 sequences |

### Audit results for DDG and Boltz-2 (2026-05-28)

**DDG (`run_DDG.py`):** ThermoMPNN-based per-cluster sequence
scoring (NOT coevolution-DDG). Code already samples
`cluster_sample_n=10` sequences per cluster, scores each against
fold-1 + fold-2 backbones, takes the difference. No deep-MSA
requirement → **fine resolution.**

**Boltz-2 (`run_Boltz2.py`):** Per-cluster fold prediction (AF3-style).
Reads `output_msa_cluster/ShallowMsa_*.a3m` and passes the cluster a3m
to Boltz as the chain's MSA. Whatever depth the cluster has, Boltz
uses. At fine resolution clusters are 10-50 seqs → effectively
shallow-MSA mode → **fine resolution.**

Minor consistency note: Boltz-2 currently uses the cluster
**consensus** sequence, not the medoid. AF2/AF3 use the medoid.
A one-line change in `run_Boltz2.py:_consensus_from_a3m` would align
them. Probably worth doing for cross-method consistency.

---

## 4. Pre-Registered Inclusion Criteria

To preserve BH-FDR validity, the inclusion criterion must be defined
before observing any p-values.

### 4.1 Fine resolution
- Include pair if it produces **≥ 10 fine clusters** after applying
  K_adapt(20,100) with min_output_size=10 / min_neff=5.
- Justification: K < 10 has insufficient phylogenetic resolution for
  the placement test.

### 4.2 Coarse resolution
- Include pair if it produces **≥ 5 coarse clusters** after applying
  K_adapt_coarse(10,30) with min_output_size=100 / min_neff=30.
- Justification: K < 5 yields trivial phylogenetic tests.

### 4.3 Reporting in the paper
- Methods section will state the count M_fine and M_coarse explicitly.
- Supplementary table lists each dropped pair with its MSA size and
  fail reason.

Predicted counts based on the v2 (20,100) clustering already
completed (with the old `min_output_size=30 / min_neff=20` filters):

- 61 pairs passed K ≥ 20 (with the old filters).
- Loosening to `min_output_size=10 / min_neff=5` is expected to
  add ~10–15 pairs. Final M_fine estimate: **75–85**.
- M_coarse will be lower (smaller MSAs may fail the min=100 filter):
  estimated **40–60**.

---

## 5. Statistical Analysis

### 5.1 Within-resolution per-method analysis (primary)

For each (pair, method) combination:

1. Run the method on its assigned resolution's clusters.
2. Per cluster: compute the method-specific per-cluster fold-preference
   score:
   - AF2/AF3/ESMFold variants: `TMdiff_centered = TM_F1 - TM_F2 - median(TM_F1 - TM_F2)`
   - DDG: `ΔΔG_F1 - ΔΔG_F2` (centered)
   - CCMpred/MSAT: a contact-map-based fold-preference score per cluster (TBD per method)
   - S4PRED: per-residue secondary-structure variation projected onto F1/F2 (TBD)
3. Run phylogenetic placement test (Fitch parsimony + permutation null) on the
   cluster-level labels to get a per-pair p_Fitch.
4. **Apply BH-FDR step-up correction at α=0.05 over the M_fine (or
   M_coarse) p-values for that method.**

Reported per method:
- Number of BH-significant pairs
- Discovery list (which pairs at FDR 0.05)
- Method's overall "signal strength" (e.g., median std(score), fraction of
  clusters with |score| > 0.1)

### 5.2 Cross-resolution / cross-method comparison (descriptive)

Two questions of interest:

**Q1: Do fine-resolution and coarse-resolution methods detect the
same pairs?**

For each pair, compute concordance:
- "Fine ensemble" = consensus of fine-resolution methods (mean of
  AF2/AF3/ESMFold2/etc. per-cluster centered scores, then collapse to
  coarse-cluster level using the fine→coarse mapping)
- "Coarse ensemble" = consensus of coarse-resolution methods (CCMpred
  / MSAT score per coarse cluster)
- Concordance metric: Pearson correlation of fine-ensemble (collapsed)
  vs coarse-ensemble across the pair's coarse clusters.

**Q2: Do methods agree on per-pair discovery?**

For each pair, list which methods called it BH-significant. Overlap
analysis (Venn / UpSet plot) shows method-specific vs method-shared
discoveries.

These are descriptive — no multiple-testing correction is applied to
cross-method comparisons, since each method tests a distinct hypothesis.

### 5.3 Cross-resolution analysis is not corrected

Each method's BH-FDR is computed within its own M (fine OR coarse, not
both). Cross-resolution comparisons are descriptive only.

---

## 6. Computational Cost Estimate

Based on the pilot timing (~28 GPU-hours per pair at K=99 for AF2+AF3):

### Fine resolution (per pair)
- AF2 (fine, ~50 clusters avg): ~14 GPU-hours
- AF3 (fine): ~14 GPU-hours
- ESMFold2 (fine, single-seq): ~10 min GPU
- S4PRED (fine, CPU): ~30 seconds per cluster
- DDG (fine): ~1 minute per cluster

Total per pair at fine: ~28 GPU-hours + minor CPU.

### Coarse resolution (per pair)
- CCMpred (coarse, ~10 super-clusters per pair): ~few CPU-hours per pair
- MSA Transformer (coarse): ~1 GPU-hour per pair
- Boltz-2 (coarse): TBD (depends on implementation)

Total per pair at coarse: ~1–2 GPU-hours.

### Total for 80 viable pairs
- Fine: 80 × 28 = 2,240 GPU-hours (AF dominates)
- Coarse: 80 × 2 = 160 GPU-hours

Grand total: **~2,400 GPU-hours.** Wallclock depends on slurm
parallelism; at 10 concurrent jobs, ~10 days.

---

## 7. Backup Plan

Before any v2 method run, back up all v1 method results.

```bash
BACKUP_DIR="backups/results_v1_$(date +%Y%m%d)"
mkdir -p "$BACKUP_DIR"
for D in Pipeline/FoldPairs/*/; do
    P=$(basename "$D")
    for SUB in output_AF output_DDG output_s4pred output_ccmpred \
               output_msa_transformer output_ESMFold output_esmfold \
               output_Boltz2 output_esm output_AlphaFold \
               output_AlphaFold3 output_msa_cluster; do
        SRC="${D}${SUB}"
        [[ -d "$SRC" ]] && cp -r "$SRC" "$BACKUP_DIR/${P}__${SUB}"
    done
done
du -sh "$BACKUP_DIR"
```

Estimated size: ~5–10 GB. Already done for `output_msa_cluster/`; rerun
to capture the other method directories.

Note: AF/ESMFold v2 outputs go to suffix-tagged dirs (`AF2_newmethod_top10`,
`output_esmfold2`), so they coexist with v1. For other methods, the
v2 runs will overwrite unless we explicitly route them elsewhere.

---

## 8. Implementation Plan (Code Changes Required)

### Phase 1: Pipeline support for two-resolution clustering
- Add `--cluster_resolution {fine, coarse, both}` to
  `run_foldswitch_pipeline.py`.
- Add `--cluster_outdir_coarse`, `--cluster_target_seqs_per_cluster_coarse`,
  `--cluster_k_min_coarse`, `--cluster_k_max_coarse`,
  `--cluster_min_output_size_coarse`, `--cluster_min_neff_coarse`.
- After clustering at both resolutions, generate `fine_to_coarse.csv`
  by walking the tree and finding which coarse cluster each fine cluster
  belongs to.

### Phase 2: Per-method clustering-input routing
For each deep-MSA method (CCMpred, MSATrans, Boltz-2):
- Add `--input_cluster_dir` argument so the wrapper reads from
  `output_msa_cluster_coarse/` instead of the default
  `output_msa_cluster/`.
- Audit the existing `run_*.py` for each method.

### Phase 3: Launcher updates
`scripts/run_newmethod_all_pairs.sh` gets new stages:
- `--stage cluster_fine` (current `cluster`)
- `--stage cluster_coarse`
- `--stage cluster_both`
- `--stage af_fine`, `--stage esmfold2_fine`, `--stage ddg_fine`, etc.
- `--stage ccmpred_coarse`, `--stage msat_coarse`, `--stage boltz2_coarse`

### Phase 4: Scoring + summarization
- `scripts/score_*.py` per method: each reads from its appropriate
  resolution's predictions and writes a CSV with
  `mode, pair, cluster, fold_pref_score`.
- `scripts/summarize_v2.py` (new, replacing `summarize_treek100_pilot.py`):
  reads all per-method CSVs at both resolutions, applies per-method BH,
  produces the headline cross-pair / cross-method table.

### Phase 5: Phylogenetic placement
- Extend `scripts/phylo_placement_treek100_pilot.py` (or rename to
  `scripts/phylo_placement_v2.py`) to handle both resolutions.
- Per-pair p_Fitch + p_NN1 computed from method's own cluster set.

### Phase 6: Cross-resolution concordance
- New `scripts/cross_resolution_concordance.py`: for each pair,
  collapse fine-resolution per-cluster scores to coarse level using
  the mapping, then correlate with coarse-resolution methods.

---

## 9. Decisions Made (2026-05-28)

User-confirmed decisions:

1. **Two resolutions** (not three) — fine + coarse.
2. **K_coarse:** `clip(N_seqs / 200, 10, 30)` (between 10 and 30
   coarse clusters per pair).
3. **DDG and S4PRED at fine first**, revisit if signal is weak.
4. **Phylogenetic placement per-method at the method's own resolution**:
   AF2 fine, MSAT coarse, etc.
5. **Per-method BH-FDR within-method** (M = M_fine or M_coarse for
   that method); cross-method comparisons descriptive only.
6. **Pre-register K ≥ 10 fine and K ≥ 5 coarse as inclusion criteria**.

---

## 10. Open Questions / TBD

1. ~~**DDG implementation audit:**~~ Resolved 2026-05-28 — DDG is
   ThermoMPNN-based, fine resolution.
2. ~~**Boltz-2 wrapper audit:**~~ Resolved 2026-05-28 — Boltz-2 reads
   cluster a3m directly, fine resolution. (Minor: change consensus
   → medoid for cross-method consistency.)
3. **CCMpred wrapper audit:** Confirm it reads from
   `output_msa_cluster/` and produces a per-cluster contact map.
   Need to add `--input_cluster_dir output_msa_cluster_coarse/` flag.
4. **MSA Transformer wrapper audit:** Same as CCMpred.
5. **Cross-resolution concordance metric:** Pearson correlation is the
   default; alternatives (Spearman, rank concordance) could be more
   robust.
6. **What happens if a pair has fine K = 12 but only coarse K = 3?**
   Pair is fine-eligible but coarse-ineligible. Reported separately.
7. **AF2_treek100_top10 → AF2_newmethod_top10 rename?** Pilot used
   `_treek100_top10` suffix; v2 plan uses `_newmethod_top10`. Either
   keep pilot results under old name (and document the equivalence)
   or re-run pilots under new name. Recommend: keep old name, document.

---

## 11. References

- Wayment-Steele et al. 2024. *Predicting multiple conformations via
  sequence clustering and AlphaFold2.* Nature 625, 832-839.
- Porter & Looger 2018; Porter 2022. (Fold-switch evolution framework.)
- Chakravarty et al. 2023. (Covert fold-switching.)
- Biohub 2026-05-27. ESMC + ESMFold2 + ESM Atlas release.
