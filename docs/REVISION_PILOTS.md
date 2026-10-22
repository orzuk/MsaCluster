# Revision pilots and parked runs

Work that is designed, costed and deliberately **not** run before the first
submission of the fold-switching evolution paper. Each entry says what question
it answers, what already exists in the repo, what it costs, how to launch it,
and what result would change the paper. Nothing here blocks submission.

Last updated 2026-10-21.

---

## 1. XCL1 curated-MSA pilot (Dishman-style), with true-clade clustering

**Question.** Our XCL1 MSA holds ~8,000 sequences of which only 4% are XCL1-like,
and only 12% of our clusters are monophyletic. Does the clade-shift scan localize
the gain of metamorphism on the XCL lineage when the alignment is a curated
vertebrate chemokine set and the clusters are true clades? This is the direct
test of whether our weak clade signal is a resolution artefact.

**Status.** Code complete, never run. Lives on branch
`feature/xcl1-dishman-pilot` (pushed, head `b98ff79`), **not merged into
master**. Five commits: `77e4a96`, `8064369`, `af8a632`, `068d090`, `b98ff79`.

**What the branch adds.**
- `scripts/build_curated_chemokine_msa.py` plus `data/pilots/xcl1_dishman/`
  (curated sequences, reference FASTAs) and `data/pilots/xcl1_dishman_species.tsv`
  (36 species). Dishman et al. used 457 curated chemokines; this rebuilds an
  equivalent set.
- `cluster_tree_monophyletic()` in `run_ClusterMSA.py` (`--cluster_alg tree_mono`),
  made usable on real gene-family trees by `8064369`.
- `--data_root` / `MSACLUSTER_DATA_DIR` so a pilot run resolves every pipeline
  path under its own root and cannot touch `Pipeline/FoldPairs/`.
- `scripts/shell/run_xcl1_dishman_pilot.sh` (driver), `scripts/xcl1_dishman_method.sbatch`,
  `scripts/xcl1_dishman_analysis.sbatch`.

**Launch.**
```bash
git checkout feature/xcl1-dishman-pilot
srun -A orzuk -p glacier --mem=8G -t 01:00:00 \
     bash scripts/shell/run_xcl1_dishman_pilot.sh --dry-run
```
Defaults: pair `2n54B_2hdmA`, 30 monophyletic clusters, min 5 sequences per
cluster, methods `af esmfold2 boltz2 bioemu ddg s4pred`, outputs under
`Pipeline/Pilots/xcl1_dishman/`.

**Cost.** 11 to 21 GPU-hours (L40S), plus CPU for the MSA and tree.

**Decision rule.** If the scan localizes the XCL clade on the curated, truly
monophyletic clustering and the signal survives the sequence-divergence
correction, the paper's "resolution, not blindness" reading is confirmed
directly and this becomes a main-text figure. If it does not, the limitation
is deeper than clustering resolution and the Discussion needs rewriting.

**RESULT (run 2026-10-22, 2.2 GPU-hours; outputs in `docs/pilots/xcl1_dishman/`).**
Partly confirmatory, not decisive; it does not yet earn a main-text figure.

The clustering did its job. mafft L-INS-i over the 994 curated chemokines and
`--cluster_alg tree_mono` (K=30, min clade 5) put the XCL lineage in clusters of
its own: `ShallowMsa_011` is 24 sequences, 100% XCL, mean 70% identity to human
XCL1, and holds the PDB query; `ShallowMsa_027` is the 5 remaining, more distant
XCL entries. In the production run the same lineage was spread over 4 of 100
clusters inside an 8,000-sequence MSA.

Boltz-2 and BioEmu then select exactly the 2-cluster XCL clade {011, 027}
(100% XCL, 24 mammalian + 5 non-mammalian) as their best split, which is the
Dishman lineage recovered exactly. Neither reaches significance: p = 0.18 for
both, because a clade of 2 of 30 clusters sits below what the scan can resolve
at this tree size (the power simulation predicts exactly this for f ~ 0.07).
The methods that do reach significance pick large, mostly-CCL clades instead:
ESMFold2 p = 0.0002 (13 clusters), ddG p = 0.0011 (14), combined p = 0.0019
(10), AF2 p = 0.043 (10), S4PRED p = 0.029 (10).

Moran's I stays at zero on the true-clade tree, from -0.026 to -0.035 for every
method, so a curated alignment and genuinely monophyletic clusters do NOT make
the tree-wide statistic see this gain. That strengthens the paper's existing
claim that Moran's I is the wrong instrument for a localized change, and it
shows the claim does not depend on our clustering being poor.

Cross-method concordance is higher here than in the screen: mean 0.678 over 6
methods, AF2~AF3 0.80, Boltz-2~BioEmu 0.786, 30 clusters.

**What would make it decisive.** Power, not clustering. The curated set holds
only 36 XCL sequences in 996, so the XCL clade is 1-2 clusters however it is
cut, and a 2-of-30 clade cannot reach significance. A decisive test needs a
chemokine set enriched for XCL and its closest sister groups, so that the
lineage spans enough clusters for the scan to resolve, or a scan statistic
with a size-matched null for very small clades.

---

## 2. True-clade (`tree_mono`) re-clustering of the full dataset

**Question.** The monophyly census reports that 12% of clusters are
monophyletic, with a median of 5 pieces per cluster, and that clade signal is
somewhat stronger where clusters are more tree-consistent (Spearman -0.32,
p = 0.0035 on uncorrected preferences). Would the whole screen look different
on clusters that are clades by construction?

**Status.** The clustering mode exists (`--cluster_alg tree_mono`, same branch
as pilot 1). Nothing has been run at dataset scale.

**Cost.** 50 to 100 GPU-hours for the five structure predictors over all 93
pairs, plus re-running every downstream test. Monophyletic clusters are fewer
and larger, so the per-family cluster counts drop and the permutation nulls get
coarser; expect lower power per family.

**Decision rule.** Run only if a referee asks whether the result is an artefact
of non-monophyletic clusters, or if pilot 1 shows a large gain. The honest
answer in the paper today is that non-monophyly can attenuate clade signal or
induce apparent structure, and we say so.

---

## 3. BioEmu2 as an eighth method

**Question.** BioEmu1 is deprecated ("This repository is no longer maintained.
Please check out BioEmu2", top of its README), so a referee will ask. Does
BioEmu2 agree with AF2 less than BioEmu1 does, and does it agree with Boltz-2
more?

**Status.** Reviewed 2026-10-21, not installed, not run. Full review:
Drive `Papers/BioEmu2_review_2026-10-21.md`. Announcement notes:
Drive `Papers/BioEmu2_announcement_2026-10-20.md`.

**What it is.** A new ~500M-parameter model trained from scratch in Aug-Sep
2026, not BioEmu1 with a bigger trunk: a 4-block MSA module feeding a 48-block
pairformer and a 16-layer all-atom EDM diffusion denoiser, with Boltz-derived
ops and preprocessing. MIT licence, code and weights, ungated. No paper, no
preprint, and the model card reports zero benchmark values.

**The load-bearing facts.**
- BioEmu1 ran on frozen AlphaFold2 Evoformer representations (its
  `get_embeds.py` says so and it downloads ~3.5 GB of AF2 weights), which
  explains the Spearman of 0.37 we observe between BioEmu and AF2 per cluster.
  BioEmu2 ships no AF2 code or weights, so the shared *computation* is gone.
- 8.06M of its 8.27M training systems (97.5%) are AlphaFold predictions
  (AFDB monomer 7.23M, AFDB multiconformation dimers 775k, Teddymer 50.6k),
  against 189,745 PDB entries and 19,434 MD trajectories. The shared *training
  target* is not gone, so the independence gain is smaller than the
  architecture suggests. Expect Spearman(BioEmu2, AF2) of 0.25 to 0.45, with a
  real chance that Spearman(BioEmu2, Boltz-2) is higher than either.
- `mode="equilibrium"` conditions on MD provenance rather than the AFDB branch
  (`src/bioemu2/structure/all_atom_structure.py:688`), which is a genuine
  mitigation we would be using.
- The model card forbids reading sampling frequency as a Boltzmann probability,
  so Delta stays ordinal, exactly as it is for BioEmu1 today.

**Integration.** Route A is a drop-in: per-cluster A3M is explicitly supported
(`docs/msa.md`), `EmuSession(mode="equilibrium", msa=...).sample(N)` yields
`samples.xtc` + `topology.cif`, which feeds our existing `safe_tm_vs_truth` /
`frac_toward_f1`. The only change to `run_bioemu.py` is `.cif` rather than
`.pdb` in `_split_ensemble_to_pdbs()`. It is all-atom native, so HPacker drops
out and the jax/protobuf `jaxfirst` workaround in `scripts/shell/RunBioEmu.sh`
is no longer needed.

**Install blockers (both soft, neither solved).** `requires-python >= 3.12`
while our venvs are 3.11, so it needs its own `bioemu2-venv`; and the default
`cuda` extra pulls cuEquivariance for CUDA 13 while Moriah has CUDA 12.4, so use
the `cuda12` extra with `--extra-index-url https://download.pytorch.org/whl/cu126`.
Pre-stage the 2.05 GB checkpoint on a login node into the existing `HF_HOME`,
since compute nodes are offline. Use pip, never `uv sync`: `pyproject.toml`
sets a default uv index pointing at an internal Microsoft mirror.

**Cost.** Minimum experiment: XCL1 `2n54B_2hdmA` (93 residues, 100 clusters)
plus the five Dishman ancestors, a few GPU-hours. Full sweep over ~6,700
clusters: 600 +/- 300 GPU-hours on an L40S at 3 to 7 min per cluster; ablate
`--num-steps` first, since 200 denoiser steps is probably overkill for an
ordinal signal. Time propagation is 5,000 to 11,000 GPU-hours and is not worth
running on anything but the ancestors.

**Decision rule.** Add as an eighth column, do not replace BioEmu1, since
having both answers the deprecation question with data. The three Spearman
correlations from the minimum experiment (against AF2, against Boltz-2, against
BioEmu1) decide whether the sweep is worth it. A result of
Spearman(BioEmu2, AF2) ~ Spearman(BioEmu1, AF2) ~ 0.37 is not a failure: it
would show that the AF2 correlation comes from the shared MSA input and shared
training distribution rather than from the shared Evoformer, which is a sharper
claim than the paper can make today.

**Schedule risk.** The only checkpoint is `bioemu2-preview` (tag `v0.1.0rc0`),
and `checkpoints.py` notes that a changed checkpoint gets a new name. If a 1.0
lands after we run the sweep, we re-run or explain. One more argument for
keeping the Science-published BioEmu1 column.

---

## 4. Condition-aware (holo) prediction for the triggered pairs

**Question.** 78 of 93 pairs switch on a trigger (ligand, oligomerization,
partner binding) that none of our seven methods can see, because every one of
them receives a bare sequence. Does supplying the deposited ligand, the
oligomeric copies or the partner chain change the per-cluster preference?

**Status.** Named as the main future direction in the paper
(main text Future directions; SI Sec. on folding-condition dependence). No code
beyond what AlphaFold3 and Boltz-2 already accept. BioEmu2 (pilot 3) would add
a third model that can take a ligand, through `Ligand(smiles=...)`.

**Cost.** Not estimated. Scoping it is the first task if this is picked up.

**Decision rule.** This is a second paper, not a revision item, unless a referee
demands a holo control on a handful of ligand-triggered hits.
