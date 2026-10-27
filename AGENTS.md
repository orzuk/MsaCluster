# Repository Guidelines

## Project Structure & Module Organization
- `run_foldswitch_pipeline.py` is the main entry point for the pipeline.
- `Pipeline/FoldPairs/<PAIR_ID>/` holds per-pair working directories and outputs (e.g., `Pipeline/FoldPairs/<PAIR_ID>/output_get_msa`); `Pipeline/Controls/` holds the 249 control proteins. SLURM/shell wrappers live in `scripts/shell/`.
- `Analysis/` contains post-processing and plotting utilities.
- `TableResults/` generates HTML/LaTeX tables and summary reports.
- `data/` stores input lists and reference tables; `docs/` is used for GitHub Pages artifacts.

## Build, Test, and Development Commands
- `python3 run_foldswitch_pipeline.py --run_mode msaclust_pipeline --foldpair_ids 1dzlA_5keqF` runs the full pipeline for one pair.
- `python3 run_foldswitch_pipeline.py --run_mode cluster_msa --foldpair_ids <PAIR_ID>` runs only MSA clustering.
- `python3 run_foldswitch_pipeline.py --run_mode plot --foldpair_ids <PAIR_ID> --plot_scope both` generates plots (requires PyMOL for 3D).
- `python3 run_foldswitch_pipeline.py --run_mode get_msa` batches over all pairs in `data/foldswitch_PDB_IDs_full.txt`.
- Tree post-processing (Moran's I, clade scans, coarse trees) runs locally on a slim copy of `Pipeline/FoldPairs/<PAIR_ID>/` holding only `output_phytree/DeepMsa_tree.nwk` and header-only `output_msa_cluster/ShallowMsa_*.a3m` (built on the cluster and streamed with `ssh huji-work "cat ...tgz"`). It needs `ete3` and `PyQt5` (`pip install --user --break-system-packages ete3 PyQt5`), `QT_QPA_PLATFORM=offscreen`, and `scripts/run_without_torch.py` to bypass the module-level torch import of `utils/plot_contacts.py`, e.g. `QT_QPA_PLATFORM=offscreen python3 scripts/run_without_torch.py scripts/phylo_placement.py --pairs ALL --labels corrected --survey-csv docs/fold_diversity_survey_corrected_9method.csv --tests morans_combined,clade_shift_combined --clade-methods AF2,AF3,ESM,Boltz2,BioEmu,DDG,S4PRED`.
- There is no pinned `requirements.txt`; install dependencies ad hoc (see the EXTERNAL DEPENDENCIES list in `CODE_STRUCTURE.txt`). BioEmu has its own setup in `INSTALL_bioemu.md`.

## Coding Style & Naming Conventions
- Python code uses 4-space indentation and PEP 8-style naming (`snake_case` for functions/vars, `UPPER_SNAKE` for constants).
- Keep new scripts in the most relevant folder (`scripts/`, `Analysis/`, or `TableResults/`); shell wrappers go in `scripts/shell/`.
- Configuration is centralized in `config.py`; extend it instead of hard-coding paths.

## Testing Guidelines
- There is no formal unit-test suite. Validate changes by running a small pipeline job on a single pair.
- `scripts/shell/test.sh` is a SLURM smoke-test wrapper; update it if you add tests.
- Name any new tests `test_*.py` and document how to run them in this file.
- `python3 tests/test_phylo_tests.py` and `python3 tests/test_power_negtraits.py` (no ete3/data needed; synthetic trees + toy a3m) cover `utils/phylo_tests.py`, `scripts/power_simulation_phylo.py` and `scripts/negative_traits_phylo.py`.
- `tests/test_phylo_tests.py` checks `utils/phylo_tests.py` (Moran's I / Fitch / NN / D unchanged vs the pre-registry code, vectorized clade-shift z vs `scipy.stats.mannwhitneyu`, planted-clade power, null calibration, and the pooled `morans_combined` test: null calibration under a joint leaf permutation, single-method agreement with the analytic Moran z, detection of a shared planted clade). Run with `python3 tests/test_phylo_tests.py` (or `python3 -m pytest tests/test_phylo_tests.py -q`); no Pipeline data needed, ete3 optional.
- `tests/test_plot_clade_shift.py` smoke-tests `scripts/plot_xcl1_clade_shift.py` (clade-shift tree + heatmap figure) on a synthetic 100-leaf tree and synthetic survey/clade CSVs. Run with `python3 tests/test_plot_clade_shift.py [OUT_DIR]` (OUT_DIR keeps the PNG/PDF for inspection); no Pipeline data needed, ete3 stubbed if missing.
- `tests/test_cluster_reps.py` checks `utils/cluster_reps.py` (coarse-tree representative leaf per cluster: `first` vs `largest_piece_medoid`, piece partition and linear-time medoid vs brute force on random trees). Run with `python3 tests/test_cluster_reps.py`; no Pipeline data needed, ete3 not required.
- `tests/test_contacts_fine.py` checks `utils/contact_utils.py` and `Analysis/cmap_analysis.py`: top-L binarisation (symmetric, exactly k pairs, respects `sep_min`), the per-output provenance manifest (`output_cmaps/<method>/_manifest.json`; hash match / re-clustered mismatch / legacy / output older than run), and that the scorer skips stale npz and records `top_l_factor` + `cluster_hash`. Run with `python3 tests/test_contacts_fine.py`; synthetic inputs, torch not required (stubbed).
- `tests/test_trait_null_ancestral.py` checks `scripts/trait_null_phylo.py` (fast I / I_norm / analytic z / clade T vs the registry, planted clade shift beats structured random traits while a fold vector drawn like a random trait gives ~uniform p_trait, Bonferroni + BH family table, vectorised random projections) and `utils/ancestral_tests.py` (fast Fitch identical to `ancestral_utils` Fitch, gain/loss/swap typing on a tiny tree, exact label-frequency preservation of the trait discretisation, trait control). Run with `python3 tests/test_trait_null_ancestral.py`; no ete3/data needed. Cluster run: `sbatch scripts/run_trait_null_ancestral.sbatch`.
- `tests/test_fp_calibration.py` checks the calibration and matched-dispersion analyses: `scripts/cluster_composition.py` (match-state-only composition features, `#` lines skipped, insertions ignored, `write_table` round trip), `scripts/concordance_fp_calibration.py` (`regress_out` identical to the `seq_divergence_correction.py` branch logic including the all-NaN and partial-NaN second-reference fallbacks, fold-free pseudo-methods landing on concordance 0.5, hand-counted sign concordance, BH monotonicity), `scripts/tm_subsample_matched.py` (exact statistics on a known matrix, and that matching the cluster budget removes the extreme-value gap between a 60-cluster and an 8-cluster target) and `scripts/ancestor_divergence_baseline.py` (Spearman vs scipy with ties, the exact 120-permutation p, the lineage ceiling bounding any monotone score, gap-aware PID). Run with `python3 tests/test_fp_calibration.py`; synthetic inputs only, no Pipeline data, no cluster, no GPU.
- `tests/test_ancestor_calibration.py` checks `scripts/ancestor_calibration.py` (the Dishman 2021 XCL1 ancestor calibration): the Delta sign conventions for TM / ThermoMPNN DDG / S4PRED SS (positive always = prefers F1 = the alternative fold), `rankdata`/`spearman_rho` vs scipy incl. ties, the EXACT permutation p over all 5! = 120 orderings (multiple of 1/120, null-calibrated, 2/120 floor from the Anc.0/Anc.2 tie), projecting an ancestor into the DeepMsa match-state frame, pseudo-pair staging + per-method collection on a synthetic pair, `--dry_run`, and a full `main(--collect_only)` producing the CSVs and figure. Run with `python3 tests/test_ancestor_calibration.py [KEEP_DIR]`; no Pipeline data, no GPU, heavy deps stubbed if missing. Cluster run: `sbatch scripts/run_ancestor_calibration.sbatch`.

## Commit & Pull Request Guidelines
- Commit messages are short, imperative, and sentence case (e.g., "Fixed Nans cmap analysis").
- PRs should include: a clear description, the run mode(s) executed, and any generated artifacts or figures.
- If changes affect outputs, call out which `Pipeline/<PAIR_ID>` directories were regenerated.

## Security & Configuration Tips
- Do not commit large generated outputs under `Pipeline/` or `docs/HTML/figs` unless required for publication.
- If running on SLURM, update `Pipeline/*_params.sh` scripts to match your environment and modules.
