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
- There is no pinned `requirements.txt`; install dependencies ad hoc (see the EXTERNAL DEPENDENCIES list in `CODE_STRUCTURE.txt`). BioEmu has its own setup in `INSTALL_bioemu.md`.

## Coding Style & Naming Conventions
- Python code uses 4-space indentation and PEP 8-style naming (`snake_case` for functions/vars, `UPPER_SNAKE` for constants).
- Keep new scripts in the most relevant folder (`scripts/`, `Analysis/`, or `TableResults/`); shell wrappers go in `scripts/shell/`.
- Configuration is centralized in `config.py`; extend it instead of hard-coding paths.

## Testing Guidelines
- There is no formal unit-test suite. Validate changes by running a small pipeline job on a single pair.
- `scripts/shell/test.sh` is a SLURM smoke-test wrapper; update it if you add tests.
- Name any new tests `test_*.py` and document how to run them in this file.

## Commit & Pull Request Guidelines
- Commit messages are short, imperative, and sentence case (e.g., "Fixed Nans cmap analysis").
- PRs should include: a clear description, the run mode(s) executed, and any generated artifacts or figures.
- If changes affect outputs, call out which `Pipeline/<PAIR_ID>` directories were regenerated.

## Security & Configuration Tips
- Do not commit large generated outputs under `Pipeline/` or `docs/HTML/figs` unless required for publication.
- If running on SLURM, update `Pipeline/*_params.sh` scripts to match your environment and modules.
