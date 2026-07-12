# MsaCluster

You can find the summary results table in this [link](https://orzuk.github.io/MsaCluster/protein_comparison_table.html), with full detailed results available for each protein pair in clickable links.

Global plots showing all pairs results can be found [here](https://orzuk.github.io/MsaCluster/pairs_global_analysis.html). 

A more detailed per-cluster analysis table is fond [here](https://orzuk.github.io/MsaCluster/protein_clusters_table.html). 

## Running from python
The main entry point is `run_foldswitch_pipeline.py`. It lets you run the entire pipeline or individual tasks.  
Usage (from command line):

```
python3 run_foldswitch_pipeline.py --run_mode {operation} --foldpair_ids <PAIR_ID> [options...]
```

If `--foldpair_ids` is not given, the program will loop over all pairs defined in `data/foldswitch_PDB_IDs_full.txt`.

### Example runs

```bash
# Cluster MSAs for a specific pair
python3 run_foldswitch_pipeline.py --run_mode cluster_msa --foldpair_ids 1dzlA_5keqF

# Run full pipeline for one pair (inline)
python3 run_foldswitch_pipeline.py --run_mode msaclust_pipeline --foldpair_ids 1dzlA_5keqF

# Run full pipeline for one pair, submitting AF2/MSA jobs via SLURM
python3 run_foldswitch_pipeline.py --run_mode msaclust_pipeline --foldpair_ids 1dzlA_5keqF --run_job_mode sbatch

# Run ESMFold on 3 sampled sequences per cluster with esm3 model on GPU
python3 run_foldswitch_pipeline.py --run_mode run_esmfold --foldpair_ids 1dzlA_5keqF --cluster_sample_n 3 --esm_model esm3 --esm_device cuda

# Run AF2 predictions for one pair (submit job)
python3 run_foldswitch_pipeline.py --run_mode run_AF --foldpair_ids 1dzlA_5keqF --run_job_mode sbatch

# Process all pairs listed in data/foldswitch_PDB_IDs_full.txt
python3 run_foldswitch_pipeline.py --run_mode get_msa
```

---

## Pipeline steps

Each `run_mode` corresponds to a stage in the pipeline:

- **`load`**: prepare sequences and structures (writes chain FASTAs).  
- **`get_msa`**: build a deep MSA seeded by both chains.  
- **`cluster_msa`**: cluster the deep MSA into shallow sub-MSAs (`--cluster_alg tree` (default), `ahc`, `hdbscan`, `tree_mono`).  
- **`run_AF`**: run AlphaFold2/AlphaFold3 for both sequences, on full MSA and per cluster (`--af_ver 2|3|both`).  
- **`run_cmap_msa_transformer`**: run MSA Transformer to generate attention/contact maps.  
- **`run_cmap_ccmpred`**: run CCMpred to generate contact maps (CPU-only).  
- **`run_esmfold`**: sample sequences from clusters and run ESMFold (`--esm_model esm2/esm3/both`).  
- **`run_boltz2`**: run Boltz-2 per cluster (`--boltz2_mode apo|holo`).  
- **`run_alphaflow`** / **`run_speachaf`** / **`run_cf_random`**: alternative per-cluster predictors.  
- **`run_bioemu`**: run BioEmu per-cluster equilibrium-ensemble fold preference.  
- **`run_ddg`**: ThermoMPNN per-cluster conformation-bias (ΔΔG) scoring.  
- **`run_s4pred`** / **`aggregate_s4pred`**: per-cluster secondary structure, folded into the survey.  
- **`tree`**: build a phylogenetic tree from the deep MSA.  
- **`postprocess`**: compute TM-scores and build per-pair/global result tables.  
- **`plot`**: generate plots (`--plot_scope pair|global|both`; PyMOL for 3D).  
- **`gen_pair_html`**: build per-pair HTML detail pages (`--capture_3d` for static Mol* PNGs).  
- **`compute_deltaG`**: compute ΔG metrics (requires PyRosetta).  
- **`clean`**: remove previous outputs for a pair.  
- **`msaclust_pipeline`**: run the full sequence (load → MSA → cluster → AF → cmap → esmfold → postprocess → plots).  

---

## Typical workflows

**Fresh run for one pair**
```bash
python3 run_foldswitch_pipeline.py --run_mode clean --foldpair_ids 1dzlA_5keqF
python3 run_foldswitch_pipeline.py --run_mode msaclust_pipeline --foldpair_ids 1dzlA_5keqF
```

**Update only AF2 predictions**
```bash
python3 run_foldswitch_pipeline.py --run_mode run_AF --foldpair_ids 1dzlA_5keqF --run_job_mode sbatch
```

**Generate plots for one pair**
```bash
python3 run_foldswitch_pipeline.py --run_mode plot --foldpair_ids 1dzlA_5keqF --plot_scope both
```

**Batch over all pairs**
```bash
python3 run_foldswitch_pipeline.py --run_mode get_msa
python3 run_foldswitch_pipeline.py --run_mode cluster_msa
python3 run_foldswitch_pipeline.py --run_mode run_AF --run_job_mode sbatch
python3 run_foldswitch_pipeline.py --run_mode run_cmap_msa_transformer
python3 run_foldswitch_pipeline.py --run_mode run_esmfold --cluster_sample_n 2
python3 run_foldswitch_pipeline.py --run_mode plot
```

---

## HURCS cluster

On HPC systems (e.g., SLURM-based clusters), add `--run_job_mode sbatch`.  
- `get_msa` submits via `scripts/shell/get_msa_params.sh` if present.  
- `run_AF` submits via `scripts/shell/RunAF_params.sh`.  
Edit these scripts to match your site’s environment and modules.

**Example:**  
```bash
python3 run_foldswitch_pipeline.py --run_mode msaclust_pipeline --foldpair_ids 1dzlA_5keqF --run_job_mode sbatch
```

---

## Other Scripts

The repository includes additional scripts for analysis and table generation.

### Analysis Folder (`Analysis/`)
- `AF_analysis.py` – compute TM-scores between AlphaFold (AF2+AF3) predictions and ground truth.
- `cmap_analysis.py` – evaluate contact maps from MSA Transformer / CCMpred vs truth.
- `esmfold_analysis.py` – compute TM-scores for ESMFold/ESM3 predictions.
- `postprocess_unified.py` – build per-pair, per-cluster, and global summary tables (absorbed the old `postprocess_global.py`).
- `class_level/` – class-level (trigger-type) analysis: does the fold-switch signal differ by trigger? (`run_all.py` orchestrates Fisher + LMM + permutation tests).
- `NotebookGen/generate_notebooks.py` – create HTML notebooks for exploration.

The bulk of the paper's quantitative analysis lives in `scripts/` (fold-diversity
survey, sequence-divergence correction, permutation nulls, BH FDR, phylogenetic
placement, ancestral reconstruction, trigger classification, PDB mining). See
`CODE_STRUCTURE.txt` for the full catalogue and the analysis-pipeline order.

### TableResults Folder (`TableResults/`)
- `gen_html_table.py` – generate interactive HTML summary table.
- `gen_latex_table.py` – generate LaTeX tables for paper.
- `gen_pair_html.py` – generate per-pair HTML pages with 3D structures, contact maps, and trees.

### Shell Wrappers (`scripts/shell/`)
Contains `.sh`/`.sbatch` wrappers for cluster jobs (`RunAF_params.sh`, `get_msa_params.sh`,
`RunBoltz2.sh`, `RunBioEmu.sh`, etc.).
Adapt them to your environment when using `--run_job_mode sbatch`.

### Documentation
- `CODE_STRUCTURE.txt` – comprehensive codebase summary, directory tree, run modes, and the full analysis-pipeline order.
- `AGENTS.md` – repository guidelines and conventions.
- `NEXT_TASKS.txt` – pending work and simplification plan.
- `INSTALL_bioemu.md` – BioEmu environment setup notes.
