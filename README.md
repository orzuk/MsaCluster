# MsaCluster

You can find the summary results table in this [link](https://steveabecassis.github.io/MsaCluster/protein_comparison_table.html), with full detailed results available for each protein pair in clickable links.

A more detailed per-cluster analysis table is fond [here](https://steveabecassis.github.io/MsaCluster/protein_clusters_table.html). 

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
- **`cluster_msa`**: cluster the deep MSA into shallow sub-MSAs.  
- **`run_AF`**: run AlphaFold2 for both sequences, on full MSA and per cluster.  
- **`run_cmap_esm`**: run MSA Transformer to generate attention/contact maps.  
- **`run_esmfold`**: sample sequences from clusters and run ESMFold (esm2/esm3).  
- **`tree`**: build a phylogenetic tree from the deep MSA.  
- **`plot`**: generate pair-specific plots (requires PyMOL).  
- **`compute_deltaG`**: compute ΔG metrics (requires PyRosetta).  
- **`clean`**: remove previous outputs for a pair.  
- **`msaclust_pipeline`**: run the full sequence (load → MSA → cluster → AF2 → cmap → esmfold → plots).  

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
python3 run_foldswitch_pipeline.py --run_mode plot --foldpair_ids 1dzlA_5keqF --global_plots
```

**Batch over all pairs**
```bash
python3 run_foldswitch_pipeline.py --run_mode get_msa
python3 run_foldswitch_pipeline.py --run_mode cluster_msa
python3 run_foldswitch_pipeline.py --run_mode run_AF --run_job_mode sbatch
python3 run_foldswitch_pipeline.py --run_mode run_cmap_esm
python3 run_foldswitch_pipeline.py --run_mode run_esmfold --cluster_sample_n 2
python3 run_foldswitch_pipeline.py --run_mode plot
```

---

## HURCS cluster

On HPC systems (e.g., SLURM-based clusters), add `--run_job_mode sbatch`.  
- `get_msa` submits via `Pipeline/get_msa_params.sh` if present.  
- `run_AF` submits via `Pipeline/RunAF_params.sh`.  
Edit these scripts to match your site’s environment and modules.

**Example:**  
```bash
python3 run_foldswitch_pipeline.py --run_mode msaclust_pipeline --foldpair_ids 1dzlA_5keqF --run_job_mode sbatch
```

---

## Other Scripts

The repository includes additional scripts for analysis and table generation.

### Analysis Folder
- `cmap_analysis.py` – analyze contact maps from MSA Transformer.  
- `cmap_viz.py` – produce visualization plots.  
- `esmfold_analysis.py` – analyze ESMFold outputs.  
- `generate_notebooks.py` – create HTML notebooks for exploration.  

### TableResults Folder
- `gen_html_table.py` – generate HTML summary table.  
- `gen_latex_table.py` – generate LaTeX tables.  
- `get_raw_data.py` – extract raw data at cluster level.  
- `get_tm_align_score.py` – compute TM-scores.  
- `summary_table.py` – produce summary tables for reporting.  

### Pipeline Folder
Contains `.sh` wrappers for cluster jobs (`RunAF_params.sh`, `get_msa_params.sh`, etc.).  
Adapt them to your environment when using `--run_job_mode sbatch`.
