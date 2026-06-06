# BioEmu (method 9) — install & run

BioEmu (Lewis et al., *Science* 2025; MIT, https://github.com/microsoft/bioemu)
is a generative diffusion model that samples the **approximate equilibrium
ensemble** of a single protein chain. We use it as the 9th per-cluster method:
for each ShallowMsa cluster we sample N structures, TM-score each against both
truth folds, and record the **population balance** between the two folds
(`frac_toward_f1`) plus the ensemble-mean preference — a direct, physical read
on bistability that the point-estimate predictors (AF2/AF3/Boltz-2) can't give.

## Why a dedicated `bioemu-venv`

BioEmu pins its own torch + an AlphaFold2-evoformer stack (for the sequence
embedding) that conflicts with the existing envs (`my-python-venv`, `af-venv`,
`esmfold2-venv`, `torch2-venv`). Install it isolated:

```bash
python3 -m venv /sci/labs/orzuk/orzuk/venvs/bioemu-venv
. /sci/labs/orzuk/orzuk/venvs/bioemu-venv/bin/activate
pip install --upgrade pip
pip install bioemu            # pulls torch + the sampling stack
pip install mdtraj            # used to split the sampled trajectory into PDBs
# (first run auto-downloads model weights to ~/.cache or $BIOEMU_CACHE)
python -m bioemu.sample --help   # VERIFY the flag names below match your version
```

If `bioemu` needs ColabFold for its MSA step, follow its README; but we **pass
our own per-cluster MSA** (see below) so the embedding uses the same alignment
AF2/AF3 use — no generic MSA fetch.

## Wiring (already in the repo)

- `run_bioemu.py` — per-cluster runner (mirrors `run_Boltz2.py`): reads
  `Pipeline/<pair>/output_msa_cluster/ShallowMsa_*.a3m`, calls the wrapper,
  splits the ensemble into per-frame PDBs (mdtraj), TM-scores each vs both
  folds via `utils.align_utils`, writes `Analysis/df_bioemu{,_raw}.csv`.
- `scripts/shell/RunBioEmu.sh` — activates `bioemu-venv` and runs
  `python -m bioemu.sample --sequence <a3m> --num_samples N --output_dir <out>`.
  **Edit the CONFIG vars + verify the CLI flags** against your installed version.
- `utils/method_labels.py` — `BioEmu` registered (order: after Boltz-2).
- `scripts/fold_diversity_survey.py` — `load_bioemu_diversity()` reads
  `df_bioemu.csv` and feeds `BioEmu` into the survey + concordance exactly like
  the other 8 methods (ensemble-mean preference as the per-cluster signal,
  `frac_toward_f1` as the vote fraction).
- `run_foldswitch_pipeline.py` — `--run_mode run_bioemu` (GPU; HEAVY_PAIR_MODES;
  sbatch fan-out; `--bioemu_num_samples`).

## Run

```bash
# one pair (inline)
python3 run_foldswitch_pipeline.py --run_mode run_bioemu --foldpair_ids 2n54B_2hdmA

# more samples = better population-balance estimate
python3 run_foldswitch_pipeline.py --run_mode run_bioemu --foldpair_ids 2n54B_2hdmA --bioemu_num_samples 100

# fan out over pairs on the cluster (one sbatch job/pair, GPU)
python3 run_foldswitch_pipeline.py --run_mode run_bioemu --foldpair_ids ALL --run_job_mode sbatch

# pilot first on the candidate pairs rather than all 93:
#   2n54B_2hdmA (XCL1), 5jytA_2qkeE (KaiB), 1zk9A_3jv6A (RelB),
#   1eboE_5fhcJ, 1nrjB_2gedB (SR-beta)

# then re-aggregate the survey so BioEmu enters concordance/Moran's/plots:
bash scripts/run_v2_analysis.sbatch       # (or however the survey is regenerated)
```

## Output / interpretation

`Analysis/df_bioemu.csv` (one row per cluster): `n_samples`, `mean_tm_fold1`,
`mean_tm_fold2`, `mean_TMdiff` (signed: >0 leans fold1), `best_tm_fold1/2`,
`frac_toward_f1` (fraction of the ensemble closer to fold 1 — the bistability
signal). In the survey it appears as method `BioEmu`, picked up with no
special-casing.

**Pilot question:** does BioEmu populate BOTH folds for the known switchers
(XCL1, KaiB)? If a cluster comes out ~50/50 / clearly bimodal, that's the new
ensemble signal AF's single prediction can't show. Caveat: BioEmu is
AF2-evoformer-based, so it is **not fully independent** of AF2/AF3 — interpret
as an equilibrium/population refinement, not an orthogonal method.
