#!/usr/bin/env python3
import argparse
import subprocess
import shlex
import sys, os, warnings, re
from glob import glob
import json
import gzip

from pathlib import Path
from typing import List, Tuple
from copy import deepcopy
import shutil
import numpy as np
import pandas as pd

from config import *
from utils.utils import pair_str_to_tuple, ensure_dir, list_protein_pairs, write_pair_pipeline_script
from utils.protein_utils import read_msa, greedy_select, extract_protein_sequence, load_seq_and_struct, process_sequence
from utils.msa_utils import write_fasta, load_fasta, build_pair_seed_a3m_from_pair  # your existing writer
from utils.phytree_utils import phytree_from_msa
from utils.protein_plot_utils import make_foldswitch_all_plots

from Analysis.postprocess_unified import post_processing_analysis
from TableResults.gen_html_table import gen_html_from_summary_table, gen_html_from_cluster_detailed_table
from TableResults.summary_table import collect_summary_tables
from Bio import Align  # PairwiseAligner (modern replacement)


RUN_MODE_DESCRIPTIONS = {
    "load":             "Load the pair from PDB.",
    "get_msa":          "Download/build deep MSAs for the pair and write DeepMsa.a3m.",
    "cluster_msa":      "Cluster the pairwise MSA into shallow clusters (ShallowMsa_XXX.a3m).",
    "run_AF":           "Run AlphaFold (AF2/AF3/both) per chain × cluster. Use --af_ver {2,3,both}.",
    "run_esmfold":      "Run ESMFold on the pair. Use --esm_model {esm2,esm3,both}.",
    "run_cmap_msa_transformer":      "Run MSA-transformer on the pair to get contact maps.",
    "run_cmap_ccmpred": "Run CCMpred on DeepMsa and all ShallowMsa_XXX clusters to get contact maps.",
    "compute_deltaG":   "Compute ΔG stability metrics (requires PyRosetta).",
    "postprocess": "Compute TM/cmap metrics and build summary/detailed tables. Use --force_rerun_postprocess TRUE to recompute." ,
    "plot":             "Generate pair-specific plots (requires PyMOL).",
    "report":           "Generate html tables and reports.",
    "clean":            "Remove previous outputs for the pair.",
    "msaclust_pipeline":"Full pipeline: get_msa → cluster_msa → AF/ESM (as configured).",
    "help":             "Print this list of run modes with one-line explanations.",
}


# ------------------------- helpers -------------------------

def _has_deltaG(pair_id: str) -> bool:
    """Return True if both chain-level ΔG files exist for this pair."""
    pA, pB = pair_str_to_tuple(pair_id)     # e.g., '1wp8C','5ejbC'
    out_dir = Path("Pipeline/output_deltaG")
    a = out_dir / f"deltaG_{pA[:-1]}.txt"   # 'deltaG_1wp8.txt'
    b = out_dir / f"deltaG_{pB[:-1]}.txt"   # 'deltaG_5ejb.txt'
    return a.exists() and b.exists()



import getpass
from pathlib import Path

def _jupyter_env_for_scratch(base: str | None = None) -> dict[str, str]:
    """
    Return an env mapping that forces Jupyter/nbconvert/IPython/matplotlib caches
    into node-local or tmp storage, avoiding ~/.jupyter and NFS quotas.
    """
    env = os.environ.copy()

    user = getpass.getuser()
    base = base or os.environ.get("TMPDIR") or f"/tmp/{user}"
    scratch = Path(base) / "msaclust_jupyter"
    # Create all dirs up front so nbconvert won’t try to mkdir in HOME
    dirs = {
        "JUPYTER_CONFIG_DIR": scratch / "jupyter_config",
        "JUPYTER_DATA_DIR":   scratch / "jupyter_data",
        "JUPYTER_RUNTIME_DIR":scratch / "jupyter_runtime",
        "IPYTHONDIR":         scratch / "ipython",
        "MPLCONFIGDIR":       scratch / "mpl",
        "XDG_CACHE_HOME":     scratch / "xdg_cache",
    }
    for k, p in dirs.items():
        p.mkdir(parents=True, exist_ok=True)
        env[k] = str(p)

    # Last resort: also point HOME itself at scratch so *any* library that
    # insists on ~/.something won’t touch your real home.
    home_sandbox = scratch / "home"
    home_sandbox.mkdir(parents=True, exist_ok=True)
    env["HOME"] = str(home_sandbox)

    # (Optional) kernel hint—nbconvert will pick the one in your venv anyway
    env.setdefault("JUPYTER_KERNEL_NAME", "python3")

    return env


# --- helper: APC (Average Product Correction) ---
def _apc(x: np.ndarray) -> np.ndarray:
    r = x.mean(axis=1, keepdims=True)
    c = x.mean(axis=0, keepdims=True)
    t = x.mean()
    return x - r @ c / max(t, 1e-8)

def _load_a3m_strip_lower(a3m_path: str) -> list[str]:
    # reuse your utils to read & de-insert A3M rows
    entries = read_msa(a3m_path)     # -> [(id, seq), ...] (already strips insertions/lowercase)
    return [s for _, s in entries]


def _submit_msaclust_pair_job(pair_id: str, args: argparse.Namespace) -> None:
    """
    Submit ONE Slurm job that runs the full pipeline for a single pair INLINE.
    This avoids nested sbatch and lets the top-level launcher submit many jobs at once.
    """

    # Logs per pair
    jobs_dir = Path(f"Pipeline/{pair_id}/jobs")
    jobs_dir.mkdir(parents=True, exist_ok=True)
    log = jobs_dir / f"msaclust_{pair_id}.out"

    # build the wrapped python command that runs the full pipeline INLINE
    py = shlex.quote(sys.executable)  # same interpreter
    script = shlex.quote(Path(__file__).resolve().as_posix())
    # Propagate the knobs you care about; force INLINE inside the job
    cmd = (
        f"{py} {script} "
        f"--run_mode msaclust_pipeline "
        f"--foldpair_ids {shlex.quote(pair_id)} "
        f"--run_job_mode inline "
        f"--af_ver {shlex.quote(getattr(args, 'af_ver','both'))} "
        f"--esm_model {shlex.quote(getattr(args, 'esm_model','both') or 'both')} "
        f"--force_rerun {shlex.quote(getattr(args, 'force_rerun','FALSE'))} "
    )
    # Optional toggles you may want to carry through:
    if getattr(args, "global_plots", False):    cmd += " --global_plots"
    if getattr(args, "plot_trees",  False):     cmd += " --plot_trees"

    # sbatch resources
    gres = getattr(args, "sbatch_gres", "gpu:1")
    cpus = int(getattr(args, "sbatch_cpus", 8))
    mem = getattr(args, "sbatch_mem", "40G")
    time = getattr(args, "sbatch_time", "24:00:00")
    part = getattr(args, "sbatch_partition", None)
    cons = getattr(args, "sbatch_constraint", None)
    acct = getattr(args, "sbatch_account", None)
    qos = getattr(args, "sbatch_qos", None)
    mail = getattr(args, "sbatch_mail", None)
    mtyp = getattr(args, "sbatch_mail_type", None)

    sbatch_opts = f"--gres={gres} --cpus-per-task={cpus} --mem={mem} --time={time}"
    if part: sbatch_opts += f" -p {shlex.quote(part)}"
    if cons: sbatch_opts += f" --constraint={shlex.quote(cons)}"
    if acct: sbatch_opts += f" -A {shlex.quote(acct)}"
    if qos:  sbatch_opts += f" --qos={shlex.quote(qos)}"
    if mail: sbatch_opts += f" --mail-user={shlex.quote(mail)}"
    if mtyp: sbatch_opts += f" --mail-type={shlex.quote(mtyp)}"

    jobname = f"msaclust_{pair_id}"
    sb = f"sbatch {sbatch_opts} -J {shlex.quote(jobname)} -o {shlex.quote(str(log))} --wrap {shlex.quote(cmd)}"
    _run(sb, "sbatch")


def _find_best_af2_pdb(out_dir: str) -> str | None:
    """
    Heuristics for ColabFold outputs:
    - prefer symlink best_model.pdb (we create it)
    - then unrelaxed_rank_001
    - then ranked_0/rank_001
    - else pick any *.pdb deterministically
    """
    p = Path(out_dir)
    cand = p / "best_model.pdb"
    if cand.exists():
        return str(cand.resolve())
    # common CF patterns
    hits = list(p.glob("*unrelaxed_rank_001_*.pdb"))
    if hits: return str(sorted(hits)[0])
    hits = list(p.glob("ranked_0*.pdb")) + list(p.glob("ranked_1*.pdb"))
    if hits: return str(sorted(hits)[0])
    hits = list(p.glob("rank_001*.pdb")) + list(p.glob("rank_000*.pdb"))
    if hits: return str(sorted(hits)[0])
    hits = list(p.glob("*.pdb"))
    return str(sorted(hits)[0]) if hits else None

def _find_best_af3_pdb(out_dir: str) -> str | None:
    """
    After cif_to_pdb.sh, AF3 run dir usually contains rank1 PDB.
    - look for *rank1*.pdb first
    - else any *.pdb
    """
    p = Path(out_dir)
    hits = list(p.rglob("*rank1*.pdb"))
    if hits: return str(sorted(hits)[0])
    hits = list(p.rglob("*.pdb"))
    return str(sorted(hits)[0]) if hits else None

def _export_canonical_best_pdbs(pair_id: str, ver: str) -> None:
    """
    Walk .../output_AF/AF{ver}/{DeepMsa|ShallowMsa_xxx}/{chain}/
    and copy the best PDB to .../output_AF/AF{ver}/<ClusterLabel>__<chain>.pdb
    """
    root = Path(f"Pipeline/{pair_id}/output_AF/AF{ver}")
    if not root.exists():
        return
    for cluster_dir in sorted(root.iterdir()):
        if not cluster_dir.is_dir():
            continue
        label = cluster_dir.name  # DeepMsa or ShallowMsa_###
        for chain_dir in sorted(cluster_dir.iterdir()):
            if not chain_dir.is_dir():
                continue
            chain_tag = chain_dir.name  # e.g., 1fzpD or 2frhA
            # Resolve best PDB
            best = _find_best_af2_pdb(chain_dir) if ver == "2" else _find_best_af3_pdb(chain_dir)
            if not best:
                print(f"[export] no PDB in {chain_dir}")
                continue
            dst = root / f"{label}__{chain_tag}.pdb"
            try:
                # copy; overwrite if exists
                shutil.copy2(best, dst)
                print(f"[export] {dst.name}  <=  {Path(best).name}")
            except Exception as e:
                print(f"[export] failed {dst}: {e}")


def _bool_from_tf(s: str) -> bool:
    return str(s).strip().upper() == "TRUE"

def _af2_has_outputs(out_dir: str) -> bool:
    p = Path(out_dir)
    if not p.exists(): return False
    # Typical ColabFold/AF2 artifacts:
    if list(p.rglob("ranked_*.pdb")): return True
    if list(p.rglob("rank_*.pdb")): return True
    if list(p.rglob("ranking_debug.json")): return True
    if (p / ".af2_ok").exists(): return True
    # last resort: any pdb in out_dir
    return bool(list(p.rglob("*.pdb")))

def _af3_has_outputs(out_dir: str) -> bool:
    p = Path(out_dir)
    if not p.exists(): return False
    # AF3 writes a run subdir with *_ranking_scores.csv and *_model.cif
    if list(p.rglob("*_ranking_scores.csv")): return True
    if list(p.rglob("*_model.cif")): return True
    if (p / ".af3_ok").exists(): return True
    return False

def _convert_existing_af3(out_dir: str, mode: str = "all") -> None:
    """
    For each AF3 run directory under out_dir (directory that contains *_ranking_scores.csv),
    run the CIF→PDB helper. Always inline; it’s quick.
    """
    base = Path(out_dir)
    run_dirs = sorted({csv.parent for csv in base.rglob("*_ranking_scores.csv")})
    if not run_dirs:
        return
    helper = Path("Pipeline") / "cif_to_pdb.sh"
    for rd in run_dirs:
        cmd = f"bash {shlex.quote(str(helper))} {shlex.quote(str(rd))} --mode {shlex.quote(mode)}"
        subprocess.run(cmd, shell=True, check=False)


def _modes_epilog() -> str:
    lines = ["Run modes:"]
    for k, v in RUN_MODE_DESCRIPTIONS.items():
        lines.append(f"  {k:<18} {v}")
    return "\n".join(lines)



def ensure_chain_fastas(pair_dir: str, pdbids: list[str], pdbchains: list[str]) -> None:
    """
    Make ./Pipeline/<pair>/fasta_chain_files/<pdb><chain>.fasta for each chain.
    If a single-record root FASTA (<pdb>.fasta) exists, reuse it; otherwise
    extract the chain sequence from the PDB.
    """
    fasta_dir = os.path.join(pair_dir, "fasta_chain_files")
    os.makedirs(fasta_dir, exist_ok=True)

    for pdb, ch in zip(pdbids, pdbchains):
        tag = f"{pdb}{ch}"
        out_fa = os.path.join(fasta_dir, f"{tag}.fasta")
        if os.path.isfile(out_fa) and os.path.getsize(out_fa) > 0:
            continue

        # Try to reuse the root FASTA if it's single-record
        root_fa = os.path.join(pair_dir, f"{pdb}.fasta")
        seq = None
        if os.path.isfile(root_fa):
            ids, seqs = load_fasta(root_fa)
            if len(seqs) == 1:
                seq = seqs[0]

        # Otherwise extract chain from the PDB (sequence-based folding anyway)
        if not seq:
            pdb_path = os.path.join(pair_dir, f"{pdb}.pdb")
            seq = extract_protein_sequence(pdb_path, chain=ch)

        write_fasta([tag], [seq], out_fa)

def _in_slurm_session() -> bool:
    # true inside srun/salloc/sbatch job contexts
    env = os.environ
    return bool(env.get("SLURM_JOB_ID") or env.get("SLURM_STEP_ID") or env.get("SLURM_NODEID"))

def _is_windows() -> bool:
    return platform.system().lower().startswith("win")


def _run_inline(cmd: str) -> None:
    print(f"[inline] {cmd}", flush=True)
    subprocess.run(cmd, shell=True, check=True)

def _run(cmd: str, mode: str) -> None:
    if mode == "inline":
        return _run_inline(cmd)
    elif mode == "sbatch":
        # You can adapt per-site wrappers here if needed
        print(f"[sbatch] {cmd}", flush=True)
        subprocess.run(cmd, shell=True, check=True)
    else:
        raise ValueError(f"Unknown run_job_mode: {mode}")


def _cluster_files(pair_id: str) -> tuple[str, list[str]]:
    deep = f"Pipeline/{pair_id}/output_get_msa/DeepMsa.a3m"   # <-- FIXED
    clusters = sorted(glob(f"Pipeline/{pair_id}/output_msa_cluster/ShallowMsa_*.a3m"))
    return deep, clusters

def _has_deep_msa(pair_id: str) -> bool:
    f = Path(f"Pipeline/{pair_id}/output_get_msa/DeepMsa.a3m")
    return f.is_file() and f.stat().st_size > 0

def _has_cluster_msas(pair_id: str) -> bool:
    return bool(list(Path(f"Pipeline/{pair_id}/output_msa_cluster").glob("ShallowMsa_*.a3m")))

def _has_cmaps(pair_id: str) -> bool:
    return bool(list(Path(f"Pipeline/{pair_id}/output_cmaps/msa_transformer").glob("*.npy")))

def _has_esm_model(pair_id: str, model: str) -> bool:
    d = Path(f"Pipeline/{pair_id}/output_esm_fold/{model}")
    if not d.is_dir():
        return False
    tsv = d / "samples_index.tsv"
    if tsv.is_file() and tsv.stat().st_size > 0:
        try:
            return len(pd.read_csv(tsv, sep="\t")) > 0
        except Exception:
            pass
    # fallback: any PDBs written by ESM
    return bool(list(d.glob("*.pdb")))


def _sample_to_tmp_fastas(pair_id: str, sample_n: int, include_deep: bool = False) -> List[str]:
    """
    For each cluster, sample up to sample_n sequences into:
      Pipeline/<pair>/_tmp_ShallowMsa_XXX.fasta

    If include_deep=True, also writes:
      Pipeline/<pair>/_tmp_DeepMsa.fasta

    Returns list of created fasta paths.
    """
    created = []
    deep, clusters = _cluster_files(pair_id)
    outdir = f"Pipeline/{pair_id}/tmp_esmfold"
    ensure_dir(outdir)


    def _write_sample(src_a3m: str, dst_fa: str) -> None:
        # If we already sampled this cluster, keep those choices
        if os.path.isfile(dst_fa) and os.path.getsize(dst_fa) > 0:
            created.append(dst_fa)
            return
        entries = read_msa(src_a3m)  # [(id, seq), ...] (may still contain '-' gaps)
        picked = greedy_select(entries, min(sample_n, len(entries)))
        names, seqs = zip(*picked) if picked else ([], [])
        # sanitize: drop '-' and any non-AA so downstream ESMFold never sees gaps
        seqs = [process_sequence(s) for s in seqs]
        if not names:
            return
        write_fasta(names, seqs, dst_fa)
        created.append(dst_fa)

    for a3m in clusters:
        cluster_id = Path(a3m).stem.split("_")[-1]  # e.g., '003'
        dst = f"{outdir}/_tmp_ShallowMsa_{cluster_id}.fasta"
        _write_sample(a3m, dst)

    if include_deep and os.path.isfile(deep):
        _write_sample(deep, f"{outdir}/_tmp_DeepMsa.fasta")

    if not created:
        print(f"[warn] No clusters or empty MSAs for {pair_id}.", flush=True)
    else:
        print(f"[ok] wrote {len(created)} sampled FASTAs for {pair_id}", flush=True)
    return created


def _postprocess_af2_run(out_dir: str):
    # link best model
    rj = os.path.join(out_dir, "ranking_debug.json")
    if os.path.isfile(rj):
        try:
            with open(rj) as f: ranking = json.load(f)
            top = (ranking.get("order") or [None])[0]
            if top:
                candidates = glob(os.path.join(out_dir, f"*{top}*.pdb"))
                if candidates:
                    best = sorted(candidates)[0]
                    link = os.path.join(out_dir, "best_model.pdb")
                    if os.path.islink(link) or os.path.exists(link):
                        os.remove(link)
                    os.symlink(os.path.basename(best), link)
        except Exception:
            pass

    # drop mmseqs env folder
    for d in glob(os.path.join(out_dir, "*_env")):
        subprocess.run(f"rm -rf {shlex.quote(d)}", shell=True, check=False)

    # remove huge AF2 "all" pickles if present (defensive cleanup)
    for pkl in glob(os.path.join(out_dir, "*_all_*alphafold2_*.pickle")):
        try:
            os.remove(pkl)
            print(f"[cleanup] removed {os.path.basename(pkl)}")
        except Exception as e:
            print(f"[cleanup] warn: {pkl}: {e}")

    # compress large AF2 JSONs (scores + PAE) to save space
    for jf in glob(os.path.join(out_dir, "*_scores_rank_*.json")) + \
              glob(os.path.join(out_dir, "*_predicted_aligned_error_v1.json")):
        if not os.path.exists(jf + ".gz"):
            try:
                with open(jf, "rb") as fin, gzip.open(jf + ".gz", "wb", compresslevel=9) as fout:
                    fout.writelines(fin)
                os.remove(jf)
                print(f"[cleanup] compressed {os.path.basename(jf)} → .json.gz")
            except Exception as e:
                print(f"[cleanup] warn: could not compress {jf}: {e}")

def _write_pair_a3m_for_chain(cluster_a3m: str, deep_a3m: str, chain_tag: str, out_path: str) -> bool:
    """
    From a single base alignment (cluster_a3m if given else deep_a3m),
    synthesize a per-chain A3M with the chain as FIRST row in the SAME column space.
    If the chain isn't present in the base A3M, align it to the base query row and project.
    """
    # If your Biopython is older, you can import pairwise2 and use format_alignment instead.

    def _ungap_upper(s: str) -> str:
        return "".join(ch for ch in (s or "") if ch.isalpha()).upper()

    # Load base alignment entries
    base_entries = read_msa(cluster_a3m) if cluster_a3m else read_msa(deep_a3m)
    if not base_entries:
        print(f"[error] Base alignment empty for {out_path}")
        return False

    # Load chain FASTA
    pair_dir = os.path.dirname(os.path.dirname(deep_a3m))  # .../Pipeline/<pair>
    chain_fa = os.path.join(pair_dir, "fasta_chain_files", f"{chain_tag}.fasta")
    ids, seqs = load_fasta(chain_fa)
    if not seqs or not seqs[0]:
        print(f"[error] No sequence in {chain_fa} for {chain_tag}")
        return False
    chain_seq = "".join(ch for ch in seqs[0].strip() if ch.isalpha())
    chain_key = chain_seq.upper()

    # Try to find an existing aligned row for this chain (by ungapped sequence)
    idx = {}
    for nm, aln in base_entries:
        idx.setdefault(_ungap_upper(aln), []).append((nm, aln))
    if chain_key in idx:
        # perfect: chain already present in base A3M
        _, chain_aln = idx[chain_key][0]
    else:
        # Need to synthesize an aligned row:
        # 1) get base "query" row (row 0) from A3M
        base_q_aln = base_entries[0][1]                # aligned with gaps
        base_q_seq = _ungap_upper(base_q_aln)          # ungapped
        # 2) align chain_seq <-> base_q_seq (global)
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        # Reasonable scores for protein global alignment (tweak if needed)
        aligner.match_score = 1.0
        aligner.mismatch_score = -1.0
        aligner.open_gap_score = -2.0
        aligner.extend_gap_score = -0.5
        aln = max(aligner.align(base_q_seq, chain_seq), key=lambda a: a.score)
        a_base, a_chain = str(aln).split("\n")[0:2]  # aligned strings without gaps marking?
        # The PairwiseAligner string format is not fixed; safer to use aligned coordinates:
        # Build projected chain row column-by-column against base_q_aln
        chain_aln_list = []
        # Map positions in base_q_aln (with gaps) to positions in base_q_seq (ungapped)
        bpos = 0  # position in base_q_seq
        # Prepare an iterator over (ops) from the alignment to know when chain advances
        # We'll reconstruct a step-function over base_q_seq indices indicating gap/match
        # Using aligned coordinates API:
        b_blocks, c_blocks = aln.aligned  # arrays of (start, end) blocks for base and chain
        # Create an array over base_q_seq length telling for each bpos whether chain has a residue or gap
        cover = np.zeros(len(base_q_seq), dtype=np.int8)  # 1 = has residue, 0 = gap in chain
        c_map = {}  # map bpos -> chain residue index (for residue retrieval)
        cpos = 0
        for (bs, be), (cs, ce) in zip(b_blocks, c_blocks):
            # gap in chain before this block:
            # bs - bpos bases correspond to gaps
            # now filled block:
            for k in range(bs, be):
                cover[k] = 1
            # advance pointers (we'll compute c indices on the fly)
        # Now reconstruct chain_aln by scanning base_q_aln columns:
        # We need the chain residue sequence in order along covered positions
        # Extract chain-aligned residues from alignment object
        chain_aligned_residues = []
        for (bs, be), (cs, ce) in zip(b_blocks, c_blocks):
            chain_aligned_residues.extend(list(chain_seq[cs:ce]))
        cair = iter(chain_aligned_residues)

        for ch in base_q_aln:
            if ch == '-':
                chain_aln_list.append('-')
            else:
                # this column corresponds to base_q_seq[bpos]
                if cover[bpos]:
                    chain_aln_list.append(next(cair, '-'))
                else:
                    chain_aln_list.append('-')
                bpos += 1
        chain_aln = "".join(chain_aln_list)


    # Write new A3M: chain first, then the base entries (skip duplicate of same ungapped seq)
    # --- NEW: drop columns where the chain has a gap so the query (first row) is ungapped ---
    keep = [ch.isalpha() for ch in chain_aln]  # keep only letters in the query row

    def _filter_cols(s: str) -> str:
        # also drop '.' if present; uppercase everything for sanity
        return "".join(ch.upper() for ch, k in zip(s, keep) if k and ch != '.')

    chain_aln_nogap = _filter_cols(chain_aln)
    base_entries_f = [(nm, _filter_cols(aln)) for (nm, aln) in base_entries]

    # Write new A3M: chain first, then filtered base entries (skip exact duplicate)
    ensure_dir(os.path.dirname(out_path))

    with open(out_path, "w") as fh:
        fh.write(f">{chain_tag}\n{chain_aln_nogap}\n")
        chain_key_up = _ungap_upper(chain_aln_nogap)
        for nm, aln in base_entries_f:
            if _ungap_upper(aln) == chain_key_up:
                continue
            fh.write(f">{nm}\n{aln}\n")
    return True


# ------------------------- tasks -------------------------
def task_clean(pair_id: str, args: argparse.Namespace) -> None:

    dry   = _bool_from_tf(getattr(args, "clean_dry_run", "TRUE"))
    level = getattr(args, "clean_level", "derived")
    pair_dir = Path(f"Pipeline/{pair_id}")

    if not pair_dir.exists():
        print(f"[clean] skip (missing): {pair_id}")
        return

    if level == "all":
        print(f"[clean:ALL] rm -rf {pair_dir}")
        if not dry:
            shutil.rmtree(pair_dir, ignore_errors=True)
        return

    # ===== derived-only (keep base inputs) =====
    # Current outputs
    rm_dirs_current = [
        "output_get_msa",
        "output_msa_cluster",
        "output_AF",
        "output_cmaps",
        "output_esm_fold",
        "output_phytree",
        "Analysis",
        # temp / logs
        "tmp_msa_files",
        "tmp_esmfold",
        "jobs",
        "logs",
    ]
    # Legacy outputs from old runs
    rm_dirs_legacy = [
        "output_cmap_esm",
        "esm_cmap_output",
        "AF_preds",
        "chain_pdb_files",
        "fasta_chain_files",
    ]
    # Root-level junk logs (old & new)
    rm_globs = [
        "run_pipeline_for_*.out",
        "RunAF.out",
        "CmapESM.out",
        "*.out",
        "*.log",
        "*.err",
    ]

    # remove known directories
    for d in rm_dirs_current + rm_dirs_legacy:
        p = pair_dir / d
        if p.exists():
            print(f"[clean] rm -rf {p}")
            if not dry:
                shutil.rmtree(p, ignore_errors=True)

    # remove any other tmp_* dirs at the top level
    for p in pair_dir.glob("tmp_*"):
        if p.is_dir():
            print(f"[clean] rm -rf {p}")
            if not dry:
                shutil.rmtree(p, ignore_errors=True)

    # remove root-level logs
    for pat in rm_globs:
        for f in pair_dir.glob(pat):
            if f.is_file():
                print(f"[clean] rm {f}")
                if not dry:
                    try: f.unlink()
                    except: pass

    # prune stale FASTA files:
    # keep only chain-specific FASTA (e.g., 2pbkB.fasta, 3njqA.fasta)
    # remove unsuffixed FASTA (e.g., 2pbk.fasta) if chain FASTA exist
    m = re.match(r"^([0-9a-zA-Z]{4})([A-Za-z])_([0-9a-zA-Z]{4})([A-Za-z])$", pair_id)
    expected = set()
    if m:
        pdb1, ch1, pdb2, ch2 = m.groups()
        expected = {f"{pdb1}{ch1}.fasta", f"{pdb2}{ch2}.fasta"}
    fasta_files = list(pair_dir.glob("*.fasta"))
    have_chain_fastas = expected and all((pair_dir / e).exists() for e in expected)
    for f in fasta_files:
        name = f.name
        if have_chain_fastas and name not in expected:
            print(f"[clean] rm stale FASTA: {f}")
            if not dry:
                try: f.unlink()
                except: pass

    # keep base inputs:
    #   Pipeline/<pair>/*.pdb, *_cif.pdb, *.fasta (chain-specific), *_pdb_contacts.npy
    print(f"[clean] done: {pair_id}")


def task_load(pair_id: str, args: argparse.Namespace) -> None:
    foldA, foldB = pair_str_to_tuple(pair_id)
    cur_family_dir = f"Pipeline/{pair_id}"
    ensure_dir(cur_family_dir)
    force = _bool_from_tf(getattr(args, "force_rerun", "FALSE"))
    load_seq_and_struct(
        cur_family_dir,
        [foldA[:-1], foldB[:-1]],
        [foldA[-1],  foldB[-1]],
        force_rerun=force)

def task_get_msa(pair_id: str, run_job_mode: str) -> None:
    """
    Build a 2-sequence seed A3M from BOTH chains of the pair, then run get_msa to
    produce Pipeline/<pair>/output_get_msa/DeepMsa.a3m.
    """
    # Ensure the two chain FASTAs exist (reuses root FASTAs or extracts from PDBs)
    pair_dir = f"Pipeline/{pair_id}"
    foldA, foldB = pair_str_to_tuple(pair_id)        # e.g. '1dzlA', '5keqF'
    pdbids     = [foldA[:-1], foldB[:-1]]            # ['1dzl','5keq']
    pdbchains  = [foldA[-1],  foldB[-1]]             # ['A','F']
    ensure_dir(os.path.join(pair_dir, "fasta_chain_files"))
    ensure_chain_fastas(pair_dir, pdbids, pdbchains)

    # Build BOTH-chains seed alignment once
    seed_a3m = build_pair_seed_a3m_from_pair(pair_id, data_dir="Pipeline")
    out_dir  = f"{pair_dir}/output_get_msa"
    ensure_dir(out_dir)

    # Prefer sbatch wrapper if it exists; else run python inline
    sbatch_script = "./Pipeline/get_msa_params.sh"
    if run_job_mode == "sbatch" and os.path.exists(sbatch_script):
        log = f"{pair_dir}/get_msa_for_{pair_id}.out"
        cmd = f"sbatch -o '{log}' {sbatch_script} '{seed_a3m}' {pair_id}"
        _run(cmd, "sbatch")
    else:
        cmd = f"python3 ./get_msa.py '{seed_a3m}' '{out_dir}' --pair {pair_id}"

        _run(cmd, "inline")


def task_cluster_msa(pair_id: str, run_job_mode: str) -> None:
    # Use the *current* Python interpreter so the same venv is used downstream
    py  = shlex.quote(sys.executable)
    pid = shlex.quote(pair_id)
    cmd = (
        f"bash -lc 'cd Pipeline/{pid} && "
        f"{py} ../../ClusterMSA_moriah.py "
        f"--keyword ShallowMsa "
        f"-i output_get_msa/DeepMsa.a3m "
        f"-o output_msa_cluster'"
    )
    _run(cmd, run_job_mode)


def task_cmap_msa_transformer(pair_id: str, run_job_mode: str) -> None:
    outdir = f"Pipeline/{pair_id}/output_cmaps/msa_transformer"
    ensure_dir(outdir)
    cmd = (
        f"python3 ./runMSATrans.py "
        f"--input_msas Pipeline/{pair_id}/output_msa_cluster "
        f"-o {outdir}"
    )
    _run(cmd, run_job_mode)


def task_cmap_ccmpred(pair_id: str, run_job_mode: str) -> None:
    """
    Run CCMpred on DeepMsa and on every ShallowMsa_XXX in output_msa_cluster.
    Outputs: Pipeline/<pair>/output_cmaps/ccmpred/<tag>.ccmpred.npy (APC-corrected)
    """

    pair_dir = Path(f"Pipeline/{pair_id}")
    out_dir = pair_dir / "output_cmaps" / "ccmpred"
    tmp_dir = pair_dir / "tmp_ccmpred"
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    ccmpred_bin = getattr(args, "ccmpred_bin", "/sci/labs/orzuk/orzuk/github/CCMpred/bin/ccmpred")
    threads = int(getattr(args, "ccmpred_threads", 8))

    def a3m_to_fa(a3m_path: Path, fasta_out: Path) -> bool:
        # read_msa returns [(id, seq)], where seq may contain lowercase inserts for A3M.
        entries = read_msa(str(a3m_path))
        if not entries:
            return False
        names, seqs = zip(*entries)
        clean = []
        for s in seqs:
            s = re.sub(r"[a-z.]", "", s)  # strip inserts; keep '-' gaps
            s = s.replace("U", "C")       # CCMpred doesn't know 'U'
            clean.append(s)
        Ls = {len(s) for s in clean}
        if len(Ls) != 1:
            # ragged columns? pad with gaps just in case
            L = max(Ls)
            clean = [s.ljust(L, "-") for s in clean]
        with open(fasta_out, "w") as f:
            for i, (n, s) in enumerate(zip(names, clean), 1):
                f.write(f">{n or f'seq{i}'}\n{s}\n")
        return True

    def apc(m: np.ndarray) -> np.ndarray:
        ri = m.mean(axis=1, keepdims=True)
        rj = m.mean(axis=0, keepdims=True)
        mu = float(m.mean()) or 1e-8
        return m - (ri @ rj) / mu

    def run_one(a3m: Path, tag: str):
        fa = tmp_dir / f"{tag}.fa"
        npy = out_dir / f"{tag}.ccmpred.npy"
        mat = out_dir / f"{tag}.ccmpred.mat"
        if npy.exists() and npy.stat().st_size > 0:
            return
        if not a3m_to_fa(a3m, fa):
            print(f"[ccmpred] skip (empty): {a3m}")
            return
        cmd = f"{shlex.quote(ccmpred_bin)} -t {threads} {shlex.quote(str(fa))} {shlex.quote(str(mat))}"
        if run_job_mode == "inline":
            subprocess.run(cmd, shell=True, check=True)
        else:
            # sbatch wrapper uses the same helper you use elsewhere
            _run(cmd, run_job_mode)
        arr = np.loadtxt(str(mat))
        arr = 0.5 * (arr + arr.T)     # symmetrize, just in case
        arr = apc(arr)                # APC correction
        np.fill_diagonal(arr, 0.0)
        np.save(str(npy), arr)

    # Deep
    deep = pair_dir / "output_get_msa" / "DeepMsa.a3m"
    if deep.exists():
        run_one(deep, "DeepMsa")

    # All clusters
    for a3m in sorted((pair_dir / "output_msa_cluster").glob("ShallowMsa_*.a3m")):
        run_one(a3m, a3m.stem)



def task_esmfold(pair_id: str, args: argparse.Namespace) -> None:
    """
    - Sample up to N seqs per cluster -> _tmp_ShallowMsa_XXX.fasta
    - Call ESMFoldHF once with correct CLI:
        python3 ./ESMFoldHF.py -input <PAIR_ID> --model {esm2|esm3} [--device ...]
    - Outputs go to: Pipeline/<pair>/output_esm_fold/<esm_model>/
    """
    if _is_windows():
        raise SystemExit("ESMFold can’t run on Windows. Run on Moriah/Linux.")

    # Prepare sampled FASTAs for clusters (and skip deep by default here)
    _sample_to_tmp_fastas(pair_id, args.cluster_sample_n, include_deep=False)

    # Ensure output dir exists (the HF script will place files under the model subdir)
    model_dir = f"Pipeline/{pair_id}/output_esm_fold/{args.esm_model}"
    ensure_dir(model_dir)

    # Correct CLI: no -i/-o
    device = args.esm_device or "auto"
    cmd = f"python3 ./ESMFoldHF.py -input {pair_id} --model {args.esm_model} --device {device}"
    _run(cmd, args.run_job_mode)

    # No cleanup!
    # tmpdir = f"Pipeline/{pair_id}/tmp_esmfold"
    # subprocess.run(f"rm -rf {shlex.quote(tmpdir)}", shell=True, check=False)


def task_af(pair_id: str, args: argparse.Namespace) -> None:
    if _is_windows():
        raise SystemExit("AlphaFold must run on Linux.")

    af_ver = str(getattr(args, "af_ver", "both")).lower()  # default: do both

    pair_dir = f"Pipeline/{pair_id}"
    foldA, foldB = pair_str_to_tuple(pair_id)
    pdbids = [foldA[:-1], foldB[:-1]]
    pdbchains = [foldA[-1],  foldB[-1]]
    chains = [foldA, foldB]
    ensure_chain_fastas(pair_dir, pdbids, pdbchains)

    deep_a3m = os.path.join(pair_dir, "output_get_msa", "DeepMsa.a3m")
    if not os.path.isfile(deep_a3m):
        print(f"[err] Missing DeepMsa: {deep_a3m}")
        return

    cluster_dir = os.path.join(pair_dir, "output_msa_cluster")
    cluster_a3ms = sorted(glob(os.path.join(cluster_dir, "ShallowMsa_*.a3m")))

    tmp_pairs_dir = os.path.join(pair_dir, "tmp_msa_files")
    ensure_dir(tmp_pairs_dir)

    # Build pair-specific A3Ms per chain for DeepMsa + each cluster (same as before)
    jobs = []  # (a3m_path, out_dir_base_name) — out_dir_base_name = "DeepMsa/<chain>" or "ShallowMsa_XXX/<chain>"
    for ch in chains:
        pair_a3m = os.path.join(tmp_pairs_dir, f"tmp_DeepMsa__{ch}.a3m")
        if _write_pair_a3m_for_chain(None, deep_a3m, ch, pair_a3m):
            jobs.append((pair_a3m, f"DeepMsa/{ch}"))

    for a3m in cluster_a3ms:
        cl_stem = Path(a3m).stem  # e.g. ShallowMsa_007
        for ch in chains:
            pair_a3m = os.path.join(tmp_pairs_dir, f"tmp_{cl_stem}__{ch}.a3m")
            if _write_pair_a3m_for_chain(a3m, deep_a3m, ch, pair_a3m):
                jobs.append((pair_a3m, f"{cl_stem}/{ch}"))

    if not jobs:
        print(f"[warn] No AF jobs to run for {pair_id}")
        return

    inside_slurm = _in_slurm_session()

    # New (respect CLI defaults)
    gres = getattr(args, "sbatch_gres", "gpu:1")
    cpus = int(getattr(args, "sbatch_cpus", 8))
    mem = getattr(args, "sbatch_mem", "40G")
    time = getattr(args, "sbatch_time", "24:00:00")
    part = getattr(args, "sbatch_partition", None)
    cons = getattr(args, "sbatch_constraint", None)
    sbatch_opts = f"--gres={gres} --cpus-per-task={cpus} --mem={mem} --time={time}"
    if part: sbatch_opts += f" -p {shlex.quote(part)}"
    if cons: sbatch_opts += f" --constraint={shlex.quote(cons)}"


    def _cmd_for(ver: str, a3m_path: str, out_dir: str) -> str:
        if ver == "2":
            return (
                f"bash ./Pipeline/RunAF2_Colabfold.sh "
                f"{shlex.quote(a3m_path)} {shlex.quote(out_dir)} "
                f"--num-models 5 --num-recycle 1 --model-type alphafold2_ptm"
            )  # do not save large pickle files
        elif ver == "3":
            # AF3 runner converts A3M->JSON and runs AF3; also export top PDB
            return (
                f"bash ./Pipeline/RunAF3_Colabfold.sh "
                f"{shlex.quote(a3m_path)} {shlex.quote(out_dir)} "
                f"--pdb=rank1"
            )
        else:
            raise ValueError(ver)

    # Run AF2, AF3, or both
    force = _bool_from_tf(getattr(args, "force_rerun_AF", "FALSE"))

    versions = ["2"] if af_ver == "2" else (["3"] if af_ver == "3" else ["2", "3"])
    for ver in versions:
        out_root = os.path.join(pair_dir, f"output_AF/AF{ver}")
        ensure_dir(out_root)
        log_dir = os.path.join(out_root, "logs")
        ensure_dir(log_dir)

        for a3m_path, base in jobs:
            out_dir = os.path.join(out_root, base)
            ensure_dir(out_dir)

            # ---------- SKIP PATH (already computed) ----------
            if not force:
                if ver == "2" and _af2_has_outputs(out_dir):
                    print(f"[skip] AF2 exists → {out_dir}")
                    # NEW: ensure best_model.pdb symlink exists so canonical export can find it
                    try:
                        _postprocess_af2_run(out_dir)  # creates/refreshes best_model.pdb from ranking_debug.json
                    except Exception as e:
                        print(f"[warn] _postprocess_af2_run failed in skip-path: {e}")
                    continue
                if ver == "3" and _af3_has_outputs(out_dir):
                    print(f"[skip] AF3 exists → {out_dir}")
                    # NEW: ensure there is a PDB for existing CIF results
                    try:
                        _convert_existing_af3(out_dir, mode="rank1")  # "all" if you want every ranked PDB
                    except Exception as e:
                        print(f"[warn] _convert_existing_af3 failed in skip-path: {e}")
                    continue

            # ---------- RUN PATH ----------
            cmd = _cmd_for(ver, a3m_path, out_dir)
            if args.run_job_mode == "sbatch" or (not inside_slurm and not getattr(args, "allow_inline_af", False)):
                stem = f"{Path(a3m_path).stem}"
                log_path = os.path.join(log_dir, f"run_AF{ver}_{pair_id}__{stem}.out")
                _run(f"sbatch {sbatch_opts} -o '{log_path}' --wrap {shlex.quote(cmd)}", "sbatch")
            else:
                _run(cmd, "inline")
                # NEW: after an inline run, immediately prepare artifacts needed for canonical export
                if ver == "2":
                    try:
                        _postprocess_af2_run(out_dir)
                    except Exception as e:
                        print(f"[warn] _postprocess_af2_run failed after inline AF2: {e}")
                elif ver == "3":
                    try:
                        _convert_existing_af3(out_dir, mode="rank1")
                    except Exception as e:
                        print(f"[warn] _convert_existing_af3 failed after inline AF3: {e}")

    # ---------- CANONICAL EXPORT (flat files) ----------
    # After running (or skipping) all jobs per version, create flat canonical copies.
    for ver in versions:
        _export_canonical_best_pdbs(pair_id, ver)


def task_tree(pair_id: str, run_job_mode: str) -> None:
    msa_file = f"Pipeline/{pair_id}/output_get_msa/DeepMsa.a3m"
    out = f"Pipeline/{pair_id}/output_phytree/DeepMsa_tree.nwk"
    ensure_dir(os.path.dirname(out))
    phytree_from_msa(msa_file, output_tree_file=out)

def task_plot(pair_id: str, args: argparse.Namespace) -> None:
    # Import PyMOL-consuming code only here
    foldA, foldB = pair_str_to_tuple(pair_id)
    pdbids = [foldA[:-1], foldB[:-1]]
    pdbchains = [foldA[-1], foldB[-1]]
    make_foldswitch_all_plots(pdbids, "Pipeline", pair_id, pdbchains,
                              plot_tree_clusters=args.plot_trees,
                              plot_contacts=True, global_plots=args.global_plots)

def task_deltaG(pair_id: str) -> None:
    # Import PyRosetta-consuming code only here
    try:
        from utils.energy_utils import compute_global_and_residue_energies
    except ImportError:
        raise SystemExit("PyRosetta utilities are unavailable in this env.")
    out_dir = "Pipeline/output_deltaG"
    ensure_dir(out_dir)
    pA, pB = pair_str_to_tuple(pair_id)
    pdb_pair = [(f"Pipeline/{pair_id}/{pA[:-1]}.pdb", f"Pipeline/{pair_id}/{pB[:-1]}.pdb")]
    compute_global_and_residue_energies(pdb_pair, [pair_id], out_dir)


def task_postprocess(foldpairs: list[str], args: argparse.Namespace) -> None:
    """
    Unified post-processing/reports:
      1) per-pair postprocess_unified (df_af/df_cmap/df_esm etc.)
      2) global CSVs (summary + detailed) and their HTML tables
      3) per-pair HTML pages
    Controlled by --reports: none | tables | html | all
    Safe to run incrementally.
    """
    # Normalize pairs to strings like 1wp8C_5ejbC
    norm_pairs = [
        p if isinstance(p, str) else f"{p[0]}_{p[1]}"
        for p in (foldpairs if isinstance(foldpairs, list) else [foldpairs])
    ]

    # 1) Per-pair metrics
    try:
        force = _bool_from_tf(getattr(args, "force_rerun_postprocess", "FALSE"))
        post_processing_analysis(force_rerun=force, pairs=norm_pairs)  # pass list, or None for discover-all
    except Exception as e:
        print(f"[postprocess] WARN post_processing_analysis: {e}")

    # 2) Global CSVs + HTML tables
    if args.reports in ("tables", "all"):
        try:
            collect_summary_tables(DATA_DIR, DETAILED_RESULTS_TABLE, SUMMARY_RESULTS_TABLE)
        except Exception as e:
            print(f"[reports] WARN collect_summary_tables: {e}")

        try:
            gen_html_from_summary_table()
        except Exception as e:
            print(f"[reports] WARN gen_html_from_summary_table: {e}")

        try:
            gen_html_from_cluster_detailed_table()
        except Exception as e:
            print(f"[reports] NOTE cluster-detailed HTML skipped: {e}")

        # After building the HTML tables, mirror the main table to repo root for GitHub Pages
        try:
            src = os.path.join(TABLES_RES, "table.html")
            dst = os.path.join(MAIN_DIR, "table.html")
            if os.path.isfile(src):
                shutil.copy2(src, dst)
                print(f"[reports] copied {src} -> {dst}")
        except Exception as e:
            print(f"[reports] WARN copying table.html to repo root: {e}")

    # 3) Per-pair HTML notebook pages
    if args.reports in ("html", "all"):
        try:
            # Keep only pairs that already have per-pair postprocess outputs
            ready = []
            for p in norm_pairs:
                if os.path.isfile(f"Pipeline/{p}/Analysis/df_af.csv"):
                    ready.append(p)
                else:
                    print(f"[html] skip {p}: missing Pipeline/{p}/Analysis/df_af.csv")

            if not ready:
                print("[html] no ready pairs; skipping notebook generation")
            else:
                pairs_arg = " ".join(ready)
                cmd = f"{shlex.quote(sys.executable)} Analysis/NotebookGen/generate_notebooks.py {pairs_arg} --kernel python3"
                subprocess.run(cmd, shell=True, check=True, env=_jupyter_env_for_scratch())
        except Exception as e:
            print(f"[reports] WARN per-pair HTML generation: {e}")


    print("[postprocess] done.", flush=True)


# All Pipeline
def task_msaclust_pipeline(pair_id: str, args: argparse.Namespace) -> None:
    """
    Full per-pair pipeline:
      load → get_msa → cluster_msa → AF (AF2/AF3) → cmaps (msa-transformer, ccmpred)
      → ESMFold → tree → ΔG → plots → per-pair postprocess (df_*.csv) → (optional) per-pair HTML.
    If --force_rerun TRUE: force all recomputations; otherwise each step skips when outputs exist.
    """
    force_all = _bool_from_tf(getattr(args, "force_rerun", "FALSE"))

    # 1) load PDB/FASTA (cheap; idempotent)
    task_load(pair_id, args)

    # 2) get_msa
    if force_all or not _has_deep_msa(pair_id):
        print("[pipeline] get_msa → running")
        task_get_msa(pair_id, "inline")
    else:
        print("[pipeline] get_msa → skip (DeepMsa exists)")

    # 3) cluster_msa
    if force_all or not _has_cluster_msas(pair_id):
        print("[pipeline] cluster_msa → running")
        task_cluster_msa(pair_id, "inline")
    else:
        print("[pipeline] cluster_msa → skip (clusters exist)")

    # 4) AF (AF2/AF3/both). Always invoke; the task itself skips completed outdirs unless forced.
    if force_all:
        setattr(args, "force_rerun_AF", "TRUE")  # only for AF stage
    print("[pipeline] AF → running (task will skip per-outdir if already complete)")
    # Default to both versions inside the full pipeline unless user explicitly asked AF2-only.
    if getattr(args, "af_ver", None) not in ("2", "3", "both"):
        args.af_ver = "both"
    task_af(pair_id, args)

    # 5) cmaps (MSA-Transformer)
    if force_all or not _has_cmaps(pair_id):
        print("[pipeline] cmap_msa_transformer → running")
        task_cmap_msa_transformer(pair_id, "inline")
    else:
        print("[pipeline] cmap_msa_transformer → skip (cmaps exist)")

    # 5b) cmaps (CCMpred) — best-effort
    try:
        print("[pipeline] cmap_ccmpred → running")
        task_cmap_ccmpred(pair_id, "inline")
    except Exception as e:
        print(f"[pipeline] cmap_ccmpred → skipped: {e}")

    # 6) ESMFold (esm2/esm3 or user-specified)
    wanted_models = ["esm2", "esm3"] if getattr(args, "esm_model", None) in (None, "both") else [args.esm_model]
    for model in wanted_models:
        a2 = deepcopy(args)
        a2.esm_model = model
        if force_all or not _has_esm_model(pair_id, model):
            print(f"[pipeline] esmfold({model}) → running")
            task_esmfold(pair_id, a2)
        else:
            print(f"[pipeline] esmfold({model}) → skip (outputs exist)")

    # 7) phylogenetic tree (needed for tree plots)
    tree_path = Path(f"Pipeline/{pair_id}/output_phytree/DeepMsa_tree.nwk")
    if force_all or not tree_path.exists():
        print("[pipeline] tree → running")
        task_tree(pair_id, "inline")
    else:
        print("[pipeline] tree → skip (exists)")

    # 8) ΔG energies (PyRosetta) — best-effort
    try:
        if force_all or not _has_deltaG(pair_id):
            print("[pipeline] deltaG → running")
            task_deltaG(pair_id)
        else:
            print("[pipeline] deltaG → skip (exists)")
    except Exception as e:
        print(f"[deltaG] skipped: {e}")

    # 9) plots (includes tree clusters if available)
    try:
        ap = deepcopy(args)
        ap.global_plots = False
        ap.plot_trees = True
        print("[pipeline] plot → running")
        task_plot(pair_id, ap)
    except Exception as e:
        print(f"[plot] skipped: {e}")

    # 10) per-pair post-process (writes Pipeline/<pair>/Analysis/df_*.csv)
    try:
        print("[pipeline] postprocess → running (per-pair)")
        # Always compute for this pair; internal readers still skip unchanged unless forced.
        post_processing_analysis(force_rerun=force_all or _bool_from_tf(getattr(args, "force_rerun_postprocess", "FALSE")),
                                 pairs=[pair_id])
    except Exception as e:
        print(f"[postprocess-inline] WARN: {e}")


    # 11) OPTIONAL: per-pair HTML (execute notebook for this pair inside the same job)
    if _bool_from_tf(getattr(args, "per_pair_html", "FALSE")):
        try:
            print("[pipeline] html → running (per-pair)")
            env = _jupyter_env_for_scratch()  # ← use scratch env
            cmd = (
                f"{shlex.quote(sys.executable)} Analysis/NotebookGen/generate_notebooks.py "
                f"{shlex.quote(pair_id)} --kernel {shlex.quote(getattr(args, 'per_pair_kernel', 'python3'))}"
            )
            # make HTML failures non-fatal so the job completes
            rc = subprocess.run(cmd, shell=True, check=False, env=env).returncode
            if rc != 0:
                print(f"[html] nbconvert failed (rc={rc}) — will be retried later.")
        except Exception as e:
            print(f"[html] skipped: {e}")
    else:
        print("[pipeline] html → skip (per-pair HTML disabled)")


# ------------------------- CLI / main -------------------------

def main():
    p = argparse.ArgumentParser(
#        description="Fold-switching pipeline runner",
        formatter_class=argparse.RawTextHelpFormatter,  # keeps newlines
        epilog=_modes_epilog(),
    )
    p.add_argument("--run_mode",
                   required=True,
                   choices=["load", "get_msa", "cluster_msa", "run_cmap_msa_transformer", "run_cmap_ccmpred",
                            "run_esmfold", "run_AF", "tree", "plot", "compute_deltaG", "clean",
                            "postprocess", "msaclust_pipeline", "help"])  # Last one is the full pipeline for a pair
    p.add_argument("--foldpair_ids", nargs="+", required=True,
                   help="List of pair IDs (e.g. 1dzlA_5keqF), or the literal token ALL")

    p.add_argument("--run_job_mode", default="inline", choices=["inline", "sbatch"])

    # AlphaFold options
    p.add_argument("--allow_inline_af", action="store_true",
                    help="Allow AF2 to run inline even if not in a Slurm session (expert only).")
    p.add_argument("--af_ver", default="both", choices=["2", "3", "both"],
                    help="Which AlphaFold to run for --run_mode run_AF")  # default do both AF2 and AF3

    p.add_argument("--force_rerun_AF", default="FALSE", choices=["TRUE", "FALSE"],
                    help="Run new AF2/AF3 predictions even if outputs exist. Default FALSE (skip if found).")

    # ESMFold options
    p.add_argument("--cluster_sample_n", type=int, default=10)
    p.add_argument("--esm_model", default=None, choices=["esm2", "esm3", "both"])
    p.add_argument("--esm_device", default="auto", choices=["auto", "cpu", "cuda", "mps"])

    # ---- CCMpred options ----
    p.add_argument("--ccmpred_bin",
                   default="/sci/labs/orzuk/orzuk/github/CCMpred/bin/ccmpred",
                   help="Path to CCMpred binary")
    p.add_argument("--ccmpred_threads", type=int, default=8,
                   help="Threads for CCMpred (-t)")

    # Post-processing options, computing metrics
    p.add_argument("--postprocess", default="FALSE", help="If TRUE, run post-processing after the selected task(s).")
    p.add_argument("--force_rerun_postprocess", default="FALSE", help="If TRUE, recompute per-pair caches.")
    p.add_argument("--per_pair_html", default="TRUE")
    p.add_argument("--per_pair_kernel", default="python3")

    # Plotting
    p.add_argument("--global_plots", action="store_true")
    p.add_argument("--plot_trees", action="store_true")

    # Output html tables and pages
    p.add_argument("--reports", default="none", choices=["none", "tables", "html", "all"],
                   help="Post-run reporting: 'tables' builds CSV+HTML tables, 'html' builds per-pair pages, 'all' does both.")
    p.add_argument("--html_pairs", nargs="+", default=["ALL"],
                   help="Pairs to render per-pair HTML for (defaults to ALL).")


    # Pipeline-wide force flag
    p.add_argument("--force_rerun", default="FALSE", choices=["TRUE", "FALSE"],
        help="If TRUE (only for --run_mode msaclust_pipeline), run every step regardless of outputs.")

    p.add_argument("--clean_level", default="derived", choices=["derived", "all"],
        help="derived: remove computed outputs but keep base inputs; all: remove the entire pair folder")
    p.add_argument("--clean_dry_run", default="TRUE", choices=["TRUE", "FALSE"],
        help="TRUE: only print what would be removed; FALSE: actually delete")

    # Slurm resource knobs for per-pair jobs (defaults safe for your cluster)
    p.add_argument("--sbatch_gres", default="gpu:1",
                   help="Slurm --gres for per-pair job (e.g., gpu:1 or gpu:a100:1)")
    p.add_argument("--sbatch_cpus", type=int, default=8,
                   help="Slurm --cpus-per-task for per-pair job")
    p.add_argument("--sbatch_mem", default="40G",
                   help="Slurm --mem for per-pair job")
    p.add_argument("--sbatch_time", default="24:00:00",
                   help="Slurm --time for per-pair job")
    p.add_argument("--sbatch_partition", default=None,
                   help="Optional Slurm -p/--partition")
    p.add_argument("--sbatch_constraint", default=None,
                   help="Optional Slurm --constraint (e.g., a100, v100)")
    p.add_argument("--sbatch_account", default=None,
                   help="Optional Slurm -A/--account")
    p.add_argument("--sbatch_qos", default=None,
                   help="Optional Slurm --qos")
    p.add_argument("--sbatch_mail", default=None,
                   help="Optional --mail-user email for notifications")
    p.add_argument("--sbatch_mail_type", default=None,
                   help="Optional --mail-type (e.g., END,FAIL,ALL)")


    args = p.parse_args()
    # allow: python run_foldswitch_pipeline.py --run_mode help
    if args.run_mode == "help":
        p.print_help()
        sys.exit(0)
    if args.run_mode == "run_esmfold" and args.esm_model is None:
        args.esm_model = "esm2"  # sensible default

    # resolve list of pairs
    if args.foldpair_ids == ["ALL"]:
        # read the same list your current script reads
#        with open("data/foldswitch_PDB_IDs_full.txt", "r") as f:
#            raw = [line.rstrip("\n") for line in f if line.strip()]
        foldpairs = list_protein_pairs() # [s.replace("\t", "_") for s in raw]
    else:
        foldpairs = args.foldpair_ids

    # --- Before the loop ---
    if args.run_mode == "postprocess":
        # Pass the full list so the function can filter to pairs that are ready
        task_postprocess(foldpairs, args)
        return  # or sys.exit(0)


    for pair_id in foldpairs:
        print(f"=== {args.run_mode} :: {pair_id} ===", flush=True)

        if type(pair_id) in [tuple, list]:  # convert format
            pair_id = f"{pair_id[0]}_{pair_id[1]}"

        if args.run_mode == "load":
            task_load(pair_id, args)

        elif args.run_mode == "get_msa":
            task_get_msa(pair_id, args.run_job_mode)

        elif args.run_mode == "cluster_msa":
            task_cluster_msa(pair_id, args.run_job_mode)

        elif args.run_mode == "run_cmap_msa_transformer":
            task_cmap_msa_transformer(pair_id, args.run_job_mode)

        elif args.run_mode == "run_cmap_ccmpred":
            task_cmap_ccmpred(pair_id, args.run_job_mode)

        elif args.run_mode == "run_esmfold":
            if args.esm_model == "both":
                for m in ("esm2", "esm3"):
                    a = deepcopy(args)
                    a.esm_model = m
                    task_esmfold(pair_id, a)
            else:
                task_esmfold(pair_id, args)

        elif args.run_mode == "run_AF":
            task_af(pair_id, args)

        elif args.run_mode == "tree":
            task_tree(pair_id, args.run_job_mode)

        elif args.run_mode == "compute_deltaG":
            task_deltaG(pair_id)

        elif args.run_mode == "plot":
            task_plot(pair_id, args)


        elif args.run_mode == "clean":  # Remove existing files to run the pipeline clean
            task_clean(pair_id, args)

        elif args.run_mode == "msaclust_pipeline":
            if args.run_job_mode == "sbatch":
                # Submit one Slurm job per pair; each job runs the whole pipeline inline.
                _submit_msaclust_pair_job(pair_id, args)
            else:
                # Run the full pipeline inline for this pair (exactly as before)
                task_msaclust_pipeline(pair_id, args)

        else:
            raise ValueError(args.run_mode)



    # If we only submitted per-pair jobs, don’t build reports now
    if args.run_mode == "msaclust_pipeline" and args.run_job_mode == "sbatch":
        print("[submit-only] Per-pair jobs have been submitted. Run reports later.", flush=True)
        return

    # Otherwise build whatever the user asked for
    if args.reports != "none":
        task_postprocess(foldpairs, args)
        print("[done html outputs]", flush=True)

if __name__ == "__main__":
    main()
