#!/usr/bin/env python3
import argparse
import subprocess
import shlex
import sys, os, warnings, re, math
from glob import glob
import json, shutil, gzip, getpass

from pathlib import Path
from typing import List, Tuple
from copy import deepcopy
import numpy as np
import pandas as pd
import time

import signal
from contextlib import contextmanager

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

from config import *
from utils.utils import pair_str_to_tuple, ensure_dir, list_protein_pairs, write_pair_pipeline_script, str2bool
from utils.msa_utils import *
from utils.align_utils import get_or_compute_true_tm
try:
    from utils.phytree_utils import phytree_from_msa
except Exception:
    phytree_from_msa = None
    print("[phytree ERROR] Can't import phytree!", flush=True)
from utils.protein_utils import *
from utils.protein_plot_utils import make_foldswitch_all_plots, global_pairs_statistics_plots
try:
    # NB: in repo this module sits under utils/. Keep your original import path.
    from utils.energy_utils import compute_global_and_residue_energies
except ImportError:
    raise SystemExit("[RosettaFold ERROR] PyRosetta utilities are unavailable in this env.")

from Analysis.postprocess_global import build_or_load_global_tables
from Analysis.postprocess_unified import post_processing_analysis, build_unified_tables_from_cluster_dfs
from TableResults.gen_html_table import *

from Bio import Align  # PairwiseAligner (modern replacement)
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1


RUN_MODE_DESCRIPTIONS = {
    "load":             "Load the pair from PDB.",
    "get_msa":          "Download/build deep MSAs for the pair and write DeepMsa.a3m.",
    "cluster_msa":      "Cluster the pairwise MSA into shallow clusters (ShallowMsa_XXX.a3m).",
    "tree":             "Build phylogenetic tree from MSA.",
    "run_AF":           "Run AlphaFold (AF2/AF3/both) per chain × cluster. Use --af_ver {2,3,both}.",
    "run_esmfold":      "Run ESMFold on the pair. Use --esm_model {esm2,esm3,both}.",
    "run_cmap_msa_transformer":      "Run MSA-transformer on the pair to get contact maps.",
    "run_cmap_ccmpred": "Run CCMpred on DeepMsa and all ShallowMsa_XXX clusters to get contact maps.",
    "compute_deltaG":   "Compute ΔG stability metrics (requires PyRosetta).",
    "postprocess": "Compute TM/cmap metrics and build summary/detailed tables. Use --force_rerun_postprocess TRUE to recompute." ,
    "plot":             "Generate pair-specific plots (requires PyMOL).",
    "clean":            "Remove previous outputs for the pair.",
    "msaclust_pipeline":"Full pipeline: get_msa → cluster_msa → AF/ESM (as configured).",
    "help":             "Print this list of run modes with one-line explanations.",
}

# Modes that are time/compute-heavy and parallelize well across pairs
HEAVY_PAIR_MODES = {
    "cluster_msa",            # new: HDBSCAN on MSA
    "run_cmap_msa_transformer",
    "run_cmap_ccmpred",
    "run_esmfold",            # your ESMFold step
    "run_AF",         # if you have this mode
    "compute_deltaG",         # if you have a Rosetta/PyRosetta ΔG step
    "plot",                   # rendering/alignments can be slow
    "gen_pair_html",          # if HTML generation per pair is non-trivial
    "tree",                   # build phylogenetic tree
    # add more as needed
}

# Suggested sbatch resources per mode (tweak to your cluster)
SBATCH_HINTS = {
    "cluster_msa":           {"time":"02:00:00", "mem":"8G",   "gpus":0, "cpus":4},
    "run_cmap_msa_transformer": {"time":"04:00:00", "mem":"8G",   "gpus":0, "cpus":4},
    "run_cmap_ccmpred":      {"time":"04:00:00", "mem":"16G",  "gpus":0, "cpus":8},
    "run_esmfold":           {"time":"06:00:00", "mem":"24G",  "gpus":1, "cpus":6},
    "run_AF":                {"time":"24:00:00","mem":"64G",  "gpus":1, "cpus":8},
    "compute_deltaG":        {"time":"06:00:00", "mem":"16G",  "gpus":0, "cpus":8},
    "plot":                  {"time":"02:00:00", "mem":"8G",   "gpus":0, "cpus":4},
    "gen_pair_html":         {"time":"01:00:00", "mem":"4G",   "gpus":0, "cpus":2},
}



# ------------------------- helpers -------------------------
@contextmanager
def _alarm_timeout(seconds: int, what: str):
    """Unix-only timeout; raises TimeoutError if block exceeds `seconds`."""
    def _handler(signum, frame):
        raise TimeoutError(f"timeout after {seconds}s in {what}")
    old = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(int(seconds))
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)

def _csv_bool(x) -> bool:
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    s = str(x).strip().lower()
    return s in ("true", "1", "yes", "y", "t")

def _renumber_kept_clusters(pair_id: str) -> bool:
    """
    After cluster_msa, keep only clusters with kept==True and renumber them
    to contiguous ShallowMsa_000.. ShallowMsa_{K-1} in:
       - Pipeline/<pair>/output_msa_cluster/
       - Pipeline/<pair>/output_msa_cluster_pairtrim/   (if present)
    Also rewrites ShallowMsa_clusters.csv accordingly.

    Returns True if any changes were made, False if no action was needed.
    """
    base = Path(f"Pipeline/{pair_id}")
    cl_csv = base / "output_msa_cluster" / "ShallowMsa_clusters.csv"
    if not cl_csv.is_file():
        return False

    df = pd.read_csv(cl_csv)
    if "cluster" not in df.columns:
        return False

    # Determine kept set
    if "kept" in df.columns:
        kept_df = df[df["kept"].apply(_csv_bool)].copy()
    else:
        kept_df = df.copy()

    if kept_df.empty:
        # Nothing to keep → remove all shallow .a3m files (both dirs) and write empty CSV
        for dname in ("output_msa_cluster", "output_msa_cluster_pairtrim"):
            d = base / dname
            if d.is_dir():
                for p in d.glob("ShallowMsa_*.a3m"):
                    try: p.unlink()
                    except Exception: pass
        kept_df[["cluster","n","neff","path"]] = [[]]  # ensure columns exist
        kept_df.to_csv(cl_csv, index=False)
        return True

    # If already contiguous 0..K-1 AND filenames are aligned, skip
    try:
        old_ids = [int(c) for c in kept_df["cluster"].tolist()]
    except Exception:
        old_ids = kept_df["cluster"].tolist()  # fallback (strings)

    sorted_old = sorted(old_ids)
    already_contiguous = (sorted_old == list(range(len(sorted_old))))
    if already_contiguous:
        # Check that all files exist with the same numbering; if so, skip
        ok = True
        for dname in ("output_msa_cluster", "output_msa_cluster_pairtrim"):
            d = base / dname
            if not d.is_dir():
                continue
            for k in range(len(sorted_old)):
                if not (d / f"ShallowMsa_{k:03d}.a3m").is_file():
                    ok = False
                    break
            if not ok: break
        if ok:
            return False  # nothing to do

    # Build mapping: old -> new (0..K-1 in cluster-sorted order)
    # We keep the original ordering by 'cluster'
    kept_df = kept_df.sort_values("cluster").reset_index(drop=True)
    mapping = {}
    for new_idx, (_, row) in enumerate(kept_df.iterrows()):
        try:
            old_idx = int(row["cluster"])
        except Exception:
            old_idx = row["cluster"]
        mapping[old_idx] = new_idx

    def _rename_in_dir(d: Path) -> None:
        if not d.is_dir():
            return

        # Phase 1: move kept → temp names (avoid collisions)
        temp_moves = []
        for old_idx, new_idx in mapping.items():
            src = d / f"ShallowMsa_{int(old_idx):03d}.a3m"
            if src.exists():
                tmp = d / f"__TMP_RENUM_{new_idx:03d}.a3m"
                try:
                    src.rename(tmp)
                    temp_moves.append((tmp, d / f"ShallowMsa_{new_idx:03d}.a3m"))
                except Exception:
                    pass  # continue best-effort

        # Remove all non-kept shallow MSAs
        for p in d.glob("ShallowMsa_*.a3m"):
            name = p.name
            m = re.search(r"ShallowMsa_(\d+)\.a3m$", name)
            if not m:
                continue
            old = int(m.group(1))
            if old not in mapping:
                try: p.unlink()
                except Exception: pass

        # Phase 2: finalize kept tmp → final
        for tmp, dst in temp_moves:
            try:
                if dst.exists():
                    dst.unlink()
                tmp.rename(dst)
            except Exception:
                pass

    # Apply to both dirs
    _rename_in_dir(base / "output_msa_cluster")
    _rename_in_dir(base / "output_msa_cluster_pairtrim")

    # Rewrite the CSV with only kept clusters and the new numbering/paths
    out_rows = []
    for _, row in kept_df.iterrows():
        old_idx = int(row["cluster"])
        new_idx = mapping[old_idx]
        n_val    = row["n"]    if "n" in row else None
        neff_val = row["neff"] if "neff" in row else None
        out_rows.append({
            "cluster": new_idx,
            "n": n_val,
            "neff": neff_val,
            "path": f"output_msa_cluster/ShallowMsa_{new_idx:03d}.a3m",
            "kept": True
        })
    pd.DataFrame(out_rows, columns=["cluster","n","neff","path","kept"]).to_csv(cl_csv, index=False)

    return True

def _truth_pdb_paths_for_pair(pair_id: str) -> tuple[str, str]:
    # Prefer chain-sliced files if present
    from utils.utils import pair_str_to_tuple
    pair_dir = Path("Pipeline") / pair_id
    a, b = pair_str_to_tuple(pair_id)  # e.g., ("2qqjA","4qdsA")
    c1, c2 = a[-1], b[-1]
    cand1 = pair_dir / "chain_pdb_files" / f"{a}.pdb"
    cand2 = pair_dir / "chain_pdb_files" / f"{b}.pdb"
    pdb1 = str(cand1 if cand1.is_file() else (pair_dir / f"{a[:-1]}.pdb"))
    pdb2 = str(cand2 if cand2.is_file() else (pair_dir / f"{b[:-1]}.pdb"))
    return pdb1, pdb2


def _prepare_pairtrim_msas(pair_id: str, cap_1024: bool = True) -> str:
    """
    Create pair-trimmed A3Ms (keep columns where S1 OR S2 has a residue) to avoid
    gap-gap columns and length explosions for MSA-Transformer/CCMpred.
    Returns the directory that contains the pair-trimmed cluster A3Ms.
    """
    deep_src = f"Pipeline/{pair_id}/output_get_msa/DeepMsa.a3m"
    clus_src = f"Pipeline/{pair_id}/output_msa_cluster"
    clus_dst = f"Pipeline/{pair_id}/output_msa_cluster_pairtrim"
    ensure_dir(clus_dst)

    pA, pB = pair_str_to_tuple(pair_id)
    s1_tokens = [pA[:-1], pA[-1]]
    s2_tokens = [pB[:-1], pB[-1]]
    max_len = 1024 if cap_1024 else None

    # Deep (also keeps a .colmap.txt alongside)
    deep_dst = f"Pipeline/{pair_id}/output_get_msa/DeepMsa_pairtrim.a3m"
    keep = trim_a3m_for_pair_union(deep_src, deep_dst, s1_tokens, s2_tokens, max_len=max_len)

    # Also make a trimmed seed for analysis (so indices match predictions)
    try:
        seed_src = f"Pipeline/{pair_id}/_seed_both.a3m"
        if os.path.isfile(seed_src):
            # project first two rows by the same 'keep'
            lines = Path(seed_src).read_text().splitlines(True)
            headers, seqs = parse_a3m_records(lines, strip_inserts=True)
            seeds2 = project_columns(seqs[:2], keep)
            with open(f"Pipeline/{pair_id}/_seed_both_pairtrim.a3m", "w") as fh:
                fh.write(">S1\n" + seeds2[0] + "\n>S2\n" + seeds2[1] + "\n")
    except Exception as e:
        print(f"[seed_pairtrim] WARN: {e}")

    # Clusters
    for a3m in sorted(glob(os.path.join(clus_src, "ShallowMsa_*.a3m"))):
        dst = os.path.join(clus_dst, os.path.basename(a3m))
        trim_a3m_for_pair_union(a3m, dst, s1_tokens, s2_tokens, max_len=max_len)

    return clus_dst


def _extract_seq_ca_only(pdb_path: str, chain: str) -> str:
    """
    Return the CA-based sequence for the given chain, using the same codepath
    that constructs the true CA contact map. This guarantees:
        len(seq) == true_cmap.shape[0]
    and preserves non-standard residues via your existing mapping tables
    (e.g., MSE→M through NONSTD_TO_STD; unknowns become 'X').
    """
    atom_array = load_structure_to_atom_array(pdb_path)  # Biopython → Biotite AtomArray
    _, _, pdb_seq, _, _ = read_seq_coord_contacts_from_pdb(atom_array, chain=chain)
    return pdb_seq


def _has_deltaG(pair_id: str) -> bool:
    """Return True if both chain-level ΔG files exist for this pair."""
    pA, pB = pair_str_to_tuple(pair_id)     # e.g., '1wp8C','5ejbC'
    out_dir = Path("Pipeline/output_deltaG")
    a = out_dir / f"deltaG_{pA[:-1]}.txt"   # 'deltaG_1wp8.txt'
    b = out_dir / f"deltaG_{pB[:-1]}.txt"   # 'deltaG_5ejb.txt'
    return a.exists() and b.exists()




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


def _repo_root() -> str:
    return str(Path(__file__).resolve().parent)

def _submit_pair_job(run_mode: str, pair_id: str, args: argparse.Namespace, extra_cli: str = "") -> None:
    gres = getattr(args, "sbatch_gres", "gpu:1")
    cpus = int(getattr(args, "sbatch_cpus", 4))
    mem  = str(getattr(args, "sbatch_mem", "32G"))
    time = str(getattr(args, "sbatch_time", "24:00:00"))
    part = getattr(args, "sbatch_partition", None)
    cons = getattr(args, "sbatch_constraint", None)
    nlist= getattr(args, "sbatch_nodelist", None)
    acct = getattr(args, "sbatch_account", None)
    qos  = getattr(args, "sbatch_qos", None)
    mail = getattr(args, "sbatch_mail", None)
    mtyp = getattr(args, "sbatch_mail_type", None)

    # --- NEW: auto-upgrade GPU for long chains when running esm3/both ---
    try:
        L = pair_max_len_from_truth(pair_id)  # or _pair_max_len_from_truth(pair_id) if you have it
    except Exception:
        L = 0
    max_len = int(getattr(args, "esm_gpu_len_threshold", ESM_MAX_LIGHT_SEQ_LEN))
    if run_mode == "run_esmfold" and wants_esm3 and L >= max_len:
        gres = getattr(args, "sbatch_gres_heavy", "gpu:l40s:1")
        print(f"[esm-sbatch] {pair_id}: max_len={L} ≥ {thr} ⇒ using {gres}", flush=True)

    inner = [
        shlex.quote(sys.executable),
        "run_foldswitch_pipeline.py",
        "--run_mode", run_mode,
        "--foldpair_ids", pair_id,
        "--run_job_mode", "inline",
    ]
    if extra_cli:
        inner.append(extra_cli)
    wrap_str = " ".join(inner)

    job_dir = Path(f"Pipeline/{pair_id}/jobs")
    job_dir.mkdir(parents=True, exist_ok=True)
    out_path = str(job_dir / f"{run_mode}_{pair_id}.out")

    def build_sbatch(with_partition: bool):
        sb = [
            "sbatch",
            "-J", f"{run_mode}_{pair_id}",
            f"--gres={gres}",
            "-c", str(cpus),
            f"--mem={mem}",
            "-t", time,
            "--parsable",
            "--chdir", _repo_root(),
            "--output", out_path,
            "--wrap", wrap_str,
        ]
        # If user specified a nodelist, prefer it and SKIP -p to avoid conflicts
        if nlist:
            sb += ["--nodelist", nlist]
        elif with_partition and part:
            sb += ["-p", part]
        if cons:
            sb += [f"--constraint={cons}"]
        if acct: sb += ["-A", acct]
        if qos:  sb += ["--qos", qos]
        if mail: sb += ["--mail-user", mail]
        if mtyp: sb += ["--mail-type", mtyp]
        return sb

    # Prefer: if nodelist is set, don’t pass partition at all
    sb = build_sbatch(with_partition=not bool(nlist))
    print("[sbatch]", " ".join(shlex.quote(x) for x in sb), flush=True)
    try:
        res = subprocess.run(sb, check=True, capture_output=True, text=True)
        jobid = res.stdout.strip()
        print(f"[submitted] {run_mode} {pair_id} → job {jobid}", flush=True)
    except subprocess.CalledProcessError as e:
        msg = (e.stderr or "").strip().lower()
        print(f"[sbatch ERROR] exit={e.returncode}", (e.stderr or "").strip(), sep="\n", flush=True)
        # Auto-retry once without partition if Slurm says partition is invalid/inactive
        if "required partition not available" in msg or "invalid partition" in msg:
            sb2 = build_sbatch(with_partition=False)
            print("[sbatch RETRY]", " ".join(shlex.quote(x) for x in sb2), flush=True)
            res2 = subprocess.run(sb2, check=True, capture_output=True, text=True)
            jobid = res2.stdout.strip()
            print(f"[submitted] {run_mode} {pair_id} → job {jobid}", flush=True)
        else:
            raise




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

    # Optional toggles to carry through:
    ps = getattr(args, "plot_scope", None)
    if ps:
        cmd += f" --plot_scope {shlex.quote(ps)}"
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

    # --- NEW: auto-upgrade GPU for long chains when running esm3/both ---
    try:
        L = pair_max_len_from_truth(pair_id)  # or _pair_max_len_from_truth(pair_id) if you have it
    except Exception:
        L = 0
    max_len = int(getattr(args, "esm_gpu_len_threshold", ESM_MAX_LIGHT_SEQ_LEN))
    if run_mode == "run_esmfold" and wants_esm3 and L >= max_len:
        gres = getattr(args, "sbatch_gres_heavy", "gpu:l40s:1")
        print(f"[esm-sbatch] {pair_id}: max_len={L} ≥ {thr} ⇒ using {gres}", flush=True)


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

def _prune_stale_fastas(pair_id: str, pair_dir: Path, dry: bool):
    # Keep chain FASTAs (e.g., 2qqjA.fasta, 4qdsA.fasta); drop stray FASTAs if chain ones exist
    m = re.match(r"^([0-9A-Za-z]{4})([A-Za-z])_([0-9A-Za-z]{4})([A-Za-z])$", pair_id)
    expected = set()
    if m:
        pdb1, ch1, pdb2, ch2 = m.groups()
        expected = {f"{pdb1}{ch1}.fasta", f"{pdb2}{ch2}.fasta"}
    have_chain = expected and all((pair_dir / e).exists() for e in expected)
    for f in pair_dir.glob("*.fasta"):
        if have_chain and f.name not in expected:
            print(f"[clean] rm stale FASTA: {f}")
            if not dry:
                try: f.unlink()
                except: pass


def _update_basic_cache(pair_id: str) -> None:
    """
    Update Pipeline/<pair>/Analysis/cache.json with:
      - n_residues_A, n_residues_B (from chain FASTAs)
      - msa_depth, msa_width (DeepMsa.a3m; width after stripping inserts)
      - seq_id (pairwise identity of the two chains; matches / non-gap positions)
      - clusters: { "ShallowMsa_000": {"n": ..., "neff": ...}, ... }
    Safe: only writes keys it can compute.
    """
    pair_dir = Path(f"Pipeline/{pair_id}")
    analysis = pair_dir / "Analysis"
    analysis.mkdir(parents=True, exist_ok=True)
    cache_path = analysis / "cache.json"
    cache = {}
    if cache_path.exists():
        try:
            cache = json.loads(cache_path.read_text())
        except Exception:
            cache = {}

    # A) residues & seq-id from chain FASTAs
    m = re.match(r"^([0-9A-Za-z]{4})([A-Za-z])_([0-9A-Za-z]{4})([A-Za-z])$", pair_id)
    if m:
        pdb1, ch1, pdb2, ch2 = m.groups()
        fa1 = pair_dir / "fasta_chain_files" / f"{pdb1}{ch1}.fasta"
        fa2 = pair_dir / "fasta_chain_files" / f"{pdb2}{ch2}.fasta"
        if fa1.exists() and fa2.exists():
            from Bio import SeqIO, Align
            s1 = str(next(SeqIO.parse(str(fa1), "fasta")).seq)
            s2 = str(next(SeqIO.parse(str(fa2), "fasta")).seq)
            cache["n_residues_A"] = len(s1)
            cache["n_residues_B"] = len(s2)
            # global alignment with identity-like scoring
            aligner = Align.PairwiseAligner()
            aligner.mode = "global"
            aligner.match_score = 1.0
            aligner.mismatch_score = 0.0
            aligner.open_gap_score = -1.0
            aligner.extend_gap_score = -0.5
            aln = aligner.align(s1, s2)[0]
            # compute identity over both-residue positions
            a, b = str(aln.aligned[0]), str(aln.aligned[1])  # not directly the strings; build from coordinates
            # Reconstruct aligned strings (Aligner doesn't expose directly; use format)
            aligned = aln.format().splitlines()
            # second and third lines are aligned strings with gaps
            if len(aligned) >= 3:
                Aline = aligned[1].strip()
                Bline = aligned[2].strip()
                both = [(x != '-') and (y != '-') for x,y in zip(Aline, Bline)]
                denom = sum(both)
                matches = sum((x == y) and m for x,y,m in zip(Aline, Bline, both))
                cache["seq_id"] = float(matches/denom) if denom>0 else 0.0

    # B) Deep MSA depth & width
    deep = pair_dir / "output_get_msa" / "DeepMsa.a3m"
    if deep.exists():
        lines = deep.read_text().splitlines(True)
        headers, seqs = parse_a3m_records(lines, strip_inserts=True)
        if seqs:
            cache["msa_depth"] = len(seqs)
            cache["msa_width"] = len(seqs[0])

    st = StageTimer("Computing Neff for clusters for cache", verbose)

    # C) Per-cluster size & Neff
    clus_dir = pair_dir / "output_msa_cluster"
    clusters = {}
    if clus_dir.exists():
        for a3m in sorted(clus_dir.glob("ShallowMsa_*.a3m")):
            name = a3m.stem  # e.g., ShallowMsa_000
            lines = a3m.read_text().splitlines(True)
            hdrs, seqs = parse_a3m_records(lines, strip_inserts=True)
            n = len(seqs)
            neff = compute_neff_from_a3m(str(a3m), id_thresh=0.8, use_query_columns=True)
            clusters[name] = {"n": int(n), "neff": int(neff)}
    if clusters:
        cache["clusters"] = clusters

    cache["_updated"] = int(__import__("time").time())
    cache_path.write_text(json.dumps(cache, indent=2))
    print(f"[cache] updated: {cache_path}")



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

        pdb_path = os.path.join(pair_dir, f"{pdb}.pdb")
        seq = _extract_seq_ca_only(pdb_path, ch)
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

def _has_loaded(pair_id: str) -> bool:
    """
    Consider 'loaded' if for BOTH folds we already have either the chain-sliced PDB
    Pipeline/<pair>/chain_pdb_files/<PDBID+CHAIN>.pdb or the base PDB
    Pipeline/<pair>/<PDBID>.pdb present and non-empty.
    """
    a, b = pair_str_to_tuple(pair_id)  # e.g., ("2qqjA", "4qdsA")
    base = Path("Pipeline") / pair_id

    # A chain or base
    a_chain = base / "chain_pdb_files" / f"{a}.pdb"
    a_base  = base / f"{a[:-1]}.pdb"

    # B chain or base
    b_chain = base / "chain_pdb_files" / f"{b}.pdb"
    b_base  = base / f"{b[:-1]}.pdb"

    def _ok(p: Path) -> bool:
        return p.is_file() and p.stat().st_size > 0

    return (_ok(a_chain) or _ok(a_base)) and (_ok(b_chain) or _ok(b_base))


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
    Build a per-chain A3M for AF:
      • For clustered A3M: try to take the mask from a seed/aligned row for this chain.
        - Accept headers: S1/S2, chain names (e.g., 1abcA/2xyzB), or any row whose ungapped == chain FASTA.
      • If no suitable row is found in the cluster A3M, fall back to using the aligned row from DeepMsa.a3m.
      • For DeepMsa (cluster_a3m is None): require an aligned row for the chain in DeepMsa and use it.
    The query row is written ungapped; all other rows are column-masked (drop columns where the query has gaps).
    """
    def _ungap_upper(s: str) -> str:
        return "".join(ch for ch in (s or "") if ch.isalpha()).upper()

    # load chain query seq
    pair_dir = os.path.dirname(os.path.dirname(deep_a3m))  # .../Pipeline/<pair>
    chain_fa = os.path.join(pair_dir, "fasta_chain_files", f"{chain_tag}.fasta")
    ids, seqs = load_fasta(chain_fa)
    if not seqs or not seqs[0]:
        print(f"[a3m] ERR no sequence in {chain_fa} for {chain_tag}")
        return False
    chain_seq = "".join(ch for ch in seqs[0].strip() if ch.isalpha())
    chain_key = chain_seq.upper()

    def read_entries(a3m_path: str):
        return read_msa(a3m_path) if a3m_path else []

    def build_index(entries):
        idx = {}
        for nm, aln in entries:
            idx.setdefault(_ungap_upper(aln), []).append((nm, aln))
        return idx

    def write_masked(outp, query_seq_ungapped, entries, mask_keep, skip_ungap):
        ensure_dir(os.path.dirname(outp))
        with open(outp, "w", encoding="utf-8") as fh:
            fh.write(f">{chain_tag}\n{query_seq_ungapped}\n")
            for nm, aln in entries:
                if _ungap_upper(aln) == skip_ungap:
                    continue
                masked = "".join(ch.upper() for ch, k in zip(aln, mask_keep) if k and ch != ".")
                fh.write(f">{nm}\n{masked}\n")

    # -------- helper to get a seed/mask from a clustered A3M --------
    def seed_from_cluster(a3m_path: str):
        if not a3m_path:
            return None
        try:
            with open(a3m_path, "r", encoding="utf-8", errors="ignore") as fin:
                rows = []
                hdr = None
                for line in fin:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    if line.startswith(">"):
                        hdr = line[1:].split()[0]  # header w/o '>'
                    else:
                        rows.append((hdr, line))
        except FileNotFoundError:
            print(f"[a3m] ERR missing {a3m_path}")
            return None
        if not rows:
            print(f"[a3m] ERR empty {a3m_path}")
            return None

        # try header-based match first (S1/S2 or chain_tag)
        L = len(rows[0][1])
        candidates = {}
        for hdr, aln in rows:
            if hdr is None or len(aln) != L:
                continue
            if hdr in {"S1", "S2", chain_tag}:  # exact header match
                candidates[hdr] = aln

        # also look for the *other* chain name if present
        # e.g., if chain_tag is 1miqB and the other is 1qs8B, either could be the “seed”
        # extract pair names from file system if available
        other_hdrs = set()
        try:
            foldA, foldB = pair_str_to_tuple(os.path.basename(pair_dir))
            other_hdrs.update([foldA, foldB])
        except Exception:
            pass
        for hdr, aln in rows:
            if hdr in other_hdrs:
                candidates[hdr] = aln

        # final fallback: any row whose ungapped == chain query
        for hdr, aln in rows:
            if _ungap_upper(aln) == chain_key:
                candidates[hdr or "MATCH_BY_SEQ"] = aln
                break

        if not candidates:
            return None

        # pick a candidate (prefer exact chain header; else any)
        chosen_hdr = chain_tag if chain_tag in candidates else next(iter(candidates.keys()))
        seed_aln = candidates[chosen_hdr]
        return seed_aln, rows, L

    # ========== clustered A3M path ==========
    if cluster_a3m:
        seed = seed_from_cluster(cluster_a3m)
        if seed is not None:
            seed_aln, all_rows, L = seed
            keep = [ch != "-" for ch in seed_aln]
            kept = sum(keep)
            if kept != len(chain_seq):
                print(f"[a3m] WARN {chain_tag}: seed mask (keep={kept}) "
                      f"!= query length ({len(chain_seq)}); using DeepMsa mask instead")
                base_entries = read_msa(deep_a3m)
                idx = {"".join(ch for ch in aln if ch.isalpha()).upper(): aln for _, aln in base_entries}
                if "".join(ch for ch in chain_seq if ch.isalpha()).upper() not in idx:
                    print(f"[a3m] ERR {chain_tag}: cannot find aligned row in DeepMsa for fallback")
                    return False
                chain_aln = idx["".join(ch for ch in chain_seq if ch.isalpha()).upper()]
                keep = [ch.isalpha() for ch in chain_aln]
                # and proceed to write with this 'keep'
            # sanity: the number of kept columns should match len(query) in the ideal case
            kept = sum(keep)
            if kept != len(chain_seq):
                print(f"[a3m] WARN {chain_tag}: keep={kept} vs |query|={len(chain_seq)} in {os.path.basename(cluster_a3m)}")
            # stream from the cluster file to preserve order
            ensure_dir(os.path.dirname(out_path))
            with open(out_path, "w", encoding="utf-8") as fh:
                fh.write(f">{chain_tag}\n{chain_seq}\n")
                for hdr, aln in all_rows:
                    if aln is None or len(aln) != L:
                        continue
                    masked = "".join(ch for ch, k in zip(aln, keep) if k and ch != ".")
                    fh.write(f">{hdr}\n{masked}\n")
            print(f"[a3m] wrote {out_path}  (mask L={L}, keep={kept}, drop={L-kept})")
            return True

        # no seed row found in cluster; fall back to DeepMsa
        print(f"[a3m] WARN {chain_tag}: no seed row found in {os.path.basename(cluster_a3m)}; falling back to DeepMsa mask")

    # ========== DeepMsa path (or fallback) ==========
    base_entries = read_entries(deep_a3m)
    if not base_entries:
        print(f"[a3m] ERR base alignment empty for {out_path}")
        return False
    idx = build_index(base_entries)
    if chain_key not in idx:
        print(f"[a3m] ERR {chain_tag}: no aligned row for chain in DeepMsa; skipping.")
        return False

    _, chain_aln = idx[chain_key][0]
    keep = [ch.isalpha() for ch in chain_aln]
    write_masked(out_path, _ungap_upper(chain_aln), base_entries, keep, skip_ungap=_ungap_upper(chain_aln))
    print(f"[a3m] wrote {out_path}  (DeepMsa mask keep={sum(keep)})")
    return True

def _submit_notebook_job(pair_id: str, kernel: str, args: argparse.Namespace) -> None:
    out = Path(f"Pipeline/{pair_id}/jobs")
    out.mkdir(parents=True, exist_ok=True)
    log = out / f"nb_{pair_id}.out"
    cpus = getattr(args, "html_cpus", 2)
    mem  = getattr(args, "html_mem",  "4G")
    time = getattr(args, "html_time", "02:00:00")
    part = getattr(args, "html_partition", None)

    cmd = f"{shlex.quote(sys.executable)} Analysis/NotebookGen/generate_notebooks.py {shlex.quote(pair_id)} --kernel {shlex.quote(kernel)}"
    sbatch = [
        "sbatch",
        f"--job-name=nb_{pair_id}",
        f"--cpus-per-task={cpus}",
        f"--mem={mem}",
        f"--time={time}",
        f"--output={log}",
        "--parsable",
        "--wrap", shlex.quote(cmd),
    ]
    if part: sbatch.insert(1, f"--partition={part}")

    subprocess.run(" ".join(sbatch), shell=True, check=True, env=_jupyter_env_for_scratch())


# ------------------------- tasks -------------------------
# add imports at file top if missing
from pathlib import Path
import shutil, re
from config import PAIR_TARGET_DIRS

def task_clean(pair_id: str, args) -> None:
    dry     = _bool_from_tf(getattr(args, "clean_dry_run", "TRUE"))
    level   = getattr(args, "clean_level", "derived")
    targets = getattr(args, "clean_targets", None)
    pair_dir = Path(f"Pipeline/{pair_id}")

    if not pair_dir.exists():
        print(f"[clean] skip (missing): {pair_id}")
        return

    # ----- TARGETED CLEAN (overrides level) -----
    if targets:
        # expand targets -> relative paths
        rel_paths = []
        for t in targets:
            dirs = PAIR_TARGET_DIRS.get(t)
            if not dirs:
                print(f"[clean] unknown target '{t}', skipping.")
                continue
            rel_paths.extend(dirs)

        # de-dupe while keeping order
        seen, uniq = set(), []
        for rel in rel_paths:
            if rel not in seen:
                seen.add(rel)
                uniq.append(rel)

        for rel in uniq:
            p = pair_dir / rel
            if p.exists():
                print(f"[clean:{','.join(targets)}] rm -rf {p}")
                if not dry:
                    print("[WET] Running actual clean")
                    shutil.rmtree(p, ignore_errors=True)
                else:
                    print("[DRY] Skipping actual clean")

        # special handling for logs: also remove root-level files
        if any(t in ("logs",) for t in targets):
            for pat in ["*.out", "*.log", "*.err", "RunAF.out", "CmapESM.out", "run_pipeline_for_*.out"]:
                for f in pair_dir.glob(pat):
                    if f.is_file():
                        print(f"[clean:logs] rm {f}")
                        if not dry:
                            print("[WET] Running actual clean")
                            try: f.unlink()
                            except: pass
                        else:
                            print("[DRY] Skipping actual clean")

        # special handling for tmp: remove any tmp_* dirs
        if any(t in ("tmp",) for t in targets):
            for p in pair_dir.glob("tmp_*"):
                if p.is_dir():
                    print(f"[clean:tmp] rm -rf {p}")
                    if not dry:
                        print("[WET] Running actual clean")
                        shutil.rmtree(p, ignore_errors=True)
                    else:
                        print("[DRY] Skipping actual clean")

        _prune_stale_fastas(pair_id, pair_dir, dry)
        return

    # ----- EXISTING LEVEL-BASED CLEAN (unchanged semantics) -----
    if level == "all":
        print(f"[clean:ALL] rm -rf {pair_dir}")
        if not dry:
            print("[WET] Running actual clean:ALL")
            shutil.rmtree(pair_dir, ignore_errors=True)
        else:
            print("[DRY] skipping clean:ALL (dry run)")
        return

    # derived-only clean (keep base inputs)
    rm_dirs_current = [
        "output_get_msa", "output_msa_cluster", "output_AF", "output_deltaG", "output_figs",
        "output_cmaps", "output_esm_fold", "output_phytree",
        "Analysis", "tmp_msa_files", "tmp_esmfold", "jobs", "logs",
    ]
    rm_dirs_legacy = ["output_cmap_esm", "esm_cmap_output", "AF_preds", "chain_pdb_files", "fasta_chain_files"]
    rm_globs = ["run_pipeline_for_*.out", "RunAF.out", "CmapESM.out", "*.out", "*.log", "*.err", "*hdbscan*.pkl"]

    for d in rm_dirs_current + rm_dirs_legacy:
        p = pair_dir / d
        if p.exists():
            print(f"[clean] rm -rf {p}")
            if not dry:
                print("[WET] Running actual clean")
                shutil.rmtree(p, ignore_errors=True)
            else:
                print("[DRY] Skipping actual clean")

    for p in pair_dir.glob("tmp_*"):
        if p.is_dir():
            print(f"[clean] rm -rf {p}")
            if not dry:
                print("[WET] Running actual clean")
                shutil.rmtree(p, ignore_errors=True)
            else:
                print("[DRY] Skipping actual clean")

    for pat in rm_globs:
        for f in pair_dir.glob(pat):
            if f.is_file():
                print(f"[clean] rm {f}")
                if not dry:
                    print("[WET] Running actual clean")
                    try: f.unlink()
                    except: pass
                else:
                    print("[DRY] Skipping actual clean")

    for pat in pair_dir.glob("*hdbscan*.pkl"):
        print("removing hdbscan pickle:", pat)
        pat.unlink()

    _prune_stale_fastas(pair_id, pair_dir, dry)


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
    pair_dir  = f"Pipeline/{pair_id}"
    foldA, foldB = pair_str_to_tuple(pair_id)
    pdbids    = [foldA[:-1], foldB[:-1]]
    pdbchains = [foldA[-1],  foldB[-1]]
    ensure_dir(os.path.join(pair_dir, "fasta_chain_files"))
    print("Making fasta chain files: ")
    ensure_chain_fastas(pair_dir, pdbids, pdbchains)

    seed_a3m = build_pair_seed_a3m_from_pair(pair_id, data_dir="Pipeline")
    out_dir  = f"{pair_dir}/output_get_msa"
    ensure_dir(out_dir)

    # Resolve repo root and sbatch script absolutely
    repo_root = Path(__file__).resolve().parent
    sbatch_script = repo_root / "Pipeline" / "get_msa_params.sh"
    sbatch_exists = sbatch_script.is_file()
    print(f"[get_msa] run_job_mode={run_job_mode} sbatch_script={sbatch_script} exists={sbatch_exists}")

    # Interpreter to use for INLINE fallback (if sbatch not used/available)
    python_for_msa = os.environ.get("GET_MSA_PYTHON", sys.executable)
    print(f"[get_msa] using interpreter: {python_for_msa}")

    if run_job_mode == "sbatch" and sbatch_exists:
        log = f"{pair_dir}/get_msa_for_{pair_id}.out"
        cmd = (
            f"sbatch -o {shlex.quote(log)} "
            f"{shlex.quote(str(sbatch_script))} "
            f"{shlex.quote(seed_a3m)} {shlex.quote(pair_id)}"
        )
        _run(cmd, "sbatch")
    else:
        # Safety net: if we’re here, run inline but with the ColabFold python if provided
        os.environ.setdefault("COLABFOLD_USE_LOCAL_MMSEQS", "0")
        os.environ.setdefault("COLABFOLD_USER_AGENT", "MsaCluster/0.1 or.zuk@mail.huji.ac.il")
        cmd = (
            f"{shlex.quote(python_for_msa)} ./run_getMSA.py "
            f"{shlex.quote(seed_a3m)} {shlex.quote(out_dir)} --pair {shlex.quote(pair_id)}"
        )
        _run(cmd, "inline")

    # After DeepMsa exists, update basic cache (depth/width/seq-id)
    try:
        _update_basic_cache(pair_id)
    except Exception as e:
        print(f"[cache] WARN get_msa: {e}")


def task_cluster_msa(pair_id: str, run_job_mode: str, args) -> None:
    py  = shlex.quote(sys.executable)
    pid = shlex.quote(pair_id)

    alg = args.cluster_alg
    # If tree mode and no explicit tree is provided, try a sensible default location
    tree_arg = ""
    if alg == "tree":
        pair_dir = os.path.join("Pipeline", pair_id)

        # 1) If user provided a path, use it as-is (convert to pair-relative if it redundantly includes the pair prefix)
        tree_path = args.cluster_tree_path
        if tree_path:
            if not os.path.isabs(tree_path):
                # If user passed "Pipeline/<pair_id>/...", strip that prefix because we cd into pair_dir
                pref = os.path.join("Pipeline", pair_id) + os.sep
                if tree_path.startswith(pref):
                    tree_path = tree_path[len(pref):]
        else:
            # 2) No path provided: choose a sensible default.
            # Preferred default: relative to the pair dir (we cd into it before calling run_ClusterMSA.py)
            tree_path = os.path.join(pair_dir, "output_phytree", "DeepMsa_tree.nwk")
            # we already `cd` into pair dir; default tree is relative to it
            tree_path = "output_phytree/DeepMsa_tree.nwk"

        tree_arg = f"--tree_path {shlex.quote(tree_path)} "


    cmd = (
        f"bash -lc 'cd Pipeline/{pid} && "
        f"{py} ../../run_ClusterMSA.py "
        f"--keyword ShallowMsa "
        f"--a3m output_get_msa/DeepMsa.a3m "
        f"-o output_msa_cluster "
        f"--cluster_alg {shlex.quote(alg)} {tree_arg} "
        f"--max_clusters {int(args.cluster_max_clusters)} "
        f"--min_output_size {int(args.cluster_min_output_size)} "
        f"--min_neff {int(args.cluster_min_neff)} "
        f"--neff_id_thresh {float(args.cluster_neff_id_thresh)} "
        f"--neff_mode approx "
        f"--neff_use_query_columns 1 "
        f"--neff_strip_inserts 1 "
        f"--neff_approx_n_hashes 3 "
        f"--neff_approx_sig_len 32 "
        f"--neff_approx_bucket_cap 5000 "
        f"--neff_approx_subsample_cap 0 "        
        f"--frac_gaps_cutoff {float(args.cluster_frac_gaps_cutoff)} "
        f"--sample_cap {int(args.cluster_sample_cap)} "
        f"--sample_seed {int(args.cluster_sample_seed)} "
        f"--min_cluster_size {int(args.hdbscan_min_cluster_size)} "
        f"--min_samples {shlex.quote(str(args.hdbscan_min_samples))} "
        f"--cluster_selection {shlex.quote(args.hdbscan_cluster_selection)} "
        f"'"
    )
    _run(cmd, run_job_mode)

    print(f"[cache] Finished clustering, updating cache")
    # After ClusterMsa exists, update basic cache (depth/width/seq-id)
    try:
        _update_basic_cache(pair_id)
    except Exception as e:
        print(f"[cache] WARN cluster_msa: {e}")


def task_cmap_msa_transformer(pair_id: str, run_job_mode: str) -> None:
    """
    Always run MSA-Transformer on the *trimmed* MSAs for both clusters and deep,
    but save the deep output as 'msa_t__DeepMsa.npy' (no 'pairtrim' in the name).
    """

    outdir = f"Pipeline/{pair_id}/output_cmaps/msa_transformer"
    ensure_dir(outdir)

    # 1) Make trimmed Deep + trimmed clusters (≤1024)
    clus_dir_pairtrim = _prepare_pairtrim_msas(pair_id, cap_1024=True)

    # 2) Clusters (trimmed)
    # NOTE the leading space before --model is important.
    cmd = (
        f"python3 ./run_MSATrans.py "
        f"--input_msas {clus_dir_pairtrim} "
        f"-o {outdir} "
        f"--model msa_t --keyword clusters --clean matched"
    )
    _run(cmd, run_job_mode)

    # 3) Deep (trimmed) — but force basename 'DeepMsa' for output naming
    deep_trim = f"Pipeline/{pair_id}/output_get_msa/DeepMsa_pairtrim.a3m"
    if os.path.isfile(deep_trim):
        # Make a tiny temp input folder and place a copy named exactly 'DeepMsa.a3m'
        # (Copy instead of symlink to avoid FS/permission surprises.)
        tmp_in = os.path.join(outdir, "_tmp_msa_inputs")
        if os.path.isdir(tmp_in):
            # Clean the temp dir so --clean matched doesn't accidentally hit old junk
            try:
                shutil.rmtree(tmp_in)
            except OSError:
                pass
        ensure_dir(tmp_in)

        deep_alias = os.path.join(tmp_in, "DeepMsa.a3m")
        shutil.copyfile(deep_trim, deep_alias)

        # Empty keyword ⇒ output file will be 'msa_t__DeepMsa.npy'
        cmd_deep = (
            f"python3 ./run_MSATrans.py "
            f"--input_msas {deep_alias} "
            f"-o {outdir} "
            f"--model msa_t --keyword DeepMsa --clean matched"
        )
        _run(cmd_deep, run_job_mode)

        # Optional: clean temp input
        try:
            shutil.rmtree(tmp_in)
        except OSError:
            pass


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
        cmd = f"{shlex.quote(CCMPRED_EXE)} -t {threads} {shlex.quote(str(fa))} {shlex.quote(str(mat))}"
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
    - Sample up to N seqs per cluster -> tmp_esmfold/*.fasta
    - Call run_ESMFoldHF once:
        python3 ./run_ESMFoldHF.py -input <PAIR_ID> --model {esm2|esm3} --device {auto|cuda|cpu}
    - Outputs: Pipeline/<pair>/output_esm_fold/<esm_model>/
    """
    if _is_windows():
        raise SystemExit("ESMFold can’t run on Windows. Run on Moriah/Linux.")

    # Prepare sampled FASTAs (your helper)
    _sample_to_tmp_fastas(pair_id, args.cluster_sample_n, include_deep=False)

    # Decide device: if sequences are very long, prefer CPU to avoid CUDA OOM
    tmp_dir = Path(f"Pipeline/{pair_id}/tmp_esmfold")
    max_len = 0
    for fa in sorted(tmp_dir.glob("*.fasta")):
        with open(fa) as f:
            for line in f:
                if line.startswith(">"): 
                    continue
                max_len = max(max_len, len(line.strip()))

    # Prefer GPU by default; only use 'cuda' unless the user explicitly says otherwise
    device = args.esm_device or "cuda"

    ensure_dir(f"Pipeline/{pair_id}/output_esm_fold/{args.esm_model}")

    # Always require a GPU when device is 'cuda' or 'auto'
    require_cuda_flag = " --require_cuda" if device in ("cuda", "auto") else ""

    cmd = (f"python3 ./run_ESMFoldHF.py -input {pair_id} "
        f"--model {args.esm_model} --device {device}{require_cuda_flag}")

    try:
        _run(cmd, args.run_job_mode)
    except Exception as e:
        print(f"[esm] ERROR on GPU: {e}")



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

    print(f"[af-plan] pair={pair_id} clusters-a3m={len(cluster_a3ms)}", flush=True)
    if not cluster_a3ms:
        print(f"[af-plan]   (no ShallowMsa_*.a3m under {cluster_dir})", flush=True)

    tmp_pairs_dir = os.path.join(pair_dir, "tmp_msa_files")
    ensure_dir(tmp_pairs_dir)

    # Build pair-specific A3Ms per chain for DeepMsa + each cluster (same as before)
    jobs = []  # (a3m_path, out_dir_base_name) — out_dir_base_name = "DeepMsa/<chain>" or "ShallowMsa_XXX/<chain>"
    for ch in chains:
        pair_a3m = os.path.join(tmp_pairs_dir, f"tmp_DeepMsa__{ch}.a3m")
        print(f"[af-plan] build pair A3M (DeepMsa): src={deep_a3m} → {pair_a3m}  chain={ch}", flush=True)
        try:
            with _alarm_timeout(180, f"_write_pair_a3m_for_chain(DeepMsa,{ch})"):
                ok = _write_pair_a3m_for_chain(None, deep_a3m, ch, pair_a3m)
        except TimeoutError as te:
            print(f"[af-plan] DeepMsa→{ch}: TIMEOUT {te}; skipping", flush=True)
            ok = False
        except Exception as e:
            print(f"[af-plan] DeepMsa→{ch}: EXC {e}; skipping", flush=True)
            ok = False
        print(f"[af-plan] DeepMsa→{ch}: {'OK' if ok else 'FAIL'}  src={deep_a3m}", flush=True)
        if ok:
            jobs.append((pair_a3m, f"DeepMsa/{ch}"))

    for a3m in cluster_a3ms:
        cl_stem = Path(a3m).stem  # e.g., ShallowMsa_007
        try:
            fsize = os.path.getsize(a3m)
        except Exception:
            fsize = -1
        print(f"[af-plan] scanning {cl_stem}  size={fsize / 1e6:.1f} MB  path={a3m}", flush=True)

        for ch in chains:
            pair_a3m = os.path.join(tmp_pairs_dir, f"tmp_{cl_stem}__{ch}.a3m")
            print(f"[af-plan] build pair A3M ({cl_stem}): src={a3m} → {pair_a3m}  chain={ch}", flush=True)
            try:
                with _alarm_timeout(180, f"_write_pair_a3m_for_chain({cl_stem},{ch})"):
                    ok = _write_pair_a3m_for_chain(a3m, deep_a3m, ch, pair_a3m)
            except TimeoutError as te:
                print(f"[af-plan] {cl_stem}→{ch}: TIMEOUT {te}; skipping", flush=True)
                ok = False
            except Exception as e:
                print(f"[af-plan] {cl_stem}→{ch}: EXC {e}; skipping", flush=True)
                ok = False

            print(f"[af-plan] {cl_stem}→{ch}: {'OK' if ok else 'FAIL'}  src={a3m}", flush=True)
            if ok:
                jobs.append((pair_a3m, f"{cl_stem}/{ch}"))

    print(f"[af-plan] total AF jobs planned: {len(jobs)}", flush=True)
    if not jobs:
        print(f"[warn] No AF jobs to run for {pair_id}", flush=True)
        print(f"[diag] DeepMsa exists? {os.path.isfile(deep_a3m)}  path={deep_a3m}", flush=True)
        print(f"[diag] Shallow a3m files found: {len(cluster_a3ms)} in {cluster_dir}", flush=True)
        return

    inside_slurm = _in_slurm_session()
    mode_auto = ("sbatch" if (args.run_job_mode == "sbatch" or
                              (not inside_slurm and not getattr(args, "allow_inline_af", False)))
                 else "inline")
    print(f"[af-plan] execution={mode_auto}  total_jobs={len(jobs)}  af_ver={af_ver}", flush=True)

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
#                _run(f"sbatch {sbatch_opts} -o '{log_path}' --wrap {shlex.quote(cmd)}", "sbatch")
                sb = f"sbatch {sbatch_opts} -o '{log_path}' --wrap {shlex.quote(cmd)}"
                print(f"[af-sbatch] {sb}", flush=True)
                _run(sb, "sbatch")

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

def task_tree(pair_id: str, args: argparse.Namespace) -> None:
    msa_file = f"Pipeline/{pair_id}/output_get_msa/DeepMsa.a3m"
    out = f"Pipeline/{pair_id}/output_phytree/DeepMsa_tree.nwk"
    ensure_dir(os.path.dirname(out))

    # Effective cap: default 5000; <=0 → no cap (use all).
    cap = getattr(args, "tree_max_seqs", 5000)
    if cap is not None and int(cap) <= 0:
        cap = None  # all sequences -> no stratification

    # Use stratified sampling when capped: ensure ≥1 per ShallowMsa_###, then fill random.
    cluster_dir = f"Pipeline/{pair_id}/output_msa_cluster" if cap is not None else None

    # If there are more clusters than the cap, bump the cap to #clusters
    eff_cap = cap
    if cluster_dir and os.path.isdir(cluster_dir):
        n_clusters = len([fn for fn in os.listdir(cluster_dir)
                          if fn.startswith("ShallowMsa_") and fn.endswith(".a3m")])
        if eff_cap is None:
            pass
        else:
            eff_cap = max(int(eff_cap), n_clusters)

    # Build tree (phytree_utils.phytree_from_msa has the stratified logic)
    t0 = time.time()
    print(f"[tree] building from {msa_file} with cap={eff_cap if eff_cap is not None else 'ALL'} ...", flush=True)

    phytree_from_msa(
        msa_file,
        output_tree_file=out,
        max_seqs=eff_cap,                 # None → all; int → capped (stratified)
        cluster_msa_dir=cluster_dir,      # used only when capped
        seed=getattr(args, "tree_seed", 123),
    )

    dt = time.time() - t0
    print(f"[tree] done → {out} in {dt:.1f}s", flush=True)


def task_plot(pair_id: str | None, args: argparse.Namespace) -> None:
    """
    Plot mode with scope:
      - pair:   only per-pair plots
      - global: only global plots (once)
      - both:   per-pair plots, then global plots
    """
    scope = getattr(args, "plot_scope", "both")
    plot_trees = getattr(args, "plot_trees", False)

    # === 1) GLOBAL-ONLY: run ONCE and RETURN ===
    if scope in ["global", "both"]:
        os.makedirs(FIGURE_RES_DIR, exist_ok=True)
        print("[plot] Generating global plots…")
        # write all global scatter plots (AF2, AF3, ESM2, ESM3, MSA-Trans, etc.)
        global_pairs_statistics_plots(output_dir=FIGURE_RES_DIR)
        print("[plot] Global plots written to", FIGURE_RES_DIR)
        if scope == "global":  # global only
            return

    # === 2) Per-pair plots ===
    pairs = list_protein_pairs() if (pair_id in (None, "ALL")) else pair_id
    if isinstance(pairs, str):
        pairs = [pairs]
    for pid in pairs:
        pA, pB = pair_str_to_tuple(pid)
        pdbids    = [pA[:-1], pB[:-1]]
        pdbchains = [pA[-1],  pB[-1]]
        print(f"=== plot :: {tuple(pdbids)} ===", flush=True)
        print("pid is: ", pA + "_" + pB)
        make_foldswitch_all_plots(
            pdbids=pdbids,
            fasta_dir="Pipeline",
            foldpair_id=pA + "_" + pB,
            pdbchains=pdbchains,
            plot_tree_clusters=bool(plot_trees),
            plot_contacts=True,
            global_plots=False,   # global is handled outside this loop
            plot3dformat=args.plot3dformat  # allow interactive 3D plots
        )

def task_deltaG(pair_id: str) -> None:
    """
    Compute global ΔG per true structure and an aligned per-residue ΔΔG profile.
    Writes:
      Pipeline/<pair>/Analysis/df_energy_<token>.csv           (per-structure per-res energies)
      Pipeline/<pair>/Analysis/df_energy_global.csv            (2 rows: ΔG1/ΔG2)
      Pipeline/<pair>/Analysis/df_ddg_aligned.csv              (aligned ΔΔG table)
      Pipeline/<pair>/output_deltaG/deltaG_diff_*.jpg          (plot)
    """

    pair_dir = Path("Pipeline") / pair_id
    out_dir  = pair_dir / "output_deltaG"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Prefer chain-sliced true PDBs when available
    a, b = pair_str_to_tuple(pair_id)        # e.g., ("2qqjA","4qdsA")
    p_chain1 = pair_dir / "chain_pdb_files" / f"{a}.pdb"
    p_chain2 = pair_dir / "chain_pdb_files" / f"{b}.pdb"
    p_raw1   = pair_dir / f"{a[:-1]}.pdb"
    p_raw2   = pair_dir / f"{b[:-1]}.pdb"
    pdb1 = str(p_chain1 if p_chain1.is_file() else p_raw1)
    pdb2 = str(p_chain2 if p_chain2.is_file() else p_raw2)

    compute_global_and_residue_energies(pdb_pairs=[(pdb1, pdb2)], foldpair_ids=[pair_id], output_dir=str(out_dir), plot_dir=None)

def task_postprocess(foldpairs: list[str], args: argparse.Namespace) -> None:
    """
    Unified post-processing / reports:
      1) (optional) per-pair postprocess_unified (df_af/df_esm/df_cmap, etc.)
      2) build unified CSVs (summary + detailed) from cached per-pair Analysis CSVs
      3) generate global HTML tables (summary + clusters) and global plots page
      4) (optional) per-pair HTML pages
    Controlled by --reports: none | tables | html | all
    Safe to run incrementally.
    """

    # Normalize pairs to strings like "1wp8C_5ejbC"
    norm_pairs = [
        p if isinstance(p, str) else f"{p[0]}_{p[1]}"
        for p in (foldpairs if isinstance(foldpairs, list) else [foldpairs])
    ]
    # For builders that auto-discover when None:
    pairs_arg = None if norm_pairs == ["ALL"] else norm_pairs

    # Map CLI flags
    cached_only = _bool_from_tf(getattr(args, "cached_only", "FALSE"))
    force_rerun = not cached_only or _bool_from_tf(getattr(args, "force_rerun_postprocess", "FALSE"))
    run_job_mode = (getattr(args, "run_job_mode", "local") or "local").lower()

    # =========================
    # PER-PAIR POSTPROCESS STEP
    # =========================
    if not cached_only:
        # Figure out which pairs to process
        if pairs_arg is None:
            # auto-discover all pair dirs in Pipeline/
            pair_dirs = [d.name for d in Path(DATA_DIR).iterdir() if d.is_dir()]
        else:
            pair_dirs = pairs_arg

        if run_job_mode in ("local", "inline", ""):
            # Local/in-process: recompute per-pair CSVs
            try:
                from Analysis.postprocess_unified import post_processing_analysis
                post_processing_analysis(force_rerun=force_rerun, pairs=pair_dirs)
                print(f"[postprocess] per-pair Analysis refreshed (local) for {len(pair_dirs)} pairs.")
            except Exception as e:
                print(f"[postprocess] WARN post_processing_analysis(local): {e}")
        elif run_job_mode == "sbatch":
            # Submit one job per pair; then a single unify job with dependency
            jobs_dir = Path("Pipeline/jobs"); jobs_dir.mkdir(parents=True, exist_ok=True)

            # sbatch resource knobs (with sane defaults)
            pp_cpus = getattr(args, "postproc_cpus", 2)
            pp_mem  = getattr(args, "postproc_mem",  "8G")
            pp_time = getattr(args, "postproc_time", "02:00:00")
            pp_part = getattr(args, "postproc_partition", None)

            job_ids = []
            for pid in pair_dirs:
                cmd = (
                    f"{shlex.quote(sys.executable)} Analysis/postprocess_unified.py "
                    f"--pairs {shlex.quote(pid)} --mode full "
                    + ("--force_rerun " if force_rerun else "")
                )
                log = jobs_dir / f"postproc_{pid}.out"
                sbatch = [
                    "sbatch",
                    "--parsable",
                    f"--job-name=pp_{pid}",
                    f"--cpus-per-task={pp_cpus}",
                    f"--mem={pp_mem}",
                    f"--time={pp_time}",
                    f"--output={log}",
                    "--wrap", shlex.quote(cmd),
                ]
                if pp_part:
                    sbatch.insert(1, f"--partition={pp_part}")
                jid = subprocess.check_output(" ".join(sbatch), shell=True, text=True).strip()
                job_ids.append(jid)
                print(f"[postprocess] submitted per-pair job {jid} for {pid}")

            # Submit a "unify" job that depends on all per-pair jobs
            dep = ":".join(job_ids) if job_ids else ""
            unify_cmd = f"{shlex.quote(sys.executable)} Analysis/postprocess_unified.py --mode unify"
            # Allow this unify to also regenerate global CSV/HTML as requested
            # We pass env var to tell postprocess_unified which reports to emit (optional).
            unify_log = jobs_dir / "postproc_unify.out"
            sbatch_unify = [
                "sbatch",
                "--parsable",
                f"--job-name=pp_unify",
                f"--cpus-per-task=2",
                f"--mem=4G",
                f"--time=01:00:00",
                f"--output={unify_log}",
                f"--dependency=afterok:{dep}" if dep else "",
                "--wrap", shlex.quote(unify_cmd),
            ]
            # Remove empty arg if dep is empty
            sbatch_unify = [x for x in sbatch_unify if x]
            unify_jid = subprocess.check_output(" ".join(sbatch_unify), shell=True, text=True).strip()
            print(f"[postprocess] submitted UNIFY job {unify_jid} depending on {len(job_ids)} jobs.")
            # When using sbatch, we stop here and let unify produce tables/html.
            print("[postprocess] sbatch mode: global tables/HTML will be generated by the unify job.")
            return
        else:
            print(f"[postprocess] WARN unknown run_job_mode={run_job_mode}; running local fallback.")
            try:
                from Analysis.postprocess_unified import post_processing_analysis
                post_processing_analysis(force_rerun=force_rerun, pairs=pair_dirs)
            except Exception as e:
                print(f"[postprocess] WARN post_processing_analysis(fallback): {e}")

    # ================================
    # GLOBAL CSV + HTML (LOCAL ONLY)
    # ================================
    # Only build immediately if we didn’t go through sbatch path.
    try:
        force_global = _bool_from_tf(getattr(args, "force_rerun_postprocess", "FALSE"))
        summary_csv, detail_csv = build_or_load_global_tables(force=force_global)
        print(f"[global] summary: {summary_csv}")
        print(f"[global] detailed: {detail_csv}")
    except Exception as e:
        print(f"[postprocess] WARN global tables: {e}")

    if args.reports in ("tables", "all"):
        try:
            build_unified_tables_from_cluster_dfs(pairs=(None if pairs_arg is None else pairs_arg), write_out=True)
            print(f"[reports] unified CSVs written:\n  {SUMMARY_RESULTS_TABLE}\n  {DETAILED_RESULTS_TABLE}")
        except Exception as e:
            print(f"[reports] WARN building unified CSVs: {e}")

        try:
            os.makedirs(TABLES_RES, exist_ok=True)
            write_global_tables(force_rerun_csv=False, fade_min_clusters=2)
        except Exception as e:
            print(f"[reports] WARN write_global_tables: {e}")

        try:
            gen_html_for_global_plots(
                images_dir=FIGURE_RES_DIR,
                output_html=os.path.join("docs", "pairs_global_analysis.html"))
            src = os.path.join("docs", "pairs_global_analysis.html")
            dst = os.path.join(MAIN_DIR, "pairs_global_analysis.html")
            if os.path.isfile(src):
                shutil.copy2(src, dst)
                print(f"[reports] copied {src} -> {dst}")
        except Exception as e:
            print(f"[reports] WARN building global plots page: {e}")

    # Per-pair HTML generation unchanged (runs local/sbatch already below)
    if args.reports in ("html", "all"):
        try:
            ready = []
            for p in (norm_pairs if pairs_arg is not None else [d.name for d in Path(DATA_DIR).iterdir() if d.is_dir()]):
                if os.path.isfile(f"Pipeline/{p}/Analysis/df_af.csv"):
                    ready.append(p)
                else:
                    print(f"[html] skip {p}: missing Pipeline/{p}/Analysis/df_af.csv")

            if not ready:
                print("[html] no ready pairs; skipping per-pair HTML generation")
            else:
                mode = getattr(args, "html_mode", "inline")
                outdir = OUTPUT_PATH_NOTEBOOKS
                script = Path("TableResults/gen_pair_html.py")
                pairs_arg_str = "ALL" if pairs_arg is None else ",".join(ready)
                cmd = (
                    f"{shlex.quote(sys.executable)} {shlex.quote(str(script))} "
                    f"--pairs {shlex.quote(pairs_arg_str)} "
                    f"--mode {shlex.quote(mode)} "
                    f"--output_dir {shlex.quote(str(outdir))}"
                )
                if run_job_mode == "sbatch":
                    jobs_dir = Path("Pipeline/jobs"); jobs_dir.mkdir(parents=True, exist_ok=True)
                    log = jobs_dir / f"genhtml_{mode}.out"
                    cpus = getattr(args, "html_cpus", 2)
                    mem = getattr(args, "html_mem", "4G")
                    time = getattr(args, "html_time", "02:00:00")
                    part = getattr(args, "html_partition", None)
                    sbatch = [
                        "sbatch", "--parsable",
                        f"--job-name=genhtml",
                        f"--cpus-per-task={cpus}",
                        f"--mem={mem}",
                        f"--time={time}",
                        f"--output={log}",
                        "--wrap", shlex.quote(cmd),
                    ]
                    if part: sbatch.insert(1, f"--partition={part}")
                    subprocess.run(" ".join(sbatch), shell=True, check=True)
                    print(f"[html] submitted gen_pair_html job for {pairs_arg_str} (mode={mode})")
                else:
                    subprocess.run(cmd, shell=True, check=True)
                    print(f"[html] generated pages for {pairs_arg_str} (mode={mode})")
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

    def _step_hdr(n: int, title: str) -> None:
        line = f"Step {n}: {title}"
        print("\n" + line + "\n" + "=" * len(line))

    def _has_cmaps_shallow_local(pid: str) -> bool:
        """Any cluster npy produced by run_MSATrans.py first pass."""
        out = Path(f"Pipeline/{pid}/output_cmaps/msa_transformer")
        # Typical files: msa_t_clusters_ShallowMsa_000.npy, ...
        return any(out.glob("msa_t_clusters_ShallowMsa_*.npy"))

    def _has_cmaps_deep_local(pid: str) -> bool:
        """Any deep npy produced by run_MSATrans.py (either pass / naming)."""
        out = Path(f"Pipeline/{pid}/output_cmaps/msa_transformer")
        # We tolerate several historical names:
        patterns = [
            "msa_t__DeepMsa*.npy",          # old naming
            "msa_t_DeepMsa_*.npy",          # new naming with keyword DeepMsa
            "msa_t_clusters_MSA_deep.npy",  # deep was written in first pass
            "msa_t_MSA_deep*.npy",
        ]
        for pat in patterns:
            if any(out.glob(pat)):
                return True
        return False

    # 1) load PDB/FASTA
    _step_hdr(1, "Load Seqs")
    if force_all or not _has_loaded(pair_id):
        print("Running load …")
        task_load(pair_id, args)
    else:
        print("Seqs already here → skipped")

    # 2) get_msa
    _step_hdr(2, "Get MSA")
    if force_all or not _has_deep_msa(pair_id):
        print("Running get_msa …")
        task_get_msa(pair_id, args.run_job_mode)
    else:
        print("DeepMsa exists → skipped")

    # 3) phylogenetic tree
    _step_hdr(3, "Phylogenetic Tree")
    try:
        tree_path = Path(f"Pipeline/{pair_id}/output_phytree/DeepMsa_tree.nwk")
        if force_all or not tree_path.exists():
            print("Running tree …")
            task_tree(pair_id, args)
        else:
            print("Tree exists → skipped")
    except Exception as e:
        print(f"Tree step skipped: {e}")

    # 4) cluster_msa
    _step_hdr(4, "Cluster MSA")
    if force_all or not _has_cluster_msas(pair_id):
        print("Running cluster_msa …")
        task_cluster_msa(pair_id, "inline", args)
    else:
        print("Clusters exist → skipped")

    # 4b) Renumber kept clusters to 000..K-1 (always attempt; no-op if already contiguous)
    try:
        changed = _renumber_kept_clusters(pair_id)
        if changed:
            print("Kept clusters renumbered to 000..K-1 and CSV updated.")
        else:
            print("Clusters already contiguous / nothing to renumber.")
    except Exception as e:
        print(f"Renumber kept clusters skipped: {e}")

    # 5) AF predictions
    _step_hdr(5, "Run AF Predictions")
    if force_all:
        setattr(args, "force_rerun_AF", "TRUE")  # only for AF stage
    # Default to both versions inside the pipeline unless user explicitly asked AF2-only.
    if getattr(args, "af_ver", None) not in ("2", "3", "both"):
        args.af_ver = "both"
    try:
        print("Invoking AF task (per-outdir auto-skip inside) …")
        task_af(pair_id, args)
    except Exception as e:
        print(f"AF step skipped: {e}")

    # 6) CMAPs (MSA-Transformer)
    _step_hdr(6, "MSA-Transformer CMAPs")
    try:
        shallow_ok = _has_cmaps_shallow_local(pair_id)
        deep_ok    = _has_cmaps_deep_local(pair_id)
        if force_all or not (shallow_ok and deep_ok):
            need_sh = "missing" if not shallow_ok else "ok"
            need_dp = "missing" if not deep_ok else "ok"
            print(f"Running MSA-Transformer … (shallow: {need_sh}, deep: {need_dp})")
            task_cmap_msa_transformer(pair_id, "inline")
        else:
            print("Shallow & deep CMAPs present → skipped")
    except Exception as e:
        print(f"MSA-Transformer step skipped: {e}")

    # 6b) CMAPs (CCMpred) — best-effort
    _step_hdr(7, "CCMpred CMAPs")
    try:
        print("Invoking CCMpred (per-tag auto-skip inside) …")
        task_cmap_ccmpred(pair_id, "inline")
    except Exception as e:
        print(f"CCMpred step skipped: {e}")

    # 7) ESMFold (esm2/esm3 or user-specified)
    _step_hdr(8, "ESMFold Predictions")
    wanted_models = ["esm2", "esm3"] if getattr(args, "esm_model", None) in (None, "both") else [args.esm_model]
    for model in wanted_models:
        try:
            a2 = deepcopy(args); a2.esm_model = model
            if force_all or not _has_esm_model(pair_id, model):
                print(f"Running esmfold({model}) …")
                task_esmfold(pair_id, a2)
            else:
                print(f"esmfold({model}) outputs exist → skipped")
        except Exception as e:
            print(f"esmfold({model}) step skipped: {e}")

    # 8) ΔG energies (PyRosetta) — best-effort
    _step_hdr(9, "ΔG Energies")
    try:
        if force_all or not _has_deltaG(pair_id):
            print("Running ΔG …")
            task_deltaG(pair_id)
        else:
            print("ΔG exists → skipped")
    except Exception as e:
        print(f"ΔG step skipped: {e}")

    # 9) plots (includes tree clusters if available)
    _step_hdr(10, "Per-Pair Plots")
    try:
        ap = deepcopy(args)
        ap.plot_scope = "pair"
        ap.plot_trees = True
        print("Generating plots …")
        task_plot(pair_id, ap)
    except Exception as e:
        print(f"Plot step skipped: {e}")

    # 10) per-pair post-process (writes Pipeline/<pair>/Analysis/df_*.csv)
    _step_hdr(11, "Per-Pair Postprocess")
    try:
        print("Postprocessing (per-pair) …")
        post_processing_analysis(
            force_rerun=force_all or _bool_from_tf(getattr(args, "force_rerun_postprocess", "FALSE")),
            pairs=[pair_id]
        )
    except Exception as e:
        print(f"Postprocess step skipped: {e}")

    # 11) OPTIONAL: per-pair HTML
    _step_hdr(12, "Per-Pair HTML")
    if getattr(args, "reports", "none") in ("html", "all"):
        try:
            print("Generating per-pair HTML …")
            env = _jupyter_env_for_scratch()  # use scratch env
            cmd = (
                f"{shlex.quote(sys.executable)} Analysis/NotebookGen/generate_notebooks.py "
                f"{shlex.quote(pair_id)} --kernel python3"
            )
            rc = subprocess.run(cmd, shell=True, check=False, env=env).returncode
            if rc != 0:
                print(f"[html] nbconvert failed (rc={rc}) — will be retried later.")
        except Exception as e:
            print(f"HTML step skipped: {e}")
    else:
        print("Per-pair HTML disabled → skipped")


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

    # --- Clustering options ---
    p.add_argument("--cluster_alg", default="tree", choices=["hdbscan", "tree", "ahc"],
                   help="Clustering algorithm for --run_mode cluster_msa (default: tree).")
    p.add_argument("--cluster_tree_path", default=None,
                   help="Path to Newick tree for --cluster_alg tree (e.g., Pipeline/<pair>/output_phytree/DeepMsa_tree.nwk).")
    p.add_argument("--cluster_max_clusters", type=int, default=100,
                   help="Maximum clusters to output (default: 100).")
    p.add_argument("--cluster_min_output_size", type=int, default=200,
                   help="Minimum sequences per cluster (default: 200).")
    p.add_argument("--cluster_min_neff", type=int, default=50,
                   help="Minimum Neff per cluster (default: 50).")
    p.add_argument("--cluster_neff_id_thresh", type=float, default=0.8,
                   help="Identity threshold for Neff (default: 0.8).")
    p.add_argument("--cluster_frac_gaps_cutoff", type=float, default=0.7,
                   help="Drop sequences with ≥ this gap fraction before clustering (default: 0.7).")
    p.add_argument("--cluster_sample_cap", type=int, default=5000,
                   help="Sample size used by HDBSCAN/AHC (default: 12000).")
    p.add_argument("--cluster_sample_seed", type=int, default=12345,
                   help="Random seed for sampling (default: 12345).")
    # HDBSCAN knobs (kept for compatibility)
    p.add_argument("--hdbscan_min_cluster_size", type=int, default=200)
    p.add_argument("--hdbscan_min_samples", default="auto")
    p.add_argument("--hdbscan_cluster_selection", choices=["eom", "leaf"], default="eom")

    # AlphaFold options
    p.add_argument("--allow_inline_af", action="store_true",
                    help="Allow AF2 to run inline even if not in a Slurm session (expert only).")
    p.add_argument("--af_ver", default="both", choices=["2", "3", "both"],
                    help="Which AlphaFold to run for --run_mode run_AF")  # default do both AF2 and AF3

    p.add_argument("--force_rerun_AF", default="FALSE", choices=["TRUE", "FALSE"],
                    help="Run new AF2/AF3 predictions even if outputs exist. Default FALSE (skip if found).")

    # ESMFold options
    p.add_argument("--cluster_sample_n", type=int, default=10)
    p.add_argument("--esm_model", default="both", choices=["esm2", "esm3", "both"])
    p.add_argument("--esm_device", default="auto", choices=["auto", "cpu", "cuda", "mps"])
    p.add_argument("--esm_gpu_len_threshold", type=int, default=ESM_MAX_LIGHT_SEQ_LEN,
               help="If max chain length ≥ this, run ESM on CPU to avoid CUDA OOM (default: 800).")

    # ---- CCMpred options ----
    p.add_argument("--ccmpred_bin", default=CCMPRED_EXE, help="Path to CCMpred binary")
    p.add_argument("--ccmpred_threads", type=int, default=8, help="Threads for CCMpred (-t)")

    # --- Tree options ---
    p.add_argument("--tree_max_seqs", type=int, default=5000,
                   help="Max sequences to include in the Deep tree. If stratified, "
                        "it’s bumped up to the number of represented clusters when needed.")
    p.add_argument("--tree_use_clusters", type=str2bool, nargs="?", const=True, default=True,
                   help="If TRUE, do stratified sampling: ≥1 seq per ShallowMsa_### cluster "
                        "that overlaps the Deep MSA, then fill the rest at random.")
    p.add_argument("--tree_seed", type=int, default=123, help="Random seed for tree sampling.")

    # Post-processing options, computing metrics
    p.add_argument("--postprocess", default="FALSE", help="If TRUE, run post-processing after the selected task(s).")
    p.add_argument("--force_rerun_postprocess", default="FALSE", help="If TRUE, recompute per-pair caches.")
    p.add_argument("--cached_only", default="FALSE", choices=["TRUE", "FALSE"],
               help="Skip per-pair metric recomputation and only build reports/HTML.")

    # Plotting
    p.add_argument("--plot_trees",   type=str2bool, nargs="?", const=True, default=False)
    p.add_argument("--plot_scope", choices=["pair", "global", "both"], default="both",
                help="In --run_mode plot: generate pair-specific plots only, global plots only, or both.")
    p.add_argument('--plot3dformat', default='static', choices = ['static', 'interactive', 'both'],
                help = '3D figure format: PNG (static), HTML (interactive), or both')

    # Output html tables and pages
    p.add_argument("--reports", default="none", choices=["none", "tables", "html", "all"],
                   help="Post-run reporting: 'tables' builds CSV+HTML tables, 'html' builds per-pair pages, 'all' does both.")
    p.add_argument("--html_pairs", nargs="+", default=["ALL"],
                   help="Pairs to render per-pair HTML for (defaults to ALL).")
    p.add_argument("--html_mode", choices=["inline", "copy", "link"], default="inline",
                   help=(
                       "How per-pair HTML includes images: "
                       "'inline' = base64 embed (single self-contained file), "
                       "'copy' = copy PNGs to publish dir and link them, "
                       "'link' = assume figures already published (no copying)."))


    # Pipeline-wide force flag
    p.add_argument("--force_rerun", default="FALSE", choices=["TRUE", "FALSE"],
        help="If TRUE (only for --run_mode msaclust_pipeline), run every step regardless of outputs.")

    p.add_argument("--clean_level", default="derived", choices=["derived", "all"],
        help="derived: remove computed outputs but keep base inputs; all: remove the entire pair folder")
    # Targeted clean (optional): choose which operation outputs to remove
    # If provided, this overrides --clean_level and only removes the selected targets.
    from config import PAIR_TARGET_DIRS  # add near other imports

    # ... inside build_argparser / parser setup ...
    p.add_argument(
        "--clean_targets",
        nargs="+",
        default=None,
        choices=sorted(PAIR_TARGET_DIRS.keys()),
        help=("Targeted clean of specific steps. "
              "Examples: --clean_targets msa_clusters | cmap | esmfold | alphafold | get_msa | plots | logs | tmp")
    )


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
    p.add_argument("--sbatch_nodelist", default=None,
               help="Comma-separated Slurm nodelist, e.g. salmon-[01-10],dogfish-[01-02]")
    p.add_argument("--sbatch_gres_heavy", default="gpu:l40s:1",
                   help="Slurm --gres to use automatically when chain length >= threshold (default: gpu:l40s:1).")

    # Slurm options for html output
    p.add_argument("--html_cpus", type=int, default=2)
    p.add_argument("--html_mem",  default="4G")
    p.add_argument("--html_time", default="02:00:00")
    p.add_argument("--html_partition", default=None)

    args = p.parse_args()
    # allow: python run_foldswitch_pipeline.py --run_mode help
    if args.run_mode == "help":
        p.print_help()
        sys.exit(0)
    if args.run_mode == "run_esmfold" and args.esm_model is None:
        args.esm_model = "esm2"  # sensible default

    # resolve list of pairs
    if args.foldpair_ids == ["ALL"]:
        foldpairs = list_protein_pairs() # [s.replace("\t", "_") for s in raw]
    else:
        foldpairs = args.foldpair_ids

    # --- Before the loop ---
    if args.run_mode == "postprocess":
        # Pass the full list so the function can filter to pairs that are ready
        task_postprocess(foldpairs, args)
        return  # or sys.exit(0)

    # Run once for all pairs
    scope = getattr(args, "plot_scope", "both")
    if args.run_mode == "plot":
        task_plot(pair_id=foldpairs, args=args)        # Run global plots ONCE
        return

    for pair_id in foldpairs:
        print(f"=== {args.run_mode} :: {pair_id} ===", flush=True)

        if type(pair_id) in [tuple, list]:  # convert format
            pair_id = f"{pair_id[0]}_{pair_id[1]}"

        # NEW: generic per-pair submit for single-step run modes
        if args.run_job_mode == "sbatch" and args.run_mode in HEAVY_PAIR_MODES:
            extras = []

            # ----- Slurm resource hints -> use *pipeline* args (sbatch_*) NOT raw sbatch flags -----
            try:
                hint = SBATCH_HINTS.get(args.run_mode, {"time": "04:00:00", "mem": "8G", "gpus": 0, "cpus": 4})
            except NameError:
                hint = {"time": "04:00:00", "mem": "8G", "gpus": 0, "cpus": 4}

            if hint.get("time"):
                extras += [f"--sbatch_time {hint['time']}"]
            if hint.get("mem"):
                extras += [f"--sbatch_mem {hint['mem']}"]
            if hint.get("cpus"):
                extras += [f"--sbatch_cpus {hint['cpus']}"]
            # your argparse exposes --sbatch_gres (not --gpus); map gpus -> gres=gpu:<n>
            if hint.get("gpus", 0):
                extras += [f"--sbatch_gres gpu:{hint['gpus']}"]
            # (optional) partition/qos/account/constraint if you set them in SBATCH_HINTS
            for k_src, k_dst in [
                ("partition", "sbatch_partition"),
                ("qos", "sbatch_qos"),
                ("account", "sbatch_account"),
                ("constraint", "sbatch_constraint"),
            ]:
                if hint.get(k_src):
                    extras += [f"--{k_dst} {hint[k_src]}"]

            # ----- Mode-specific passthrough flags you already had -----
            if args.run_mode == "run_esmfold":
                if getattr(args, "esm_model", None):
                    extras += [f"--esm_model {args.esm_model}"]
                if getattr(args, "esm_device", None):
                    extras += [f"--esm_device {args.esm_device}"]

            if args.run_mode == "tree":
                # forward tree knobs into the inner inline run
                extras += [f"--tree_max_seqs {int(args.tree_max_seqs)}"]
                extras += [f"--tree_seed {int(args.tree_seed)}"]
                # (optional) only if/when tree_use_clusters starts being honored in task_tree
                # extras += [f"--tree_use_clusters {str(args.tree_use_clusters)}"]

            if args.run_mode == "cluster_msa":
                # pass clustering knobs through to the inner inline run
                extras += [
                    f"--cluster_alg {args.cluster_alg}",
                    # tree path is optional; task_cluster_msa will auto-pick default if not set
                ]
                if getattr(args, "cluster_tree_path", None):
                    extras += [f"--cluster_tree_path {args.cluster_tree_path}"]

                extras += [
                    f"--cluster_max_clusters {int(args.cluster_max_clusters)}",
                    f"--cluster_min_output_size {int(args.cluster_min_output_size)}",
                    f"--cluster_min_neff {int(args.cluster_min_neff)}",
                    f"--cluster_neff_id_thresh {float(args.cluster_neff_id_thresh)}",
                    f"--cluster_frac_gaps_cutoff {float(args.cluster_frac_gaps_cutoff)}",
                    f"--cluster_sample_cap {int(args.cluster_sample_cap)}",
                    f"--cluster_sample_seed {int(args.cluster_sample_seed)}",
                    # keep HDBSCAN knobs for compatibility (AHC ignores them; tree path ignores them)
                    f"--hdbscan_min_cluster_size {int(args.hdbscan_min_cluster_size)}",
                    f"--hdbscan_min_samples {str(args.hdbscan_min_samples)}",
                    f"--hdbscan_cluster_selection {str(args.hdbscan_cluster_selection)}",
                ]

            if args.run_mode == "postprocess":
                if getattr(args, "reports", "none") != "none":
                    extras += [f"--reports {args.reports}"]

            _submit_pair_job(args.run_mode, pair_id, args, " ".join(extras))
            continue


        if args.run_mode == "load":
            task_load(pair_id, args)

        elif args.run_mode == "get_msa":
            task_get_msa(pair_id, args.run_job_mode)

        elif args.run_mode == "cluster_msa":
            task_cluster_msa(pair_id, args.run_job_mode, args)

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
            task_tree(pair_id, args)
#            task_tree(pair_id, args.run_job_mode)

        elif args.run_mode == "compute_deltaG":
            task_deltaG(pair_id)

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
    if args.run_job_mode == "sbatch" and args.run_mode in {
        "msaclust_pipeline", "run_esmfold", "run_cmap_msa_transformer", "run_cmap_ccmpred", "postprocess"
    }:
        print("[submit-only] Per-pair jobs have been submitted. Run reports later.", flush=True)
        return


    # Otherwise build whatever the user asked for
    if args.reports != "none":
        task_postprocess(foldpairs, args)
        print("[done html outputs]", flush=True)

if __name__ == "__main__":
    main()
