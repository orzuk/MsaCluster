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
                    shutil.rmtree(p, ignore_errors=True)

        # special handling for logs: also remove root-level files
        if any(t in ("logs",) for t in targets):
            for pat in ["*.out", "*.log", "*.err", "RunAF.out", "CmapESM.out", "run_pipeline_for_*.out"]:
                for f in pair_dir.glob(pat):
                    if f.is_file():
                        print(f"[clean:logs] rm {f}")
                        if not dry:
                            try: f.unlink()
                            except: pass

        # special handling for tmp: remove any tmp_* dirs
        if any(t in ("tmp",) for t in targets):
            for p in pair_dir.glob("tmp_*"):
                if p.is_dir():
                    print(f"[clean:tmp] rm -rf {p}")
                    if not dry:
                        shutil.rmtree(p, ignore_errors=True)

        _prune_stale_fastas(pair_id, pair_dir, dry)
        return

    # ----- EXISTING LEVEL-BASED CLEAN (unchanged semantics) -----
    if level == "all":
        print(f"[clean:ALL] rm -rf {pair_dir}")
        if not dry:
            shutil.rmtree(pair_dir, ignore_errors=True)
        return

    # derived-only clean (keep base inputs)
    rm_dirs_current = [
        "output_get_msa", "output_msa_cluster", "output_AF",
        "output_cmaps", "output_esm_fold", "output_phytree",
        "Analysis", "tmp_msa_files", "tmp_esmfold", "jobs", "logs",
    ]
    rm_dirs_legacy = ["output_cmap_esm", "esm_cmap_output", "AF_preds", "chain_pdb_files", "fasta_chain_files"]
    rm_globs = ["run_pipeline_for_*.out", "RunAF.out", "CmapESM.out", "*.out", "*.log", "*.err"]

    for d in rm_dirs_current + rm_dirs_legacy:
        p = pair_dir / d
        if p.exists():
            print(f"[clean] rm -rf {p}")
            if not dry:
                shutil.rmtree(p, ignore_errors=True)

    for p in pair_dir.glob("tmp_*"):
        if p.is_dir():
            print(f"[clean] rm -rf {p}")
            if not dry:
                shutil.rmtree(p, ignore_errors=True)

    for pat in rm_globs:
        for f in pair_dir.glob(pat):
            if f.is_file():
                print(f"[clean] rm {f}")
                if not dry:
                    try: f.unlink()
                    except: pass

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
            f"--model msa_t --keyword '' --clean matched"
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

    # 1) Per-pair metrics (heavy) — skip if --cached_only TRUE
    cached_only = _bool_from_tf(getattr(args, "cached_only", "FALSE"))
    if not cached_only:
        try:
            force = _bool_from_tf(getattr(args, "force_rerun_postprocess", "FALSE"))
            post_processing_analysis(force_rerun=force, pairs=pairs_arg)
        except Exception as e:
            print(f"[postprocess] WARN post_processing_analysis: {e}")

    # === Build/refresh global cached tables (no heavy compute here) ===
    try:
        from Analysis.postprocess_global import build_or_load_global_tables
        force_global = _bool_from_tf(getattr(args, "force_rerun_postprocess", "FALSE"))
        summary_csv, detail_csv = build_or_load_global_tables(force=force_global)
        print(f"[global] summary: {summary_csv}")
        print(f"[global] detailed: {detail_csv}")
    except Exception as e:
        print(f"[postprocess] WARN global tables: {e}")

    # 2) + 3) Global CSVs + HTML tables
    if args.reports in ("tables", "all"):
        # 2) Build unified CSVs ...
        try:
            build_unified_tables_from_cluster_dfs(pairs=pairs_arg, write_out=True)
            print(f"[reports] unified CSVs written:\n  {SUMMARY_RESULTS_TABLE}\n  {DETAILED_RESULTS_TABLE}")
        except Exception as e:
            print(f"[reports] WARN building unified CSVs: {e}")

        # 3.1) Generate HTML tables (single source-of-truth)
        try:
            os.makedirs(TABLES_RES, exist_ok=True)
            write_global_tables(force_rerun_csv=False, fade_min_clusters=2)
        except Exception as e:
            print(f"[reports] WARN write_global_tables: {e}")

        # 3.2) Global plots page (unchanged)
        try:
            gen_html_for_global_plots(
                images_dir=FIGURE_RES_DIR,
                output_html=os.path.join("docs", "pairs_global_analysis.html"))
            # Mirror to repo root for GitHub Pages (optional)
            src = os.path.join("docs", "pairs_global_analysis.html")
            dst = os.path.join(MAIN_DIR, "pairs_global_analysis.html")
            if os.path.isfile(src):
                shutil.copy2(src, dst)
                print(f"[reports] copied {src} -> {dst}")
        except Exception as e:
            print(f"[reports] WARN building global plots page: {e}")


    # 4) Per-pair HTML pages (new generator)
    if args.reports in ("html", "all"):
        try:
            # keep only pairs that already have per-pair postprocess outputs
            ready = []
            for p in (
            norm_pairs if pairs_arg is not None else [d.name for d in Path(DATA_DIR).iterdir() if d.is_dir()]):
                if os.path.isfile(f"Pipeline/{p}/Analysis/df_af.csv"):
                    ready.append(p)
                else:
                    print(f"[html] skip {p}: missing Pipeline/{p}/Analysis/df_af.csv")

            if not ready:
                print("[html] no ready pairs; skipping per-pair HTML generation")
            else:
                mode = getattr(args, "html_mode", "inline")
                outdir = OUTPUT_PATH_NOTEBOOKS  # from config.py
                script = Path("TableResults/gen_pair_html.py")

                # If user said ALL at CLI we respect that; otherwise pass explicit list
                pairs_arg_str = "ALL" if pairs_arg is None else ",".join(ready)

                cmd = (
                    f"{shlex.quote(sys.executable)} {shlex.quote(str(script))} "
                    f"--pairs {shlex.quote(pairs_arg_str)} "
                    f"--mode {shlex.quote(mode)} "
                    f"--output_dir {shlex.quote(str(outdir))}"
                )

                if args.run_job_mode == "sbatch":
                    jobs_dir = Path("Pipeline/jobs");
                    jobs_dir.mkdir(parents=True, exist_ok=True)
                    log = jobs_dir / f"genhtml_{mode}.out"
                    cpus = getattr(args, "html_cpus", 2)
                    mem = getattr(args, "html_mem", "4G")
                    time = getattr(args, "html_time", "02:00:00")
                    part = getattr(args, "html_partition", None)

                    sbatch = [
                        "sbatch",
                        f"--job-name=genhtml",
                        f"--cpus-per-task={cpus}",
                        f"--mem={mem}",
                        f"--time={time}",
                        f"--output={log}",
                        "--parsable",
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

    # 1) load PDB/FASTA (cheap; idempotent)
    task_load(pair_id, args)

    # 2) get_msa
    if force_all or not _has_deep_msa(pair_id):
        print("[pipeline] get_msa → running")
#        task_get_msa(pair_id, "inline")
        task_get_msa(pair_id, args.run_job_mode)

    else:
        print("[pipeline] get_msa → skip (DeepMsa exists)")

    # 3) cluster_msa
    if force_all or not _has_cluster_msas(pair_id):
        print("[pipeline] cluster_msa → running")
        task_cluster_msa(pair_id, "inline", args)
    else:
        print("[pipeline] cluster_msa → skip (clusters exist)")

    # 4) AF (AF2/AF3/both). Always invoke; the task itself skips completed outdirs unless forced.
    if force_all:
        setattr(args, "force_rerun_AF", "TRUE")  # only for AF stage
    print("[pipeline] AF → running (task will skip per-outdir if already complete)")
    # Default to both versions inside the full pipeline unless user explicitly asked AF2-only.
    if getattr(args, "af_ver", None) not in ("2", "3", "both"):
        args.af_ver = "both"
    try:
        task_af(pair_id, args)
    except Exception as e:
        print(f"[pipeline] AF → skipped: {e}")        

    # 5) cmaps (MSA-Transformer)
    try:
        if force_all or not _has_cmaps(pair_id):
            print("[pipeline] cmap_msa_transformer → running")
            task_cmap_msa_transformer(pair_id, "inline")
        else:
            print("[pipeline] cmap_msa_transformer → skip (cmaps exist)")
    except Exception as e:
        print(f"[pipeline] cmap_msa_transformer → skipped: {e}")

    # 5b) cmaps (CCMpred) — best-effort
    try:
        print("[pipeline] cmap_ccmpred → running")
        task_cmap_ccmpred(pair_id, "inline")
    except Exception as e:
        print(f"[pipeline] cmap_ccmpred → skipped: {e}")

    # 6) ESMFold (esm2/esm3 or user-specified)    
    wanted_models = ["esm2", "esm3"] if getattr(args, "esm_model", None) in (None, "both") else [args.esm_model]
    for model in wanted_models:
        try:
            a2 = deepcopy(args); a2.esm_model = model
            if force_all or not _has_esm_model(pair_id, model):
                print(f"[pipeline] esmfold({model}) → running")
                task_esmfold(pair_id, a2)
            else:
                print(f"[pipeline] esmfold({model}) → skip (outputs exist)")
        except Exception as e:
            print(f"[pipeline] esmfold({model}) → skipped: {e}")

    # 7) phylogenetic tree (needed for tree plots)
    try:
        tree_path = Path(f"Pipeline/{pair_id}/output_phytree/DeepMsa_tree.nwk")
        if force_all or not tree_path.exists():
            print("[pipeline] tree → running")
            task_tree(pair_id, args)
        else:
            print("[pipeline] tree → skip (exists)")
    except Exception as e:
        print(f"[pipeline] tree → skipped: {e}")

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
        ap.plot_scope = "pair"
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
    if getattr(args, "reports", "none") in ("html", "all"):
        try:
            print("[pipeline] html → running (per-pair)")
            env = _jupyter_env_for_scratch()  # use scratch env
            cmd = (
                f"{shlex.quote(sys.executable)} Analysis/NotebookGen/generate_notebooks.py "
                f"{shlex.quote(pair_id)} --kernel python3"
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

    # --- Clustering options ---
    p.add_argument("--cluster_alg", default="ahc", choices=["hdbscan", "tree", "ahc"],
                   help="Clustering algorithm for --run_mode cluster_msa (default: ahc).")
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
    p.add_argument("--cluster_frac_gaps_cutoff", type=float, default=0.6,
                   help="Drop sequences with ≥ this gap fraction before clustering (default: 0.6).")
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
    p.add_argument("--esm_model", default=None, choices=["esm2", "esm3", "both"])
    p.add_argument("--esm_device", default="auto", choices=["auto", "cpu", "cuda", "mps"])
    p.add_argument("--esm_gpu_len_threshold", type=int, default=800,
               help="If max chain length ≥ this, run ESM on CPU to avoid CUDA OOM (default: 800).")

    # ---- CCMpred options ----
    p.add_argument("--ccmpred_bin",
                   default="/sci/labs/orzuk/orzuk/github/CCMpred/bin/ccmpred",
                   help="Path to CCMpred binary")
    p.add_argument("--ccmpred_threads", type=int, default=8,
                   help="Threads for CCMpred (-t)")

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
