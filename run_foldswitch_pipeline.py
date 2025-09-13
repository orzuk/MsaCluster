#!/usr/bin/env python3
import argparse
import subprocess
import shlex
import sys, warnings
from glob import glob
import json
from pathlib import Path
from typing import List, Tuple
from copy import deepcopy

from config import *
from utils.utils import pair_str_to_tuple, ensure_dir, write_pair_pipeline_script
from utils.protein_utils import read_msa, greedy_select, extract_protein_sequence, load_seq_and_struct
from utils.msa_utils import write_fasta, load_fasta, build_pair_seed_a3m_from_pair  # your existing writer
from utils.phytree_utils import phytree_from_msa

from textwrap import dedent

RUN_MODE_DESCRIPTIONS = {
    "get_msa":          "Download/build deep MSAs for the pair and write DeepMsa.a3m.",
    "cluster_msa":      "Cluster the pairwise MSA into shallow clusters (ShallowMsa_XXX.a3m).",
    "run_AF":           "Run AlphaFold (AF2/AF3/both) per chain × cluster. Use --af_ver {2,3,both}.",
    "run_esmfold":      "Run ESMFold on the pair. Use --esm_model {esm2,esm3,both}.",
    "run_cmap_msa_transformer":      "Run MSA-transformer on the pair to get contact maps.",
    "compute_deltaG":   "Compute ΔG stability metrics (requires PyRosetta).",
    "clean":            "Remove previous outputs for the pair.",
    "msaclust_pipeline":"Full pipeline: get_msa → cluster_msa → AF/ESM (as configured).",
    "help":             "Print this list of run modes with one-line explanations.",
}


# ------------------------- helpers -------------------------

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
        entries = read_msa(src_a3m)  # [(id, seq), ...]
        picked = greedy_select(entries, min(sample_n, len(entries)))
        names, seqs = zip(*picked) if picked else ([], [])
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


def _write_pair_a3m_for_chain(cluster_a3m: str, deep_a3m: str, chain_tag: str, out_path: str) -> bool:
    """
    From a single base alignment (cluster_a3m if given else deep_a3m),
    synthesize a per-chain A3M with the chain as FIRST row in the SAME column space.
    If the chain isn't present in the base A3M, align it to the base query row and project.
    """
    from Bio import Align  # PairwiseAligner (modern replacement)
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
        import numpy as np
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
    ensure_dir(os.path.dirname(out_path))
    with open(out_path, "w") as fh:
        fh.write(f">{chain_tag}\n{chain_aln}\n")
        chain_key_up = _ungap_upper(chain_aln)
        for nm, aln in base_entries:
            if _ungap_upper(aln) == chain_key_up:
                continue
            fh.write(f">{nm}\n{aln}\n")
    return True


# ------------------------- tasks -------------------------
def task_clean(pair_id: str, _args) -> None:
    base = f"Pipeline/{pair_id}"
    to_rm_dirs = [
        "output_get_msa",
        "output_msa_cluster",
        "output_cmaps",
        "output_cmap_esm",   # legacy
        "output_esm_fold",
        "output_AF",
        "output_phytree",
        "fasta_chain_files",
        "AF_preds",          # legacy
    ]
    to_rm_files = [
        "_seed_both.a3m",
        "run_AF_for_{p}.out".format(p=pair_id),
        "get_msa_for_{p}.out".format(p=pair_id),
        "cluster_msa_for_{p}.out".format(p=pair_id),
        "_tmp_ShallowMsa_*.fasta",
    ]
    for d in to_rm_dirs:
        path = os.path.join(base, d)
        if os.path.exists(path):
            print(f"[clean] rm -rf {path}")
            subprocess.run(f"rm -rf {shlex.quote(path)}", shell=True, check=False)
    for pattern in to_rm_files:
        for f in glob(os.path.join(base, pattern)):
            print(f"[clean] rm -f {f}")
            os.remove(f)


def task_load(pair_id: str, run_job_mode: str) -> None:
    # Keep your existing loader through protein_utils; left as-is
    foldA, foldB = pair_str_to_tuple(pair_id)
    cur_family_dir = f"Pipeline/{pair_id}"
    ensure_dir(cur_family_dir)
    load_seq_and_struct(cur_family_dir, [foldA[:-1], foldB[:-1]], [foldA[-1], foldB[-1]])

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
    cmd = (
        f"bash -lc 'cd Pipeline/{pair_id} && "
        f"python3 ../../ClusterMSA_moriah.py "
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

    # cleanup
    tmpdir = f"Pipeline/{pair_id}/tmp_esmfold"
    subprocess.run(f"rm -rf {shlex.quote(tmpdir)}", shell=True, check=False)


def task_af(pair_id: str, args: argparse.Namespace) -> None:
    if _is_windows():
        raise SystemExit("AlphaFold must run on Linux.")

    af_ver = str(getattr(args, "af_ver", "2")).lower()  # "2", "3", or "both"

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
    jobs = []  # (a3m_path, out_dir_base_name) — out_dir_base_name = "DeepMsa" or "ShallowMsa_XXX"
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
    sbatch_opts = "--gres=gpu:a100:1 --cpus-per-task=8 --mem=40G --time=24:00:00"

    def _cmd_for(ver: str, a3m_path: str, out_dir: str) -> str:
        if ver == "2":
            return (
                f"bash ./Pipeline/RunAF2_Colabfold.sh "
                f"{shlex.quote(a3m_path)} {shlex.quote(out_dir)} "
                f"--num-models 5 --num-recycle 1 --model-type alphafold2_ptm --save-all"
            )
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
    versions = ["2"] if af_ver == "2" else (["3"] if af_ver == "3" else ["2", "3"])
    for ver in versions:
        out_root = os.path.join(pair_dir, f"output_AF/AF{ver}")
        ensure_dir(out_root)
        log_dir = os.path.join(out_root, "logs"); ensure_dir(log_dir)

        for a3m_path, base in jobs:
            out_dir = os.path.join(out_root, base); ensure_dir(out_dir)
            cmd = _cmd_for(ver, a3m_path, out_dir)
            if args.run_job_mode == "sbatch" or (not inside_slurm and not args.allow_inline_af):
                stem = f"{Path(a3m_path).stem}"
                log_path = os.path.join(log_dir, f"run_AF{ver}_{pair_id}__{stem}.out")
                _run(f"sbatch {sbatch_opts} -o '{log_path}' --wrap {shlex.quote(cmd)}", "sbatch")
            else:
                _run(cmd, "inline")


def task_tree(pair_id: str, run_job_mode: str) -> None:
    msa_file = f"Pipeline/{pair_id}/output_get_msa/DeepMsa.a3m"
    out = f"Pipeline/{pair_id}/output_phytree/DeepMsa_tree.nwk"
    ensure_dir(os.path.dirname(out))
    phytree_from_msa(msa_file, output_tree_file=out)

def task_plot(pair_id: str, args: argparse.Namespace) -> None:
    # Import PyMOL-consuming code only here
    from utils.protein_plot_utils import make_foldswitch_all_plots
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


# All Pipeline
def task_msaclust_pipeline(pair_id: str, args: argparse.Namespace) -> None:
    """
    If --run_job_mode inline: run steps sequentially in-process.
    If --run_job_mode sbatch: create a single SBATCH job that runs all steps sequentially.
    """
    if args.run_job_mode == "inline":
        # run sequentially here
        task_load(pair_id, "inline")
        task_get_msa(pair_id, "inline")
        task_cluster_msa(pair_id, "inline")
        if args.run_mode == "msaclust_pipeline" and getattr(args, "af_ver", "2") == "2":
            args.af_ver = "both" # Run AF2 and AF3
        task_af(pair_id, args)  # this will run inline inside current shell
        task_cmap_msa_transformer(pair_id, "inline")
        for model in ("esm2", "esm3"):
            a2 = deepcopy(args)
            if hasattr(a2, "esm_model"):
                a2.esm_model = model
            task_esmfold(pair_id, a2)  # inline
        try:
            ap = deepcopy(args)
            ap.global_plots = False
            if hasattr(ap, "plot_trees"): ap.plot_trees = False
            task_plot(pair_id, ap)
        except Exception as e:
            print(f"[plot] skipped: {e}")
        return

    # === sbatch: one big job per pair ===
    script_path = write_pair_pipeline_script(pair_id, args)
    # submit it; we rely on SBATCH header inside the script; we only pass -o already embedded
    cmd = f"sbatch {shlex.quote(script_path)}"
    print(f"[sbatch] {cmd}", flush=True)
    subprocess.run(cmd, shell=True, check=True)



# ------------------------- CLI / main -------------------------

def main():
    p = argparse.ArgumentParser(
#        description="Fold-switching pipeline runner",
        formatter_class=argparse.RawTextHelpFormatter,  # keeps newlines
        epilog=_modes_epilog(),
    )
    p.add_argument("--run_mode",
                   required=True,
                   choices=["load", "get_msa", "cluster_msa", "run_cmap_msa_transformer",
                            "run_esmfold", "run_AF", "tree", "plot", "compute_deltaG", "clean",
                            "msaclust_pipeline", "help"])  # Last one is the full pipeline for a pair
    p.add_argument("--foldpair_ids", nargs="+", default=["ALL"],
                   help="e.g. 1dzlA_5keqF (default: ALL in data list)")
    p.add_argument("--run_job_mode", default="inline", choices=["inline", "sbatch"])

    # AlphaFold options
    p.add_argument("--allow_inline_af", action="store_true",
                   help="Allow AF2 to run inline even if not in a Slurm session (expert only).")
    p.add_argument("--af_ver", default="2", choices=["2", "3", "both"],
                   help="Which AlphaFold to run for --run_mode run_AF")

    # ESMFold options
    p.add_argument("--cluster_sample_n", type=int, default=10)
    p.add_argument("--esm_model", default=None, choices=["esm2", "esm3", "both"])
    p.add_argument("--esm_device", default="auto", choices=["auto", "cpu", "cuda", "mps"])

    # plotting
    p.add_argument("--global_plots", action="store_true")
    p.add_argument("--plot_trees", action="store_true")

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
        with open("data/foldswitch_PDB_IDs_full.txt", "r") as f:
            raw = [line.rstrip("\n") for line in f if line.strip()]
        foldpairs = [s.replace("\t", "_") for s in raw]
    else:
        foldpairs = args.foldpair_ids

    for pair_id in foldpairs:
        print(f"=== {args.run_mode} :: {pair_id} ===", flush=True)

        if args.run_mode == "load":
            task_load(pair_id, args.run_job_mode)

        elif args.run_mode == "get_msa":
            task_get_msa(pair_id, args.run_job_mode)

        elif args.run_mode == "cluster_msa":
            task_cluster_msa(pair_id, args.run_job_mode)

        elif args.run_mode == "run_cmap_msa_transformer":
            task_cmap_msa_transformer(pair_id, args.run_job_mode)


        elif args.run_mode == "run_esmfold":
            if args.esm_model == "both":
                from copy import deepcopy
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

        elif args.run_mode == "plot":
            task_plot(pair_id, args)

        elif args.run_mode == "compute_deltaG":
            task_deltaG(pair_id)

        elif args.run_mode == "clean":  # Remove existing files to run the pipeline clean
            task_clean(pair_id, args)

        elif args.run_mode == "msaclust_pipeline":
            task_msaclust_pipeline(pair_id, args)

        else:
            raise ValueError(args.run_mode)


    print("[done]", flush=True)


if __name__ == "__main__":
    main()
