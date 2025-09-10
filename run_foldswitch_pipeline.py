#!/usr/bin/env python3
import argparse
import subprocess
import shlex
import sys
from glob import glob
from pathlib import Path
from typing import List, Tuple
from copy import deepcopy

from config import *
from utils.utils import pair_str_to_tuple, ensure_dir, write_pair_pipeline_script
from utils.protein_utils import read_msa, greedy_select, extract_protein_sequence, load_seq_and_struct
from utils.msa_utils import write_fasta, load_fasta, build_pair_seed_a3m_from_pair  # your existing writer
from utils.phytree_utils import phytree_from_msa


# ------------------------- helpers -------------------------

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

def _cluster_files(pair_id: str) -> Tuple[str, List[str]]:
    """Return (deep_msa, list_of_cluster_msas) for a pair."""
    base = f"Pipeline/{pair_id}/output_msa_cluster"
    deep = f"{base}/DeepMsa.a3m"
    clusters = sorted(glob(f"{base}/ShallowMsa_*.a3m"))
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
    outdir = f"Pipeline/{pair_id}"
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
        cmd = (f"bash ./Pipeline/RunAF2_Colabfold.sh --python "
            f"./get_msa.py '{seed_a3m}' '{out_dir}' --name DeepMsa")
        _run(cmd, "inline")


def task_cluster_msa(pair_id: str, run_job_mode: str) -> None:
    # Your existing script; unchanged except paths
    cmd = (
        f"python3 ./ClusterMSA_moriah.py "
        f"--keyword ShallowMsa "
        f"-i Pipeline/{pair_id}/output_get_msa/DeepMsa.a3m "
        f"-o Pipeline/{pair_id}/output_msa_cluster"
    )
    _run(cmd, run_job_mode)

def task_cmap_esm(pair_id: str, run_job_mode: str) -> None:
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

    # Correct CLI: no -i/-o, no --esm_version
    device = args.esm_device or "auto"
    cmd = f"python3 ./ESMFoldHF.py -input {pair_id} --model {args.esm_model} --device {device}"
    _run(cmd, args.run_job_mode)


def task_af2(pair_id: str, args: argparse.Namespace) -> None:
    """
    Run AF2 (colabfold_batch via wrapper) on:
      - the two chain FASTAs, and
      - EVERY MSA: DeepMsa.a3m + all ShallowMsa_*.a3m
    No sampling; MSAs are used as-is.
    """
    if _is_windows():
        raise SystemExit("AlphaFold2 must run on Linux.")

    pair_dir = f"Pipeline/{pair_id}"
    out_root = f"{pair_dir}/output_AF/AF2"
    ensure_dir(out_root)

    # (1) Make sure chain FASTAs exist
    foldA, foldB = pair_str_to_tuple(pair_id)
    pdbids     = [foldA[:-1], foldB[:-1]]
    pdbchains  = [foldA[-1],  foldB[-1]]
    ensure_dir(os.path.join(pair_dir, "fasta_chain_files"))
    ensure_chain_fastas(pair_dir, pdbids, pdbchains)

    chain_fastas = [
        os.path.join(pair_dir, "fasta_chain_files", f"{pdbids[0]}{pdbchains[0]}.fasta"),
        os.path.join(pair_dir, "fasta_chain_files", f"{pdbids[1]}{pdbchains[1]}.fasta"),
    ]

    # (2) Collect ALL MSAs: Deep + clusters
    deep, clusters = _cluster_files(pair_id)  # DeepMsa.a3m, list of ShallowMsa_*.a3m
    a3ms = []
    if os.path.isfile(deep):
        a3ms.append(deep)
    a3ms.extend([c for c in clusters if os.path.isfile(c)])

    # (3) Build commands: one run per input; clean output structure
    def _out_dir_for(inp_path: str) -> str:
        stem = Path(inp_path).stem  # e.g., "1fzpA", "DeepMsa", "ShallowMsa_003"
        d = os.path.join(out_root, stem)
        ensure_dir(d)
        return d

    def _cmd_for(inp_path: str) -> str:
        out_dir = _out_dir_for(inp_path)
        # Wrapper activates af2-venv and runs colabfold_batch
        # You can add/tune flags (models/recycles/seeds) as needed.
        return (
            f"bash ./Pipeline/RunAF2_Colabfold.sh "
            f"--num-models 1 --num-recycle 1 "
            f"{shlex.quote(inp_path)} {shlex.quote(out_dir)}"
        )

    inputs = chain_fastas + a3ms
    if not inputs:
        print(f"[warn] No AF2 inputs for {pair_id}", flush=True)
        return

    inside_slurm = _in_slurm_session()
    for ip in inputs:
        cmd = _cmd_for(ip)
        if args.run_job_mode == "sbatch" or (not inside_slurm and not args.allow_inline_af):
            log_path = os.path.join(pair_dir, f"run_AF_for_{pair_id}_{Path(ip).stem}.out")
            _run(f"sbatch -o '{log_path}' --wrap {shlex.quote(cmd)}", "sbatch")
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
        task_af2(pair_id, args)  # this will run inline inside current shell
        task_cmap_esm(pair_id, "inline")
        for model in ("esm2", "esm3"):
            a2 = deepcopy(args)
            if hasattr(a2, "esm_model"):
                a2.esm_model = model
            elif hasattr(a2, "esm_version"):
                a2.esm_version = model
            else:
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
    p = argparse.ArgumentParser("FoldSwitch pipeline (clean)")
    p.add_argument("--run_mode",
                   required=True,
                   choices=["load", "get_msa", "cluster_msa", "run_cmap_esm",
                            "run_esmfold", "run_AF", "tree", "plot", "compute_deltaG", "clean",
                            "msaclust_pipeline"])  # Last one is the full pipeline for a pair
    p.add_argument("--foldpair_ids", nargs="+", default=["ALL"],
                   help="e.g. 1dzlA_5keqF (default: ALL in data list)")
    p.add_argument("--run_job_mode", default="inline", choices=["inline", "sbatch"])

    # AlphaFold options
    p.add_argument("--allow_inline_af", action="store_true",
                   help="Allow AF2 to run inline even if not in a Slurm session (expert only).")


    # ESMFold options
    p.add_argument("--cluster_sample_n", type=int, default=10)
    p.add_argument("--esm_model", default=None, choices=["esm2", "esm3"])
    p.add_argument("--esm_version", default=None, choices=["esm2", "esm3"],  # back-compat alias
                   help="Deprecated alias of --esm_model")
    p.add_argument("--esm_device", default="auto", choices=["auto", "cpu", "cuda", "mps"])

    # plotting
    p.add_argument("--global_plots", action="store_true")
    p.add_argument("--plot_trees", action="store_true")

    args = p.parse_args()
    if args.esm_model is None and args.esm_version is not None:
        args.esm_model = args.esm_version
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

        elif args.run_mode == "run_cmap_esm":
            task_cmap_esm(pair_id, args.run_job_mode)

        elif args.run_mode == "run_esmfold":
            task_esmfold(pair_id, args)

        elif args.run_mode == "run_AF":
            task_af2(pair_id, args)

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
