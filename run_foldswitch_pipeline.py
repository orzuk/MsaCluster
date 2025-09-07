#!/usr/bin/env python3
import argparse
import os
import sys
import subprocess
import platform
from glob import glob
from pathlib import Path
from typing import List, Tuple

from config import *
from utils.utils import pair_str_to_tuple
from utils.protein_utils import read_msa, greedy_select, extract_protein_sequence
from utils.msa_utils import write_fasta, load_fasta  # your existing writer


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

def _is_windows() -> bool:
    return platform.system().lower().startswith("win")

def _ensure_dir(p: str) -> None:
    Path(p).mkdir(parents=True, exist_ok=True)

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
    _ensure_dir(outdir)

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

def task_load(pair_id: str, run_job_mode: str) -> None:
    # Keep your existing loader through protein_utils; left as-is
    from utils.protein_utils import load_seq_and_struct
    foldA, foldB = pair_str_to_tuple(pair_id)
    cur_family_dir = f"Pipeline/{pair_id}"
    _ensure_dir(cur_family_dir)
    load_seq_and_struct(cur_family_dir, [foldA[:-1], foldB[:-1]], [foldA[-1], foldB[-1]])

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
    _ensure_dir(outdir)
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
    _ensure_dir(model_dir)

    # Correct CLI: no -i/-o, no --esm_version
    device = args.esm_device or "auto"
    cmd = f"python3 ./ESMFoldHF.py -input {pair_id} --model {args.esm_model} --device {device}"
    _run(cmd, args.run_job_mode)


def task_af2(pair_id: str, args: argparse.Namespace) -> None:
    """
    Submit AF2 for this pair via Slurm, after ensuring the two chain FASTAs exist:
      Pipeline/<pair>/fasta_chain_files/<pdb><chain>.fasta
    The sbatch script is your existing: ./Pipeline/RunAF_params.sh
    """
    if _is_windows():
        raise SystemExit("AlphaFold2 must run on the Linux cluster. Use --run_job_mode sbatch.")

    if args.run_job_mode != "sbatch":
        raise SystemExit("AlphaFold2 must be submitted via Slurm. Re-run with: --run_job_mode sbatch")

    # ensure chain FASTAs (reuses root FASTAs if single-record; else extracts from PDB)
    foldA, foldB = pair_str_to_tuple(pair_id)        # e.g. '1fzpD','2frhA'
    pdbids     = [foldA[:-1], foldB[:-1]]            # ['1fzp','2frh']
    pdbchains  = [foldA[-1],  foldB[-1]]             # ['D','A']
    pair_dir   = f"Pipeline/{pair_id}"
    _ensure_dir(os.path.join(pair_dir, "fasta_chain_files"))
    ensure_chain_fastas(pair_dir, pdbids, pdbchains)

    # output root (AF2 predictions live here)
    out_root = f"{pair_dir}/output_AF/AF2"
    _ensure_dir(out_root)

    # submit your existing AF2 params script
    log_path = f"{pair_dir}/run_AF_for_{pair_id}.out"
    cmd = f"sbatch -o '{log_path}' ./Pipeline/RunAF_params.sh {pair_id}"
    _run(cmd, "sbatch")


def task_tree(pair_id: str, run_job_mode: str) -> None:
    from utils.phytree_utils import phytree_from_msa
    msa_file = f"Pipeline/{pair_id}/output_get_msa/DeepMsa.a3m"
    out = f"Pipeline/{pair_id}/output_phytree/DeepMsa_tree.nwk"
    _ensure_dir(os.path.dirname(out))
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
    _ensure_dir(out_dir)
    pA, pB = pair_str_to_tuple(pair_id)
    pdb_pair = [(f"Pipeline/{pair_id}/{pA[:-1]}.pdb", f"Pipeline/{pair_id}/{pB[:-1]}.pdb")]
    compute_global_and_residue_energies(pdb_pair, [pair_id], out_dir)


# ------------------------- CLI / main -------------------------

def main():
    p = argparse.ArgumentParser("FoldSwitch pipeline (clean)")
    p.add_argument("--run_mode",
                   required=True,
                   choices=["load", "get_msa", "cluster_msa",
                            "run_cmap_esm", "run_esmfold",
                            "run_AF", "tree", "plot", "compute_deltaG"])
    p.add_argument("--foldpair_ids", nargs="+", default=["ALL"],
                   help="e.g. 1dzlA_5keqF (default: ALL in data list)")
    p.add_argument("--run_job_mode", default="inline", choices=["inline", "sbatch"])

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

        else:
            raise ValueError(args.run_mode)

    print("[done]", flush=True)


if __name__ == "__main__":
    main()
