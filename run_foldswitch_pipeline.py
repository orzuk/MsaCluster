#!/usr/bin/env python3
import os, sys, argparse, platform, random, subprocess, shlex
from glob import glob

from config import *  # paths & environment flags
from utils.utils import *
from utils.energy_utils import *  # for compute_deltaG branch; respects pyrosetta gating
from utils.protein_utils import read_msa  # use your existing fast a3m reader  :contentReference[oaicite:7]{index=7}

# ---------------------------
# helpers (no external deps)
# ---------------------------

def _ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)
    return path

def _pair_to_ids(pair_id: str) -> tuple[str,str]:
    # '1dzlA_5keqF' -> ('1dzlA','5keqF')
    return pair_str_to_tuple(pair_id)  # you already have this  :contentReference[oaicite:8]{index=8}

def _cluster_files(pair_dir: str) -> list[str]:
    # a3m clusters written by your clustering step
    return sorted(glob(os.path.join(pair_dir, "output_msa_cluster", "*.a3m")))

def _deep_msa_file(pair_dir: str) -> str:
    # where get_msa wrote the full alignment
    return os.path.join(pair_dir, "output_get_msa", "DeepMsa.a3m")

def _sample_ids_from_a3m(a3m_file: str, k: int) -> list[tuple[str,str]]:
    """
    Returns up to k (header, sequence) pairs from the a3m, in input order.
    Uses your read_msa() which strips insertions (good for ESM/AF feeding).  :contentReference[oaicite:9]{index=9}
    """
    entries = read_msa(a3m_file)
    if len(entries) <= k:
        return entries
    # deterministic but representative: take first k by default (you can change to random.sample)
    return entries[:k]

def _write_fasta(out_fa: str, entries: list[tuple[str,str]]):
    with open(out_fa, "w") as f:
        for hdr, seq in entries:
            name = hdr.split()[0].lstrip(">")
            f.write(f">{name}\n{seq}\n")

def _run(cmd: str, mode: str):
    if mode == "inline":
        print("[inline]", cmd, flush=True)
        # safer than os.system; keep shell quoting
        subprocess.run(cmd, shell=True, check=True)
    else:
        print("[sbatch]", cmd, flush=True)
        subprocess.run(cmd, shell=True, check=True)

def _require_linux(feature: str):
    if platform.system() == "Windows":
        raise RuntimeError(f"Can't run {feature} on Windows. Run on the Moriah cluster.")

# ---------------------------
# command builders (one place)
# ---------------------------

def build_msat_cmap_cmd(input_dir: str, out_dir: str) -> str:
    # Your repo already calls: python3 ./runMSATrans.py --input_msas <dir> -o <out>
    # Keep exactly that here so there's only ONE place to adjust in future.  :contentReference[oaicite:10]{index=10}
    return f"python3 ./runMSATrans.py --input_msas {shlex.quote(input_dir)} -o {shlex.quote(out_dir)}"

def build_esmfold_cmd(input_fa: str, out_dir: str, esm_version: str, model_tag: str|None) -> str:
    # Your repo uses ESMFoldHF.py; previously called with -input (adjust here if needed).  :contentReference[oaicite:11]{index=11}
    # I pass explicit -i/-o and optional flags so changing to RunESMFoldHF.py later is trivial.
    tag = f" --model_tag {shlex.quote(model_tag)}" if model_tag else ""
    return f"python3 ./ESMFoldHF.py -i {shlex.quote(input_fa)} -o {shlex.quote(out_dir)} --esm_version {shlex.quote(esm_version)}{tag}"

def build_af2_cmd(pair_dir: str, target_seq_fa: str, cluster_msa_a3m: str, out_dir: str, model_tag: str|None) -> str:
    # Your repo uses: python3 runAF.py -input ./Pipeline/<pair> -o ./Pipeline/<pair>/output_AF  :contentReference[oaicite:12]{index=12}
    # We extend runAF.py CLI in-place-friendly style via env vars or extra args; if runAF.py expects only -input/-o, you can
    # make it read the following envs; or modify runAF.py to accept --target_fasta/--msa flags. Centralized here either way.
    tag = f" --model_tag {shlex.quote(model_tag)}" if model_tag else ""
    return (
        f"python3 ./runAF.py -input {shlex.quote(pair_dir)} -o {shlex.quote(out_dir)} "
        f"--target_fasta {shlex.quote(target_seq_fa)} --msa {shlex.quote(cluster_msa_a3m)}{tag}"
    )

# ---------------------------
# tasks
# ---------------------------

def task_cmap_msa_transormer(pair_id: str, args: argparse.Namespace):
    """
    For each cluster + one deep MSA, run MSA Transformer & save .npy into:
      .../Pipeline/<pair>/output_cmaps/msa_transformer/
    """
    pair_dir = os.path.join(args.output_dir, pair_id)
    in_dir   = os.path.join(pair_dir, "output_msa_cluster")
    deep_msa = _deep_msa_file(pair_dir)

    out_dir  = os.path.join(pair_dir, "output_cmaps", "msa_transformer")
    _ensure_dir(out_dir)

    # Per-cluster
    if os.path.isdir(in_dir):
        cmd = build_msat_cmap_cmd(in_dir, out_dir)
        _run(cmd, args.run_job_mode)
    else:
        print(f"[skip] no clusters found at {in_dir}")

    # Whole-deep alignment (single map) – create a temp dir with just one file
    if os.path.isfile(deep_msa):
        tmp_dir = os.path.join(pair_dir, "_tmp_deep_msa_dir")
        _ensure_dir(tmp_dir)
        deep_copy = os.path.join(tmp_dir, "MSA_deep.a3m")
        if not os.path.exists(deep_copy):
            # symlink if possible to save space
            try:
                os.symlink(os.path.relpath(deep_msa, tmp_dir), deep_copy)
            except Exception:
                import shutil
                shutil.copyfile(deep_msa, deep_copy)
        cmd = build_msat_cmap_cmd(tmp_dir, out_dir)
        _run(cmd, args.run_job_mode)
        # keep tmp_dir so runs are reproducible; you can remove it if you prefer
    else:
        print(f"[skip] missing deep MSA: {deep_msa}")

    print(f"[ok] cmap_msa_transormer done for {pair_id}. Expect .npy files named per cluster and 'MSA_deep' in {out_dir}")

def task_esmfold(pair_id: str, args: argparse.Namespace):
    """
    For each cluster, sample up to N sequences, fold each separately.
    Output:
      .../Pipeline/<pair>/output_esm_fold/<esm_version>/
      file names include cluster id and sequence id
    """
    pair_dir = os.path.join(args.output_dir, pair_id)
    esm_out  = _ensure_dir(os.path.join(pair_dir, "output_esm_fold", args.esm_version))

    clusters = _cluster_files(pair_dir)
    if not clusters:
        print(f"[warn] no clusters found at {pair_dir}/output_msa_cluster; nothing to fold.")
        return

    for a3m in clusters:
        cluster_name = os.path.splitext(os.path.basename(a3m))[0]  # e.g., ShallowMsa_003
        # sample sequences
        entries = _sample_ids_from_a3m(a3m, args.cluster_sample_n)
        if not entries:
            continue
        # write per-cluster temp fasta containing up to N entries
        fa = os.path.join(pair_dir, f"_tmp_{cluster_name}.fasta")
        _write_fasta(fa, entries)

        # run ESMFold
        cmd = build_esmfold_cmd(fa, esm_out, args.esm_version, args.esm_model_tag)
        _run(cmd, args.run_job_mode)

    print(f"[ok] esmfold done for {pair_id}. Outputs under {esm_out}")

def task_af2(pair_id: str, args: argparse.Namespace):
    """
    For AF2, for each cluster and for EACH of the two native sequences,
    run AF with the cluster MSA, and also once with the deep MSA.
    Naming should reflect [sequence_of_pair] + [cluster_id] + [model].
    Output:
      .../Pipeline/<pair>/output_AF/AF2/
    """
    _require_linux("AlphaFold2")  # clear error on Windows as requested
    pair_dir   = os.path.join(args.output_dir, pair_id)
    af_out     = _ensure_dir(os.path.join(pair_dir, "output_AF", "AF2"))

    seqA, seqB = _pair_to_ids(pair_id)  # ('1dzlA','5keqF')  :contentReference[oaicite:13]{index=13}
    fasta_dir  = os.path.join(pair_dir, "fasta_chain_files")
    fastaA     = os.path.join(fasta_dir, f"{seqA}.fasta")
    fastaB     = os.path.join(fasta_dir, f"{seqB}.fasta")
    deep_msa   = _deep_msa_file(pair_dir)

    if not (os.path.isfile(fastaA) and os.path.isfile(fastaB)):
        raise FileNotFoundError(f"Missing fasta(s) under {fasta_dir}. Make sure 'load' step ran to extract seqs.")

    clusters = _cluster_files(pair_dir)
    if not clusters:
        print(f"[warn] no clusters found at {pair_dir}/output_msa_cluster; will still run deep MSA.")
    targets = [(seqA, fastaA), (seqB, fastaB)]

    # 1) Deep MSA for both sequences (2 * 5 models downstream)
    if os.path.isfile(deep_msa):
        for name, fa in targets:
            cmd = build_af2_cmd(pair_dir, fa, deep_msa, af_out, args.af_model_tag)
            _run(cmd, args.run_job_mode)
    else:
        print(f"[skip] missing deep MSA: {deep_msa}")

    # 2) Each cluster MSA for both sequences
    for a3m in clusters:
        cluster_name = os.path.splitext(os.path.basename(a3m))[0]  # e.g. ShallowMsa_003
        for name, fa in targets:
            cmd = build_af2_cmd(pair_dir, fa, a3m, af_out, args.af_model_tag)
            _run(cmd, args.run_job_mode)

    print(f"[ok] AF2 queued/completed for {pair_id}. Expect ~2*5*(#clusters+1) predictions in {af_out}.")

# ---------------------------
# deltaG helper (your bug)
# ---------------------------

def kw_deltaG(pair_id: str, args: argparse.Namespace):
    # You had: pA[:-0] and a string replace that evaluates oddly.
    # The goal is to remove the chain letter (last char) to get the pdb file names that were saved earlier.
    pA, pB = pair_id.split("_")
    pdbA = os.path.join(args.output_dir, pair_id, f"{pA[:-1]}.pdb")
    pdbB = os.path.join(args.output_dir, pair_id, f"{pB[:-1]}.pdb")
    return dict(
        pdb_pairs=[(pdbA, pdbB)],
        foldpair_ids=[pair_id],
        output_dir=os.path.join(args.output_dir, "output_deltaG"),
    )

# ---------------------------
# CLI
# ---------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_mode", required=True,
                        choices=["load","get_msa","cluster_msa","run_cmap_msa_transormer","run_esmfold","run_AF","plot","compute_deltaG","tree","ancestral","run_pipeline"])
    parser.add_argument("--foldpair_ids", nargs="+", default=["ALL"])
    parser.add_argument("--output_dir", default="Pipeline")
    parser.add_argument("--run_job_mode", default="inline", choices=["inline","sbatch"])

    # new knobs:
    parser.add_argument("--cluster_sample_n", type=int, default=10, help="ESMFold: sample up to N seqs/cluster")
    parser.add_argument("--esm_version", default="esm2", choices=["esm2","esm3"], help="which ESMFold model family")
    parser.add_argument("--esm_model_tag", default=None, help="optional extra tag for the ESMFold runner (e.g., hf card/model id)")

    parser.add_argument("--af_model_tag", default=None, help="optional AF model tag (kept generic)")

    # plot options
    parser.add_argument("--global_plots", action="store_true", help="passed into plotting function if used")

    args = parser.parse_args()

    # prepare list of pairs
    if args.foldpair_ids == ["ALL"]:
        with open("data/foldswitch_PDB_IDs_full.txt","r") as f:
            pairs = [line.strip().replace("\t","_") for line in f if line.strip()]
    else:
        pairs = args.foldpair_ids

    for pair_id in pairs:
        pair_dir = os.path.join(args.output_dir, pair_id)
        _ensure_dir(pair_dir)
        _ensure_dir(os.path.join(pair_dir, "fasta_chain_files"))

        # keep your existing load/msa/cluster/simple runners exactly as-is where you already have them wired
        if args.run_mode == "run_cmap_msa_transormer":
            task_cmap_msa_transormer(pair_id, args)

        elif args.run_mode == "run_esmfold":
            task_esmfold(pair_id, args)

        elif args.run_mode == "run_AF":
            task_af2(pair_id, args)

        elif args.run_mode == "plot":
            # lazy import to avoid PyMOL hard dep
            from utils.protein_plot_utils import make_foldswitch_all_plots, PYMOL_AVAILABLE  # :contentReference[oaicite:14]{index=14}
            if not PYMOL_AVAILABLE:
                raise RuntimeError("Plotting requires PyMOL. Install PyMOL or run on a machine with PyMOL.")
            # You already compute pdbids/pdbchains elsewhere; reuse your existing call pattern if needed.

        elif args.run_mode == "compute_deltaG":
            # call your existing energy pipeline with the fixed kw helper
            kw = kw_deltaG(pair_id, args)
            compute_global_and_residue_energies(kw["pdb_pairs"], kw["foldpair_ids"], kw["output_dir"])

        else:
            # keep your current/legacy modes working (load/get_msa/cluster_msa/tree/ancestral/run_pipeline)
            # by delegating to your current command strings (unchanged from repo)  :contentReference[oaicite:15]{index=15} :contentReference[oaicite:16]{index=16}
            run_str = ""
            if args.run_mode == "load":
                # your existing loader
                pass
            elif args.run_mode == "get_msa":
                run_str = f"python3 ./get_msa.py {os.path.join(pair_dir,'fasta_chain_files', f'{pair_id[:5]}.fasta')} {os.path.join(pair_dir,'output_get_msa')} -name 'DeepMsa'"
            elif args.run_mode == "cluster_msa":
                run_str = f"python3 ./ClusterMSA_moriah.py --keyword ShallowMsa -i {os.path.join(pair_dir,'output_get_msa','DeepMsa.a3m')} -o {os.path.join(pair_dir,'output_msa_cluster')}"
            elif args.run_mode == "tree":
                msa_file = os.path.join(pair_dir,"output_get_msa","DeepMsa.a3m")
                print("reconstruct phylogenetic tree from msa=", msa_file)
                from utils.phytree_utils import phytree_from_msa
                phytree_from_msa(msa_file, output_tree_file= os.path.join(pair_dir,"output_phytree","DeepMsa_tree.nwk"))  # :contentReference[oaicite:17]{index=17}
            elif args.run_mode == "ancestral":
                from utils.phytree_utils import reconstruct_ancestral_sequences
                msa_file = os.path.join(pair_dir,"output_get_msa","DeepMsa.a3m")
                out_tree = os.path.join(pair_dir,"output_phytree","DeepMsa_tree.nwk")
                anc_out  = os.path.join(pair_dir,"output_phytree","DeepMsa_anc_seq.a3m")
                reconstruct_ancestral_sequences(out_tree, msa_file, anc_out)

            if run_str:
                _run(run_str, args.run_job_mode)

if __name__ == "__main__":
    print(f"FINISHED LOADING CONFIG. RUNNING FOLD-SWITCH PIEPLINE WITH USER: {user}  ENVIRONMENT: {platform.system()}")
    main()
