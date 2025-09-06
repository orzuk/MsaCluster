#!/usr/bin/env python3
# Rub pipeline from a list of PDBs

import os
import sys
import platform
import argparse
from typing import Callable, Optional, Dict, Any

# --- use your modules (do not reinvent) ---
from config import *  # paths + env
from utils.utils import *
from utils.protein_utils import *
from utils.protein_plot_utils import *
from utils.phytree_utils import *
from utils.energy_utils import *  # ΔG

# Notes:
# - Inline/sbatch command templates are taken from your current script:
#   get_msa/cluster_msa/MSA-Transformer/ESMFold/AF/tree (sbatch + inline)  :contentReference[oaicite:0]{index=0} :contentReference[oaicite:1]{index=1} :contentReference[oaicite:2]{index=2}
# - Path conventions for PDB filenames per pair follow your code (store as <pair>/<4-letter>.pdb)  :contentReference[oaicite:3]{index=3}
# - `pair_str_to_tuple("1dzlA_5keqF") -> ("1dzlA","5keqF")` is in your utils  :contentReference[oaicite:4]{index=4}

PIPELINE_DIR = "Pipeline"  # keep your relative path used throughout scripts
PAIRS_FILE_DEFAULT = "data/foldswitch_PDB_IDs_full.txt"


def read_pairs_file(pairs_file: str):
    """
    Parse your TSV-like list to:
      - foldpair_ids: ["1dzlA_5keqF", ...]
      - pdbids:       [["1dzl","5keq"], ...]
      - pdbchains:    [["A","F"], ...]
    This mirrors your original logic.  :contentReference[oaicite:5]{index=5}
    """
    with open(pairs_file, "r") as fh:
        rows = [line.rstrip() for line in fh]
    foldpair_ids = [s.replace("\t", "_") for s in rows]
    pdbids = [s.split("\t") for s in rows]
    pdbchains = [[s[0][-1], s[1][-1]] for s in pdbids]
    pdbids = [[s[0][:-1], s[1][:-1]] for s in pdbids]
    return foldpair_ids, pdbids, pdbchains


def ensure_pair_dirs(pair_id: str):
    """Create per-pair dirs once (same logic you had).  :contentReference[oaicite:6]{index=6}"""
    pair_dir = f"{PIPELINE_DIR}/{pair_id}"
    for d in [pair_dir, f"{pair_dir}/chain_pdb_files", f"{pair_dir}/fasta_chain_files"]:
        if not os.path.exists(d):
            os.makedirs(d, exist_ok=True)
    return pair_dir


def pair_ctx(pair_id: str) -> Dict[str, Any]:
    """
    Context dict for templates and inline runners:
      pdb1_id/pdb2_id: 4-letter ids
      chain1/chain2:   chain letters
      pdb1_path/pdb2_path: Pipeline/<pair>/<id>.pdb
      fasta_path: fasta for the *first* fold (kept as in your script)  :contentReference[oaicite:7]{index=7}
    """
    pA, pB = pair_str_to_tuple(pair_id)  # e.g. "1dzlA","5keqF"
    pdb1_id, chain1 = pA[:-1], pA[-1]
    pdb2_id, chain2 = pB[:-1], pB[-1]
    pair_dir = ensure_pair_dirs(pair_id)
    pdb1_path = f"{pair_dir}/{pdb1_id}.pdb"
    pdb2_path = f"{pair_dir}/{pdb2_id}.pdb"
    fasta_path = f"{pair_dir}/fasta_chain_files/{pdb1_id+chain1}.fasta"
    return dict(
        pair_id=pair_id, pair_dir=pair_dir,
        pdb1_id=pdb1_id, pdb2_id=pdb2_id, chain1=chain1, chain2=chain2,
        pdb1_path=pdb1_path, pdb2_path=pdb2_path, fasta_path=fasta_path
    )


def ensure_inputs_for_fasta(ctx: Dict[str, Any]):
    """
    Your pipeline calls get_fasta_chain_seq on the first PDB file even when run_mode != 'load'.
    If the PDB doesn't exist yet, fetch via load_seq_and_struct first to avoid a crash.  :contentReference[oaicite:8]{index=8}
    """
    if not os.path.isfile(ctx["pdb1_path"]) or not os.path.isfile(ctx["pdb2_path"]):
        # fetch/save PDBs using your existing loader
        load_seq_and_struct(ctx["pair_dir"], [ctx["pdb1_id"], ctx["pdb2_id"]], [ctx["chain1"], ctx["chain2"]])
    if not os.path.isfile(ctx["fasta_path"]):
        get_fasta_chain_seq(ctx["pdb1_path"], ctx["pdb1_id"] + ctx["chain1"], ctx["pair_dir"])


def run_shell(cmd: str):
    print(">>", cmd)
    rc = os.system(cmd)
    if rc != 0:
        raise RuntimeError(f"Command failed (exit={rc}): {cmd}")


def _err(msg: str):
    raise RuntimeError(msg)


# ---- Option registry: define each option ONCE here ----
# Each entry provides:
#   - 'inline': either a callable(ctx,args) or a format-string shell command
#   - 'sbatch': same as above (or None if not supported)
#   - 'requires': list of capabilities to validate upfront
#   - 'help': short description
TASKS: Dict[str, Dict[str, Any]] = {
    "load": {
        "help": "Fetch PDBs/Seqs for the pair and write local files.",
        "requires": [],
        "inline": lambda ctx, args: load_seq_and_struct(
            ctx["pair_dir"], [ctx["pdb1_id"], ctx["pdb2_id"]], [ctx["chain1"], ctx["chain2"]]
        ),
        "sbatch": None,  # Not a cluster job
    },
    "get_msa": {
        "help": "Create MSAs (DeepMSA) for the first chain fasta.",
        "requires": [],
        # Inline command template from your current script  :contentReference[oaicite:9]{index=9}
        "inline": "python3 ./get_msa.py {fasta} ./Pipeline/{pair}/output_get_msa -name 'DeepMsa'",
        # sbatch template from your current script  :contentReference[oaicite:10]{index=10}
        "sbatch": "sbatch -o './Pipeline/{pair}/get_msa_for_{pair}.out' ./Pipeline/get_msa_params.sh {fasta} {pair}",
        "pre": ensure_inputs_for_fasta,
    },
    "cluster_msa": {
        "help": "Cluster MSAs.",
        "requires": [],
        # Inline command  :contentReference[oaicite:11]{index=11}
        "inline": "python3 ./ClusterMSA_moriah.py --keyword ShallowMsa -i ./Pipeline/{pair}/output_get_msa/DeepMsa.a3m -o ./Pipeline/{pair}/output_msa_cluster",
        # sbatch command  :contentReference[oaicite:12]{index=12}
        "sbatch": "sbatch -o './Pipeline/{pair}/cluster_msa_for_{pair}.out' ./Pipeline/ClusterMSA_params.sh {pair}",
    },
    "run_cmap_esm": {
        "help": "Run MSA Transformer cmap from clustered MSAs.",
        "requires": [],
        # Inline  :contentReference[oaicite:13]{index=13}
        "inline": "python3 ./runMSATrans.py --input_msas ./Pipeline/{pair}/output_msa_cluster -o ./Pipeline/{pair}/output_cmap_esm",
        # sbatch  :contentReference[oaicite:14]{index=14}
        "sbatch": "sbatch ./Pipeline/CmapESM_params.sh Pipeline/{pair}",
    },
    "run_esmfold": {
        "help": "Run ESMFold.",
        "requires": [],
        # Inline (your script supports `-input <pair_id>`)  :contentReference[oaicite:15]{index=15}
        "inline": "python3 ./ESMFoldHF.py -input {pair}",
        # sbatch (same as your current hardcoded path)  :contentReference[oaicite:16]{index=16}
        "sbatch": "sbatch /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/RunEsmFold_params.sh {pair}",
    },
    "run_AF": {
        "help": "Run AlphaFold2 models.",
        "requires": [],
        # Inline  :contentReference[oaicite:17]{index=17}
        "inline": "python3 runAF.py -input ./Pipeline/{pair} -o ./Pipeline/{pair}/output_AF",
        # sbatch  :contentReference[oaicite:18]{index=18}
        "sbatch": "sbatch -o './Pipeline/{pair}/run_AF_for_{pair}.out' ./Pipeline/RunAF_params.sh {pair}",
    },
    "tree": {
        "help": "Reconstruct phylogenetic tree from DeepMSA.",
        "requires": [],
        # Inline: call your function directly (you already do this)  :contentReference[oaicite:19]{index=19}
        "inline": lambda ctx, args: phytree_from_msa(
            f"Pipeline/{ctx['pair_id']}/output_get_msa/DeepMsa.a3m",
            output_tree_file=f"Pipeline/{ctx['pair_id']}/output_phytree/DeepMsa_tree.nwk",
        ),
        # sbatch template  :contentReference[oaicite:20]{index=20}
        "sbatch": "sbatch -o './Pipeline/{pair}/tree_reconstruct_for_{pair}.out' ./Pipeline/tree_reconstruct_params.sh {pair}",
    },
    "compute_deltaG": {
        "help": "Compute PyRosetta ΔG for both PDBs and compare residue energies.",
        "requires": ["linux", "pyrosetta"],
        # Inline: use your energy_utils API  :contentReference[oaicite:21]{index=21}
        "inline": lambda ctx, args: compute_global_and_residue_energies(
            [(ctx["pdb1_path"], ctx["pdb2_path"])],
            [ctx["pair_id"]],
            f"{PIPELINE_DIR}/output_deltaG",
        ),
        "sbatch": None,  # keep it inline-only for now
        "pre": lambda ctx: (os.path.isfile(ctx["pdb1_path"]) and os.path.isfile(ctx["pdb2_path"])) or
                           load_seq_and_struct(ctx["pair_dir"], [ctx["pdb1_id"], ctx["pdb2_id"]], [ctx["chain1"], ctx["chain2"]])
    },
    "plot": {
        "help": "Make all plots for the pair.",
        "requires": [],
        # Inline: call your plotting entrypoint  :contentReference[oaicite:22]{index=22}
        "inline": lambda ctx, args: make_foldswitch_all_plots(
            [ctx["pdb1_id"], ctx["pdb2_id"]], PIPELINE_DIR, ctx["pair_id"], [ctx["chain1"], ctx["chain2"]],
            plot_tree_clusters=bool(getattr(args, "plot_trees", False)),
            plot_contacts=True,
            global_plots=bool(getattr(args, "global_plots", False))
        ),
        "sbatch": None,
    },
    # Examples to add later:
    # "alphafold3": { ... },
    # "asr": { ... },
}


def validate_requires(task_name: str, spec: Dict[str, Any], run_job_mode: str):
    sys_os = platform.system()
    if "linux" in spec.get("requires", []):
        if sys_os != "Linux":
            _err(f"{task_name}: Linux-only. You're on {sys_os}.")
    if "pyrosetta" in spec.get("requires", []):
        if not PYROSETTA_AVAILABLE:
            _err("Can't run PyRosetta in your current environment. Re-run on Linux cluster.")  # :contentReference[oaicite:23]{index=23}
    if run_job_mode == "job" and spec.get("sbatch") is None:
        _err(f"{task_name}: sbatch mode not supported. Use --run_job_mode inline.")


def run_one_pair(run_mode: str, pair_id: str, args: argparse.Namespace):
    spec = TASKS.get(run_mode)
    if not spec:
        raise KeyError(f"Unknown run_mode: {run_mode}. Try --run_mode help")

    ctx = pair_ctx(pair_id)

    # Any per-option precondition (e.g., generate fasta or ensure PDBs exist)
    pre = spec.get("pre")
    if pre:
        # support both pre(ctx) and pre(ctx,args)
        try:
            pre(ctx, args)  # type: ignore
        except TypeError:
            pre(ctx)  # type: ignore

    validate_requires(run_mode, spec, args.run_job_mode)

    # Execute as inline python or shell / sbatch command
    if args.run_job_mode == "job":
        cmd = spec["sbatch"]
        if callable(cmd):
            return cmd(ctx, args)
        else:
            # fill minimal placeholders
            filled = cmd.format(
                pair=ctx["pair_id"], fasta=ctx["fasta_path"], pair_dir=ctx["pair_dir"]
            )
            return run_shell(filled)
    else:
        inline = spec["inline"]
        if callable(inline):
            return inline(ctx, args)
        else:
            filled = inline.format(
                pair=ctx["pair_id"], fasta=ctx["fasta_path"], pair_dir=ctx["pair_dir"]
            )
            return run_shell(filled)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Fold-switch pipeline runner")
    p.add_argument("--run_mode", required=True,
                   help=f"One of: {', '.join(sorted(TASKS.keys()))} (or 'help')")
    p.add_argument("--foldpair_ids", default="ALL",
                   help="Comma-separated list of pairs (e.g., 1dzlA_5keqF), or 'ALL'")
    p.add_argument("--run_job_mode", default="inline", choices=["inline", "job"],
                   help="inline (local python/shell) or job (sbatch)")
    p.add_argument("--pdbids_file", default=PAIRS_FILE_DEFAULT)
    # plotting flags (used by 'plot')
    p.add_argument("--global_plots", action="store_true", help="Enable global plots")
    p.add_argument("--plot_trees", action="store_true", help="Draw tree clusters")
    return p.parse_args()


def main():
    args = parse_args()
    if args.run_mode == "help":
        print("Available run modes:")
        for k in sorted(TASKS.keys()):
            print(f"  - {k:16s} {TASKS[k]['help']}")
        return

    # Gather pairs list
    foldpair_ids, pdbids, pdbchains = read_pairs_file(args.pdbids_file)

    # Determine target pairs
    if args.foldpair_ids == "ALL":
        targets = foldpair_ids
    else:
        targets = [s.strip() for s in args.foldpair_ids.split(",") if s.strip()]

    # Basic Windows guard for sbatch mode
    if platform.system() == "Windows" and args.run_job_mode == "job":
        _err("sbatch jobs are not supported on Windows. Use --run_job_mode inline.")

    for pair in targets:
        print(f"\n=== [{args.run_mode}] on pair: {pair} ===")
        run_one_pair(args.run_mode, pair, args)


if __name__ == "__main__":
    main()
