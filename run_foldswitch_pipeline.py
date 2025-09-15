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
from utils.protein_utils import read_msa, greedy_select, extract_protein_sequence, load_seq_and_struct, process_sequence
from utils.msa_utils import write_fasta, load_fasta, build_pair_seed_a3m_from_pair  # your existing writer
from utils.phytree_utils import phytree_from_msa

from Analysis.postprocess_unified import post_processing_analysis
from TableResults.gen_html_table import gen_html_from_summary_table, gen_html_from_cluster_detailed_table,


RUN_MODE_DESCRIPTIONS = {
    "get_msa":          "Download/build deep MSAs for the pair and write DeepMsa.a3m.",
    "cluster_msa":      "Cluster the pairwise MSA into shallow clusters (ShallowMsa_XXX.a3m).",
    "run_AF":           "Run AlphaFold (AF2/AF3/both) per chain × cluster. Use --af_ver {2,3,both}.",
    "run_esmfold":      "Run ESMFold on the pair. Use --esm_model {esm2,esm3,both}.",
    "run_cmap_msa_transformer":      "Run MSA-transformer on the pair to get contact maps.",
    "compute_deltaG":   "Compute ΔG stability metrics (requires PyRosetta).",
    "postprocess": "Compute TM/cmap metrics and build summary/detailed tables. Use --force_rerun_postprocess TRUE to recompute." ,
    "plot":             "Generate pair-specific plots (requires PyMOL).",
    "clean":            "Remove previous outputs for the pair.",
    "msaclust_pipeline":"Full pipeline: get_msa → cluster_msa → AF/ESM (as configured).",
    "help":             "Print this list of run modes with one-line explanations.",
}


# ------------------------- helpers -------------------------

import re, shutil

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
            import pandas as pd
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


def task_postprocess(pair_id: str | None, args) -> None:
    pairs = None
    if pair_id:
        pairs = [pair_id]
    force = _bool_from_tf(getattr(args, "force_rerun_postprocess", "FALSE"))
    print(f"[postprocess] pairs={pairs or 'ALL'} | force_rerun={force}")
    post_processing_analysis(force_rerun=force, pairs=pairs)

# All Pipeline
def task_msaclust_pipeline(pair_id: str, args: argparse.Namespace) -> None:
    """
    If --force_rerun TRUE: run all steps unconditionally.
    Else: run step-by-step and only execute steps whose outputs are missing.
    AF step is always invoked (it internally skips already-complete jobs).
    """
    force_all = _bool_from_tf(getattr(args, "force_rerun", "FALSE"))

    # 1) load PDBs/FASTA (cheap: always safe)
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

    # 4) AF (AF2/AF3/both): always call; the task itself skips completed out_dirs unless forced
    if force_all:
        # bubble a TRUE into AF-only force without changing your AF default elsewhere
        setattr(args, "force_rerun_AF", "TRUE")
    print("[pipeline] AF → running (task will skip per-outdir if already complete)")
    # If user asked specifically for AF2 only inside pipeline, keep it; otherwise default to both inside pipeline
    if getattr(args, "af_ver", "2") == "2" and not force_all:
        pass  # respect user's choice
    else:
        # default to both in the full pipeline unless user overrode
        args.af_ver = "both"
    task_af(pair_id, args)

    # 5) cmap (MSA-Transformer)
    if force_all or not _has_cmaps(pair_id):
        print("[pipeline] cmap_msa_transformer → running")
        task_cmap_msa_transformer(pair_id, "inline")
    else:
        print("[pipeline] cmap_msa_transformer → skip (cmaps exist)")

    # 6) ESMFold (both models if esm_model==both; otherwise the selected one)
    wanted_models = ["esm2", "esm3"] if getattr(args, "esm_model", None) in (None, "both") else [args.esm_model]
    for model in wanted_models:
        a2 = deepcopy(args)
        a2.esm_model = model
        need = force_all or not _has_esm_model(pair_id, model)
        if need:
            print(f"[pipeline] esmfold({model}) → running")
            task_esmfold(pair_id, a2)
        else:
            print(f"[pipeline] esmfold({model}) → skip (outputs exist)")

    # 7) plots (optional to gate; keep as-is or skip if you detect finished artifacts)
    try:
        ap = deepcopy(args)
        ap.global_plots = False
        if hasattr(ap, "plot_trees"):
            ap.plot_trees = False
        print("[pipeline] plot → running")
        task_plot(pair_id, ap)
    except Exception as e:
        print(f"[plot] skipped: {e}")


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
                            "postprocess", "msaclust_pipeline", "help"])  # Last one is the full pipeline for a pair
    p.add_argument("--foldpair_ids", nargs="+", default=["ALL"],
                   help="e.g. 1dzlA_5keqF (default: ALL in data list)")
    p.add_argument("--run_job_mode", default="inline", choices=["inline", "sbatch"])

    # AlphaFold options
    p.add_argument("--allow_inline_af", action="store_true",
                   help="Allow AF2 to run inline even if not in a Slurm session (expert only).")
    p.add_argument("--af_ver", default="both", choices=["2", "3", "both"],
                   help="Which AlphaFold to run for --run_mode run_AF")  # default do both AF2 and AF3

    p.add_argument(
        "--force_rerun_AF",
        default="FALSE",
        choices=["TRUE", "FALSE"],
        help="Run new AF2/AF3 predictions even if outputs exist. Default FALSE (skip if found).")

    # ESMFold options
    p.add_argument("--cluster_sample_n", type=int, default=10)
    p.add_argument("--esm_model", default=None, choices=["esm2", "esm3", "both"])
    p.add_argument("--esm_device", default="auto", choices=["auto", "cpu", "cuda", "mps"])

    # Post-processing options, computing metrics
    p.add_argument("--postprocess", default="FALSE", help="If TRUE, run post-processing after the selected task(s).")
    p.add_argument("--force_rerun_postprocess", default="FALSE", help="If TRUE, recompute per-pair caches.")

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
            task_load(pair_id, args)

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

        elif args.run_mode == "compute_deltaG":
            task_deltaG(pair_id)

        elif args.run_mode == "postprocess":
            task_postprocess(pair_id, args)             # Allow optional single-pair positional arg the same way other modes do

        elif args.run_mode == "plot":
            task_plot(pair_id, args)


        elif args.run_mode == "clean":  # Remove existing files to run the pipeline clean
            task_clean(pair_id, args)

        elif args.run_mode == "msaclust_pipeline":
            task_msaclust_pipeline(pair_id, args)

        else:
            raise ValueError(args.run_mode)


    print("[done run+processing]", flush=True)

    # ---- After all pairs processed, optionally build reports ----
    if args.reports in ("tables", "all"):
        # build global CSVs (summary + detailed)
        # This will read existing per-pair CSVs (and only recompute if needed)
        post_processing_analysis(force_rerun=False, pairs=None)

        # turn CSVs into HTML tables
        gen_html_from_summary_table()
        gen_html_from_cluster_detailed_table()

    if args.reports in ("html", "all"):
        # per-pair HTML pages from the notebook
        import subprocess, shlex
        pairs_arg = " ".join(args.html_pairs)
        cmd = f"{shlex.quote(sys.executable)} Analysis/NotebookGen/generate_notebooks.py {pairs_arg} --kernel python3"
        subprocess.run(cmd, shell=True, check=True)

    print("[done html outputs]", flush=True)

if __name__ == "__main__":
    main()
