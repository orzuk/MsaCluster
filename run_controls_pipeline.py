#!/usr/bin/env python3
"""
Controls pipeline runner -- single-protein (non-pair) negative controls.

Runs the same MSA-clustering + structure-prediction workflow on single proteins
(e.g., Porter 2022 IDP controls) to measure cluster structural diversity.
If a protein is *not* a fold-switcher, all cluster predictions should look alike
(high pairwise TM-scores), giving a baseline for fold-switching detection.

Usage examples
--------------
    # Run MSA for one control protein
    python run_controls_pipeline.py --run_mode get_msa --protein_ids Q9NZC2 \\
           --control_type PorterIDP --run_job_mode inline

    # Full pipeline for all Porter IDP controls (sbatch on cluster)
    python run_controls_pipeline.py --run_mode full_pipeline --protein_ids ALL \\
           --control_type PorterIDP --run_job_mode sbatch

    # Postprocess: compute cluster diversity for all controls
    python run_controls_pipeline.py --run_mode postprocess --protein_ids ALL

Stages
------
    load         Write FASTA (already done by setup_porter_idp.py; placeholder)
    get_msa      Build deep MSA via ColabFold for a single sequence
    tree         Build phylogenetic tree from deep MSA
    cluster_msa  Cluster the deep MSA into shallow sub-alignments
    run_AF       Run AlphaFold2/3 on each cluster MSA (monomer mode)
    run_esmfold  Run ESMFold on sampled sequences per cluster
    postprocess  Compute cluster structural diversity (pairwise TM-scores)
    full_pipeline  Run all stages sequentially
"""

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from copy import deepcopy
from glob import glob
from pathlib import Path
from typing import List, Optional

# ---------------------------------------------------------------------------
# Project imports
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(__file__))

from config import (
    CONTROLS_BASE_DIR, PORTER_IDP_DIR, PIPELINE_BASE_DIR,
    CONTACT_CUTOFF, ESM_MAX_LIGHT_SEQ_LEN,
)
from utils.utils import list_control_proteins, ensure_dir, str2bool
from utils.msa_utils import load_fasta, write_fasta

# Optional heavy imports (may be unavailable in light environments)
try:
    from utils.phytree_utils import phytree_from_msa
except Exception:
    phytree_from_msa = None

try:
    from utils.protein_utils import read_msa, greedy_select
except Exception:
    read_msa = greedy_select = None

try:
    from utils.align_utils import compute_tmscore_align
except Exception:
    compute_tmscore_align = None

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _run(cmd: str, mode: str = "inline") -> None:
    """Execute *cmd* inline (blocking) or via sbatch."""
    print(f"[run:{mode}] {cmd}", flush=True)
    if mode == "sbatch":
        subprocess.run(cmd, shell=True, check=True)
    else:
        subprocess.run(cmd, shell=True, check=True)


def _protein_dir(protein_id: str, control_type: str = "PorterIDP") -> str:
    return os.path.join(CONTROLS_BASE_DIR, control_type, protein_id)


def _find_fasta(protein_id: str, control_type: str = "PorterIDP") -> Optional[str]:
    """Locate the FASTA file for a control protein."""
    pdir = _protein_dir(protein_id, control_type)
    for name in [f"{protein_id}.fasta", "sequence.fasta"]:
        p = os.path.join(pdir, name)
        if os.path.isfile(p):
            return p
    # fallback: first .fasta in dir
    fastas = glob(os.path.join(pdir, "*.fasta"))
    return fastas[0] if fastas else None


def _read_sequence(fasta_path: str) -> str:
    """Read the first sequence from a FASTA file (gap-free)."""
    ids, seqs = load_fasta(fasta_path)
    if not seqs:
        raise FileNotFoundError(f"No sequences in {fasta_path}")
    return "".join(ch for ch in seqs[0] if ch.isalpha())


# ---------------------------------------------------------------------------
# Pipeline stages
# ---------------------------------------------------------------------------

def task_load(protein_id: str, control_type: str, args: argparse.Namespace) -> None:
    """Ensure FASTA exists.  (Already created by setup_porter_idp.py.)"""
    pdir = _protein_dir(protein_id, control_type)
    ensure_dir(pdir)
    fa = _find_fasta(protein_id, control_type)
    if fa and os.path.isfile(fa):
        seq = _read_sequence(fa)
        print(f"[load] {protein_id}: {len(seq)} residues  ({fa})")
    else:
        print(f"[load] WARN no FASTA for {protein_id} in {pdir}")


def task_get_msa(protein_id: str, control_type: str, args: argparse.Namespace) -> None:
    """
    Build a deep MSA for a single protein via ColabFold.

    Unlike the fold-pair pipeline (which uses run_getMSA.py with a 2-sequence
    seed A3M + seed-locking), controls call ColabFold directly with one
    sequence and write the resulting MSA as DeepMsa.a3m.

    NOTE: run_getMSA.py requires a 2-seq seed and cannot be reused here.
    """
    pdir = _protein_dir(protein_id, control_type)
    fa = _find_fasta(protein_id, control_type)
    if not fa:
        print(f"[get_msa] ERR no FASTA for {protein_id}")
        return
    seq = _read_sequence(fa)

    out_dir = os.path.join(pdir, "output_get_msa")
    ensure_dir(out_dir)

    # Write 1-sequence seed A3M (for reference/compatibility)
    seed_a3m = os.path.join(out_dir, "seed.a3m")
    with open(seed_a3m, "w") as fh:
        fh.write(f">{protein_id}\n{seq}\n")

    deep_msa = os.path.join(out_dir, "DeepMsa.a3m")
    if os.path.isfile(deep_msa) and os.path.getsize(deep_msa) > 100:
        print(f"[get_msa] DeepMsa.a3m already exists for {protein_id}, skipping")
        return

    # Write a small helper script that calls ColabFold for one sequence.
    # We use the same _call_get_msa_and_templates wrapper from run_getMSA.py.
    helper_script = os.path.join(pdir, "_get_msa_single.py")
    repo_root = os.path.abspath(os.path.dirname(__file__))
    with open(helper_script, "w") as fh:
        fh.write(f"""\
#!/usr/bin/env python3
\"\"\"Build a monomer deep MSA for a single protein via ColabFold.\"\"\"
import sys, os
sys.path.insert(0, {repo_root !r})
os.environ.setdefault("COLABFOLD_USE_LOCAL_MMSEQS", "0")
os.environ.setdefault("COLABFOLD_USER_AGENT", "MsaCluster/controls")

from pathlib import Path
from run_getMSA import _call_get_msa_and_templates
from colabfold.batch import msa_to_str

protein_id = {protein_id !r}
seq = {seq !r}
result_dir = Path({out_dir !r})
result_dir.mkdir(parents=True, exist_ok=True)

print(f"[get_msa] Building monomer MSA for {{protein_id}} ({{len(seq)}} aa)")

tup = _call_get_msa_and_templates(
    sequences=[seq], result_dir=result_dir, jobname="DeepMsa"
)
unpaired, paired, q_uniq, q_card, _ = tup
a3m_text = msa_to_str(unpaired, None, q_uniq, q_card)

# Write DeepMsa.a3m with the query as first row
out_path = result_dir / "DeepMsa.a3m"
with open(out_path, "w") as fh:
    # Ensure query row is first (ColabFold usually includes it)
    fh.write(f">{{protein_id}}\\n{{seq}}\\n")
    # Append the rest (skip query if ColabFold already included it)
    for line in a3m_text.splitlines():
        fh.write(line + "\\n")

# Count rows
n_rows = sum(1 for line in open(out_path) if line.startswith(">"))
print(f"[get_msa] DeepMsa.a3m written: {{out_path}} (rows={{n_rows}})")
""")

    python_for_msa = os.environ.get("GET_MSA_PYTHON", sys.executable)

    if args.run_job_mode == "sbatch":
        log = os.path.join(pdir, f"get_msa_{protein_id}.out")
        cmd = (
            f"sbatch -o {shlex.quote(log)} "
            f"--wrap {shlex.quote(f'{shlex.quote(python_for_msa)} {shlex.quote(helper_script)}')}"
        )
        _run(cmd, "sbatch")
    else:
        cmd = f"{shlex.quote(python_for_msa)} {shlex.quote(helper_script)}"
        _run(cmd, "inline")


def task_tree(protein_id: str, control_type: str, args: argparse.Namespace) -> None:
    """Build a phylogenetic tree from the deep MSA."""
    if phytree_from_msa is None:
        raise RuntimeError("phytree_from_msa unavailable (missing deps).")

    pdir = _protein_dir(protein_id, control_type)
    msa_file = os.path.join(pdir, "output_get_msa", "DeepMsa.a3m")
    out = os.path.join(pdir, "output_phytree", "DeepMsa_tree.nwk")
    ensure_dir(os.path.dirname(out))

    cap = getattr(args, "tree_max_seqs", 5000)
    if cap is not None and int(cap) <= 0:
        cap = None

    t0 = time.time()
    print(f"[tree] {protein_id}: building from {msa_file} cap={cap or 'ALL'}", flush=True)
    phytree_from_msa(msa_file, output_tree_file=out, max_seqs=cap,
                     seed=getattr(args, "tree_seed", 123))
    print(f"[tree] {protein_id}: done in {time.time() - t0:.1f}s", flush=True)


def task_cluster_msa(protein_id: str, control_type: str, args: argparse.Namespace) -> None:
    """Cluster the deep MSA into shallow sub-alignments."""
    pdir = _protein_dir(protein_id, control_type)
    py = shlex.quote(sys.executable)
    pid_q = shlex.quote(protein_id)

    alg = args.cluster_alg
    tree_arg = ""
    if alg == "tree":
        tree_path = "output_phytree/DeepMsa_tree.nwk"
        tree_arg = f"--tree_path {shlex.quote(tree_path)} "

    cluster_script = os.path.abspath("run_ClusterMSA.py")
    cmd = (
        f"bash -lc 'cd {shlex.quote(pdir)} && "
        f"{py} {shlex.quote(cluster_script)} "
        f"--keyword ShallowMsa "
        f"--a3m output_get_msa/DeepMsa.a3m "
        f"-o output_msa_cluster "
        f"--cluster_alg {shlex.quote(alg)} {tree_arg}"
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
    _run(cmd, args.run_job_mode)


def task_af(protein_id: str, control_type: str, args: argparse.Namespace) -> None:
    """
    Run AlphaFold on each cluster MSA in monomer mode.

    Unlike fold-pair mode (which builds pair-specific A3Ms per chain),
    controls simply feed each cluster A3M directly to AF -- the MSA already
    contains only one protein's sequences with a consistent column frame.
    """
    pdir = _protein_dir(protein_id, control_type)
    af_ver = str(getattr(args, "af_ver", "2")).lower()

    deep_a3m = os.path.join(pdir, "output_get_msa", "DeepMsa.a3m")
    if not os.path.isfile(deep_a3m):
        print(f"[af] ERR missing DeepMsa: {deep_a3m}")
        return

    cluster_dir = os.path.join(pdir, "output_msa_cluster")
    cluster_a3ms = sorted(glob(os.path.join(cluster_dir, "ShallowMsa_*.a3m")))

    # Build list of (a3m_path, label) jobs
    jobs = [(deep_a3m, "DeepMsa")]
    for a3m in cluster_a3ms:
        jobs.append((a3m, Path(a3m).stem))

    print(f"[af] {protein_id}: {len(jobs)} jobs (1 deep + {len(cluster_a3ms)} clusters)")
    if not jobs:
        return

    versions = ["2"] if af_ver == "2" else (["3"] if af_ver == "3" else ["2", "3"])
    force = str(getattr(args, "force_rerun_AF", "FALSE")).upper() == "TRUE"

    # Slurm options
    gres = getattr(args, "sbatch_gres", "gpu:1")
    cpus = int(getattr(args, "sbatch_cpus", 8))
    mem = getattr(args, "sbatch_mem", "40G")
    stime = getattr(args, "sbatch_time", "24:00:00")
    part = getattr(args, "sbatch_partition", None)
    cons = getattr(args, "sbatch_constraint", None)
    acct = getattr(args, "sbatch_account", None)
    sbatch_opts = f"--gres={gres} --cpus-per-task={cpus} --mem={mem} --time={stime}"
    if part:
        sbatch_opts += f" -p {shlex.quote(part)}"
    if cons:
        sbatch_opts += f" --constraint={shlex.quote(cons)}"
    if acct:
        sbatch_opts += f" -A {shlex.quote(acct)}"

    for ver in versions:
        out_root = os.path.join(pdir, f"output_AF/AF{ver}")
        log_dir = os.path.join(out_root, "logs")
        ensure_dir(out_root)
        ensure_dir(log_dir)

        for a3m_path, label in jobs:
            out_dir = os.path.join(out_root, label)
            ensure_dir(out_dir)

            # Skip if outputs already exist
            if not force and _af_has_outputs(out_dir, ver):
                print(f"[af] skip AF{ver} {label} (outputs exist)")
                continue

            if ver == "2":
                cmd = (
                    f"bash ./scripts/shell/RunAF2_Colabfold.sh "
                    f"{shlex.quote(a3m_path)} {shlex.quote(out_dir)} "
                    f"--num-models 5 --num-recycle 1 --model-type alphafold2_ptm"
                )
            else:
                cmd = (
                    f"bash ./scripts/shell/RunAF3_Colabfold.sh "
                    f"{shlex.quote(a3m_path)} {shlex.quote(out_dir)} "
                    f"--pdb=rank1"
                )

            if args.run_job_mode == "sbatch":
                log_path = os.path.join(log_dir, f"AF{ver}_{protein_id}__{label}.out")
                sb = f"sbatch {sbatch_opts} -o '{log_path}' --wrap {shlex.quote(cmd)}"
                _run(sb, "sbatch")
            else:
                _run(cmd, "inline")

    # Export canonical best PDBs
    for ver in versions:
        _export_canonical_best_pdbs(protein_id, ver, control_type)


def _af_has_outputs(out_dir: str, ver: str) -> bool:
    p = Path(out_dir)
    if not p.exists():
        return False
    if ver == "2":
        return bool(list(p.rglob("*.pdb")))
    else:
        return bool(list(p.rglob("*_model.cif")) or list(p.rglob("*.pdb")))


def _export_canonical_best_pdbs(protein_id: str, ver: str,
                                control_type: str = "PorterIDP") -> None:
    """Copy the best PDB from each AF cluster dir to a flat canonical name."""
    root = Path(_protein_dir(protein_id, control_type)) / f"output_AF/AF{ver}"
    if not root.exists():
        return
    for cluster_dir in sorted(root.iterdir()):
        if not cluster_dir.is_dir() or cluster_dir.name == "logs":
            continue
        label = cluster_dir.name
        best = _find_best_pdb(cluster_dir)
        if not best:
            continue
        dst = root / f"{label}.pdb"
        try:
            shutil.copy2(best, dst)
            print(f"[export] {dst.name} <= {Path(best).name}")
        except Exception as e:
            print(f"[export] failed {dst}: {e}")


def _find_best_pdb(d: Path) -> Optional[str]:
    """Heuristic: prefer ranked_0/rank_001/best_model, then any PDB."""
    for pattern in ["ranked_0*.pdb", "rank_001*.pdb", "best_model.pdb"]:
        hits = list(d.rglob(pattern))
        if hits:
            return str(sorted(hits)[0])
    pdbs = list(d.rglob("*.pdb"))
    return str(sorted(pdbs)[0]) if pdbs else None


def task_esmfold(protein_id: str, control_type: str, args: argparse.Namespace) -> None:
    """
    Sample sequences from each cluster MSA and run ESMFold.

    run_ESMFoldHF.py is pair-oriented (uses DATA_DIR), so for controls
    we prepare sampled FASTAs and call the ESMFold Python API directly
    in a small subprocess script.
    """
    pdir = _protein_dir(protein_id, control_type)
    cluster_dir = os.path.join(pdir, "output_msa_cluster")
    tmp_dir = os.path.join(pdir, "tmp_esmfold")
    ensure_dir(tmp_dir)

    if read_msa is None or greedy_select is None:
        print("[esmfold] ERR: read_msa / greedy_select not available")
        return

    sample_n = int(getattr(args, "cluster_sample_n", 10))
    cluster_a3ms = sorted(glob(os.path.join(cluster_dir, "ShallowMsa_*.a3m")))

    if not cluster_a3ms:
        print(f"[esmfold] no cluster MSAs for {protein_id}")
        return

    # Sample sequences per cluster -> individual FASTA files
    for a3m in cluster_a3ms:
        cl_id = Path(a3m).stem.split("_")[-1]
        dst = os.path.join(tmp_dir, f"_tmp_ShallowMsa_{cl_id}.fasta")
        if os.path.isfile(dst) and os.path.getsize(dst) > 0:
            continue
        entries = read_msa(a3m)
        picked = greedy_select(entries, min(sample_n, len(entries)))
        if not picked:
            continue
        names, seqs = zip(*picked)
        seqs = ["".join(ch for ch in s if ch.isalpha()) for s in seqs]
        write_fasta(list(names), list(seqs), dst)

    esm_model = getattr(args, "esm_model", "esm2")
    device = getattr(args, "esm_device", "cuda")

    out_esm = os.path.join(pdir, "output_esm_fold", esm_model)
    ensure_dir(out_esm)

    # Build a small inline Python command that folds each sampled FASTA
    # using ESMFold (HuggingFace) directly, writing PDBs to out_esm.
    fastas = sorted(glob(os.path.join(tmp_dir, "_tmp_ShallowMsa_*.fasta")))
    if not fastas:
        print(f"[esmfold] no sampled FASTAs for {protein_id}")
        return

    # Write a helper script that folds FASTAs using run_ESMFoldHF's internals.
    # We invoke it as a subprocess to keep GPU memory isolated.
    helper_script = os.path.join(pdir, "_run_esm_control.py")
    with open(helper_script, "w") as fh:
        fh.write(f"""\
#!/usr/bin/env python3
import sys, os
sys.path.insert(0, {os.path.abspath('.') !r})
from pathlib import Path
from utils.msa_utils import load_fasta

fasta_dir = {tmp_dir !r}
outdir = Path({out_esm !r})
outdir.mkdir(parents=True, exist_ok=True)
model_tag = {esm_model !r}
device = {device !r}

# Collect all sequences from sampled FASTAs
sequences = []
for fa in sorted(Path(fasta_dir).glob("_tmp_ShallowMsa_*.fasta")):
    ids, seqs = load_fasta(str(fa))
    cluster_label = fa.stem.replace("_tmp_", "")  # e.g., ShallowMsa_001
    for name, seq in zip(ids, seqs):
        clean = "".join(ch for ch in seq if ch.isalpha())
        sequences.append((f"{{cluster_label}}_{{name}}", clean))

if not sequences:
    print("[esmfold] no sequences to fold")
    sys.exit(0)

print(f"[esmfold] folding {{len(sequences)}} sequences with {{model_tag}} on {{device}}")

if model_tag == "esm2":
    from run_ESMFoldHF import run_esm2_fold, write_normalized_outputs
    result = run_esm2_fold(sequences, device, outdir=outdir)
else:
    from run_ESMFoldHF import run_esm3_fold_streaming, write_normalized_outputs
    result = run_esm3_fold_streaming(sequences, device, outdir, model_tag)

write_normalized_outputs(result, outdir, {protein_id !r}, model_tag, sequences, device)
print(f"[esmfold] done, outputs in {{outdir}}")
""")

    cmd = f"{shlex.quote(sys.executable)} {shlex.quote(helper_script)}"
    try:
        _run(cmd, args.run_job_mode)
    except Exception as e:
        print(f"[esmfold] ERROR: {e}")
    finally:
        # Clean up helper script
        try:
            os.remove(helper_script)
        except OSError:
            pass


def task_postprocess(protein_ids: List[str], control_type: str,
                     args: argparse.Namespace) -> None:
    """
    Compute cluster structural diversity for each control protein.

    For each protein, find all predicted PDBs across AF clusters,
    compute pairwise TM-scores, and write Analysis/df_diversity.csv.
    """
    import numpy as np
    import pandas as pd

    if compute_tmscore_align is None:
        print("[postprocess] ERR: compute_tmscore_align unavailable")
        return

    for protein_id in protein_ids:
        pdir = _protein_dir(protein_id, control_type)
        analysis_dir = os.path.join(pdir, "Analysis")
        ensure_dir(analysis_dir)

        # Collect best PDBs per cluster from AF2 (prefer AF2 for consistency)
        cluster_pdbs = {}
        for ver in ["2", "3"]:
            af_root = Path(pdir) / f"output_AF/AF{ver}"
            if not af_root.exists():
                continue
            for pdb in sorted(af_root.glob("*.pdb")):
                label = pdb.stem  # e.g. "ShallowMsa_001" or "DeepMsa"
                if label not in cluster_pdbs:
                    cluster_pdbs[label] = str(pdb)

        # Also check ESMFold outputs
        for esm_model in ["esm2", "esm3"]:
            esm_dir = Path(pdir) / "output_esm_fold" / esm_model
            if not esm_dir.exists():
                continue
            for pdb in sorted(esm_dir.glob("*.pdb")):
                label = f"esm_{pdb.stem}"
                if label not in cluster_pdbs:
                    cluster_pdbs[label] = str(pdb)

        if len(cluster_pdbs) < 2:
            print(f"[postprocess] {protein_id}: <2 cluster PDBs found, skipping diversity")
            continue

        # Compute pairwise TM-scores
        labels = sorted(cluster_pdbs.keys())
        n = len(labels)
        tm_matrix = np.full((n, n), np.nan)
        rows = []

        for i in range(n):
            tm_matrix[i, i] = 1.0
            for j in range(i + 1, n):
                try:
                    tm = compute_tmscore_align(
                        cluster_pdbs[labels[i]],
                        cluster_pdbs[labels[j]]
                    )
                    if isinstance(tm, (list, tuple)):
                        tm = float(tm[0])
                    else:
                        tm = float(tm)
                except Exception as e:
                    print(f"[postprocess] TM-score failed {labels[i]} vs {labels[j]}: {e}")
                    tm = np.nan
                tm_matrix[i, j] = tm
                tm_matrix[j, i] = tm
                rows.append({
                    "protein_id": protein_id,
                    "cluster_A": labels[i],
                    "cluster_B": labels[j],
                    "TM_score": tm,
                })

        # Summary statistics
        # Exclude diagonal and DeepMsa for cluster-only diversity
        shallow_labels = [l for l in labels if l.startswith("ShallowMsa")]
        shallow_tms = []
        for i, li in enumerate(labels):
            if li not in shallow_labels:
                continue
            for j, lj in enumerate(labels):
                if lj not in shallow_labels or j <= i:
                    continue
                if not np.isnan(tm_matrix[i, j]):
                    shallow_tms.append(tm_matrix[i, j])

        all_tms = [r["TM_score"] for r in rows if not np.isnan(r["TM_score"])]
        summary = {
            "protein_id": protein_id,
            "control_type": control_type,
            "n_clusters": len(shallow_labels),
            "n_total_pdbs": n,
            "mean_pairwise_TM": float(np.mean(all_tms)) if all_tms else np.nan,
            "median_pairwise_TM": float(np.median(all_tms)) if all_tms else np.nan,
            "min_pairwise_TM": float(np.min(all_tms)) if all_tms else np.nan,
            "max_pairwise_TM": float(np.max(all_tms)) if all_tms else np.nan,
            "std_pairwise_TM": float(np.std(all_tms)) if all_tms else np.nan,
            "cluster_mean_TM": float(np.mean(shallow_tms)) if shallow_tms else np.nan,
            "cluster_median_TM": float(np.median(shallow_tms)) if shallow_tms else np.nan,
        }

        # Write per-pair CSV
        df_pairs = pd.DataFrame(rows)
        df_pairs.to_csv(os.path.join(analysis_dir, "df_diversity_pairs.csv"), index=False)

        # Write summary as single-row CSV
        df_summary = pd.DataFrame([summary])
        df_summary.to_csv(os.path.join(analysis_dir, "df_diversity.csv"), index=False)

        print(f"[postprocess] {protein_id}: {len(rows)} pairs, "
              f"mean_TM={summary['mean_pairwise_TM']:.3f}, "
              f"n_clusters={summary['n_clusters']}")

    # Build global summary across all controls
    _build_global_controls_summary(protein_ids, control_type)


def _build_global_controls_summary(protein_ids: List[str], control_type: str) -> None:
    """Merge per-protein df_diversity.csv into a global controls summary."""
    import pandas as pd

    frames = []
    for pid in protein_ids:
        csv = os.path.join(_protein_dir(pid, control_type), "Analysis", "df_diversity.csv")
        if os.path.isfile(csv):
            frames.append(pd.read_csv(csv))

    if not frames:
        print("[global] no diversity CSVs found")
        return

    df = pd.concat(frames, ignore_index=True)
    out = os.path.join(CONTROLS_BASE_DIR, control_type, "controls_diversity_summary.csv")
    df.to_csv(out, index=False)
    print(f"[global] wrote {out} ({len(df)} proteins)")

    # Print summary statistics
    valid = df.dropna(subset=["mean_pairwise_TM"])
    if len(valid):
        print(f"[global] mean of means: {valid['mean_pairwise_TM'].mean():.3f}")
        print(f"[global] median of means: {valid['mean_pairwise_TM'].median():.3f}")


def task_full_pipeline(protein_id: str, control_type: str,
                       args: argparse.Namespace) -> None:
    """Run all stages sequentially for a single control protein."""
    print(f"=== full_pipeline :: {protein_id} ===", flush=True)

    # Force inline for the inner stages when running the full pipeline
    inner_args = deepcopy(args)
    inner_args.run_job_mode = "inline"

    task_load(protein_id, control_type, inner_args)
    task_get_msa(protein_id, control_type, inner_args)

    if phytree_from_msa is not None and args.cluster_alg == "tree":
        task_tree(protein_id, control_type, inner_args)

    task_cluster_msa(protein_id, control_type, inner_args)
    task_af(protein_id, control_type, inner_args)

    # ESMFold is optional (needs GPU)
    try:
        task_esmfold(protein_id, control_type, inner_args)
    except Exception as e:
        print(f"[full_pipeline] ESMFold skipped: {e}")

    task_postprocess([protein_id], control_type, inner_args)
    print(f"=== full_pipeline :: {protein_id} DONE ===", flush=True)


def _submit_protein_job(protein_id: str, control_type: str,
                        run_mode: str, args: argparse.Namespace) -> None:
    """Submit a single-protein pipeline job via sbatch."""
    pdir = _protein_dir(protein_id, control_type)
    jobs_dir = os.path.join(CONTROLS_BASE_DIR, control_type, "jobs")
    ensure_dir(jobs_dir)

    log = os.path.join(jobs_dir, f"{run_mode}_{protein_id}.out")

    # Build the inner command (runs inline inside the sbatch job)
    inner_cmd = (
        f"{shlex.quote(sys.executable)} {shlex.quote(os.path.abspath(__file__))} "
        f"--run_mode {run_mode} "
        f"--protein_ids {shlex.quote(protein_id)} "
        f"--control_type {shlex.quote(control_type)} "
        f"--run_job_mode inline "
        f"--cluster_alg {args.cluster_alg} "
        f"--af_ver {getattr(args, 'af_ver', '2')}"
    )

    gres = getattr(args, "sbatch_gres", "gpu:1")
    cpus = int(getattr(args, "sbatch_cpus", 8))
    mem = getattr(args, "sbatch_mem", "40G")
    stime = getattr(args, "sbatch_time", "24:00:00")
    part = getattr(args, "sbatch_partition", None)
    cons = getattr(args, "sbatch_constraint", None)
    acct = getattr(args, "sbatch_account", None)

    sbatch_opts = f"--gres={gres} --cpus-per-task={cpus} --mem={mem} --time={stime}"
    if part:
        sbatch_opts += f" -p {shlex.quote(part)}"
    if cons:
        sbatch_opts += f" --constraint={shlex.quote(cons)}"
    if acct:
        sbatch_opts += f" -A {shlex.quote(acct)}"

    sb = f"sbatch {sbatch_opts} -o {shlex.quote(log)} --wrap {shlex.quote(inner_cmd)}"
    print(f"[submit] {sb}", flush=True)
    _run(sb, "sbatch")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="Controls pipeline: single-protein MSA clustering + structure prediction",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--run_mode", required=True,
                   choices=["load", "get_msa", "tree", "cluster_msa",
                            "run_AF", "run_esmfold", "postprocess",
                            "full_pipeline"])
    p.add_argument("--protein_ids", nargs="+", required=True,
                   help="UniProt IDs (e.g., Q9NZC2 P12345) or ALL")
    p.add_argument("--control_type", default="PorterIDP",
                   help="Control set name (subdirectory under Pipeline/Controls/)")
    p.add_argument("--run_job_mode", default="inline", choices=["inline", "sbatch"])

    # Clustering options (same as fold-pair pipeline)
    p.add_argument("--cluster_alg", default="tree", choices=["hdbscan", "tree", "ahc"])
    p.add_argument("--cluster_max_clusters", type=int, default=100)
    p.add_argument("--cluster_min_output_size", type=int, default=200)
    p.add_argument("--cluster_min_neff", type=int, default=50)
    p.add_argument("--cluster_neff_id_thresh", type=float, default=0.8)
    p.add_argument("--cluster_frac_gaps_cutoff", type=float, default=0.7)
    p.add_argument("--cluster_sample_cap", type=int, default=5000)
    p.add_argument("--cluster_sample_seed", type=int, default=12345)
    p.add_argument("--hdbscan_min_cluster_size", type=int, default=200)
    p.add_argument("--hdbscan_min_samples", default="auto")
    p.add_argument("--hdbscan_cluster_selection", choices=["eom", "leaf"], default="eom")

    # AF options
    p.add_argument("--af_ver", default="2", choices=["2", "3", "both"])
    p.add_argument("--force_rerun_AF", default="FALSE", choices=["TRUE", "FALSE"])

    # ESMFold options
    p.add_argument("--cluster_sample_n", type=int, default=10)
    p.add_argument("--esm_model", default="esm2", choices=["esm2", "esm3", "both"])
    p.add_argument("--esm_device", default="cuda", choices=["auto", "cpu", "cuda"])

    # Tree options
    p.add_argument("--tree_max_seqs", type=int, default=5000)
    p.add_argument("--tree_seed", type=int, default=123)

    # Slurm options
    p.add_argument("--sbatch_gres", default="gpu:1")
    p.add_argument("--sbatch_cpus", type=int, default=8)
    p.add_argument("--sbatch_mem", default="40G")
    p.add_argument("--sbatch_time", default="24:00:00")
    p.add_argument("--sbatch_partition", default=None)
    p.add_argument("--sbatch_constraint", default=None)
    p.add_argument("--sbatch_account", default=None,
                   help="Slurm account (e.g., course-52017-25 to use higher fair-share).")

    args = p.parse_args()

    # Resolve protein list
    if args.protein_ids == ["ALL"]:
        proteins = list_control_proteins(args.control_type)
        if not proteins:
            sys.exit(f"No control proteins found under {CONTROLS_BASE_DIR}/{args.control_type}")
        print(f"[controls] Found {len(proteins)} proteins for {args.control_type}")
    else:
        proteins = args.protein_ids

    # Postprocess handles the full list at once
    if args.run_mode == "postprocess":
        task_postprocess(proteins, args.control_type, args)
        return

    # Per-protein dispatch
    for protein_id in proteins:
        print(f"=== {args.run_mode} :: {protein_id} ({args.control_type}) ===", flush=True)

        # sbatch mode for full_pipeline: submit one job per protein
        if args.run_job_mode == "sbatch" and args.run_mode == "full_pipeline":
            _submit_protein_job(protein_id, args.control_type, "full_pipeline", args)
            continue

        if args.run_mode == "load":
            task_load(protein_id, args.control_type, args)
        elif args.run_mode == "get_msa":
            task_get_msa(protein_id, args.control_type, args)
        elif args.run_mode == "tree":
            task_tree(protein_id, args.control_type, args)
        elif args.run_mode == "cluster_msa":
            task_cluster_msa(protein_id, args.control_type, args)
        elif args.run_mode == "run_AF":
            task_af(protein_id, args.control_type, args)
        elif args.run_mode == "run_esmfold":
            if args.esm_model == "both":
                for m in ("esm2", "esm3"):
                    a = deepcopy(args)
                    a.esm_model = m
                    task_esmfold(protein_id, args.control_type, a)
            else:
                task_esmfold(protein_id, args.control_type, args)
        elif args.run_mode == "full_pipeline":
            task_full_pipeline(protein_id, args.control_type, args)


if __name__ == "__main__":
    main()
