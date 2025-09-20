#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ESMFoldHF.py — unified runner for ESMFold (ESM2-backed) and ESM3

What this does
- `--model {esm2, esm3}` (default: esm2)
- Auto-detects device (CUDA/MPS/CPU); override with env `FORCE_CPU=1` or `--device`
- **Uses ONLY your existing utils** to get sequences (no new FASTA/PDB readers):
    * utils.msa_utils.load_fasta
    * utils.protein_utils.load_pdb_structure
    * utils.protein_utils.extract_protein_sequence
- Writes normalized outputs to: Pipeline/<pair>/output_esm_fold/<model>/
    structure.pdb, prediction.json, logs.txt
- Compatible with your repo layout (uses config.py MAIN_DIR/DATA_DIR and utils/*)

Add to config.py if you plan to run ESM3 via a local repo script:
    ESM_PATH = "/mnt/c/Code/Github/esm"         # path where you cloned Meta's ESM repo
    ESM3_INFER_SCRIPT = "scripts/infer_esm3.py"   # script (relative to ESM_PATH) that prints PDB to stdout

Run examples
    python3 ESMFoldHF.py -input 1fzpD_2frhA                 # esm2 by default
    python3 ESMFoldHF.py -input 1fzpD_2frhA --model esm3    # run esm3 via subprocess

"""
from __future__ import annotations
import os, sys
import json, time
import argparse
import subprocess
from typing import Dict, List, Tuple, Optional
from transformers import EsmForProteinFolding, AutoTokenizer
from pathlib import Path


# Quiet TensorFlow chatter if present
os.environ.setdefault("TRANSFORMERS_NO_TF", "1")
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")

# --------------------------------------------------------------------------------------
# Project config & utils (DO NOT reinvent)
# --------------------------------------------------------------------------------------
from config import *  # assume present
ESM_PATH = globals().get("ESM_PATH", os.environ.get("ESM_PATH", ""))  # <= add this
try:
    from config import ESM3_INFER_SCRIPT
except Exception:
    ESM3_INFER_SCRIPT = Path(MAIN_DIR) / "scripts/infer_esm3.py"

# Your utils — rely ONLY on these for I/O
from utils.msa_utils import load_fasta
from utils.protein_utils import extract_protein_sequence, process_sequence

# Core deps
import torch
import numpy as np  # noqa: F401


def _cluster_fastas_dir(pair_id: str) -> Path:
    return Path(DATA_DIR) / pair_id / "tmp_esmfold"

def _sequences_from_cluster_fastas(pair_id: str) -> List[Tuple[str, str]]:
    """
    Read sequences from Pipeline/<pair>/tmp_esmfold/_tmp_ShallowMsa_*.fasta
    Each record becomes one folding job. Names include the cluster id and
    a running index so you can trace back later.
    """
    tmp = _cluster_fastas_dir(pair_id)
    if not tmp.exists():
        return []

    seqs: List[Tuple[str, str]] = []
    # Shallow clusters
    for fa in sorted(tmp.glob("_tmp_ShallowMsa_*.fasta")):
        cl = fa.stem.split("_")[-1]  # e.g., "003"
        ids, sseqs = load_fasta(str(fa))  # returns (names, seqs)
        for i, (rid, s) in enumerate(zip(ids, sseqs), start=1):
            # Compact, sortable id. Keep original id in filename via rid if you prefer.
            name = f"ShallowMsa_{cl}__sample_{i:03d}"
            seqs.append((name, s))

    # (optional) Deep set if you ever call _sample_to_tmp_fastas(..., include_deep=True)
    deep = tmp / "_tmp_DeepMsa.fasta"
    if deep.exists():
        ids, sseqs = load_fasta(str(deep))
        for i, (rid, s) in enumerate(zip(ids, sseqs), start=1):
            name = f"DeepMsa__sample_{i:03d}"
            seqs.append((name, s))

    return seqs


# --------------------------------------------------------------------------------------
# Device selection
# --------------------------------------------------------------------------------------

def pick_device() -> str:
    if os.environ.get("FORCE_CPU", "0") == "1":
        return "cpu"
    if torch.cuda.is_available():
        return "cuda"
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return "mps"
    return "cpu"

# --------------------------------------------------------------------------------------
# Paths & sequence loading (using existing utils only)
# --------------------------------------------------------------------------------------
def ensure_outdir(pair_id: str, model_tag: str) -> Path:
    out = Path(DATA_DIR) / pair_id / "output_esm_fold" / model_tag
    out.mkdir(parents=True, exist_ok=True)
    return out

def get_sequences_for_pair(pair_id: str) -> List[Tuple[str, str]]:
    """
    NEW order:
    1) If tmp_esmfold samples exist -> fold ALL of them (this is your 10×K).
    2) Else fall back to the original pair FASTA (if present).
    3) Else extract the two reference chains from the PDBs.
    """
    sampled = _sequences_from_cluster_fastas(pair_id)
    if sampled:
        print(f"[info] Using {len(sampled)} cluster-sampled sequences from tmp_esmfold/")
        return sampled

    # --- original fallbacks ---
    seqs = maybe_read_pair_fasta(Path(DATA_DIR) / pair_id)
    if seqs:
        return seqs
    return sequences_from_pdbs(Path(DATA_DIR) / pair_id, pair_id)


def parse_pair_id(pair_id: str) -> Tuple[Tuple[str, str], Tuple[str, str]]:
    """Split like 1fzpD_2frhA -> (("1fzp","D"),("2frh","A"))"""
    left, right = pair_id.split("_", 1)
    def split(pc: str) -> Tuple[str, str]:
        pdbid, chain = pc[:-1], pc[-1]
        return pdbid.lower(), chain
    return split(left), split(right)


def maybe_read_pair_fasta(pdir: Path) -> Optional[List[Tuple[str, str]]]:
    """Try typical FASTA names via existing utils.msa_utils.load_fasta.
    load_fasta returns (IDs, seqs) -> convert to List[(id, seq)].
    """
    for rel in ("seqs.fasta", "sequences.fasta", "input.fasta", "fasta/seqs.fasta"):
        fp = pdir / rel
        if fp.exists():
            ids, seqs = load_fasta(str(fp))
            records = list(zip(ids, seqs))
            if records:
                print(f"[info] Using FASTA: {fp}")
                return records
    return None


def sequences_from_pdbs(pdir: Path, pair_id: str) -> List[Tuple[str, str]]:
    """Use existing utils.protein_utils helpers to load PDBs and extract sequences.
    extract_protein_sequence(pdb_file: str, chain: Optional[str], ca_only: bool=False) -> str
    """
    (pdb1, ch1), (pdb2, ch2) = parse_pair_id(pair_id)
    pdb_paths = [pdir / f"{pdb1}.pdb", pdir / f"{pdb2}.pdb"]
    names = [f"{pdb1}{ch1}", f"{pdb2}{ch2}"]
    chains = [ch1, ch2]

    seqs: List[Tuple[str, str]] = []
    for name, chain, pdb_path in zip(names, chains, pdb_paths):
        if not pdb_path.exists():
            hits = list(pdir.glob(f"{pdb_path.stem}*.pdb"))
            if not hits:
                raise FileNotFoundError(f"Expected PDB file not found: {pdb_path}")
            pdb_path = hits[0]
        sequence = extract_protein_sequence(str(pdb_path), chain=chain)
        if not sequence:
            raise RuntimeError(f"Empty sequence for {name} extracted from {pdb_path}")
        seqs.append((name, str(sequence)))
    return seqs


def load_esmfold(device="cuda"):
    # Try Meta fair-esm path (needs OpenFold)
    try:
        import esm
        return ("meta-esmfold", esm.pretrained.esmfold_v1().eval())
    except Exception as e_meta:
        # Fallback to HF port (no system OpenFold install required)
        model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1").to(device).eval()
        return ("hf-esmfold", model)

# --------------------------------------------------------------------------------------
# ESM2 runner (ESMFold)
# --------------------------------------------------------------------------------------
# --- add this helper anywhere near your ESM2 runner (top-level, not inside a function): ---
def _detect_esm2_checkpoint() -> str | None:
    """
    Try to point to where ESMFold weights live. fair-esm uses torch.hub cache.
    """
    try:
        import torch, os
        hub = Path(os.environ.get("TORCH_HOME", Path.home() / ".cache" / "torch"))
        for d in [hub / "hub" / "checkpoints", hub / "checkpoints"]:
            if d.exists():
                hits = sorted(d.glob("*esm*fold*.*"), key=lambda p: p.stat().st_mtime, reverse=True)
                if hits:
                    return str(hits[0])
    except Exception:
        pass
    return None



def run_esm2_fold(
    seqs,
    device,
    chunk_size: int = 64,
    num_recycles: int = 1,
    amp_dtype=None,     # torch.float16 / torch.bfloat16 / None
):
    """
    Run structure prediction using ESMFold (Meta or HF port).
    Memory savers:
      - trunk.chunk_size
      - num_recycles
      - autocast (fp16/bf16) when on CUDA
    """
    backend, model = load_esmfold(device)  # returns ("meta-esmfold"| "hf-esmfold", model)
    if device != "cpu":
        model = model.to(device)
    # memory knobs
    if hasattr(model, "trunk"):
        try:
            model.trunk.chunk_size = int(chunk_size)
        except Exception:
            pass
        if hasattr(model, "set_num_recycles"):
            try: model.set_num_recycles(int(num_recycles))
            except Exception: pass
        elif hasattr(model.trunk, "num_recycles"):
            try: model.trunk.num_recycles = int(num_recycles)
            except Exception: pass

    use_amp = (device == "cuda" and amp_dtype is not None)

    outputs = []
    for name, seq in seqs:
        print(f"[{backend}] predicting {name} (len={len(seq)}) on {device} …", flush=True)
        t0 = time.time()
        # sanitize (strip gaps/non-letters)
        try:
            seq = process_sequence(seq)
        except Exception:
            seq = "".join(ch for ch in seq if ch.isalpha())

        with torch.no_grad():
            if use_amp:
                with torch.cuda.amp.autocast(dtype=amp_dtype):
                    pdb_str = model.infer_pdb(seq)
            else:
                pdb_str = model.infer_pdb(seq)

        dt = time.time() - t0
        print(f"[{backend}] done in {dt:.1f}s")
        outputs.append({
            "name": name,
            "pdb": pdb_str,
            "plddt": None,
            "pae": None,
            "residue_index": list(range(1, len(seq) + 1)),
        })
    return {"backend": backend, "chains": outputs}




# --------------------------------------------------------------------------------------
# ESM3 runner (subprocess to your local ESM repo script)
# --------------------------------------------------------------------------------------

def run_esm3_fold(seqs: List[Tuple[str, str]], device: str) -> Dict:
    """
    Call a local ESM3 inference script that prints PDB to stdout.
    Expected CLI (adjust your script accordingly):
        python <ESM_PATH>/<ESM3_INFER_SCRIPT> --sequence "ACDE..." --device <cpu|cuda|mps>
    """
    if not ESM_PATH:
        raise RuntimeError("ESM_PATH is not set in config.py or environment.")
    script = Path(ESM_PATH) / ESM3_INFER_SCRIPT
    if not script.exists():
        raise FileNotFoundError(f"ESM3 inference script not found: {script}")

    outputs = []
    for name, seq in seqs:
        print(f"[esm3] predicting {name} (len={len(seq)}) via {script} …", flush=True)

        # 1) sanitize (removes '-', whitespace, non-ACDEFGHIKLMNPQRSTVWY)
        orig_len = len(seq)
        seq = process_sequence(seq)
        if not seq:
            print(f"[esm3] skip {name}: empty after sanitization (was {orig_len} aa)")
            continue
        if len(seq) != orig_len:
            print(f"[esm3] {name}: sanitized length {orig_len} -> {len(seq)}")

        # 2) write a tiny temp FASTA (safer than pushing the sequence via argv)
        import tempfile, os
        with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tf:
            tf.write(f">{name}\n{seq}\n")
            fasta_path = tf.name

        try:
            # 3) run from the ESM repo folder; pass --fasta
            env = os.environ.copy()
            # keep your current PYTHONPATH trick (if your infer script needs MsaCluster utils)
            repo_root = Path(__file__).resolve().parent
            env["PYTHONPATH"] = f"{repo_root}:{env.get('PYTHONPATH', '')}"

            cmd = [sys.executable, str(script), "--fasta", fasta_path, "--device", device]
            res = subprocess.run(
                cmd, capture_output=True, text=True,
                cwd=str(Path(script).parent),  # <— run inside the ESM repo
                env=env
            )
            if res.returncode != 0:
                raise RuntimeError(f"ESM3 subprocess failed: {res.stderr[:500]}")

            pdb_str = res.stdout
            outputs.append({
                "name": name,
                "pdb": pdb_str,
                "plddt": None,
                "pae": None,
                "residue_index": list(range(1, len(seq) + 1)),
            })
        finally:
            try:
                os.unlink(fasta_path)
            except Exception:
                pass
    return {"backend": "esm3-subprocess", "chains": outputs}

# --------------------------------------------------------------------------------------
# Normalization to ESM2-style layout
# --------------------------------------------------------------------------------------

# New function allowing us to normalize outputs for both ESM2 and ESM3
def write_normalized_outputs(result, outdir: Path, pair_id: str, model_tag: str, sequences, device: str) -> None:
    """
    For many samples we keep it flat (single directory) but encode the cluster
    in the filename, e.g. ShallowMsa_003__sample_007_<model>.pdb
    """
    outdir.mkdir(parents=True, exist_ok=True)

    index_rows = []
    for ch in result.get("chains", []):
        name = (ch.get("name") or "unknown").strip()  # e.g., ShallowMsa_003__sample_007
        pdb_txt = (ch.get("pdb") or "").strip()
        if not pdb_txt:
            continue
        pdb_path = outdir / f"{name}_{model_tag}.pdb"
        with open(pdb_path, "w") as f:
            f.write(pdb_txt if pdb_txt.endswith("\n") else pdb_txt + "\n")
        index_rows.append((name, str(pdb_path)))

    # Optional combined multi-MODEL file (handy for quick viewing)
    combo_path = outdir / f"{pair_id}__{model_tag}_combined.pdb"
    with open(combo_path, "w") as f:
        for i, ch in enumerate(result.get("chains", []), start=1):
            pdb_txt = (ch.get("pdb") or "").strip()
            if pdb_txt:
                f.write(f"MODEL     {i}\n{pdb_txt}\nENDMDL\n")

    # Index TSV
    with open(outdir / "samples_index.tsv", "w") as tsv:
        tsv.write("name\tpdb_path\n")
        for name, path in index_rows:
            tsv.write(f"{name}\t{path}\n")

    # Metadata + log (unchanged)
    meta_path = outdir / "pair_prediction_metadata.json"
    meta = {
        "pair_id": pair_id,
        "model": model_tag,
        "backend": result.get("backend"),
        "device": device,
        "n_predictions": len(index_rows),
        "outputs": {
            "per_sample_index": str(meta_path.parent / "samples_index.tsv"),
            "combined": str(combo_path),
        },
        "timestamps": {"written": time.strftime("%Y-%m-%d %H:%M:%S")},
    }
    with open(meta_path, "w") as jf:
        json.dump(meta, jf, indent=2)

    with open(outdir / "pair_prediction.log", "a") as lf:
        lf.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] pair={pair_id} model={model_tag} "
                 f"backend={result.get('backend')} device={device} "
                 f"n={len(index_rows)}\n")



# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------

# --- in main(), add timing around each stage and print a summary: ---
def main(argv: Optional[List[str]] = None) -> int:
    print(f"Running ESM, going to read parameters!", flush=True)

    parser = argparse.ArgumentParser(description="Run ESMFold (ESM2) or ESM3 and normalize outputs.")
    parser.add_argument("-input", dest="pair_id", required=True, help="Pair id (e.g., 1fzpD_2frhA) OR an existing directory path")
    parser.add_argument("--model", choices=["esm2", "esm3"], default="esm2")
    parser.add_argument("--device", choices=["auto", "cpu", "cuda", "mps"], default="auto")
    parser.add_argument("--esm_chunk", type=int, default=64, help="Chunk size for trunk attention (smaller → less VRAM).")
    parser.add_argument("--esm_recycles", type=int, default=1, help="Number of recycles (1–3; fewer → less VRAM).")
    parser.add_argument("--esm_dtype", choices=["auto","fp32","fp16","bf16"], default="auto",
               help="Autocast dtype on CUDA (auto→bf16 if available, else fp16).")
    args = parser.parse_args(argv)
    print(f"Running ESM, finished reading args parameters!", flush=True)

    t_start = time.time()
    device = pick_device() if args.device == "auto" else args.device

    # decide AMP dtype (optional)
    use_cuda = (device == "cuda")
    amp_dtype = None
    if args.esm_dtype == "bf16" or (args.esm_dtype == "auto" and use_cuda and torch.cuda.is_bf16_supported()):
        amp_dtype = torch.bfloat16
    elif args.esm_dtype == "fp16" or (args.esm_dtype == "auto" and use_cuda):
        amp_dtype = torch.float16

    dev_str = device
    if device == "cuda":
        try: dev_str += f" ({torch.cuda.get_device_name(0)})"
        except Exception: pass
    print(f"Running ESM on device: {dev_str}", flush=True)

    t0 = time.time()
    pair_path = Path(DATA_DIR) / args.pair_id
    if not pair_path.exists():
        raise FileNotFoundError(f"Pair directory not found: {pair_path}")
    sequences = get_sequences_for_pair(args.pair_id)
    t1 = time.time()

    model_tag = args.model
    outdir = ensure_outdir(args.pair_id, model_tag)

    t2 = time.time()
    if model_tag == "esm2":
        result = run_esm2_fold(
            sequences, device,
            chunk_size=args.esm_chunk,
            num_recycles=args.esm_recycles,
            amp_dtype=amp_dtype,
        )
    else:
        result = run_esm3_fold(sequences, device)
    t3 = time.time()

    write_normalized_outputs(result, outdir, args.pair_id, model_tag, sequences, device)
    t4 = time.time()

    print(f"[timing] load_seqs: {t1-t0:.1f}s | predict: {t3-t2:.1f}s | write: {t4-t3:.1f}s | total: {t4-t_start:.1f}s", flush=True)
    print(f"[done] Outputs: {outdir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise SystemExit(1)