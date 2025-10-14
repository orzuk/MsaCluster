#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scripts/infer_esm3.py
Reads a protein sequence (via --sequence or --fasta using your utils) and prints a PDB to stdout.
Logs go to stderr.

Priority:
  1) Try EvolutionaryScale ESM3 (if available in this environment).
  2) Fallback to 🤗 Transformers ESMFold ("facebook/esmfold_v1").

Usage:
  python scripts/infer_esm3.py --sequence "ACDE..." --device cuda
  python scripts/infer_esm3.py --fasta /path/to/seqs.fasta --device cpu

Env:
  ESM3_MODEL (optional), HF_HOME / HF_HUB_OFFLINE (optional for huggingface cache/offline)
"""
from __future__ import annotations
from pathlib import Path
import os, sys, argparse, traceback, glob
from typing import Optional, List, Tuple
import contextlib
# repo root = parent of scripts/
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# --- use project utils (no re-implementations) ---
from utils.msa_utils import load_fasta   # returns (ids, seqs)


# Reduce CUDA fragmentation / improve allocation behavior
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF",
                      "expandable_segments:True,max_split_size_mb:64")

def eprint(*a, **k): print(*a, file=sys.stderr, **k)

def run_multi_fasta_dir(fasta_dir: str, device: str = "cuda", model_id: str = "facebook/esmfold_v1") -> int:
    """
    Process all *.fasta in `fasta_dir` with a single HF ESMFold model load on GPU.
    Prints one block per success to STDOUT delimited by:
        >>>PDB_START <name>\n ...PDB...\n>>>PDB_END <name>
    Logs go to STDERR. Returns 0 on full success, 2 if any sample failed.
    """
    try:
        from transformers import AutoTokenizer, EsmForProteinFolding
        import torch
    except Exception as e:
        eprint(f"[hf-esmfold] import failure: {e}")
        return 2

    pat1 = os.path.join(fasta_dir, "*.fasta")
    pat2 = os.path.join(fasta_dir, "*.fa")
    fastas = sorted(set(glob.glob(pat1) + glob.glob(pat2)))
    eprint(f"[hf-esmfold] scanning {fasta_dir} → found {len(fastas)} FASTAs")
    if not fastas:
        eprint(f"[hf-esmfold] no FASTAs in {fasta_dir} (patterns tried: {pat1}, {pat2})")
        return 2

    # Load once on GPU
    eprint(f"[hf-esmfold] single-load {model_id} on {device}")
    model = EsmForProteinFolding.from_pretrained(model_id).to(device).eval()

    # --- memory-friendly config (ADD) ---
    # Try to lower model chunk size if available (safe no-op if missing)
    for attr in ("set_chunk_size", "set_chunk_size_"):
        if hasattr(model, attr):
            try:
                getattr(model, attr)(64)  # if still OOM later, try 32
            except Exception:
                pass

    if torch.cuda.is_available():
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.set_float32_matmul_precision("medium")
        model = model.half()  # cast weights to fp16 on CUDA
        amp_ctx = torch.amp.autocast("cuda", dtype=torch.float16)
    else:
        amp_ctx = contextlib.nullcontext()
    # ------------------------------------

    tok = AutoTokenizer.from_pretrained(model_id)

    any_fail = False
    for fa in fastas:
        name = os.path.splitext(os.path.basename(fa))[0]
        # read first sequence
        try:
            with open(fa) as f:
                seq = "".join(line.strip() for line in f if not line.startswith(">"))
            if not seq:
                raise ValueError("empty sequence")
        except Exception as e:
            eprint(f"[hf-esmfold] FAIL read {name}: {e}")
            any_fail = True
            continue

        try:
            batch = tok([seq], return_tensors="pt", add_special_tokens=False)

            # NEW (only half floating tensors; keep indices as long)
            for k, v in list(batch.items()):
                if not hasattr(v, "to"):
                    continue
                v = v.to(device)
                # Only cast *floating-point* tensors on CUDA
                if device == "cuda" and v.dtype in (torch.float32, torch.bfloat16, torch.float64):
                    v = v.half()
                batch[k] = v

            with torch.no_grad(), amp_ctx:
                out = model(**batch)

            to_pdb = getattr(model, "to_pdb", None) or getattr(model, "output_to_pdb", None)
            if to_pdb is None:
                raise RuntimeError("Model lacks to_pdb/output_to_pdb")
            pdb = to_pdb(out)

            # Emit a parseable block for the parent
            sys.stdout.write(f">>>PDB_START {name}\n")
            sys.stdout.write(pdb if pdb.endswith("\n") else (pdb + "\n"))
            sys.stdout.write(f">>>PDB_END {name}\n")
            sys.stdout.flush()

        except RuntimeError as e:
            # OOM or similar — log and continue
            eprint(f"[hf-esmfold] OOM/FAIL for {name}: {e}")
            any_fail = True
            continue
        except Exception as e:
            eprint(f"[hf-esmfold] FAIL {name}: {e}")
            any_fail = True
            continue

    return 2 if any_fail else 0



def read_one_sequence(args) -> Tuple[str, str]:
    if args.sequence:
        seq = args.sequence.strip().upper().replace(" ", "")
        if not seq or any(c in seq for c in ">\t\n\r"):
            raise ValueError("Invalid --sequence.")
        return ("seq", seq)
    # FASTA path: use your utils.msa_utils.load_fasta
    ids, seqs = load_fasta(args.fasta)
    if not ids or not seqs or not seqs[0]:
        raise ValueError(f"No sequences in FASTA: {args.fasta}")
    return (ids[0] or "seq", seqs[0].strip().upper())

def resolve_device(req: str) -> str:
    if req and req.lower() in {"cpu","cuda","mps"}:
        return req.lower()
    try:
        import torch
        if torch.cuda.is_available(): return "cuda"
        if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available(): return "mps"
    except Exception:
        pass
    return "cpu"

# -------------------- Try ESM3 (EvolutionaryScale package) -------------------
def try_esm3(sequence: str, device: str) -> Optional[str]:
    try:
        import esm  # EvolutionaryScale package (not fair-esm)
    except Exception as e:
        eprint(f"[esm3] esm package not importable: {e}")
        return None

    # APIs vary by version; try a few known entry points conservatively.
    candidates = []
    # explicit small/medium/turbo open models, if exposed:
    for attr in ["esm3_sm_open_v1", "esm3_md_open_v1", "esm3_turbo_open_v1"]:
        fn = getattr(getattr(esm, "pretrained", object), attr, None)
        if callable(fn):
            candidates.append(("pretrained."+attr, lambda: fn()))
    # generic getter, if available:
    get_pretrained = getattr(getattr(esm, "pretrained", object), "get_pretrained_model", None)
    if callable(get_pretrained):
        model_name = os.environ.get("ESM3_MODEL", "esm3_sm_open_v1")
        candidates.append(("pretrained.get_pretrained_model:"+model_name, lambda: get_pretrained(model_name)))

    if not candidates:
        eprint("[esm3] No usable ESM3 loader found in this `esm` package. Skipping.")
        return None

    last_err = None
    for tag, make in candidates:
        try:
            eprint(f"[esm3] trying loader: {tag}")
            model = make()
            # move-to-device variants
            try: model = model.to(device)
            except Exception:
                try: model.to(device=device)
                except Exception: pass

            # prefer simple single-call API names if present
            for api in ["infer_pdb", "sequence_to_pdb", "predict_structure", "fold"]:
                if hasattr(model, api):
                    eprint(f"[esm3] using model.{api} on {device}")
                    out = getattr(model, api)(sequence)
                    if isinstance(out, str) and "ATOM" in out:
                        return out
                    # try converting dict/objects to PDB
                    for conv in ["to_pdb", "output_to_pdb", "as_pdb"]:
                        if hasattr(model, conv):
                            out = getattr(model, conv)(out)
                            if isinstance(out, str) and "ATOM" in out:
                                return out
                    raise RuntimeError("ESM3 returned non-PDB and no converter worked.")
            # If we get here, your local ESM3 doesn’t expose a stable structure API yet.
            raise RuntimeError("No recognized ESM3 structure API.")
        except Exception as e:
            last_err = e
            eprint(f"[esm3] loader {tag} failed: {e}")
            eprint(traceback.format_exc())

    eprint(f"[esm3] all attempts failed. Last error: {last_err}")
    return None

# -------------------- HF ESMFold fallback (no OpenFold needed) ----------------
def run_hf_esmfold(sequence: str, device: str) -> str:
    from transformers import EsmForProteinFolding, AutoTokenizer
    import torch
    model_id = "facebook/esmfold_v1"

    eprint(f"[hf-esmfold] loading {model_id} on {device}")
    model = EsmForProteinFolding.from_pretrained(model_id).to(device).eval()

    # ↓↓↓ reduce memory footprint on CUDA ↓↓↓
    try:
        import torch
        if device == "cuda":
            torch.set_float32_matmul_precision("medium")
            model = model.half()
    except Exception as _:
        pass

    if hasattr(model, "infer_pdb"):
        # Some builds support AMP internally; still try to help CUDA
        try:
            import torch
            if device == "cuda":
                with torch.cuda.amp.autocast(dtype=torch.float16):
                    return model.infer_pdb(sequence)
        except Exception:
            return model.infer_pdb(sequence)

    # manual tokenization path
    tok = AutoTokenizer.from_pretrained(model_id)
    batch = tok([sequence], return_tensors="pt", add_special_tokens=False)
    try:
        import torch
        if device == "cuda":
            batch = {k: (v.to(device).half() if isinstance(v, torch.Tensor) else v) for k, v in batch.items()}
        else:
            batch = {k: (v.to(device) if isinstance(v, torch.Tensor) else v) for k, v in batch.items()}
    except Exception:
        batch = {k: (v.to(device) if isinstance(v, torch.Tensor) else v) for k, v in batch.items()}

    with torch.no_grad():
        try:
            import torch
            if device == "cuda":
                with torch.cuda.amp.autocast(dtype=torch.float16):
                    out = model(**batch)
            else:
                out = model(**batch)
        except Exception:
            out = model(**batch)

    for conv in ["to_pdb", "output_to_pdb"]:
        if hasattr(model, conv):
            return getattr(model, conv)(out)

    raise RuntimeError("[hf-esmfold] Could not convert outputs to PDB.")

# -------------------------------------- CLI ----------------------------------
def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="ESM3/ESMFold helper")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--sequence",  type=str, help="AA sequence (single)")
    g.add_argument("--fasta",     type=str, help="FASTA file (first record)")
    g.add_argument("--fasta_dir", type=str, help="DIR with many *.fasta (batch mode)")
    ap.add_argument("--device",   type=str, default="auto", help="cpu|cuda|mps|auto")
    args = ap.parse_args(argv)

    device = resolve_device(args.device)

    # --- Batch mode: keep one model in memory on GPU ---
    if args.fasta_dir:
        rc = run_multi_fasta_dir(args.fasta_dir, device=device, model_id="facebook/esmfold_v1")
        return rc

    # --- Single sequence / single FASTA (existing behavior) ---
    name, seq = read_one_sequence(args)

    print("[infer_esm3] try_esm3 on input sequence:", name, seq)
    pdb_txt = try_esm3(seq, device)          # your existing native-ESM3 attempt
    if pdb_txt is None:
        pdb_txt = run_hf_esmfold(seq, device)  # your existing HF fallback

    if not isinstance(pdb_txt, str) or "ATOM" not in pdb_txt:
        eprint("[infer_esm3] backend returned non-PDB.")
        return 3

    sys.stdout.write(pdb_txt if pdb_txt.endswith("\n") else (pdb_txt + "\n"))
    sys.stdout.flush()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as e:
        eprint("Fatal error in infer_esm3.py:", e)
        eprint(traceback.format_exc())
        raise SystemExit(3)
