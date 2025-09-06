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
import sys
# repo root = parent of scripts/
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import os, sys, argparse, traceback
from typing import Optional, List, Tuple

def eprint(*a, **k): print(*a, file=sys.stderr, **k)

# --- use your project utils (no re-implementations) ---
# Expect this script to run from repo root or that PYTHONPATH includes repo root.
try:
    from utils.msa_utils import load_fasta   # returns (ids, seqs)
except Exception as ex:
    eprint("[infer_esm3] Failed to import project utils. Ensure you run from repo root or set PYTHONPATH.")
    eprint("Import error:", ex)
    sys.exit(3)

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
    if hasattr(model, "infer_pdb"):
        return model.infer_pdb(sequence)
    # manual tokenization if infer_pdb is absent
    tok = AutoTokenizer.from_pretrained(model_id)
    batch = tok([sequence], return_tensors="pt", add_special_tokens=False)
    batch = {k: (v.to(device) if isinstance(v, torch.Tensor) else v) for k,v in batch.items()}
    with torch.no_grad():
        out = model(**batch)
    for conv in ["to_pdb", "output_to_pdb"]:
        if hasattr(model, conv):
            return getattr(model, conv)(out)
    raise RuntimeError("[hf-esmfold] Could not convert outputs to PDB.")

# -------------------------------------- CLI ----------------------------------
def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="ESM3 inference helper (prints PDB to stdout).")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--sequence", type=str, help="AA sequence")
    g.add_argument("--fasta",    type=str, help="FASTA path (first record used)")
    ap.add_argument("--device",  type=str, default="auto", help="cpu|cuda|mps|auto")
    args = ap.parse_args(argv)

    name, seq = read_one_sequence(args)
    device = resolve_device(args.device)

    # Try ESM3 first, then HF fallback
    pdb_txt = try_esm3(seq, device)
    if pdb_txt is None:
        pdb_txt = run_hf_esmfold(seq, device)

    if not isinstance(pdb_txt, str) or "ATOM" not in pdb_txt:
        eprint("[infer_esm3] backend returned non-PDB.")
        return 3

    # print only PDB to stdout
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
