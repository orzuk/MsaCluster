# import itertools
# from pathlib import Path
import time
from utils.protein_utils import *

import numpy as np
import torch
# from scipy.spatial.distance import squareform, pdist, cdist
# import pandas as pd
import matplotlib as mpl
import argparse
from glob import glob
import esm
mpl.use("agg")

torch.set_grad_enabled(False)

# === MSA helper utilities (seed-frame operations) ===
def _a3m_cols(seq: str):
    """Return A3M 'columns' (uppercase AA or '-') for a raw a3m row."""
    return [c for c in seq if c.isupper() or c == '-']

def _ungapped_seed_indices(seq: str):
    """Indices (0..L_seed-1) of columns where this row has a residue (uppercase)."""
    cols = _a3m_cols(seq)
    return [i for i, c in enumerate(cols) if c.isupper()]

def _least_gappy_first(rows):
    """
    rows: list[(header, sequence_str)]
    Return rows sorted descending by #uppercase in seed-frame columns;
    this makes row-1 (query) the least-gappy sequence.
    """
    return sorted(rows, key=lambda hs: sum(ch.isupper() for ch in _a3m_cols(hs[1])), reverse=True)


if __name__ == '__main__':

    p = argparse.ArgumentParser(description="""
    Run esm-1b and MSA transformer model.
    """)
    p.add_argument("--input_msas", nargs='*', action='store',
                   help='Path to MSAs to use in prediction (dir with *.a3m or list of *.a3m files).')
    p.add_argument("-o", action="store", help='Output directory for contact maps.')
    p.add_argument("-saveformat", action="store", help='Output format: text or pickle (default: .npy).')
    p.add_argument("--model", action='store', default='msa_t',
                   help="Model: `esm1b` or `msa_t` (default: 'msa_t').")
    p.add_argument('--keyword', action='store', default='', help="Keyword tag for this prediction.")
    p.add_argument('--clean', choices=['none', 'matched', 'all'], default='matched',
                    help = "What to delete in output dir before writing. "
                    "'none' = keep everything; 'matched' = delete only files that match "
                    "this run's model/keyword prefix; 'all' = delete all *.npy (legacy).")

    p.add_argument("--test", action='store_true', help='Test mode: only first 3 inputs.')
    p.add_argument("--parallel", action='store_true', help='(unused) Runs in parallel using Pandarallel.')
    args = p.parse_args()

    # --- I/O setup ---
    os.makedirs(args.o, exist_ok=True)
    print('####'); print(args.input_msas); print('####')
    print("Is running with GPU? " + str(torch.cuda.is_available()))
    start_time = time.time()

    if not args.input_msas:
        raise ValueError("No --input_msas provided.")

    # Accept either a single directory or explicit .a3m files
    if len(args.input_msas) == 1 and os.path.isdir(args.input_msas[0]):
        input_dir = args.input_msas[0]
        input_a3ms = sorted(glob(os.path.join(input_dir, "*.a3m")))
    else:
        input_a3ms = [p for p in args.input_msas if p.endswith(".a3m")]
        # also allow a directory accidentally mixed in the list
        extra_dirs = [d for d in args.input_msas if os.path.isdir(d)]
        for d in extra_dirs:
            input_a3ms.extend(sorted(glob(os.path.join(d, "*.a3m"))))
        input_a3ms = sorted(set(input_a3ms))

    if args.test:
        input_a3ms = input_a3ms[:3]

    # Read cluster MSAs
    msas = {
        os.path.splitext(os.path.basename(msa_fil))[0]: read_msa(msa_fil)
        for msa_fil in input_a3ms
    }

    # Derive pair root to locate DeepMsa.* only when the user explicitly passed a DeepMsa file
    if len(args.input_msas) == 1 and os.path.isdir(args.input_msas[0]):
        base = args.input_msas[0]
    else:
        base = os.path.dirname(input_a3ms[0]) if input_a3ms else "."
    pair_root = (base
                 .replace("/output_msa_cluster_pairtrim", "")
                 .replace("\\output_msa_cluster_pairtrim", "")
                 .replace("/output_msa_cluster", "")
                 .replace("\\output_msa_cluster", ""))

    # Only auto-add DeepMsa if the input list itself contains a DeepMsa file name
    is_input_explicit_deep = any("DeepMsa" in os.path.basename(p) for p in input_a3ms)
    if is_input_explicit_deep:
        deep_candidates = [
            os.path.join(pair_root, "output_get_msa", "DeepMsa_pairtrim.a3m"),  # preferred
            os.path.join(pair_root, "output_get_msa", "DeepMsa.a3m_pairtrim"),  # legacy alt
            os.path.join(pair_root, "output_get_msa", "DeepMsa.a3m"),  # original
        ]
        for deep_path in deep_candidates:
            if os.path.isfile(deep_path):
                print(f"[MSAT] Using Deep MSA: {deep_path}")
                msas["MSA_deep"] = read_msa(deep_path)
                break

    # For esm1b we need single sequences (first row of each MSA)
    sequences = {name: msa[0] for name, msa in msas.items()}

    # --- Load model ---
    if args.model == 'esm1b':
        mdl, mdl_alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
    elif args.model == 'msa_t':
        mdl, mdl_alphabet = esm.pretrained.esm_msa1b_t12_100M_UR50S()
    else:
        raise ValueError(f"Unknown --model '{args.model}' (expected 'esm1b' or 'msa_t').")

    mdl = mdl.eval()
    batch_converter = mdl_alphabet.get_batch_converter()

    # Apply cleaning policy (default: keep everything)
    if args.clean == 'all':
        patterns = ['*.npy', '*.npz']
    elif args.clean == 'matched':
        patterns = []
        if args.keyword:
            prefix = f"{args.model}_{args.keyword}_"
            patterns = [f"{prefix}*.npy", f"{prefix}*.npz"]
    else:
        patterns = []

    if patterns:
        for pattern in patterns:
            print(f"Removing old files: {args.o}/{pattern}")
            for f_old in glob(os.path.join(args.o, pattern)):
                try:
                    os.remove(f_old)
                except OSError:
                    pass


    # --- Predict ---
    if args.model == 'esm1b':
        for name, inputs in sequences.items():
            batch_labels, batch_strs, batch_tokens = batch_converter([inputs])
            batch_tokens = batch_tokens.to(next(mdl.parameters()).device)
            print('ESM1b predicting...')
            pred = mdl.predict_contacts(batch_tokens)[0].detach().cpu().numpy()
            if args.saveformat == "text":
                np.savetxt(f"{args.o}/{args.model}_{args.keyword}_{name}.npy", pred)
            else:
                np.save(f"{args.o}/{args.model}_{args.keyword}_{name}.npy", pred)
            print(f"wrote {args.o}/{args.model}_{args.keyword}_{name}.npy")

    elif args.model == 'msa_t':
        for name, msa_rows in msas.items():
            # 1) Put the least-gappy sequence first; this will be the query (row 1)
            msa_rows = _least_gappy_first(msa_rows)

            # 2) Subsample depth (keeps row 1)
            inputs = greedy_select(msa_rows, num_seqs=128)

            # 3) Seed-frame column bookkeeping BEFORE tokenization
            C = [_a3m_cols(s) for (_, s) in inputs]
            L_seed = len(C[0])
            assert all(len(x) == L_seed for x in C), "Selected rows differ in seed-frame length"

            # Columns to keep after removing columns that are all '-' across the selected rows
            keep_cols = [j for j in range(L_seed) if any(row[j].isupper() for row in C)]

            # Row-1 (query) non-gap indices in SEED frame
            q_hdr, q_seq = inputs[0]
            q_mask = set(_ungapped_seed_indices(q_seq))

            # Final idx in SEED frame = query non-gaps ∩ kept columns
            idx_seed = np.asarray([j for j in keep_cols if j in q_mask], dtype=np.int32)

            # 4) Build the sequences actually fed to the tokenizer (kept columns only)
            inputs_kept = []
            for (h, s) in inputs:
                cols = _a3m_cols(s)
                seq_kept = "".join(cols[j] for j in keep_cols)
                inputs_kept.append((h, seq_kept))

            # 5) Tokenize the kept-frame inputs
            batch_labels, batch_strs, batch_tokens = batch_converter([inputs_kept])
            batch_tokens = batch_tokens.to(next(mdl.parameters()).device)

            # ---- SAFETY CAP on sequence length ----
            seqlen = batch_tokens.size(-1)
            if seqlen > 1024:
                print(f"[MSAT] Warning: MSA length {seqlen} > 1024; cropping to 1024.")
                batch_tokens = batch_tokens[:, :, :1024]
                # If you ever hit this branch, you'd also want to trim idx_seed accordingly.
                # For your current lengths (~198), this won't trigger.

            # 6) Predict and save
            print('MSA-Transformer predicting...')
            pred = mdl.predict_contacts(batch_tokens)[0].detach().cpu().numpy()  # (L_pred, L_pred)

            out_base = f"{args.o}/{args.model}_{args.keyword}_{name}"
            np.savez(
                out_base + ".npz",
                cmap=pred.astype(np.float32),
                idx=idx_seed,                         # len(idx_seed) == pred.shape[0]
                keep_cols=np.asarray(keep_cols, np.int32),  # optional, handy for debugging
            )
            print(f"wrote {out_base}.npz  (cmap {pred.shape}, idx {len(idx_seed)})")


    print("Finished! Runtime for " + str(len(msas)) + " alignments = "
          + str(time.time() - start_time) + " seconds")
