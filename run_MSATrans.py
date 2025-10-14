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

    # Derive pair root to locate DeepMsa.* robustly (works for .../output_msa_cluster and .../output_msa_cluster_pairtrim)
    if len(args.input_msas) == 1 and os.path.isdir(args.input_msas[0]):
        base = args.input_msas[0]
    else:
        base = os.path.dirname(input_a3ms[0]) if input_a3ms else "."
    pair_root = (base
                 .replace("/output_msa_cluster_pairtrim", "")
                 .replace("\\output_msa_cluster_pairtrim", "")
                 .replace("/output_msa_cluster", "")
                 .replace("\\output_msa_cluster", ""))

    deep_candidates = [
        os.path.join(pair_root, "output_get_msa", "DeepMsa_pairtrim.a3m"),  # preferred
        os.path.join(pair_root, "output_get_msa", "DeepMsa.a3m_pairtrim"),  # legacy alt
        os.path.join(pair_root, "output_get_msa", "DeepMsa.a3m"),           # original
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
    pattern = None
    if args.clean == 'all':
        pattern = '*.npy'
    elif args.clean == 'matched':
        # delete only files produced by this run configuration
        if args.keyword:  # be conservative: with empty keyword, do nothing
            prefix = f"{args.model}_{args.keyword}_"
            pattern = f"{prefix}*.npy"

    if pattern:
        print(f"Removing old files: {args.o}/{pattern}")
        old = glob(os.path.join(args.o, pattern))
        print(old)
        for f_old in old:
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
            msa_rows = _least_gappy_first(msa_rows)

            # reduce MSA rows to a manageable depth (keeps first row as query)
            inputs = greedy_select(msa_rows, num_seqs=128)
            batch_labels, batch_strs, batch_tokens = batch_converter([inputs])
            batch_tokens = batch_tokens.to(next(mdl.parameters()).device)

            # ---- SAFETY CAP on sequence length ----
            seqlen = batch_tokens.size(-1)
            if seqlen > 1024:
                print(f"[MSAT] Warning: MSA length {seqlen} > 1024; cropping to 1024.")
                batch_tokens = batch_tokens[:, :, :1024]


            # --------------------------------------- SAVE TO NPZ FILES ---------------------------------------

            print('MSA-Transformer predicting...')
            pred = mdl.predict_contacts(batch_tokens)[0].detach().cpu().numpy()  # (L_pred, L_pred)

            # --- Build idx (seed-frame indices of row-1 non-gap columns) ---
            # NOTE: row-1 in the batch is msa_rows[0] after your selection/subsampling logic.
            # IMPORTANT: if you subselect depth BEFORE tokenization, ensure msa_rows[0] is still first.
            q_hdr, q_seq = inputs[0] if 'inputs' in locals() else msa_rows[0]
            idx_seed = np.asarray(_ungapped_seed_indices(q_seq), dtype=np.int32)  # in seed frame

            # If you ALSO drop all-gap columns across selected rows before tokenization, intersect here:
            # keep_cols = np.asarray([j for j in range(len(_a3m_cols(q_seq))) if any(_a3m_cols(s)[j].isupper() for _,s in inputs)], dtype=np.int32)
            # idx_seed = idx_seed[np.isin(idx_seed, keep_cols)]

            # --- Save .npz with both cmap and idx ---
            out_base = f"{args.o}/{args.model}_{args.keyword}_{name}"
            np.savez(
                out_base + ".npz",
                cmap=pred.astype(np.float32),
                idx=idx_seed,  # indices in 0..L_seed-1 (seed pairwise alignment)
                # keep_cols=keep_cols if you maintain them (optional)
            )
            print(f"wrote {out_base}.npz  (cmap {pred.shape}, idx {len(idx_seed)})")

    print("Finished! Runtime for " + str(len(msas)) + " alignments = "
          + str(time.time() - start_time) + " seconds")
