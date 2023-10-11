from typing import List, Tuple, Optional, Dict, NamedTuple, Union, Callable
import itertools
import os
import string
from pathlib import Path
import time
from protein_utils import *

import numpy as np
import torch
from scipy.spatial.distance import squareform, pdist, cdist
import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
from Bio import SeqIO
import biotite.structure as bs
import biotite.structure.io as bsio
from tqdm import tqdm
import argparse
import esm

import requests

torch.set_grad_enabled(False)

# This is an efficient way to delete lowercase characters and insertion characters from a string
deletekeys = dict.fromkeys(string.ascii_lowercase)
deletekeys["."] = None
deletekeys["*"] = None
translation = str.maketrans(deletekeys)


def read_sequence(filename: str) -> Tuple[str, str]:
    """ Reads the first (reference) sequences from a fasta or MSA file."""
    record = next(SeqIO.parse(filename, "fasta"))
    return record.description, str(record.seq)


def remove_insertions(sequence: str) -> str:
    """ Removes any insertions into the sequence. Needed to load aligned sequences in an MSA. """
    return sequence.translate(translation)


def read_msa(filename: str) -> List[Tuple[str, str]]:
    """ Reads the sequences from an MSA file, automatically removes insertions."""
    return [(record.description, remove_insertions(str(record.seq))) for record in SeqIO.parse(filename, "fasta")]


# Subsampling MSA


def plot_contacts_and_predictions(
    predictions: Union[torch.Tensor, np.ndarray],
    contacts: Union[torch.Tensor, np.ndarray],
    ax: Optional[mpl.axes.Axes] = None,
    # artists: Optional[ContactAndPredictionArtists] = None,
    cmap: str = "Blues",
    ms: float = 1,
    title: Union[bool, str, Callable[[float], str]] = True,
    animated: bool = False,
) -> None:

    if isinstance(predictions, torch.Tensor):
        predictions = predictions.detach().cpu().numpy()
    if isinstance(contacts, torch.Tensor):
        contacts = contacts.detach().cpu().numpy()
    if ax is None:
        ax = plt.gca()

    seqlen = contacts.shape[0]
    relative_distance = np.add.outer(-np.arange(seqlen), np.arange(seqlen))
    bottom_mask = relative_distance < 0
    masked_image = np.ma.masked_where(bottom_mask, predictions)
    invalid_mask = np.abs(np.add.outer(np.arange(seqlen), -np.arange(seqlen))) < 6
    predictions = predictions.copy()
    predictions[invalid_mask] = float("-inf")

    topl_val = np.sort(predictions.reshape(-1))[-seqlen]
    pred_contacts = predictions >= topl_val
    true_positives = contacts & pred_contacts & ~bottom_mask
    false_positives = ~contacts & pred_contacts & ~bottom_mask
    other_contacts = contacts & ~pred_contacts & ~bottom_mask

    if isinstance(title, str):
        title_text: Optional[str] = title
    elif title:
        long_range_pl = compute_precisions(predictions, contacts, minsep=24)[
            "P@L"
        ].item()
        if callable(title):
            title_text = title(long_range_pl)
        else:
            title_text = f"Long Range P@L: {100 * long_range_pl:0.1f}"
    else:
        title_text = None

    img = ax.imshow(masked_image, cmap=cmap, animated=animated)
    oc = ax.plot(*np.where(other_contacts), "o", c="grey", ms=ms)[0]
    fn = ax.plot(*np.where(false_positives), "o", c="r", ms=ms)[0]
    tp = ax.plot(*np.where(true_positives), "o", c="b", ms=ms)[0]
    ti = ax.set_title(title_text) if title_text is not None else None
    # artists = ContactAndPredictionArtists(img, oc, fn, tp, ti)

    ax.axis("square")
    ax.set_xlim([0, seqlen])
    ax.set_ylim([0, seqlen])
    plt.savefig('%s.pdf' % title, bbox_inches='tight')

if __name__=='__main__':

    p = argparse.ArgumentParser(description=
    """
    Run esm-1b and MSA transformer model.
    """)



    p.add_argument("--input_msas", nargs='*', action='store',help='Path to msas to use in prediction.')
    p.add_argument("-o", action="store", help='name of output directory to write contact maps to.')
    p.add_argument("--model", action='store', default='msa_t', help="Model: `esm1b` or `msa_t` (default is 'msa_t')")
    p.add_argument('--keyword', action='store', default='', help="Keyword for this prediction")
    p.add_argument("--test", action='store_true', help='Tests first 3 constructs.')
    p.add_argument("--parallel", action='store_true', help='Runs in parallel using Pandarallel.')

    args = p.parse_args()


    # class args:
    #     input_msas = '/Users/steveabecassis/Desktop/cluster_files/colabfold_new/output_pipeline_1ebo/cluster_msa_output'
    #     o = '/Users/steveabecassis/Desktop/cluster_files/colabfold_new/output_pipeline_1ebo'
    #     model = 'msa_t'
    #     keyword = 'test'
    #     test = False
    #     parallel = False




    os.makedirs(args.o, exist_ok=True)
    args.test = True
    if args.test:
        args.input_msas = args.input_msas[:3]

    # msa_cluster_path = args.input_msas
    # input_msas = [f'{msa_cluster_path}' + msa for msa in os.listdir(msa_cluster_path) if 'a3m' in str(msa)]

    #
    # msas = {
    #     os.path.basename(msa_fil).replace('.a3m',''): read_msa(msa_fil)
    #     for msa_fil in input_msas
    # }
    print('####')
    print(args.input_msas)
    print('####')
    print("Is running with GPU? " + str(torch.cuda.is_available()))

    start_time = time.time()
    msas = {
        os.path.basename(msa_fil).replace('.a3m', ''): read_msa(args.input_msas[0] + '/'+ msa_fil)
        for msa_fil in os.listdir(args.input_msas[0])
    }

    # msas = {
    #     os.path.basename(msa_fil).replace('.a3m',''): read_msa(msa_fil)
    #     for msa_fil in args.input_msas
    # }

    sequences = {
        name: msa[0] for name, msa in msas.items()
    }

    if args.model=='esm1b':
        mdl, mdl_alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
    elif args.model=='msa_t':
        mdl, mdl_alphabet = esm.pretrained.esm_msa1b_t12_100M_UR50S()
    else:
        print('model not understood.')

    mdl = mdl.eval()
    batch_converter = mdl_alphabet.get_batch_converter()

    if args.model=='esm1b':
        for name, inputs in sequences.items():
            batch_labels, batch_strs, batch_tokens = batch_converter([inputs])
            batch_tokens = batch_tokens.to(next(mdl.parameters()).device)
            print('predicting...')
            pred = mdl.predict_contacts(batch_tokens)[0].cpu()
            pred = pred.detach().cpu().numpy()
            print(np.sum(pred))
            np.savetxt("%s/%s_%s_%s.npy" %  (args.o, args.model, args.keyword, name), pred)
            print("wrote %s/%s_%s_%s.npy" % (args.o, args.model, args.keyword, name))

    elif args.model=='msa_t':
        for name, inputs in msas.items():
            inputs = greedy_select(inputs, num_seqs=128) # can change this to pass more/fewer sequences
            batch_labels, batch_strs, batch_tokens = batch_converter([inputs])
            batch_tokens = batch_tokens.to(next(mdl.parameters()).device)
            print('predicting...')
            pred = mdl.predict_contacts(batch_tokens)[0].cpu()
            pred = pred.detach().cpu().numpy()
            print(np.sum(pred))
            np.savetxt("%s/%s_%s_%s.npy" % (args.o, args.model, args.keyword, name), pred)
            print("wrote %s/%s_%s_%s.npy" % (args.o, args.model, args.keyword, name))

    print("Finished! Runtime for " + str(len(msas.items())) + " alignments = " + str(time.time()-start_time) + " seconds")
