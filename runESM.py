# import itertools
# import os
# from pathlib import Path
import time
from scripts.protein_utils import *

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


# Subsampling MSA

if __name__=='__main__':

    p = argparse.ArgumentParser(description=
    """
    Run esm-1b and MSA transformer model.
    """)

    p.add_argument("--input_msas", nargs='*', action='store',help='Path to msas to use in prediction.')
    p.add_argument("-o", action="store", help='name of output directory to write contact maps to.')
    p.add_argument("-saveformat", action="store", help='output file format (text or pickle).')
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

    # New! insert to MSA also the full alignment (calculating cmap can take long for this one)
    msas["MSA_deep"] = read_msa(args.input_msas[0].replace(
        "output_msa_cluster", "output_get_msa/DeepMsa.a3m"))

    # msas = {
    #     os.path.basename(msa_fil).replace('.a3m',''): read_msa(msa_fil)
    #     for msa_fil in args.input_msas
    # }

    sequences = {
        name: msa[0] for name, msa in msas.items()
    }

    if args.model == 'esm1b':
        mdl, mdl_alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
    elif args.model == 'msa_t':
        mdl, mdl_alphabet = esm.pretrained.esm_msa1b_t12_100M_UR50S()
    else:
        print('model not understood.')

    mdl = mdl.eval()
    batch_converter = mdl_alphabet.get_batch_converter()

    # Remove old cmaps in the same directory!
    print("Removing old files: " + args.o + '/*.npy')
    print(glob(args.o + '/*.npy'))
    for f_old in glob(args.o + '/*.npy'):
        os.remove(f_old)

    if args.model == 'esm1b':
        for name, inputs in sequences.items():
            batch_labels,batch_strs,batch_tokens = batch_converter([inputs])
            batch_tokens = batch_tokens.to(next(mdl.parameters()).device)
            print('ESM1b predicting...')
            pred = mdl.predict_contacts(batch_tokens)[0].cpu()
            pred = pred.detach().cpu().numpy()
            print(np.sum(pred))
            if args.saveformat == "text":
                np.savetxt("%s/%s_%s_%s.npy" % (args.o, args.model, args.keyword, name), pred)
            else:  # Save compact form
                np.save("%s/%s_%s_%s.npy" % (args.o, args.model, args.keyword, name), pred)
            print("wrote %s/%s_%s_%s.npy" % (args.o, args.model, args.keyword, name))

    elif args.model == 'msa_t':
        for name, inputs in msas.items():
            inputs = greedy_select(inputs, num_seqs=128) # can change this to pass more/fewer sequences
            batch_labels, batch_strs, batch_tokens = batch_converter([inputs])
            batch_tokens = batch_tokens.to(next(mdl.parameters()).device)
            print('MSA-Transformer predicting...')
            pred = mdl.predict_contacts(batch_tokens)[0].cpu()
            pred = pred.detach().cpu().numpy()
            print(np.sum(pred))
            if args.saveformat == "text":
                np.savetxt("%s/%s_%s_%s.npy" % (args.o, args.model, args.keyword, name), pred)
            else:
                np.save("%s/%s_%s_%s.npy" % (args.o, args.model, args.keyword, name), pred)
            print("wrote %s/%s_%s_%s.npy" % (args.o, args.model, args.keyword, name))

    print("Finished! Runtime for " + str(len(msas.items())) + " alignments = " + str(time.time()-start_time) + " seconds")
