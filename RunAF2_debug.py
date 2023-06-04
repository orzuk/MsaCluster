import sys
import os
import argparse
import hashlib
import jax
import jax.numpy as jnp
import numpy as np
import re
import subprocess
from glob import glob
import pandas as pd
from collections import namedtuple



class args:
    output_dir = '/Users/steveabecassis/Desktop/Kaib/alphafold/'
    input_msa= '/Users/steveabecassis/PycharmProjects/AF_Cluster/output/2QKEE_002.a3m'
    recycles = 1
    model_num = 1
    seed = 0
    verbose = True
    deterministic = True
    af2_dir = '/Users/steveabecassis/Downloads/alphafold_params_2022-03-02/'


sys.path.append(args.af2_dir)
# sys.path.append('/opt/conda/lib/python3.10/site-packages')
# sys.path.append('/content/alphafold')

from alphafold.model import model
from alphafold.model import config
from alphafold.model import data

from alphafold.data import parsers
from alphafold.data import pipeline

from alphafold.common import protein
from alphafold.common import residue_constants

"""
Create an AlphaFold model runner
name -- The name of the model to get the parameters from. Options: model_[1-5]
"""
def make_model_runner(name, recycles):
  cfg = config.model_config(name)

  cfg.data.common.num_recycle = recycles
  cfg.model.num_recycle = recycles
  cfg.data.eval.num_ensemble = 1
  if args.deterministic:
    cfg.data.eval.masked_msa_replace_fraction = 0.0
    cfg.model.global_config.deterministic = True

  print('to here')
  params = data.get_model_haiku_params(name, args.af2_dir + 'data/')
  print('read params')

  return model.RunModel(cfg, params)




import requests as r
from Bio import SeqIO
from io import StringIO
def get_seq_from_uniprotId(cID):
    # cID='A0A3N5IQ47'
    baseUrl="http://www.uniprot.org/uniprot/"
    currentUrl=baseUrl+cID+".fasta"
    response = r.post(currentUrl)
    cData=''.join(response.text)
    Seq=StringIO(cData)
    pSeq=list(SeqIO.parse(Seq,'fasta'))
    return str(pSeq[0].seq)

"""
Create a feature dictionary for input to AlphaFold
runner - The model runner being invoked. Returned from `make_model_runner`
sequence - The target sequence being predicted
templates - The template features being added to the inputs
seed - The random seed being used for data processing
"""
def make_processed_feature_dict(runner, a3m_file, name="test", templates=None, seed=0,msa_idx=1):
  feature_dict = {}

  # assuming sequence is first entry in msa


  with open(a3m_file,'r') as msa_fil:
      sequence = msa_fil.read().splitlines()[msa_idx].strip()
    # seqId = msa_fil.read().splitlines()[msa_idx].strip()[1:]
    # # seq_org_msa = msa_fil.read().splitlines()[msa_idx+1].strip()[1:]
    # sequence = get_seq_from_uniprotId(seqId)
    # gap_idxs = [i for i in range(len(sequence)) if sequence[i] == '-']
    # sequence = sequence.replace('-','')

  feature_dict.update(pipeline.make_sequence_features(sequence, name, len(sequence)))

  with open(a3m_file,'r') as msa_fil:
    msa = pipeline.parsers.parse_a3m(msa_fil.read())

  # with open(a3m_file, 'r') as msa_fil:
  #   seqId = msa_fil.read().splitlines()[msa_idx].strip()[1:]
  #   sequence = get_seq_from_uniprotId(seqId)
  #   gap_idxs = [i for i in range(len(sequence)) if sequence[i] == '-']
  #   sequence = sequence.replace('-', '')
  #
  # with open(a3m_file, 'r') as msa_fil:
  #   seq_org_msa = msa_fil.read().splitlines()[msa_idx+1].strip()[1:]
  #   start_idx = re.search(r'[A-Z]',seq_org_msa).start()
  #   start_seq = seq_org_msa[start_idx:start_idx + 5]
  #   start_org_seq_idx = sequence.find(start_seq)

        # seq_org_msa[:7] = sequence[:7]



  # for i in range(len(msa)):
  #     seq = msa.sequences[i]
  #     clean_seq = ''.join([seq[idx] for idx in range(len(seq)) if idx not in gap_idxs ])
  #     msa.sequences[i] = clean_seq

  # for i in range(len(msa.deletion_matrix)):
  #     del_mat = msa.deletion_matrix[i][:len(msa.sequences[0])]
  #     msa.deletion_matrix[i] = del_mat

  feature_dict.update(pipeline.make_msa_features([msa]))

  if templates is not None:
    feature_dict.update(templates)
  else:
    feature_dict.update(empty_placeholder_template_features(num_templates=0, num_res=len(sequence)))


  processed_feature_dict = runner.process_features(feature_dict, random_seed=seed)

  return processed_feature_dict

"""
Make a set of empty features for no-template evalurations
"""
def empty_placeholder_template_features(num_templates, num_res):
  return {
      'template_aatype': np.zeros(
          (num_templates, num_res,
           len(residue_constants.restypes_with_x_and_gap)), dtype=np.float32),
      'template_all_atom_masks': np.zeros(
          (num_templates, num_res, residue_constants.atom_type_num),
          dtype=np.float32),
      'template_all_atom_positions': np.zeros(
          (num_templates, num_res, residue_constants.atom_type_num, 3),
          dtype=np.float32),
      'template_domain_names': np.zeros([num_templates], dtype=object),
      'template_sequence': np.zeros([num_templates], dtype=object),
      'template_sum_probs': np.zeros([num_templates], dtype=np.float32),
  }

"""
Package AlphaFold's output into an easy-to-use dictionary
prediction_result - output from running AlphaFold on an input dictionary
processed_feature_dict -- The dictionary passed to AlphaFold as input. Returned by `make_processed_feature_dict`.
"""
def parse_results(prediction_result, processed_feature_dict):
  b_factors = prediction_result[0]['plddt'][:,None] * prediction_result[0]['structure_module']['final_atom_mask']
  dist_bins = jax.numpy.append(0,prediction_result[0]["distogram"]["bin_edges"])
  dist_mtx = dist_bins[prediction_result[0]["distogram"]["logits"].argmax(-1)]
  contact_mtx = jax.nn.softmax(prediction_result[0]["distogram"]["logits"])[:,:,dist_bins < 8].sum(-1)

  out = {"unrelaxed_protein": protein.from_prediction(processed_feature_dict, prediction_result[0], b_factors=b_factors),
        "plddt": prediction_result[0]['plddt'],
        "pLDDT": prediction_result[0]['plddt'].mean(),
        "dists": dist_mtx,
        "adj": contact_mtx}

  out.update({"pae": prediction_result[0]['predicted_aligned_error'],
              "pTMscore": prediction_result[0]['ptm']})
  return out

def write_results(result, pdb_out_path):
  plddt = float(result['pLDDT'])
  ptm = float(result["pTMscore"])
  print('plddt: %.3f' % plddt)
  print('ptm  : %.3f' % ptm)

  pdb_lines = protein.to_pdb(result["unrelaxed_protein"])
  with open(pdb_out_path, 'w') as f:
    f.write(pdb_lines)

#
# model_name = "model_{}_ptm".format(args.model_num)
# results_key = model_name + "_seed_{}".format(args.seed)
# runner = make_model_runner(model_name, args.recycles)
#
# print('made model runner')
#
# name=os.path.basename(args.input_msa).replace('.a3m','')
#
# features = make_processed_feature_dict(runner, args.input_msa, name=name, seed=args.seed)
# result = parse_results(runner.predict(features, random_seed=args.seed), features)
# write_results(result, args.output_dir+name+'.pdb')
#
#

