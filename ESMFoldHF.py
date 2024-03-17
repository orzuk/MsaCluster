import os

# from ESMFold import convert_outputs_to_pdb, save_string_as_pdb
# from transformers import AutoTokenizer, EsmForProteinFolding
# from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
# from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from argparse import  ArgumentParser
import torch
from random import sample
import random
random.seed(10)

from transformers import AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
import argparse
from Bio import SeqIO
from Bio import PDB
from Bio.PDB import PDBParser
import os
# os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:128'
#iminuit==1.5.4
#tmscoring



# Specific conversion for atoms
def convert_outputs_to_pdb(outputs):
    final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
    outputs = {k: v.detach().to("cpu").numpy() for k, v in outputs.items()}
    final_atom_positions = final_atom_positions.detach().cpu().numpy()
    final_atom_mask = outputs["atom37_atom_exists"]
    pdbs = []
    for i in range(outputs["aatype"].shape[0]):
        aa = outputs["aatype"][i]
        pred_pos = final_atom_positions[i]
        mask = final_atom_mask[i]
        resid = outputs["residue_index"][i] + 1
        pred = OFProtein(
            aatype=aa,
            atom_positions=pred_pos,
            atom_mask=mask,
            residue_index=resid,
            b_factors=outputs["plddt"][i],
            chain_index=outputs["chain_index"][i] if "chain_index" in outputs else None,
        )
        pdbs.append(to_pdb(pred))
    return pdbs

def extract_protein_sequence(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("protein_structure", pdb_file)

    residue_sequence = ""

    # Iterate through the structure and extract the residue sequence
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue):
                    residue_sequence += PDB.Polypeptide.three_to_one(residue.get_resname())

    return residue_sequence

# Save string to pdb
def save_string_as_pdb(pdb_string, file_path):
    with open(file_path, 'w') as pdb_file:
        pdb_file.write(pdb_string)


# Read fasta file
def get_sequence_from_fasta(fasta_file_path):
    sequence = None
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        sequence = str(record.seq)
        break
    return sequence

if __name__ == '__main__':
    print('######################## memory_summary ####################################')
    print(torch.cuda.memory_summary())
    print('############################################################################')


    parser = ArgumentParser()
    parser.add_argument("-input",help="Should be a path")
    parser.add_argument("-name", help="msa name")

    args = parser.parse_args()

    device = "cuda:0" if torch.cuda.is_available() else "cpu"
    # device = "cuda" if torch.cuda.is_available() else "cpu"
    print("Running ESM-Fold on device: " + device)
    torch.backends.cuda.matmul.allow_tf32 = True

    input_ = args.input
    print('input:',input_)
    fold_pair = input_.replace("Pipeline/", "",1)  # The '1' indicates to replace the first occurrence only
    input_path = f'./Pipeline/{fold_pair}/output_msa_cluster'

    msas_files = os.listdir(f'./Pipeline/{fold_pair}/output_msa_cluster')
    print('Load model...!')
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    model.trunk.set_chunk_size(64)
    model.esm.float()
    model = model.to(device)
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    print('Finish to load model !')

    print('######################## memory_summary ####################################')
    print(torch.cuda.memory_summary())
    print('############################################################################')


    folds = fold_pair.split("_")
    fold1 = folds[0]
    fold2 = folds[1]
    path = f'./Pipeline/{fold_pair}'
    print(f'seq fold1:{path}/chain_pdb_files/{fold1}.pdb')
    seq_fold1 = extract_protein_sequence(f'{path}/chain_pdb_files/{fold1}.pdb')
    seq_fold2 = extract_protein_sequence(f'{path}/chain_pdb_files/{fold2}.pdb')

    print('################################seq_fold1############################################')
    print(seq_fold1)
    print('############################################################################')
    print('###################################seq_fold2#########################################')
    print(seq_fold2)
    print('############################################################################')



    inputs = tokenizer([seq_fold1], return_tensors="pt", add_special_tokens=False, padding=True).to(device)
    outputs = model(**inputs)
    folded_positions = outputs.positions
    pdb = convert_outputs_to_pdb(outputs)
    print(f'./Pipeline/{fold_pair}/output_esm_fold/{fold1}_esm.pdb')
    save_string_as_pdb(pdb[0], f'./Pipeline/{fold_pair}/output_esm_fold/{fold1}_esm.pdb')


    print('######################## memory_summary ####################################')
    print(torch.cuda.memory_summary())
    print('############################################################################')


    # if len(seq_fold2)
    inputs = tokenizer([seq_fold2], return_tensors="pt", add_special_tokens=False, padding=True).to(device)
    outputs = model(**inputs)
    folded_positions = outputs.positions
    pdb = convert_outputs_to_pdb(outputs)
    save_string_as_pdb(pdb[0], f'./Pipeline/{fold_pair}/output_esm_fold/{fold2}_esm.pdb')


    print('######################## memory_summary ####################################')
    print(torch.cuda.memory_summary())
    print('############################################################################')


    for msa in msas_files:
        with open(f'{input_path}/{msa}', 'r') as msa_fil:
            seq = msa_fil.read().splitlines()
        msa_name = msa[:-4]
        seqs = [i.replace('-', '') for i in seq if '>' not in i]
        if len(seqs) > 10:
            seqs = sample(seqs,10)

        # model.esm = model.esm.half()


        for i in range(len(seqs)):
            try:
                print(f'Get ESM prediction {i}...')
                inputs = tokenizer([seqs[i]], return_tensors="pt", add_special_tokens=False,padding=True).to(device)
                outputs = model(**inputs)
                folded_positions = outputs.positions
                print(f'Finish ESM prediction {i}!')

                print(f'Write pdb output {i}...!')
                pdb = convert_outputs_to_pdb(outputs)
                print(f'./Pipeline/{fold_pair}/output_esm_fold/{msa_name}_{i}.pdb')
                print(pdb[0])
                save_string_as_pdb(pdb[0], f'./Pipeline/{fold_pair}/output_esm_fold/{msa_name}_{i}.pdb')
                print(f'Finish to write pdb output {i} !')
            except:
                continue



