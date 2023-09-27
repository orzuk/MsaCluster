import os

from transformers import AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from argparse import  ArgumentParser
import torch
from random import sample
import random
random.seed(10)
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


def save_string_as_pdb(pdb_string, file_path):
    with open(file_path, 'w') as pdb_file:
        pdb_file.write(pdb_string)



if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument("-input", default="input", help="Should be a .a3m file")
    parser.add_argument("-output", help="Directory to write the results to")
    parser.add_argument("-name", help="msa name")

    args = parser.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"

    # class args:
    #     input = './2QKEE_002.a3m'
    #     output = './'
    #     name = 'test'
    input_path = './Pipeline/output/output_msa_cluster'
    msas_files = os.listdir(args.input)
    print('Load model...!')
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    print('Finish to load model !')

    for msa in msas_files:
        with open(f'{input_path}/{msa}', 'r') as msa_fil:
            seq = msa_fil.read().splitlines()
        msa_name = msa[:-4]
        seqs = [i.replace('-', '') for i in seq if '>' not in i]
        if len(seqs) > 30:
            seqs = sample(seqs,30)

        # model.esm = model.esm.half()
        model.trunk.set_chunk_size(64)
        model.esm.float()

        for i in range(len(seqs)):
            try:
                print(f'Get ESM prediction {i}...')
                inputs = tokenizer([seqs[i]], return_tensors="pt", add_special_tokens=False,padding=True)
                outputs = model(**inputs)
                folded_positions = outputs.positions
                print(f'Finish ESM prediction {i}!')

                print(f'Write pdb output {i}...!')
                pdb = convert_outputs_to_pdb(outputs)
                print(pdb[0])
                save_string_as_pdb(pdb[0], f'./Pipeline/output/esm_fold_output/{msa_name}_{i}.pdb')
                print(f'Finish to write pdb output {i} !')
            except:
                continue