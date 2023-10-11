from transformers import AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
import tmscoring
import argparse
from Bio import SeqIO

#iminuit==1.5.4
#tmscoring


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


def get_sequence_from_fasta(fasta_file_path):
    sequence = None
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        sequence = str(record.seq)
        break
    return sequence


if __name__ == '__main__':

    p = argparse.ArgumentParser(description= '')

    p.add_argument("-i", nargs='*', action='store',help='Path to msas to use in prediction.')
    p.add_argument("-o", action="store", help='name of output directory to write contact maps to.')
    # p.add_argument("--model", action='store', default='msa_t', help="Model: `esm1b` or `msa_t` (default is 'msa_t')")
    # p.add_argument('--keyword', action='store', default='', help="Keyword for this prediction")
    # p.add_argument("--test", action='store_true', help='Tests first 3 constructs.')
    # p.add_argument("--parallel", action='store_true', help='Runs in parallel using Pandarallel.')

    args = p.parse_args()

    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1")
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")

    # fasta_file_path = '/Users/steveabecassis/Desktop/FS_fasta_files/rcsb_pdb_1EBO.fasta'
    fasta_file_path = args.i
   # sequence = get_sequence_from_fasta(fasta_file_path)
    sequence = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"

    inputs = tokenizer([sequence], return_tensors="pt", add_special_tokens=False)  # A tiny random peptide
    outputs = model(**inputs)
    folded_positions = outputs.positions

    pdb = convert_outputs_to_pdb(outputs)
    print(pdb[0])

    save_string_as_pdb(pdb[0], args.o)

    # # p_5jyt = '/Users/steveabecassis/Desktop/A3M_Fasta_hurcs_files/Kaib/PDB/5jyt.pdb'
    # # p_2qke = '/Users/steveabecassis/Desktop/A3M_Fasta_hurcs_files/Kaib/PDB/2qke.pdb'
    # # tmscoring.get_tm(p_5jyt, p_2qke)
    #
    # p_1ebo = '/Users/steveabecassis/Desktop/FS_fasta_files/1ebo.pdb'
    # p_5fhc = '/Users/steveabecassis/Desktop/FS_fasta_files/5fhc.pdb'
    #
    # p_2qke_temp = tmscoring.get_tm(p_1ebo, p_5fhc)
    # p_2qke_temp = tmscoring.get_tm('/Users/steveabecassis/Desktop/pdb_test_esm.pdb', p_5fhc)
    #
    # print(float(p_2qke_temp))
