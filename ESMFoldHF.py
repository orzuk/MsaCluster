from transformers import AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein



model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1")
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
inputs = tokenizer(["MLKNVQVQLV"], return_tensors="pt", add_special_tokens=False)  # A tiny random peptide
outputs = model(**inputs)
folded_positions = outputs.positions

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

print('Finish ESM prediction !')

pdb = convert_outputs_to_pdb(outputs)
print(pdb[0])
save_string_as_pdb(pdb[0], './TEST_ESM.pdb')

print('FINISH!!!!')
