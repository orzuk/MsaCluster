from pathlib import Path
import sys
import os

print("Current dir: ", os.getcwd())

sys.path.insert(0, '/sci/labs/orzuk/steveabecassis/MsaCluster')
from colabfold.batch import get_msa_and_templates
from colabfold.batch import msa_to_str
from colabfold.batch import get_queries

from argparse import ArgumentParser


# sequence = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
# result_dir = Path('./')
# jobname = 'test'
# (unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_features) = get_msa_and_templates(query_sequences=sequence,jobname=jobname,result_dir= result_dir ,use_templates=False,custom_template_path='',pair_mode='',msa_mode='mmseqs2_uniref_env')
# msa = msa_to_str(unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality)
# result_dir.joinpath(f"{jobname}.a3m").write_text(msa)
#

def main():
    parser = ArgumentParser()
    parser.add_argument("input", default="input",help="Should be a .fasta file")
    parser.add_argument("results", help="Directory to write the results to")
    parser.add_argument("-name", help="msa name")
    # parser.add_argument("--msa-mode",default="mmseqs2_uniref_env",choices=["mmseqs2_uniref_env","mmseqs2_uniref","single_sequence"],help="Choose the Db" )
    args = parser.parse_args()
    result_dir = Path(args.results)
    queries, is_complex = get_queries(args.input, 'random')
    print("queries:")
    print(queries)
    print("is_complex:")
    print(is_complex)
    sequence = queries[0][1]
    print("sequence:")
    print(sequence)
    (unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_features) = \
        get_msa_and_templates(query_sequences=sequence, jobname=args.name, result_dir= result_dir,
                              use_templates=False, custom_template_path='',pair_mode='', msa_mode='mmseqs2_uniref_env')
    print("unpaired msa:")
    print(unpaired_msa)
    print("paired msa:")
    print(paired_msa)
    print("query seqs unique:")
    print(query_seqs_unique)
    msa = msa_to_str(unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality)
    result_dir.joinpath(f"{args.name}.a3m").write_text(msa)


if __name__ == "__main__":
    main()


