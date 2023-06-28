import os
import requests as r
from Bio import SeqIO
from io import StringIO
from pathlib import Path
from colabfold.batch import get_msa_and_templates
from colabfold.batch import msa_to_str
from colabfold.batch import get_queries
from tqdm import tqdm
from argparse import ArgumentParser

def get_seq_from_uniprotId(cID):
    # cID='A0A3N5IQ47'
    baseUrl="http://www.uniprot.org/uniprot/"
    currentUrl=baseUrl+cID+".fasta"
    response = r.post(currentUrl)
    cData=''.join(response.text)
    Seq=StringIO(cData)
    pSeq=list(SeqIO.parse(Seq,'fasta'))
    return str(pSeq[0].seq)


def get_seqids(msa_file_path):
    with open(msa_file_path, 'r') as msa_fil:
        seqids = msa_fil.read().splitlines()#[idx].strip()[1:]
        # sequence = get_seq_from_uniprotId(seqId)
        ids = []
        for seqid in seqids:
            if '>' in seqid:
                ids.append(seqid)
        ids = ids[1:]
        return ids

def main():
    parser = ArgumentParser()
    parser.add_argument("input", default="input", help="Should be a folder of a3m files")
    parser.add_argument("output", default="output", help="Output directory")
    args = parser.parse_args()
    # msas_folder = '/Users/steveabecassis/Desktop/test_sample_msa'
    msas_folder   = args.input
    output_folder = args.output
    for msa in tqdm(os.listdir(msas_folder)):
        if ('a3m' not in msa):
            continue
        os.makedirs(f'{output_folder}/{msa[:-4]}')
        seqids = get_seqids(f'{msas_folder}/{msa}')
        count = 2
        for sample_cnt in range(count):
            for id in seqids:
                try:
                    seq = get_seq_from_uniprotId(id[1:])
                    if (len(seq) > 0):
                        seq.remove(id)
                        break
                    else:
                        continue
                except:
                    continue
            with open(f'{output_folder}/{msa[:-4]}/{id[1:]}.fasta', 'w') as f_out:
                f_out.write('>' + id + '\n' + seq)

            result_dir = Path(f'{output_folder}/{msa[:-4]}')
            queries, is_complex = get_queries(f'{output_folder}/{msa[:-4]}/{id[1:]}.fasta', 'random')
            sequence = queries[0][1]
            (unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality,
             template_features) = get_msa_and_templates(query_sequences=sequence, jobname=msa[:-4],
                                                        result_dir=result_dir, use_templates=False,
                                                        custom_template_path='', pair_mode='',
                                                        msa_mode='mmseqs2_uniref_env')
            msa_sample = msa_to_str(unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality)
            result_dir.joinpath(f"./{msa[:-4]}_FULL_{str(sample_cnt)}.a3m").write_text(msa_sample)

            with open(f'{result_dir}/{msa[:-4]}_FULL_{str(sample_cnt)}.a3m', 'r') as msa_fil:
                msa_split = msa_fil.read().splitlines()
                algn_msa_list = []
                for i in seqids:
                    for idx in range(len(msa_split)):
                        if i in msa_split[idx]:
                            algn_msa_list.append(msa_split[idx])
                            algn_msa_list.append(msa_split[idx + 1])
            if (len(algn_msa_list) == 0):
                count = count - 1
                continue
            df = open(f'{output_folder}/{msa[:-4]}/{msa[:-4]}_SAMPLE_{str(sample_cnt)}.a3m', 'w')
            flag = 1
            for i in algn_msa_list:
                if flag == 1:
                    df.write("%s\n" % '>101')
                    df.write("%s\n" % seq[:108])
                    flag = 0
                df.write("%s\n" % i)
                df.write('\n')





if __name__ == '__main__':
    main()


















