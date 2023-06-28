import os
import requests as r
from Bio import SeqIO
from io import StringIO
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

def get_seqId_indexes(a3m_file):
    seq_indexes = []
    seq_dict = []
    cnt_succ = 0
    seq_cnt = len(open(a3m_file).read().splitlines())
    for idx in range(seq_cnt):
        try:
            with open(a3m_file, 'r') as msa_fil:
                seqId = msa_fil.read().splitlines()[idx].strip()[1:]
                sequence = get_seq_from_uniprotId(seqId)
            if (len(sequence)>1):
                cnt_succ += 1
                seq_indexes.append(int(idx))
                seq_dict.append({'seqId': seqId,'sequence':sequence})
            if(cnt_succ == 2):
                return seq_indexes,seq_dict
            else:
                continue
        except:
            continue


def main():
    parser = ArgumentParser()
    parser.add_argument("input",default="input",help="Should be a folder of msa files" )
    parser.add_argument("results", help="Directory to write the fasta files to")
    parser.add_argument("name", help="Fasta files prefix name")
    args = parser.parse_args()
    msa_folder = os.listdir(args.input)
    for msa in msa_folder:
        seq_indexes,seq_dict = get_seqId_indexes(msa)
        fasta_file = open(f'{str(msa)}_sample_0.fasta','fasta')
        fasta_file.write('>'+seq_dict[0]['seqId'] + '/n' + seq_dict[0]['sequence'])

        fasta_file = open(f'{str(msa)}_sample_1.fasta', 'fasta')
        fasta_file.write('>' + seq_dict[1]['seqId'] + '/n' + seq_dict[1]['sequence'])

if __name__ == '__main__':
    main()







