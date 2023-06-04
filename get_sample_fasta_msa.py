import os
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

if __name__ == '__main__':

    # Create fasta seq
    msa_folder = os.listdir()
    for msa in msa_folder:
        seq_indexes,seq_dict = get_seqId_indexes(msa)
        fasta_file = open('fasta_file.fasta','fasta')
        fasta_file.write('>'+seq_dict[0]['seqId'] + '/n' + seq_dict[0]['sequence'])








