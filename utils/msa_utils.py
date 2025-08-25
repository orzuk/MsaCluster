from Bio import SeqIO
from Bio import pairwise2


def lprint(string, f):
    print(string)
    f.write(string + '\n')


def load_fasta(fil):
    seqs, IDs = [], []
    with open(fil) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = ''.join([x for x in record.seq])
            IDs.append(record.id)
            seqs.append(seq)
    return IDs, seqs


def write_fasta(names, seqs, outfile='tmp.fasta'):
    with open(outfile, 'w') as f:
        for nm, seq in list(zip(names, seqs)):
            f.write(">%s\n%s\n" % (nm, seq))



def get_align_indexes(seqA, seqB):
    alignments = pairwise2.align.globalxx(seqA, seqB, one_alignment_only=True)
    best_align = alignments[0]
    seqA = best_align.seqA
    seqB = best_align.seqB
    cursA = 0
    cursB = 0
    seqA_idxs = []
    seqB_idxs = []
    for aa in range(len(seqA)):
        if (seqA[aa] != '-') & (seqB[aa] != '-'):
            seqA_idxs.append(cursA)
            seqB_idxs.append(cursB)
            cursA += 1
            cursB += 1
        if (seqA[aa] == '-') & (seqB[aa] != '-'):
            cursB += 1
        if (seqA[aa] != '-') & (seqB[aa] == '-'):
            cursA += 1
    return seqA_idxs, seqB_idxs


import requests

import requests
import gzip
import shutil
import os
from Bio import AlignIO
import requests
from io import StringIO

def download_and_parse_pfam_msa(pfam_id, alignment_type="seed"):
    """
    Downloads and parses a Pfam MSA using the REST API (SEED only).

    Parameters:
        pfam_id (str): Pfam family ID (e.g., 'PF00069')
        alignment_type (str): Only 'seed' is supported currently

    Returns:
        MultipleSeqAlignment: Parsed MSA object from Biopython
    """
    if alignment_type != "seed":
        raise NotImplementedError("Only 'seed' alignment is available via the Pfam REST API")

    url = f"https://www.ebi.ac.uk/interpro/api/pfam/entry/pfam/{pfam_id}/alignment/seed/stockholm/"

    response = requests.get(url)
    if response.status_code != 200:
        raise Exception(f"Download failed with status code {response.status_code}: {url}")

    sto_text = response.text
    msa_io = StringIO(sto_text)
    alignment = AlignIO.read(msa_io, "stockholm")

    return alignment
