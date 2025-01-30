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
