from Bio import SeqIO


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


