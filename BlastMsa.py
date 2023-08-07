from Bio import Entrez, SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from biotite.sequence.align import align_multiple
from Bio.Blast import NCBIWWW, NCBIXML


def fetch_sequence(uniprot_id):
    Entrez.email = "steve.abecassis@mail.huji.ac.il"  # Replace with your email
    handle = Entrez.efetch(db="protein", id=uniprot_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record.seq

def perform_blast_search(protein_sequence):
    result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence,hitlist_size=2500)
    blast_records = NCBIXML.parse(result_handle)
    homologous_sequences = []

    for record in blast_records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                homologous_sequences.append(hsp.sbjct)

    return homologous_sequences

def save_sequences_to_fasta(protein_sequence, homologous_sequences):
    with open("input_sequences.fasta", "w") as f:
        f.write(">Query\n" + str(protein_sequence) + "\n")
        for i, seq in enumerate(homologous_sequences):
            f.write(f">Homolog_{i + 1}\n{seq}\n")


def perform_msa(sequences):
    output_alignment_file = "output_alignment.fasta"
    clustalo_cline = ClustalOmegaCommandline(infile="input_sequences.fasta", outfile=output_alignment_file, verbose=True, auto=True)

    stdout, stderr = clustalo_cline()

    from biotite.sequence.align import align_multiple

    align_multiple(sequences, penalty=-10, terminal_penalty=-10, gap_penalty=-1)

    return output_alignment_file



if __name__ == "__main__":
    # main()
    protein_sequence = 'MKQIEDKIEEILSKIYHIENEIARIKKLIGEADGLIEGLRQLANETTQALQLFLRATTELRTFSILNRKAIDFLLQRWGGTCHILGPDCRIEPHDWTKNITDKIDQIIHDFVDKTLPDQGDNDNWWTGWRQ'

    # Perform a BLAST search to find homologous sequences
    homologous_sequences = perform_blast_search(protein_sequence)

    # Save sequences to FASTA file
    save_sequences_to_fasta(protein_sequence,homologous_sequences)

    # Perform MSA using Clustal Omega
    output_alignment_file = perform_msa(homologous_sequences)

    from biotite.sequence import ProteinSequence
    from biotite.sequence.align import MultipleSequenceAligner
    from biotite.sequence.io.fasta import FastaFile

    align_multiple(sequences =homologous_sequences, penalty=-10, terminal_penalty=-10, gap_penalty=-1)

    from biotite.sequence import ProteinSequence
    from biotite.sequence.align import SubstitutionMatrix

    homologous_sequences_biotite = [ProteinSequence(i.replace('-','')) for i in homologous_sequences][:50]
    matrix = SubstitutionMatrix.std_protein_matrix()
    alignment, order, tree, distances = align_multiple(homologous_sequences_biotite, matrix)

    print(tree.to_newick(labels=[str(i) for i in range(len(homologous_sequences_biotite))], include_distance=False))
    print(alignment)
    len(alignment)
    print(alignment)





    with FastaFile("output_alignment.fasta", "w") as fasta:
        for seq in alignment:
            fasta.write(seq)
    # Print the output alignment file
    with open(output_alignment_file, "r") as f:
        print(f.read())

    print('Finish !')

