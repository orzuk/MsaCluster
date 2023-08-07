from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO

# Replace 'input_sequences.fasta' with the path to your input FASTA file.
input_file = 'input_sequences.fasta'

# Replace 'output_alignment.fasta' with the desired output file name.
output_file = 'output_alignment.fasta'

# Run mafft and create the alignment.
mafft_cline = MafftCommandline(input=input_file)
stdout, stderr = mafft_cline()

# Save the alignment to the output file.
with open(output_file, "w") as handle:
    handle.write(stdout)

# Read the alignment into a Biopython object (optional, but can be useful for further analysis).
alignment = AlignIO.read(output_file, "fasta")

# Print the alignment (optional).
print(alignment)
