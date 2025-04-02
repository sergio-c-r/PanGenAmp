import os
from Bio import SeqIO

folder_path = "genomes"
output_file = "coregenes.fasta"
list_file = "list.txt"
sequence_identifiers = set()

with open(list_file, "r") as list_file:
    for line in list_file:
        sequence_identifiers.add(line.strip())

gene_sequences = []

for filename in os.listdir(folder_path):
    if filename.endswith(".ffn"):
        fasta_file_path = os.path.join(folder_path, filename)

        for record in SeqIO.parse(fasta_file_path, "fasta"):
            # Check if the record's identifier is in the list
            if record.id in sequence_identifiers:
                gene_sequences.append(record)

with open(output_file, "w") as output:
    SeqIO.write(gene_sequences, output, "fasta")

print(f"Gene sequences listed in '{list_file}' have been extracted and saved to '{output_file}'.")

