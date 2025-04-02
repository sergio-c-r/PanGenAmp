import os
import shutil
import subprocess
from time import sleep
from Bio import SeqIO

def split_fasta(input_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    with open(input_file, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        num_records = len(records)
        for i in range(0, num_records, 1):
            sub_records = records[i:i+1]
            output_file = os.path.join(output_dir, f"{i//1 + 1}.fasta")
            with open(output_file, "w") as out_handle:
                SeqIO.write(sub_records, out_handle, "fasta")

def run_blastn(input_file, output_file):
    command = f'blastn -query {input_file} -db nt -out {output_file} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -entrez_query "not Vibrio europaeus" -max_target_seqs 1 -task blastn -remote'
    subprocess.run(command, shell=True)

folder_name = "tmp2"
os.makedirs(folder_name, exist_ok=True)

shutil.copy("coregenes.fasta", folder_name)

split_fasta("coregenes.fasta", folder_name)

output_dir = os.path.join(folder_name, "blast_results")
os.makedirs(output_dir, exist_ok=True)

for fasta_file in os.listdir(folder_name):
    if fasta_file.endswith(".fasta"):
        fasta_path = os.path.join(folder_name, fasta_file)
        output_name = os.path.splitext(fasta_file)[0] + "_blastn_results.txt"
        output_path = os.path.join(output_dir, output_name)
        run_blastn(fasta_path, output_path)
        sleep(15)  

output_concatenated = os.path.join(folder_name, "out_blast_rm.txt")
with open(output_concatenated, 'w') as outfile:
    for fname in os.listdir(output_dir):
        if fname.endswith('.txt'):
            with open(os.path.join(output_dir, fname), 'r') as infile:
                outfile.write(infile.read())

shutil.move(os.path.join(output_dir, "coregenes_blastn_results.txt"), ".")

shutil.rmtree(folder_name)

print("Process completed successfully.")

