import os
import subprocess
import gzip
from Bio import SeqIO

tmp_folder = "tmp"
concatenated_file = "concatenated.fasta"
with open(concatenated_file, "w") as concat_file:
    for file_name in os.listdir(tmp_folder):
        if file_name.endswith(".gz"):
            with gzip.open(os.path.join(tmp_folder, file_name), "rt") as gzip_file:
                concat_file.write(gzip_file.read())

for file_name in os.listdir(tmp_folder):
    if file_name.endswith(".fna"):
        os.remove(os.path.join(tmp_folder, file_name))

db_name = "process_DB"
subprocess.run(["makeblastdb", "-in", concatenated_file, "-dbtype", "nucl", "-out", db_name])

query_file = "coregenes.fasta"
output_file = "result_output_file.txt"
subprocess.run(["blastn", "-query", query_file, "-db", db_name, "-out", output_file, "-max_target_seqs", "1", "-outfmt", "6", "-task", "megablast"])

for file_name in os.listdir():
    if "process_DB" in file_name:
        os.remove(file_name)

blast_hits = []
with open(output_file, "r") as result_file:
    for line in result_file:
        blast_hits.append(line.split("\t")[0])

updated_sequences = []
with open(query_file, "r") as coregenes_file:
    for record in SeqIO.parse(coregenes_file, "fasta"):
        if record.id not in blast_hits:
            updated_sequences.append(record)

with open(query_file, "w") as coregenes_file:
    SeqIO.write(updated_sequences, coregenes_file, "fasta")

for file_name in os.listdir(tmp_folder):
    os.remove(os.path.join(tmp_folder, file_name))
os.rmdir(tmp_folder)

if os.path.exists(output_file):
    os.remove(output_file)

if os.path.exists(concatenated_file):
    os.remove(concatenated_file)

