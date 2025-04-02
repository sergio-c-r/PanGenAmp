import argparse
import pandas as pd
import re
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging
import sys

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

parser = argparse.ArgumentParser(description="Filter rows based on percentage identity (h), maximum coverage (c), and minimum length (l).")
parser.add_argument("-p", "--percentage_identity", type=float, default=100, help="Maximum percentage identity value (float, default: 100)")
parser.add_argument("-c", "--maximum_coverage", type=float, default=75, help="Maximum coverage value (float, default: 75)")
parser.add_argument("-l", "--minimum_length", type=int, default=300, help="Minimum length allowed (integer, default: 300)")
args = parser.parse_args()

input_file = "coregenes_blastn_resultsF1.txt"
try:
    df = pd.read_csv(input_file, sep='\t', header=None)
except Exception as e:
    logging.error(f"Failed to read input file {input_file}: {e}")
    sys.exit(1)

df[12] = df[12].apply(lambda x: float(re.search(r'(\d+(\.\d+)?)', str(x)).group(1)) if re.search(r'(\d+(\.\d+)?)', str(x)) else 0)

filtered_df = df[(df[2] <= args.percentage_identity) & (df[12] <= args.maximum_coverage) & (df[3] >= args.minimum_length)]

output_filename = "candidates.txt"
try:
    filtered_df.to_csv(output_filename, sep='\t', header=False, index=False)
    logging.info(f"Filtered data saved to {output_filename}")
except Exception as e:
    logging.error(f"Failed to write filtered data to {output_filename}: {e}")
    sys.exit(1)

output_seq_filename = "candidates_seq.fasta"
filtered_genes = set(filtered_df[0])
sequences = []

try:
    with open("coregenes.fasta", "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            gene_id = record.id.split()[0]  
            if gene_id in filtered_genes:
                sequences.append(record)
except Exception as e:
    logging.error(f"Failed to process coregenes.fasta: {e}")
    sys.exit(1)

try:
    SeqIO.write(sequences, output_seq_filename, "fasta")
    logging.info(f"Matching sequences saved to {output_seq_filename}")
except Exception as e:
    logging.error(f"Failed to write sequences to {output_seq_filename}: {e}")
    sys.exit(1)

output_protseq_filename = "candidates_protseq.fasta"
translated_records = []

for seq_record in sequences:
    gene_id = seq_record.id.split()[0]
    translated_sequence = seq_record.seq.translate(table=1, to_stop=True, cds=False, gap='X')
    translated_header = f"{gene_id}_protein"
    translated_record = SeqRecord(translated_sequence, id=translated_header, description="")
    translated_records.append(translated_record)

try:
    SeqIO.write(translated_records, output_protseq_filename, "fasta")
    logging.info(f"Translated protein sequences saved to {output_protseq_filename}")
except Exception as e:
    logging.error(f"Failed to write translated protein sequences to {output_protseq_filename}: {e}")
    sys.exit(1)

blastp_command = f"blastp -query {output_protseq_filename} -db refseq_select -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -remote -out out_RefSeq -num_alignments 25"
try:
    subprocess.run(blastp_command, shell=True, check=True)
    logging.info("BLASTP search completed.")
except subprocess.CalledProcessError as e:
    logging.error(f"BLASTP search failed: {e}")
    sys.exit(1)

