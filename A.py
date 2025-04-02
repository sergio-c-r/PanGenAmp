import argparse
import os
import subprocess
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def run_blast(input_file, genus, specific_name):
    output_file = os.path.join(os.path.dirname(input_file), 'blast_seq.txt')
    blast_command = [
        'blastn', '-query', input_file, '-db', 'nt', '-out', output_file,
        '-outfmt', '6 sseqid sseq',
        '-entrez_query', f'not {genus} {specific_name}[Organism]',
        '-task', 'blastn', '-remote'
    ]
    subprocess.run(blast_command, check=True)
    return output_file

def create_multifasta(blast_output):
    fasta_output = os.path.splitext(blast_output)[0] + '.fasta'
    with open(blast_output, 'r') as blast_file, open(fasta_output, 'w') as fasta_file:
        for line in blast_file:
            sseqid, sseq = line.strip().split('\t')
            fasta_file.write(f'>{sseqid}\n{sseq}\n')
    return fasta_output

def merge_fasta_files(input_dir, merged_fasta_filename):
    fasta_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.fasta')]
    with open(merged_fasta_filename, 'w') as merged_fasta_file:
        for fasta_file in fasta_files:
            with open(fasta_file, 'r') as infile:
                shutil.copyfileobj(infile, merged_fasta_file)

def main():
    parser = argparse.ArgumentParser(description="Run BLAST and create multifasta")
    parser.add_argument('-i', '--input', required=True, help="Path to the input file")
    parser.add_argument('-g', '--genus', required=True, help="Genus name to exclude")
    parser.add_argument('-s', '--specific_name', required=True, help="Specific name to exclude")
    
    args = parser.parse_args()
    
    input_file = args.input
    genus = args.genus
    specific_name = args.specific_name
    
    blast_output = run_blast(input_file, genus, specific_name)
    
    multifasta_file = create_multifasta(blast_output)
    
    input_dir = os.path.dirname(input_file)
    merged_fasta_filename = os.path.join(input_dir, 'merged_sequences.fasta')
    merge_fasta_files(input_dir, merged_fasta_filename)

if __name__ == "__main__":
    main()

<
