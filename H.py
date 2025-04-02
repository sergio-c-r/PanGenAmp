import os
import shutil
import pandas as pd
from Bio import SeqIO

os.makedirs('tmp3', exist_ok=True)
shutil.copy('candidates.txt', 'tmp3/candidates.txt')
shutil.copy('core_genes', 'tmp3/core_genes')

with open('tmp3/core_genes', 'r') as file:
    core_genes_content = file.read()

core_genes_content = core_genes_content.replace(' ', '\t')

with open('tmp3/core_genes', 'w') as file:
    file.write(core_genes_content)

candidates_df = pd.read_csv('tmp3/candidates.txt', sep='\t', header=None)
index_values = candidates_df.iloc[:, 0]
index_values.to_csv('tmp3/index.txt', index=False, header=False)

core_genes_df = pd.read_csv('tmp3/core_genes', sep='\t', header=None)

selected_genes_df = core_genes_df[core_genes_df.iloc[:, 1].isin(index_values)]
selected_genes_df.to_csv('tmp3/selected_gene_fam.txt', sep='\t', index=False, header=False)

with open('tmp3/selected_gene_fam.txt', 'r') as selected_file:
    for line in selected_file:
        gene_name, *gene_values = line.strip().split('\t')
        
                folder_name = os.path.join('tmp3', gene_name)
        os.makedirs(folder_name, exist_ok=True)
        
               txt_file_path = os.path.join(folder_name, gene_name + '.txt')
        with open(txt_file_path, 'w') as txt_file:
            for value in gene_values:
                txt_file.write(value + '\n')

with open('genes.ffn', 'wb') as outfile:
    for filename in os.listdir('./genomes'):
        if filename.endswith('.ffn'):
            with open(os.path.join('./genomes', filename), 'rb') as readfile:
                shutil.copyfileobj(readfile, outfile)
shutil.move('genes.ffn', 'tmp3/genes.ffn')

for folder_name in os.listdir('tmp3'):
    folder_path = os.path.join('tmp3', folder_name)
    if os.path.isdir(folder_path):
        txt_file_path = os.path.join(folder_path, f'{folder_name}.txt')
        with open(txt_file_path, 'r') as txt_file:
            gene_values = [line.strip() for line in txt_file]

        with open('tmp3/genes.ffn', 'r') as ffn_file:
            sequences = SeqIO.parse(ffn_file, 'fasta')
            for index, seq_record in enumerate(sequences, start=1):
                if any(gene_value in seq_record.id for gene_value in gene_values):
                    output_file_path = os.path.join(folder_path, f'{folder_name}_s{index}.fasta')
                    with open(output_file_path, 'w') as output_file:
                        SeqIO.write(seq_record, output_file, 'fasta')

        fasta_files_list = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.fasta')]
        fasta_files = " ".join(fasta_files_list)

        sketch_command = f"bindash sketch --outfname={os.path.join(folder_path, 'genes.sketch')} {fasta_files}"
        os.system(sketch_command)

        dist_command = f"bindash dist {os.path.join(folder_path, 'genes.sketch')} {os.path.join(folder_path, 'genes.sketch')} > {os.path.join(folder_path, f'{folder_name}.dist.txt')}"
        os.system(dist_command)

summary_lines = []

for folder_name in os.listdir('tmp3'):
    folder_path = os.path.join('tmp3', folder_name)
    if os.path.isdir(folder_path):
        dist_file_path = os.path.join(folder_path, f'{folder_name}.dist.txt')
        if os.path.exists(dist_file_path):
            with open(dist_file_path, 'r') as dist_file:
                min_fraction = None
                for line in dist_file:
                    parts = line.strip().split('\t')
                    if len(parts) == 5:
                        fraction = parts[4]
                        num, denom = map(int, fraction.split('/'))
                        current_fraction = num / denom
                        if min_fraction is None or current_fraction < min_fraction:
                            min_fraction = current_fraction
                if min_fraction is not None:
                    summary_lines.append(f"{folder_name}\t{min_fraction}")

summary_file_path = 'IS_JI_summary.txt'
with open(summary_file_path, 'w') as summary_file:
    summary_file.write('\n'.join(summary_lines))

os.makedirs('alignment', exist_ok=True)
for folder_name in os.listdir('tmp3'):
    shutil.move(os.path.join('tmp3', folder_name), os.path.join('alignment', folder_name))

for folder_name in os.listdir('alignment'):
    folder_path = os.path.join('alignment', folder_name)
    if os.path.isdir(folder_path):
        for file_name in os.listdir(folder_path):
            file_path = os.path.join(folder_path, file_name)
            if not file_name.endswith('.fasta'):
                os.remove(file_path)

for item_name in os.listdir('alignment'):
    item_path = os.path.join('alignment', item_name)
    if os.path.isfile(item_path):
        os.remove(item_path)

shutil.rmtree('tmp3')

