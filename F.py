import pandas as pd

input_file = "coregenes_blastn_results.txt"
output_file = "coregenes_blastn_resultsF1.txt"

df = pd.read_csv(input_file, sep='\t', header=None)

df_unique = df.drop_duplicates(subset=0, keep=False)

df_unique.to_csv(output_file, sep='\t', header=False, index=False)

print(f"Filtered data saved to {output_file}")

