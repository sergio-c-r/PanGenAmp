def reformat_multifasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        sequence = ""
        header = ""
        for line in infile:
            if line.startswith(">"):
                if sequence:
                    outfile.write(header + "\n" + sequence + "\n")
                header = line.strip()
                sequence = ""
            else:
                sequence += line.strip()
        if sequence:
            outfile.write(header + "\n" + sequence + "\n")

input_file = 'aligned_all2.fas'
output_file = 'reformatted_aligned_all2.fas'
reformat_multifasta(input_file, output_file)

