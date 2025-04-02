input_file_name = "core_genes"
output_file_name = "list.txt"

try:
    with open(input_file_name, 'r') as input_file:
        second_words = []

        for line in input_file:
            words = line.split()
            
            if len(words) > 1:
                # Append the second word to the list
                second_words.append(words[1])

    with open(output_file_name, 'w') as output_file:
        # Write the second words to the output file
        output_file.write('\n'.join(second_words))

    print(f"Codes have been extracted to '{output_file_name}'.")
except FileNotFoundError:
    print(f"File '{input_file_name}' not found.")
except Exception as e:
    print(f"An error occurred: {str(e)}")
