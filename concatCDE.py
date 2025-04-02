import subprocess
import argparse

def count_fasta_symbols(file_path):
    count = 0
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(">"):
                count += 1
    return count

def run_script(script, keyword, other):
    subprocess.run(["python", script, "-g", keyword, "-s", other])

def main():
    parser = argparse.ArgumentParser(description='Concatenate and execute sub-scripts.')
    parser.add_argument('-g', '--keyword', required=True, help='Genus for C.py and E.py i.e. Vibrio')
    parser.add_argument('-s', '--other', required=True, help='Species option for C.py and E.py i.e. cholerae')
    args = parser.parse_args()

    while True:
        run_script("C.py", args.keyword, args.other)
        run_script("D.py", args.keyword, args.other)

        symbols_count = count_fasta_symbols("coregenes.fasta")
        print("Number of genes in coregenes.fasta:", symbols_count)

        if symbols_count > 150:
            print("Count > 150. Restarting C.py.")
        else:
            run_script("E.py", args.keyword, args.other)
            break

if __name__ == "__main__":
    main()

