# PanGenAmp
A series of Python scripts to select species-specific candidate genes for primer design, based on Roary pangenome analysis.  

## Input Structure

-> **Roary output files**, including `core_genes`.
-> A folder named `genomes/` containing `.ffn` files from Prokka.

---

## Pipeline Overview

### `A.py`
-> Extracts gene codes from each gene family listed in `core_genes`.
-> Outputs a list of gene identifiers.

### `B.py`
-> Iterates through `.ffn` files in the `genomes/` folder.
-> Extracts matching core gene sequences.

### `C.py`
-> Indexes available species from [NCBI RefSeq Bacteria](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/).
-> Downloads random genomes from the same genus (but not from the same species) using options -g <Genus> -s <Species> (i.e. -g Vibrio -s europaeus).
  
### `D.py`
-> Concatenates downloaded genomes.
-> Builds a BLAST database.
-> Runs BLAST using core sequences from B.py.
-> Removes matching sequences from the core dataset.

### `E.py`
-> Performs remote BLAST of the remaining core sequences against the NCBI nt database (excluding the input species).

### `concatCDE.py`
-> Concatenate C.py and D.py
-> Runs E.py if fewer than 150 core sequences remain.

### `F.py`
-> Remove duplicated results of the online BLAST.

### `G.py`
-> Filters BLASTN results by:
- Percentage identity (-p, default: 100).
- Coverage (-c, default: 75).
- Minimum sequence length (-l, default: 300 bp).
-> Generates:
- candidates.txt: Tab-separated BLASTN results.
- candidates_seq.fasta: DNA sequences of candidates.
- candidates_protseq.fasta: Protein sequences of candidates.
- out_RefSeq: Remote BLASTP results for protein candidates.

### `H.py`
-> Uses Bindash to calculate Jaccard Index of candidate sequences
-> Output: IS_JI_summary.txt

### `I.py`
-> For each candidate, creates a file combining the sequences in the pangenome and the matches in the remote blast.

### `J.py`
-> Reformats multi-FASTA files for better compatibility with sequence alignment tools.
