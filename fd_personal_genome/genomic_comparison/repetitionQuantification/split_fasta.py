import sys

def split_fasta(input_file):
    with open(input_file, 'r') as f:
        seqs = {}
        current_seq = None
        for line in f:
            if line.startswith('>'):
                seq_name = line[1:].split()[0]  # keep only what's before the first space
                seq_name = seq_name.replace("#", "_").replace(":", "_").replace("-", "_")  # replace special chars
                current_seq = seq_name
                seqs[current_seq] = []
            else:
                seqs[current_seq].append(line.strip())
    
    for seq_name, seq_data in seqs.items():
        with open(f"{seq_name}.fasta", 'w') as outf:
            outf.write(f">{seq_name}\n")
            outf.write("\n".join(seq_data))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python split_fasta.py <input_file>")
        sys.exit(1)
    
    split_fasta(sys.argv[1])