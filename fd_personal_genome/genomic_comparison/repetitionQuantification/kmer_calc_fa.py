import sys

def compute_kmer_frequencies(fasta_file, k):
    """
    Compute kmer frequencies of a sequence in a FASTA file.

    Parameters:
    fasta_file (str): Path to the FASTA file.
    k (int): Size of the kmer.

    Returns:
    dict: Dictionary with kmer frequencies.
    """
    kmer_freqs = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq = line.strip()
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                if kmer in kmer_freqs:
                    kmer_freqs[kmer] += 1
                else:
                    kmer_freqs[kmer] = 1
    return kmer_freqs

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py <fasta_file> <k>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    k = int(sys.argv[2])

    kmer_freqs = compute_kmer_frequencies(fasta_file, k)
    for kmer, freq in kmer_freqs.items():
        print(f"{kmer}\t{freq}")