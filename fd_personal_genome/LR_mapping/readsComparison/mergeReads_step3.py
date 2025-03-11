import argparse
import gzip
import pandas as pd



def merge_files(file1, file2, output_file):
    with open(file2, 'rt') as f2:
        lines = [line.strip().split('\t') for line in f2]
    file2_data = {line[0]: line[1:] for line in lines}

    # Read file1 in chunks to avoid memory issues
    chunksize = 10 ** 6
    chunks1 = pd.read_csv(file1, sep='\t', header=None, chunksize=chunksize)

    # Merge data from file1 with file2
    with open(output_file, 'w') as outf:
        for chunk in chunks1:
            merged_chunk = chunk.apply(lambda row: pd.concat([row, pd.Series(file2_data.get(row[0], []))], ignore_index=True), axis=1)
            #merged_chunk = merged_chunk.replace("", pd.NA)
            outf.write(merged_chunk.to_csv(index=False, sep='\t', header=False, na_rep="NA"))


#output_file = 'test.csv'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge two tab-separated value files based on the first column.')
    parser.add_argument('file1', type=str, help='Path to the first file')
    parser.add_argument('file2', type=str, help='Path to the second file')
    parser.add_argument('output_file', type=str, help='Path to the output file')
    args = parser.parse_args()
    merge_files(args.file1, args.file2, args.output_file)