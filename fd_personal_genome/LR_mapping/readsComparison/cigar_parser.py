import re
import sys

def parse_cigar_file(file_path):
    print("M\tI\tD\tN\tS\tH\tP\t=\tX\ttotal\tunaligned\tcoverage")

    with open(file_path, 'r') as file:
        for line in file:
            cigar = line.strip()
            correspondance_counts = {"M": 0, "I": 0, "D": 0, "N": 0, "S": 0, "H": 0, "P": 0, "=": 0, "X": 0}
            matches = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
            for match in matches:
                count, correspondance = int(match[0]), match[1]
                correspondance_counts[correspondance] += count
            
            total = correspondance_counts["M"] + correspondance_counts["I"] + correspondance_counts["D"] + correspondance_counts["S"] + correspondance_counts["H"] + correspondance_counts["P"] 
            unaligned = correspondance_counts["S"] + correspondance_counts["H"]
            coverage = 0
            if (total > 0):
                coverage = (total - unaligned) / total

            print(f"{correspondance_counts['M']}\t{correspondance_counts['I']}\t{correspondance_counts['D']}\t{correspondance_counts['N']}\t{correspondance_counts['S']}\t{correspondance_counts['H']}\t{correspondance_counts['P']}\t{correspondance_counts['=']}\t{correspondance_counts['X']}\t{total}\t{unaligned}\t{coverage}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python cigar_parser.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    parse_cigar_file(file_path)

if __name__ == "__main__":
    main()