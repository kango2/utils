import argparse
from Bio import SeqIO

def filter_fasta_by_length(input_file, output_file, min_length=None, max_length=1000000):
    """
    Filters sequences in a FASTA file based on length.

    Args:
    input_file (str): Path to the input FASTA file.
    output_file (str): Path to the output FASTA file.
    min_length (int, optional): Minimum length of sequences to keep. If None, no lower limit is applied.
    max_length (int): Maximum length of sequences to keep.
    """
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if (min_length is None or len(record.seq) >= min_length) and len(record.seq) <= max_length:
                SeqIO.write(record, outfile, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Filter FASTA sequences by length")
    parser.add_argument("--input_file", type=str, required=True, help="Path to the input FASTA file")
    parser.add_argument("--output_file", type=str, required=True, help="Path to the output FASTA file")
    parser.add_argument("--min_length", type=int, help="Minimum length of sequences to keep. Optional.")
    parser.add_argument("--max_length", type=int, default=1000000, help="Maximum length of sequences to keep (default: 1000000)")

    args = parser.parse_args()

    filter_fasta_by_length(args.input_file, args.output_file, args.min_length, args.max_length)

if __name__ == "__main__":
    main()
