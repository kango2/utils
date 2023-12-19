import argparse
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
import sys

def remove_duplicate_sequences(input_file, output_file):
    """
    Removes duplicate sequences from a FASTA file, retaining the first occurrence.

    Args:
    input_file (str): Path to the input FASTA file.
    output_file (str): Path to the output FASTA file.
    """
    seen_sequences = {}
    unique_sequences = []
    total_input_sequences = 0

    with open(input_file, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            total_input_sequences += 1
            sequence_str = str(record.seq)
            if sequence_str not in seen_sequences:
                seen_sequences[sequence_str] = [record.id]
                unique_sequences.append(record)
            else:
                seen_sequences[sequence_str].append(record.id)

    with open(output_file, "w") as outfile:
        writer = FastaWriter(outfile, wrap=70)
        writer.write_file(unique_sequences)

    # Print information about retained and removed sequences
    print(f"Total sequences in input file: {total_input_sequences}", file=sys.stderr)
    print(f"Total unique sequences in output file: {len(unique_sequences)}", file=sys.stderr)
    print("retained\tremoved", file=sys.stderr)
    for seq, ids in seen_sequences.items():
        if len(ids) > 1:
            print(f"{ids[0]}\t{','.join(ids[1:])}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description="Remove duplicate sequences from a FASTA file")
    parser.add_argument("--input_file", type=str, required=True, help="Path to the input FASTA file")
    parser.add_argument("--output_file", type=str, required=True, help="Path to the output FASTA file")

    args = parser.parse_args()

    remove_duplicate_sequences(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
