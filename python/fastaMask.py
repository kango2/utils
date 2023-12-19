from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.SeqIO.FastaIO import FastaWriter
import argparse
import hashlib
import re

def calculate_md5(seq):
    """
    Calculate the MD5 checksum of a sequence.

    Args:
    seq (str): DNA sequence

    Returns:
    str: MD5 checksum of the sequence
    """
    md5hash = hashlib.md5(seq.encode('utf-8')).hexdigest()
    return md5hash

def update_md5_in_description(description, md5_checksum):
    """
    Update the MD5 checksum in the sequence description.

    Args:
    description (str): Original sequence description
    md5_checksum (str): MD5 checksum to be included or updated in the description

    Returns:
    str: Updated sequence description
    """
    if 'M5:' in description:
        return re.sub(r'M5:[a-fA-F0-9]{32}', f'M5:{md5_checksum}', description)
    else:
        return f"{description} M5:{md5_checksum}"

def mask_regions(fasta_file, bed_file, output_file):
    """
    Masks regions in a FASTA file specified by a BED file and updates or adds MD5 checksum to descriptions.

    Args:
    fasta_file (str): Path to the input FASTA file.
    bed_file (str): Path to the BED file with regions to mask.
    output_file (str): Path to the output FASTA file with masked sequences.
    """
    # Load sequences from the FASTA file
    sequences = list(SeqIO.parse(fasta_file, "fasta"))

    # Read and process the BED file
    with open(bed_file, 'r') as bed:
        for line in bed:
            chrom, start, end = line.strip().split()[:3]
            start, end = int(start), int(end)
            for seq in sequences:
                if seq.id == chrom:
                    mutable_seq = MutableSeq(str(seq.seq))
                    mutable_seq[start:end] = 'N' * (end - start)
                    seq.seq = mutable_seq

    # Write masked and annotated sequences to the output file with 70 characters per line
    with open(output_file, "w") as output:
        writer = FastaWriter(output, wrap=70)
        writer.write_header()
        for seq in sequences:
            # Calculate MD5 checksum and update the description
            md5_checksum = calculate_md5(str(seq.seq))
            seq.description = update_md5_in_description(seq.description, md5_checksum)
            writer.write_record(seq)
        writer.write_footer()

def main():
    parser = argparse.ArgumentParser(description="Mask regions in a FASTA file using a BED file and update/add MD5 checksum to descriptions")
    parser.add_argument("--fasta_file", type=str, required=True, help="Path to the input FASTA file")
    parser.add_argument("--bed_file", type=str, required=True, help="Path to the BED file with regions to mask")
    parser.add_argument("--output_file", type=str, required=True, help="Path to the output FASTA file")

    args = parser.parse_args()
    mask_regions(args.fasta_file, args.bed_file, args.output_file)

if __name__ == "__main__":
    main()
