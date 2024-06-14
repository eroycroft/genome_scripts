from Bio import SeqIO
import numpy as np
import argparse

# Function to calculate heterozygosity in a given sequence
def calculate_heterozygosity(sequence):
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    total_bases = 0
    for base in sequence:
        if base in counts:
            counts[base] += 1
            total_bases += 1
    heterozygous_sites = sum(count for count in counts.values() if count > 0) - max(counts.values())
    return heterozygous_sites / total_bases if total_bases > 0 else 0

# Parse command line arguments
parser = argparse.ArgumentParser(description='Calculate heterozygosity in non-overlapping windows of a genome.')
parser.add_argument('genome_file', metavar='genome_file', type=str, help='Path to the genome file in FASTA format (reference genome)')
parser.add_argument('-w', '--window_size', metavar='window_size', type=int, default=1000000, help='Size of the window in base pairs (default: 1 Mb)')
parser.add_argument('-o', '--output_file', metavar='output_file', type=str, default='heterozygosity_results.csv', help='Output file name (default: heterozygosity_results.csv)')
args = parser.parse_args()

# Read the reference genome from the FASTA file
genome_sequences = SeqIO.parse(args.genome_file, "fasta")

# Define window size
window_size = args.window_size

# Initialize lists to store heterozygosity for each sequence and corresponding header names
heterozygosity_per_sequence = []
header_names = []
bin_numbers = []

# Iterate over each sequence in the reference genome FASTA file
for record in genome_sequences:
    header_names.append(record.id)
    genome_sequence = record.seq
    # Calculate heterozygosity in non-overlapping windows for this sequence
    num_windows = len(genome_sequence) // window_size
    for i in range(num_windows):
        start = i * window_size
        end = start + window_size
        window_sequence = genome_sequence[start:end]
        heterozygosity = calculate_heterozygosity(window_sequence)
        heterozygosity_per_sequence.append(heterozygosity)
        bin_numbers.append(i + 1)  # Add bin number for this window

# Write the results to the output file
with open(args.output_file, 'w') as f:
    f.write("Header,Bin Number,Heterozygosity\n")
    for header, bin_number, heterozygosity in zip(header_names, bin_numbers, heterozygosity_per_sequence):
        f.write(f"{header},{bin_number},{heterozygosity:.6f}\n")
print("Results written to", args.output_file)
