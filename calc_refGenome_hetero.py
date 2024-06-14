from Bio import SeqIO
import numpy as np
import argparse

#calculate heterozygosity in a reference genome in non-overlapping windows 
# -r your_reference.fasta
# -w 1000000 (window size in base pairs)
# -o output.txt
# eroycroft 2024

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
parser.add_argument('genome_file', metavar='genome_file', type=str, help='Path to the genome file in FASTA format')
parser.add_argument('-r', '--reference_genome', metavar='reference_genome', type=str, default=None, help='Path to the reference genome file in FASTA format')
parser.add_argument('-w', '--window_size', metavar='window_size', type=int, default=1000000, help='Size of the window in base pairs (default: 1 Mb)')
parser.add_argument('-o', '--output_file', metavar='output_file', type=str, default='heterozygosity_results.txt', help='Output file name (default: heterozygosity_results.txt)')
args = parser.parse_args()

# Read the FASTA file
genome_sequence = SeqIO.read(args.genome_file, "fasta").seq

# Read the reference genome if provided
if args.reference_genome:
    reference_sequence = SeqIO.read(args.reference_genome, "fasta").seq
    if len(reference_sequence) != len(genome_sequence):
        raise ValueError("Length of the reference genome must match the length of the input genome.")

# Define window size
window_size = args.window_size

# Calculate heterozygosity in non-overlapping windows
num_windows = len(genome_sequence) // window_size
heterozygosity_per_window = []
for i in range(num_windows):
    start = i * window_size
    end = start + window_size
    window_sequence = genome_sequence[start:end]
    if args.reference_genome:
        reference_window = reference_sequence[start:end]
        if len(set(reference_window)) != 1:
            raise ValueError("Reference genome window contains more than one base.")
        if reference_window[0] != 'N':
            heterozygosity = calculate_heterozygosity(window_sequence)
            heterozygosity_per_window.append(heterozygosity)
        else:
            heterozygosity_per_window.append(np.nan)
    else:
        heterozygosity = calculate_heterozygosity(window_sequence)
        heterozygosity_per_window.append(heterozygosity)

# Write the results to the output file
with open(args.output_file, 'w') as f:
    for i, heterozygosity in enumerate(heterozygosity_per_window):
        if np.isnan(heterozygosity):
            f.write(f"Window {i+1}: Reference genome window contains N bases.\n")
        else:
            f.write(f"Window {i+1}: Heterozygosity = {heterozygosity:.6f}\n")
print("Results written to", args.output_file)
