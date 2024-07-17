#heteroWindows.py ERoycroft 2024
#The script counts heterozygous sites in windows of 1000 base pairs (or specified window size)
#and the takes the average count across 1000 such 1000 bp windows (or specify aggregate window size)

#e.g. average het sites/1kb each Mb

 #   Window Size: Each window is 1000 base pairs.
 #   Aggregation: aggregates the counts and other statistics over 1000 x 1000 bp windows.
 #   Output: calculates the averages and prints the results for each set of 1000 windows.

#takes a single sample gvcf file, and outputs a table out counts per window

#usage
#python3 heteroWindows.py input.g.vcf.gz output.txt


#options
# --window_size Set the base pair window size.
# --averaging_window Set the number of windows to average.
# --resume Resume from the most recent line if the script is interrupted.

#If --resume is specified and the output file exists, the script will read the last line to determine the contig and window block to resume from


import pysam
import sys
import argparse
import os

def count_heterozygous_sites(gvcf_file, output_file, window_size, averaging_window, resume):
    vcf = pysam.VariantFile(gvcf_file)
    start_from_contig = None
    start_from_window = 0

    # Check for resume option
    if resume and os.path.exists(output_file):
        with open(output_file, 'r') as out:
            last_line = out.readlines()[-1]
            if last_line.startswith("Contig"):
                last_line = None  # Handle case where file exists but is empty apart from headers
            if last_line:
                last_line = last_line.strip().split("\t")
                start_from_contig = last_line[0]
                start_from_window = int(last_line[3]) * averaging_window

    with open(output_file, 'a') as out:
        # Write headers if not resuming
        if not resume or not os.path.exists(output_file):
            out.write("Contig\tAvg_Window_Length\tAvg_Heterozygous_Count\tWindow_Block_Number\tAvg_Proportion_Heterozygous\n")

        for contig in vcf.header.contigs:
            if start_from_contig and contig != start_from_contig:
                continue

            contig_name = contig
            contig_length = vcf.header.contigs[contig].length
            window_start = start_from_window if start_from_contig else 0
            window_count = 0
            total_windows = 0

            agg_window_length = 0
            agg_heterozygous_count = 0
            agg_valid_sites_count = 0
            agg_window_count = 0

            print(f"Processing contig: {contig_name}")

            while window_start < contig_length:
                window_end = min(window_start + window_size, contig_length)
                heterozygous_count = 0
                valid_sites_count = 0

                try:
                    for rec in vcf.fetch(contig_name, window_start, window_end):
                        if 'GT' in rec.samples[0] and rec.samples[0]['GT'] == (0, 1):
                            heterozygous_count += 1
                        valid_sites_count += 1
                except ValueError as e:
                    print(f"Error fetching records for {contig_name}:{window_start}-{window_end}: {e}")
                    break

                window_count += 1
                total_windows += 1
                window_length = window_end - window_start
                proportion_heterozygous = heterozygous_count / window_length if window_length > 0 else 0

                agg_window_length += window_length
                agg_heterozygous_count += heterozygous_count
                agg_valid_sites_count += valid_sites_count
                agg_window_count += 1

                # Every averaging_window windows, output the aggregated data
                if agg_window_count == averaging_window or window_start + window_size >= contig_length:
                    avg_window_length = agg_window_length / agg_window_count
                    avg_heterozygous_count = agg_heterozygous_count / agg_window_count
                    avg_proportion_heterozygous = (agg_heterozygous_count / agg_window_length) if agg_window_length > 0 else 0

                    out.write(f"{contig_name}\t{avg_window_length:.2f}\t{avg_heterozygous_count:.2f}\t{total_windows // averaging_window}\t{avg_proportion_heterozygous:.6f}\n")

                    agg_window_length = 0
                    agg_heterozygous_count = 0
                    agg_valid_sites_count = 0
                    agg_window_count = 0

                window_start += window_size

                if window_count % 10000 == 0:
                    print(f"Processed {window_count} windows on contig {contig_name}")

            print(f"Finished processing contig: {contig_name}")

            if start_from_contig:
                start_from_contig = None  # Continue normally after resuming the specific contig

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count heterozygous sites in gVCF file.')
    parser.add_argument('input_gvcf', type=str, help='Input gVCF file')
    parser.add_argument('output_txt', type=str, help='Output text file')
    parser.add_argument('--window_size', type=int, default=1000, help='Window size in base pairs (default: 1000)')
    parser.add_argument('--averaging_window', type=int, default=1000, help='Number of windows to average (default: 1000)')
    parser.add_argument('--resume', action='store_true', help='Resume from the most recent line if interrupted')

    args = parser.parse_args()

    count_heterozygous_sites(args.input_gvcf, args.output_txt, args.window_size, args.averaging_window, args.resume)
