#!/usr/bin/env python3
"""
Extract one random (or all) biallelic SNPs from multi-sample alignments (FASTA/Phylip).
Handles alignments per locus in a directory, and concatenates SNPs across loci.

Usage:
  python extract_snps.py --input <dir> --out-prefix <prefix> [options]

Required arguments:
  --input <dir>            Directory containing alignment files (.fasta or .phy).
  --out-prefix <prefix>    Prefix for output files.

Optional arguments:
  --taxa-file <file>       File containing one taxon/sample name per line.
                           Only these taxa are used when identifying/filtering SNPs.
                           If not provided, all taxa are used.
  --all-snps               Output all SNPs from each locus instead of one at random.
  --min-maf <int>          Minimum count of minor allele (default 1).
  --max-missing <float>    Maximum fraction of missing data allowed per SNP (0–1, default 0.2).
  --replicates <int>       Number of random replicates to run (default 1).
  --seed <int>             Random seed for reproducibility (optional).

Outputs:
  1. <prefix>_repN_nucleotides.fasta  → Concatenated SNPs in nucleotide/IUPAC codes.
  2. <prefix>_repN_binary.fasta       → Same SNPs, encoded as 0/1/2/NA.
  3. <prefix>_repN_binary.csv         → CSV version of binary SNPs with loci as columns.
  4. <prefix>_repN_loci.txt           → Ordered list of locus files contributing SNPs.
"""

import os
import random
import argparse
import csv
from Bio import AlignIO
from collections import defaultdict


# IUPAC heterozygous codes
iupac_hetero = {
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'G', 'C'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'}
}

valid_bases = {'A', 'C', 'G', 'T'}


def parse_args():

    parser = argparse.ArgumentParser(
        description="Extract random biallelic SNPs from alignments"
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input directory of alignment files (.fasta/.phy)"
    )

    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Prefix for output files"
    )

    parser.add_argument(
        "--taxa-file",
        help="File containing one taxon/sample name per line. "
             "Only these taxa are used for SNP filtering. "
             "If not provided, all taxa are used."
    )

    parser.add_argument(
        "--all-snps",
        action="store_true",
        help="Output all SNPs instead of one random per locus"
    )

    parser.add_argument(
        "--min-maf",
        type=int,
        default=1,
        help="Minimum count of minor allele (default 1)"
    )

    parser.add_argument(
        "--max-missing",
        type=float,
        default=0.2,
        help="Max fraction missing per SNP (default 0.2)"
    )

    parser.add_argument(
        "--replicates",
        type=int,
        default=1,
        help="Number of random replicates (default 1)"
    )

    parser.add_argument(
        "--seed",
        type=int,
        help="Random seed (optional)"
    )

    return parser.parse_args()


def is_valid_base(char):

    return char in valid_bases or char in iupac_hetero


def get_alleles(char):

    if char in valid_bases:
        return {char}

    return iupac_hetero.get(char, set())


def encode_binary(col, alleles):

    allele1, allele2 = sorted(alleles)

    encoded = []

    for c in col:

        if c == 'N' or c in ['-', 'X']:
            encoded.append("NA")

        elif c == allele1:
            encoded.append("0")

        elif c == allele2:
            encoded.append("2")

        elif (
            c in iupac_hetero
            and iupac_hetero[c] == {allele1, allele2}
        ):
            encoded.append("1")

        else:
            encoded.append("NA")

    return encoded


def process_alignment(file, min_maf, max_missing, selected_taxa=None):

    try:
        aln = AlignIO.read(file, "fasta")

    except:

        try:
            aln = AlignIO.read(file, "phylip-relaxed")

        except:

            return None, None

    # All taxa in the alignment
    all_samples = [rec.id for rec in aln]

    # Sequences indexed by sample name
    all_seqs = {
        rec.id: str(rec.seq).upper()
        for rec in aln
    }

    # ---------------------------------------------------------
    # Select taxa used for SNP filtering
    # ---------------------------------------------------------

    if selected_taxa is None:

        samples = all_samples

    else:

        samples = [
            s for s in all_samples
            if s in selected_taxa
        ]

    # No requested taxa found in this alignment
    if not samples:

        return None, None

    # Only retain sequences for selected taxa
    seqs = {
        s: all_seqs[s]
        for s in samples
    }

    aln_len = aln.get_alignment_length()

    snp_sites = []

    for i in range(aln_len):

        col = [
            seqs[s][i]
            for s in samples
            if i < len(seqs[s])
        ]

        if not any(is_valid_base(c) for c in col):
            continue

        alleles = set()
        missing = 0

        for c in col:

            if c in valid_bases:

                alleles.add(c)

            elif c in iupac_hetero:

                alleles |= iupac_hetero[c]

            else:

                missing += 1

        # Must be exactly biallelic
        if len(alleles) != 2:
            continue

        # Missing-data filter
        if missing / len(samples) > max_missing:
            continue

        # Count alleles
        counts = defaultdict(int)

        for c in col:

            for a in get_alleles(c):

                counts[a] += 1

        # Must still have exactly two alleles
        if len(counts) != 2:
            continue

        # Minimum minor allele count
        if min(counts.values()) < min_maf:
            continue

        # -----------------------------------------------------
        # NEW FILTER:
        # Require variation in genotype state
        #
        # This excludes sites where all individuals are
        # heterozygous (e.g. 1,1,1,1,1...), which are
        # biallelic at the allele level but have zero
        # genotype variance and therefore cannot contribute
        # to PCA.
        # -----------------------------------------------------

        encoded = encode_binary(
            col,
            sorted(alleles)
        )

        nonmissing_genotypes = [
            e for e in encoded
            if e != "NA"
        ]

        # Must have at least two different genotype states
        if len(set(nonmissing_genotypes)) < 2:
            continue

        # If we get here, this is a usable SNP
        snp_sites.append((i, col))

    return samples, snp_sites


def main():

    args = parse_args()

    # ---------------------------------------------------------
    # Random seed
    # ---------------------------------------------------------

    if args.seed is not None:

        random.seed(args.seed)

    # ---------------------------------------------------------
    # Read taxa file
    # ---------------------------------------------------------

    selected_taxa = None

    if args.taxa_file:

        with open(args.taxa_file) as f:

            selected_taxa = {
                line.strip()
                for line in f
                if line.strip()
            }

        print(
            f"Using {len(selected_taxa)} taxa from "
            f"{args.taxa_file}"
        )

    else:

        print(
            "No taxa file supplied: using all taxa."
        )

    # ---------------------------------------------------------
    # Find alignment files
    # ---------------------------------------------------------

    files = [
        os.path.join(args.input, f)
        for f in os.listdir(args.input)
        if f.endswith((".fasta", ".fa", ".phy"))
    ]

    # ---------------------------------------------------------
    # Process replicates
    # ---------------------------------------------------------

    for r in range(1, args.replicates + 1):

        nuc_data = defaultdict(list)
        bin_data = defaultdict(list)
        loci_used = []

        for f in files:

            samples, snps = process_alignment(
                f,
                args.min_maf,
                args.max_missing,
                selected_taxa
            )

            if not snps:

                continue

            # One random SNP per locus, unless --all-snps
            chosen = (
                snps
                if args.all_snps
                else [random.choice(snps)]
            )

            for idx, col in chosen:

                # Identify the two alleles
                alleles = set()

                for c in col:

                    alleles |= get_alleles(c)

                alleles = sorted(list(alleles))[:2]

                # Nucleotide data
                for s, c in zip(samples, col):

                    nuc_data[s].append(c)

                # Binary data
                encoded = encode_binary(
                    col,
                    alleles
                )

                for s, e in zip(samples, encoded):

                    bin_data[s].append(e)

                loci_used.append(
                    os.path.basename(f)
                )

        # -----------------------------------------------------
        # Write nucleotide FASTA
        # -----------------------------------------------------

        out_nuc = (
            f"{args.out_prefix}_rep{r}_nucleotides.fasta"
        )

        with open(out_nuc, "w") as out:

            for s, seq in nuc_data.items():

                out.write(
                    f">{s}\n{''.join(seq)}\n"
                )

        # -----------------------------------------------------
        # Write binary FASTA
        # -----------------------------------------------------

        out_bin = (
            f"{args.out_prefix}_rep{r}_binary.fasta"
        )

        with open(out_bin, "w") as out:

            for s, seq in bin_data.items():

                out.write(
                    f">{s}\n{''.join(seq)}\n"
                )

        # -----------------------------------------------------
        # Write loci list
        # -----------------------------------------------------

        out_loci = (
            f"{args.out_prefix}_rep{r}_loci.txt"
        )

        with open(out_loci, "w") as out:

            out.write(
                "\n".join(loci_used) + "\n"
            )

        # -----------------------------------------------------
        # Write binary CSV
        # -----------------------------------------------------

        out_csv = (
            f"{args.out_prefix}_rep{r}_binary.csv"
        )

        headers = []

        for fname in loci_used:

            base = os.path.basename(fname)

            parts = base.split("_")

            if len(parts) >= 2:

                headers.append(
                    "_".join(parts[:2])
                )

            else:

                headers.append(base)

        with open(out_csv, "w", newline="") as f:

            writer = csv.writer(f)

            writer.writerow(
                [""] + headers
            )

            for s in sorted(bin_data.keys()):

                writer.writerow(
                    [s] + bin_data[s]
                )


if __name__ == "__main__":

    main()
