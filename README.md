# genome_scripts
Some scripts for exploring whole genome re-sequencing data

## heteroWindows.py
Counts heterozygous sites in non-overlapping windows of a specified length (e.g. 1000bp), in a .g.vcf file then takes the average count across agregated windows, e.g. average heterozygous sites/kb per 1000 windows (specify aggregate window size)

### Example:
python3 heteroWindows.py [input.g.vcf.gz] [output.txt]

### Optional:
 --window_size [Set the base pair window size, default 1000bp]
 
 --averaging_window [Set the number of windows to average, default 1000]
 
 --resume [Resume from the most recent line if the script is interrupted, read from the last line of file to determine the block to resume from]



## extract_snps.py
Extracts one random (or all) biallelic SNPs from multi-sample sequence alignments (FASTA/Phylip), enforcing custom rules (see below). 
Handles alignments per locus in a directory, and outputs concatenated nucleotide and binary SNPs in .fasta format, and a .csv in binary format (to be imported into R)

### Example:
 python extract_snps.py --input <dir> --out-prefix <prefix> [options]

### Required:
  --input <dir> [Directory containing alignment files in .fasta or .phy format]
  
  --out-prefix <prefix> [Prefix for output files]

### Optional:
  --taxa-file <file> [File containing one sample name per line. Only these taxa are used when identifying/filtering SNPs if specified]
  
  --all-snps [Output all SNPs from each locus instead of one at random]
  
  --min-maf <int> [Minimum count of minor allele, default = 1]
  
  --max-missing <float> [Maximum fraction of missing data allowed per SNP, i.e. between 0–1, default 0.2]
  
  --replicates <int> [Number of random replicates to run, default = 1]
  
  --seed <int> [Random seed for reproducibility]

### Outputs:
  1. <prefix>_repN_nucleotides.fasta [Concatenated SNPs in nucleotide/IUPAC codes]
  2. <prefix>_repN_binary.fasta [same SNPs, encoded as 0/1/2/NA]
  3. <prefix>_repN_binary.csv [CSV version of binary SNPs with loci as columns]
  4. <prefix>_repN_loci.txt [Ordered list of locus files contributing SNPs]
