# genome_scripts
Useful scripts for whole genome analysis

# calc_refGenome_hetero.py
Calculates heterozygosity in non-overlapping windows of specified length (default 1Mb), in a .fasta format reference genome with one or multiple contigs/chromosomes 

Example:
python3 [ref_genome.fasta] -w [window_size_bp] -output [outfile_name.txt]
