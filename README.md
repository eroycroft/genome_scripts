# genome_scripts
Some scripts for exploring whole genome re-sequencing data

### heteroWindows.py
Counts heterozygous sites in non-overlapping windows of a specified length (e.g. 1000bp), in a .g.vcf file then takes the average count across agregated windows, e.g. average heterozygous sites/kb per 1000 windows (specify aggregate window size)

### Example:
python3 heteroWindows.py [input.g.vcf.gz] [output.txt]

### Optional:
 --window_size [Set the base pair window size, default 1000bp]
 
 --averaging_window [Set the number of windows to average, default 1000]
 
 --resume [Resume from the most recent line if the script is interrupted, read from the last line of file to determine the block to resume from]

