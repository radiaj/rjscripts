# rjscripts
Scripts to run analysis on cluster.

see http://davetang.org/wiki/tiki-index.php?page=SAMTools for samtools commands

# Converting miR sequences for script use
$ module add fastx_tools/0.0.13

You can download the hairpin.fa file from http://www.mirbase.org/ftp.shtml
$ cd /RQusagers/johnsonr/miRNA_reference
$ fasta_formatter -i hairpin.fa -o hairpin2.fa -w 0
$ fasta_formatter -i mature.fa -o mature2.fa -w 0
$ tr 'U' 'T' < hairpin2.fa > hairpindna2.fa

# Convert U to T for bowtie
$ tr 'U' 'T' < hairpin.fa > hairpindna.fa
'# fasta_nucleotide_changer  -i hairpin2.fa -o hairpindna.fa -d # doesn't work with Y and other nuclueotide sequence

# Create reads fasta file for quantifier (test)
$ cd /exec5/GROUP/johnsonr/johnsonr/DUCHAINE_LAB/fl_17-92Precursor
$ fastq_to_fasta -i fl_17-92Precursor.fq -o fl_17-92Precursor.fa -v

Input: 55796966 reads.
Output: 55336058 reads.
discarded 460908 (0%) low-quality reads.
