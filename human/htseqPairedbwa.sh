#!/bin/bash

# Add the modules needed for the analysis
module add HTSeq/0.6.1p1
module add SAMtools

# Variables
FILENAME=$1
myGTF="/RQexec/johnsonr/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
count=0

while read mySAMPLE
        do
                let count++
                echo "$count $mySAMPLE"

                # Go to the STAR directory in the PARK_LAB folder
                cd $SCRATCH/PARK_LAB/RNAseq/"${mySAMPLE}"/A*/illumina_wtss/bwa_mem_aligned

                # python -m HTSeq.scripts.count -f sam "${mySAMPLE}"Aligned.out.sam "${myGTF}" >  "${mySAMPLE}".cnts
                # htseq-counts --stranded=no --mode=intersection-nonempty -r pos accepted_hits.bam annotation.gtf > output.txt
                python -m HTSeq.scripts.count --stranded=no --mode=intersection-nonempty -r pos "${mySAMPLE}"*.bam "${myGTF}" > "${mySAMPLE}"cnts.txt
                
done < $FILENAME

echo -e "\nTotal $count lines read"
