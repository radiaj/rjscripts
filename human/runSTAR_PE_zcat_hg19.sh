#!/bin/bash

# Add the modules needed for the analysis
module add STAR/2.4.2
module add SAMtools

# Set variables
FILENAME=$1
myDIR=$2
mySCRIPTDIR=$3
myGENOMEDIR="/RQexec/johnsonr/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes"

# Start the count variable from 0
count=0

#------------------------------------------------------------------------------
## Run from here for all samples 
#------------------------------------------------------------------------------
# Go to the directory with your scripts
cd "${mySCRIPTDIR}"

while read mySAMPLE
        do
                let count++
                echo "$count $mySAMPLE"

                # Go to the STAR directory where the results will be stored
                # mkdir $SCRATCH/STAR # command to create folder if not already done
                cd $SCRATCH/STAR
                
                # Create the folder to store the STAR results
                mkdir "${mySAMPLE}"_STAR
                cd "${mySAMPLE}"_STAR
                
                # Align the RNAseq reads to the genome with STAR
                STAR --genomeDir "${myGENOMEDIR}"  --quantMode GeneCounts --runThreadN 4 --readFilesIn "${myDIR}"/"${mySAMPLE}"*1.f*q.gz "${myDIR}"/"${mySAMPLE}"*2.f*q.gz --readFilesCommand zcat --outFileNamePrefix "${mySAMPLE}"

done < $FILENAME

echo -e "\nTotal $count lines read"
