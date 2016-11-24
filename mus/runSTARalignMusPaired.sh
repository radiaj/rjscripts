#!/bin/bash

# Add the modules needed for the analysis
module add STAR
module add SAMtools

# Variables
FILENAME=$1
myDIR=$2
myRESDIR=$3
count=0

cd /RQexec/johnsonr/JONES_LAB/scripts
while read mySIZE
        do
                let count++
                echo "$count $mySIZE"

                # Align the RNAseq reads to the genome with STAR
                cd "${myRESDIR}"/

                mkdir "${mySIZE}"_STAR1
                cd "${mySIZE}"_STAR1
                STAR --genomeDir /RQexec/johnsonr/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/  --runThreadN 4 --readFilesIn "${myDIR}"/"${mySIZE}"*1.f*q.gz "${myDIR}"/"${mySIZE}"*2.f*q.gz --readFilesCommand zcat --outFileNamePrefix "${mySIZE}"

done < $FILENAME

echo -e "\nTotal $count lines read"
