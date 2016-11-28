#!/bin/bash

# Add the modules needed for the analysis
module add STAR/2.4.2
module add SAMtools

# Variables
FILENAME=$1
myDIR=$2
myRESDIR=$3
count=0

# Go to the STAR directory in the PARK_LAB folder
#cd /RQexec/johnsonr/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/

# Build a STAR genome index into the Chromosomes folder
#STAR  --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles mm10.fa --runThreadN 4 --sjdbGTFfile /RQexec/johnsonr/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf  --sjdbOverhang 10

cd /RQexec/johnsonr/JONES_LAB/scripts
while read mySIZE
        do
                let count++
                echo "$count $mySIZE"

                # Align the RNAseq reads to the genome with STAR
                cd "${myRESDIR}"/

                mkdir "${mySIZE}"_STARt
                cd "${mySIZE}"_STARt
                STAR --genomeDir /RQexec/johnsonr/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/  --quantMode GeneCounts  --runThreadN 4 --readFilesIn "${myDIR}"/"${mySIZE}"*1.trim.f*q.gz "${myDIR}"/"${mySIZE}"*2.trim.f*q.gz --readFilesCommand zcat --outFileNamePrefix "${mySIZE}"

done < $FILENAME

echo -e "\nTotal $count lines read"
