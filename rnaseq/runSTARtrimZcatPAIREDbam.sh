#!/bin/bash

# Add the modules needed for the analysis
module add STAR/2.4.2
module add SAMtools

# Variables
FILENAME=$1
myDIR=$2
count=0

                # Go to the STAR directory
#               cd /RQexec/johnsonr/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/

                # Build a STAR genome index into the Chromosomes folder
#               STAR  --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles hg19.fa --runThreadN 4 --sjdbGTFfile /RQexec/johnsonr/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf  --sjdbOverhang 99

cd /RQusagers/johnsonr/PARK_LAB/scripts
while read mySIZE
        do
                let count++
                echo "$count $mySIZE"

                # Go to the STAR directory in the PARK_LAB folder
                # Align the RNAseq reads to the genome with STAR
                cd /RQexec/johnsonr/PARK_LAB/STAR
                mkdir "${mySIZE}"_STAR

                # Create a sorted BAM files
                cd "${mySIZE}"_STAR

                STAR --genomeDir /RQexec/johnsonr/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes  --quantMode GeneCounts --runThreadN 4 --outSAMtype BAM SortedByCoordinate --readFilesIn "${myDIR}"/"${mySIZE}"*1.trim.f*q.gz "${myDIR}"/"${my\
SIZE}"*2.trim.f*q.gz --readFilesCommand zcat --outFileNamePrefix "${mySIZE}"

                # Create index file to view bam in IGV
                samtools index "${mySIZE}"Aligned.sortedByCoord.out.bam.bam "${mySIZE}"Aligned.sortedByCoord.out.bam.bai

done < $FILENAME

echo -e "\nTotal $count lines read"

