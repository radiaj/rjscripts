#!/bin/bash

#Variables
myDate1=`date "+%m%d%y"`
FILENAME=$1
myDIR=$2
# If $3 is not passed, set the working directory to current directory
myBAMDIR="${3:-${PWD}}"

count=0

while read mySIZE
do
      let count++
      echo "$count $mySIZE"

                #Go to folder where the RNAseq bam files are stored
                cd "${myBAMDIR}"/

                #Get the bam file folder to convert the  bam file to fastq R1 and R2
                java -jar $HOME/picard-tools-1.119/SamToFastq.jar VALIDATION_STRINGENCY=SILENT INPUT="${mySIZE}".bam FASTQ="${myDIR}"/"${mySIZE}"_R1.fastq SECOND_END_FASTQ="${myDIR}"/"${mySIZE}"_R2.fastq

                # Compress the fasta file for downstream analysis and less space
                cd "${myDIR}"/
                gzip -c "${mySIZE}"_R1.fastq > "${mySIZE}"_R1.fq.gz &
                gzip -c "${mySIZE}"_R2.fastq > "${mySIZE}"_R2.fq.gz &

done < $FILENAME

echo -e "\nTotal $count Lines read"

wait
