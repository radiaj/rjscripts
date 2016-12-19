#!/bin/bash
# Add the modules

#Variables
myDate1=`date "+%m%d%y"`
FILENAME=$1
myDIR=$2
myRNAdir=$3
count=0

cd "${myDIR}"/

while read mySIZE
do
      let count++
      echo "$count $mySIZE"

                #Go to folder where the RNAseq bam files are stored
                cd "${myDIR}"/"${mySIZE}"/

                # concatenate the the files
                cat  *_R1_*.fastq.gz > "${mySIZE}"_R1.fq.gz &
                cat  *_R2_*.fastq.gz > "${mySIZE}"_R2.fq.gz &
                wait

                mv "${mySIZE}"_R1.fq.gz "${myRNAdir}" &
                mv "${mySIZE}"_R2.fq.gz "${myRNAdir}" &
                wait

done < $FILENAME

echo -e "\nTotal $count Lines read"

wait

