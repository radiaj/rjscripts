#!/bin/bash

#Variables
myDate1=`date "+%m%d%y"`
FILENAME=$1
myDIR="${2:-${PWD}}"
count=0

while read mySIZE
do
      let count++
      echo "$count $mySIZE"

                # Go to folder where the RNAseq bam files are stored
                cd "${myDIR}"

                # Compress the fasta file for downstream analysis and less space
#                gzip -c "${mySIZE}".f*q > "${mySIZE}".fq.gz
                 gzip -c "${mySIZE}"_1.f*q > "${mySIZE}"_1.fq.gz
                 gzip -c "${mySIZE}"_2.f*q > "${mySIZE}"_2.fq.gz
                # remove the fasta file
#                rm "${mySIZE}".f*q

done < $FILENAME

echo -e "\nTotal $count Lines read"
