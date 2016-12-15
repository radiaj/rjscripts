#!/bin/bash

# Variables
FILENAME=$1
myFILEDIR="${2:-${PWD}}"
count=0

cd "${myFILEDIR}"
cp ~/Trimmomatic-0.36/adapters/*.fa .
mkdir logs

while read pre
        do
                let count++
                echo "$count $pre"

java -jar $HOME/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog logs/$pre.log ${pre}_R1.f*q.gz ${pre}_R2.f*q.gz ${pre}_1.trim.fq.gz ${pre}_1.unpaired.fq.gz ${pre}_2.trim.fq.gz ${pre}_2.unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done < $FILENAME

echo -e "\nTotal $count lines read"





