#!/bin/bash

# Modules to add
module add python/2.7.2
module add gcc/5.4.0

# Variables
FILENAME=$1
count=0

while read pre
        do
                let count++
                echo "$count $pre"
                $HOME/trim_galore --paired --fastqc_args "--outdir fastQC/outputDir/" --trim1 ${pre}_R1.fq.gz ${pre}_R2.fq.gz

done < $FILENAME

echo -e "\nTotal $count lines read"
