#!/bin/bash

# Add the modules needed for the analysis
module add SAMtools
# Make sure samtools is installed on linux

# Variables
FILENAME=$1
# If $2 is not passed, set the working directory to current directory
myFILEDIR="${2:-${PWD}}"
count=0

while read mySIZE
        do
                let count++
                echo "$count $mySIZE"
                echo "Directory set at $myFILEDIR"
                # Go to the STAR directory in the PARK_LAB folder
                cd "${myFILEDIR}"/"${mySIZE}"_STAR/

                # Convert the sam file directly to a sorted bam file
                samtools sort  -o "${mySIZE}"_sorted.bam -O bam "${mySIZE}"Aligned.out.bam
                samtools index "${mySIZE}"_sorted.bam "${mySIZE}"_sorted.bai &

done < $FILENAME

echo -e "\nTotal $count lines read"

wait




