#!/bin/bash

# Variables
FILENAME=$1
mytmpDIR=$2
count=0

cd /user/yourDirectory

# Get index files (run once only)
xenome index -v -l index.log -T 1 -M 8 -P idx -H mouse.fa -G human.fa --tmp-dir "${mytmpDIR}"

while read mySIZE myFASTQ1 myFASTQ2
do
    let count++
    echo "$count $mySIZE"
    echo "fastq1 $myFASTQ1"

    xenome classify -T 1 -M 8 -P idx  --tmp-dir "${mytmpDIR}" --pairs --host-name mouse --graft-name human -i ./"${myFASTQ1}" -i ./"${myFASTQ2}" --output-filename-prefix "${mySIZE}"

done < $FILENAME

echo -e "\nTotal $count lines read"
