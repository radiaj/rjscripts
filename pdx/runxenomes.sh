#--------------------------
# To run xenomes
#--------------------------
#!/bin/bash


cd /media/cyberlogic/a384d9c7-6dda-4d1b-b32d-9a6b0e007931/toto/PARK_LAB_DATA_FILES/PDX_RNAseq_Files/
# xenome index -T 8 -P idx -H mm10.fa -G hg19.fa

while read mySIZE myFASTQ1 myFASTQ2
        do
                let count++
                echo "$count $mySIZE"
                echo "fastq1 $myFASTQ1"

                xenome classify -T 8 -P idx --pairs --host-name mm10 --graft-name hg19 -i ./"${myFASTQ1}" -i ./"${myFASTQ2}" --output-filename-prefix "${mySIZE}"

done < $FILENAME

echo -e "\nTotal $count lines read"
