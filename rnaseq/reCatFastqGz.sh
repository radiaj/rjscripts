
#!/bin/bash
# Add the modules

#Variables
myDate1=`date "+%m%d%y"`
FILENAME=$1
myRNAdir=$2
myDIR="${3:-${PWD}}"
count=0

cd "${myDIR}"/

while read mySIZE
do
      let count++
      echo "$count $mySIZE"

                #Go to folder where the RNAseq bam files are stored
                cd "${myDIR}"/"${mySIZE}"/

                gzip -k -d *.fastq.gz #Keeps the original files

                # concatenate the the files
                cat  *_R1_*.fastq > "${mySIZE}"_R1.fastq
                cat  *_R2_*.fastq > "${mySIZE}"_R2.fastq

                # Then gzip the file
                gzip -c "${mySIZE}"_R1.fastq > "${mySIZE}"_R1.fq.gz
                gzip -c "${mySIZE}"_R2.fastq > "${mySIZE}"_R2.fq.gz

                # Move the files into the directory needed to run trimmomatic
                mv "${mySIZE}"_R1.fq.gz "${myRNAdir}"
                mv "${mySIZE}"_R2.fq.gz "${myRNAdir}"


done < $FILENAME

echo -e "\nTotal $count Lines read"

cd "${myRNAdir}"/

# Run trimmomatic with the command below (modify accordingly)
# ~/trimmomaticReads.sh samplesDec12b.txt
