
#!/bin/bash

# Add the modules needed for the analysis
module add STAR/2.4.2
module add SAMtools

# Variables
FILENAME=$1
#myDIR=$2
count=0

# Go to the STAR directory
cd /RQexec/johnsonr/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/

# Build a STAR genome index into the Chromosomes folder
STAR  --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles mm10.fa --runThreadN 4 --sjdbGTFfile /RQexec/johnsonr/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf  --sjdbOverhang 149

cd /RQusagers/johnsonr/PARK_LAB/scripts
while read mySIZE myFASTQ1 myFASTQ2
	do
		let count++
		echo "$count $mySIZE" 
		echo "fastq1 $myFASTQ1"
	       
		# Go to the STAR directory in the PARK_LAB folder
		# Align the RNAseq reads to the genome with STAR
		cd /RQexec/johnsonr/PARK_LAB/STAR
		mkdir "${mySIZE}"_STAR_mus
		
		# Create a sorted BAM files
		cd "${mySIZE}"_STAR_mus

		STAR --genomeDir /RQexec/johnsonr/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/  --quantMode GeneCounts --runThreadN 4 --outSAMtype BAM Unsorted --readFilesIn "${myFASTQ1}" "${myFASTQ2}" --readFilesCommand zcat --outReadsUnmapped Fastx --outFileNamePrefix "${mySIZE}"				

done < $FILENAME

echo -e "\nTotal $count lines read"
