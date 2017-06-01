#!/bin/bash

# Add the modules needed for the analysis
module add BWA/0.7.7
module add bowtie/0.12.7
module add mirdeep2/2.0.0.5

# Variables
FILENAME=$1
myDIR=$2
myREF="/RQexec/johnsonr/Mus_musculus/UCSC/mm10/Sequence/BowtieIndex"
myMIR="/RQusagers/johnsonr/miRNA_reference"
count=0

# For mouse mm10 and Dr. Duchaines precursor miRNA files
while read mySIZE
	do
		let count++
		echo "$count $mySIZE" 

		# Go to the BWA directory in the PARK_LAB folder
		cd /exec5/GROUP/johnsonr/johnsonr/DUCHAINE_LAB/BWA

		# Align the reads with bwa
#		bwa aln -n 1 -o 0 -e 0 -k 1 -t 4 /RQexec/johnsonr/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/version0.6.0/genome.fa "${myDIR}"/"${mySIZE}"/"${mySIZE}"_1.fq > "${mySIZE}"_1.sai       
#		bwa aln -n 1 -o 0 -e 0 -k 1 -t 4 /RQexec/johnsonr/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/version0.6.0/genome.fa "${myDIR}"/"${mySIZE}"/"${mySIZE}"_2.fq > "${mySIZE}"_2.sai
		
		# Create a paired end sam file
#		bwa samse -f "${mySIZE}".sam /RQexec/johnsonr/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/version0.6.0/genome.fa "${mySIZE}_1".sai "${mySIZE}_2".sai "${myDIR}"/"${mySIZE}"/"${mySIZE}"_1.fa "${myDIR}"/"${mySIZE}"/"${mySIZE}"_2.fa

bwa_sam_converter.pl -i "${mySIZE}".sam -c -o "${mySIZE}"_reads_collapsed.fa -a "${mySIZE}"_reads_collapsed_vs_genome.arf

#mkdir "${mySIZE}"_mirDEEP2
cd  "${mySIZE}"_mirDEEP2
miRDeep2.pl ../"${mySIZE}"_reads_collapsed.fa "${myREF}"/genome.fa ../"${mySIZE}"_reads_collapsed_vs_genome.arf "${myMIR}"/mature.fa none "${myMIR}"/hairpin.fa -P -t Mouse 2> report.log

done < $FILENAME

echo -e "\nTotal $count lines read"
