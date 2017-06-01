#!/bin/bash

# Add the modules needed for the analysis
#module add picard_tools
module add Java/1.8.0_40
module add GATK/3.3-0
module add SAMtools

# Variables
FILENAME=$1
myDIR=$2
REFERENCE="$HOME/ucsc.hg19.fasta"
count=0

cd /RQusagers/johnsonr/PARK_LAB/scripts

# Run the variant calling pipeline on each patient
while read mySIZE
        do
                let count++
                echo "$count $mySIZE"

                # Go to the  directory with the bam files
                cd "${myDIR}"/"${mySIZE}"

                # Sort bam file
                samtools sort "${mySIZE}".bam "${mySIZE}".sort
                
                #Index 
				samtools index "${mySIZE}".sort.bam

				#mark duplicate
				java -Xmx5g -jar $HOME/picard-tools-1.119/MarkDuplicates.jar INPUT="${mySIZE}".sort.bam OUTPUT="${mySIZE}".sort.deduped.bam METRICS_FILE="${mySIZE}".duplicates REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE

				#Realignment based on known insert sites (Using Java 1.7 from now on as required by GATK)
				java -Xmx5g -jar $HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REFERENCE -I "${mySIZE}".sort.deduped.bam -known $HOME/1000G_phase1.indels.hg19.sites.vcf -known $HOME/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -o "${mySIZE}".realign.intervals -S LENIENT

				java -Xmx5g -jar $HOME/GenomeAnalysisTK.jar -T IndelRealigner -R $REFERENCE -I "${mySIZE}".sort.deduped.bam -targetIntervals "${mySIZE}".realign.intervals -known $HOME/1000G_phase1.indels.hg19.sites.vcf -known $HOME/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -o "${mySIZE}".sort.deduped.realigned.bam -S LENIENT

				java -Xmx5g -jar $HOME/GenomeAnalysisTK.jar -T BaseRecalibrator -R $REFERENCE -l INFO -I "${mySIZE}".sort.deduped.realigned.bam -knownSites $HOME/1000G_phase1.indels.hg19.sites.vcf -knownSites $HOME/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites $HOME/dbsnp_138.hg19.vcf -o "${mySIZE}".recalibration_report.grp -S LENIENT 

				java -Xmx5g -jar $HOME/GenomeAnalysisTK.jar -T PrintReads -R $REFERENCE -l INFO -I "${mySIZE}".sort.deduped.realigned.bam -BQSR "${mySIZE}".recalibration_report.grp -o "${mySIZE}".sort.deduped.realigned.recalibrated.bam -S LENIENT

				###########################################################################
				# GATK Variant Calling -  HaplotypeCaller
				# Set -nct, outmode, emit_thresh, call_threh,

				#outmode="EMIT_ALL_CONFIDENT_SITES"
				#emit_thresh=20	#Threshold for tagging possible variants
				#call_thresh=30	#Threshold for tagging _good_ variants
				#hetrate=0.03	#Popgen heterozygosity rate (that is, for any two random chrom in pop, what is rate of mismatch).               Human is ~0.01, so up maize to ~0.03
				#minBaseScore=30	#Minimum Phred base score to count a base (20 = 0.01 error, 30=0.001 error, etc)

echo "calling variants...."

   java -jar $HOME/GenomeAnalysisTK.jar \
     -R $REFERENCE \
     -T HaplotypeCaller \
     -I "${mySIZE}".sort.deduped.realigned.recalibrated.bam \
     --dbsnp $HOME/dbsnp_138.hg19.vcf \
     -stand_call_conf 30 \
     -stand_emit_conf 20 \
     -o "${mySIZE}"_output.raw.snps.indels.vcf

echo -e "\n Done $mySIZE."


done < $FILENAME

echo -e "\nTotal $count lines read"

# Test run 
# ./GATKvarPipeline.sh sidvariants.txt "/RQexec/johnsonr/SID_LAB/"
