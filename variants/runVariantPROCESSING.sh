#!/bin/bash
#PBS -A bap-052-aa
#PBS -l nodes=1:m48G:ppn=16
#PBS -l walltime=36:00:00
#PBS -o outputfileRUN1
#PBS -e errorfileM1RUN1
#PBS -N star001
#PBS -r n

# Add the modules needed for the analysis
#module add picard_tools
module add Java/1.8.0_40
module add GATK/3.3-0
module add SAMtools

# Variables
FILENAME=$1
count=0
myREFDir="/RQexec/johnsonr/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/"

# Create the reference fasta file index (only needs to be run once)
# cd /RQexec/johnsonr/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/
# samtools faidx /RQexec/johnsonr/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/hg19.fa

# Create the hg19 dict file need for SplitNCigarReads
# java -jar $HOME/picard-tools-1.119/CreateSequenceDictionary.jar R= hg19.fa O= hg19.dict


cd /RQusagers/johnsonr/PARK_LAB/scripts

# Run the variant calling pipeline on each patient
while read mySIZE
        do
                let count++
                echo "$count $mySIZE"

                # Go to the STAR directory in the PARK_LAB folder
                cd /RQexec/johnsonr/PARK_LAB/STAR/"${mySIZE}"_STAR1

                java -jar $HOME/picard-tools-1.119/AddOrReplaceReadGroups.jar I="${mySIZE}"Aligned.out.sam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample # works

                java -jar $HOME/picard-tools-1.119/MarkDuplicates.jar I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics  # works

#               java -jar $HOME/GenomeAnalysisTK.jar -T SplitNCigarReads -R "${myREFDir}"/hg19.fa -I dedupped.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS # works in terminal
#               java -jar $HOME/GenomeAnalysisTK.jar -T HaplotypeCaller -R "${myREFDir}"/hg19.fa -I split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o "${mySIZE}"output.vcf # works in terminal

 #      java -jar $HOME/GenomeAnalysisTK.jar -T VariantFiltration -R "${myREFDir}"/hg19.fa -V  "${mySIZE}"output.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o "${mySIZE}"filtered.vcf # works in the terminal

echo -e "\n Done $mySIZE."


done < $FILENAME

echo -e "\nTotal $count lines read"
