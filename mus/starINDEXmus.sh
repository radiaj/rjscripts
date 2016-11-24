#!/bin/bash

# Add the modules needed for the analysis
module add STAR
module add SAMtools

# Go to the STAR directory in the PARK_LAB folder
#cd /RQexec/johnsonr/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/

# Build a STAR genome index into the Chromosomes folder
#STAR  --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles mm10.fa --runThreadN 4 --sjdbGTFfile /RQexec/johnsonr/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf  --sjdbOverhang 99

