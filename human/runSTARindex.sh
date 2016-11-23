#!/bin/bash

# Add the modules needed for the analysis
module add STAR/2.4.2
module add SAMtools

# Go to the Chromosomes directory
cd /RQexec/johnsonr/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/

# Build a STAR genome index into the Chromosomes folder
STAR  --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles hg19.fa --runThreadN 4 --sjdbGTFfile /RQexec/johnsonr/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf  --sjdbOverhang 99

