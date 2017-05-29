#!/bin/bash

# Add the modules needed for the analysis
#module add BWA/0.7.7
module add bowtie/0.12.7
module add mirdeep2/2.0.0.5
module add SAMtools/1.1
module add fastx_tools/0.0.13

# Variables
FILENAME=$1
myDIR=$2
myREF="/RQexec/johnsonr/Mus_musculus/UCSC/mm10/Sequence/BowtieIndex"
myMIR="/RQusagers/johnsonr/miRNA_reference"
count=0

# Build the bowtie index for hairpin.fa
#cd  "${myMIR}"/
#bowtie-build "${myMIR}"/hairpindna.fa  hairpindna_bowtie
#bowtie-build "${myMIR}"/hairpindna.fa  hairpindna
#bowtie-build "${myMIR}"/mature.fa  mature

cd /RQusagers/johnsonr/PARK_LAB/scripts
# For mouse mm10 and Dr. Duchaines precursor miRNA files
while read mySIZE
        do
                let count++
                echo "$count $mySIZE"

                # Go to the sample directory in the DUCHAINE_LAB folder
                cd  "${myDIR}"/"${mySIZE}"/
                pwd
                # Prepare the combined fastq.gz file
#               cat "${mySIZE}"_1.fq "${mySIZE}"_2.fq > "${mySIZE}".fq
                # Or for the matures only
#               mv clean.fq "${mySIZE}".fq

                # Generate the fasta file for the quantifier.pl
#               fastq_to_fasta -i "${mySIZE}".fq -o  "${mySIZE}".fa -v > fasta.log

                # Align the reads with bowtie
                bowtie -n 1 -l 8 -a --best --strata --phred33-quals -p 4 --un "${mySIZE}"_unmapped.fq "${myMIR}"/hairpindna_bowtie "${mySIZE}".fq "${mySIZE}"_bowtie.sam

                # Try quantifier script from mirDeep2
                # Convert the sam file to a bam file
                cat "${mySIZE}"_bowtie.sam | grep -v ^@  | awk '{print ">"$1" "$2"\n"$6}' > "${mySIZE}".mapped.fa

                # Run the mapper on the mapped fasta file
                mapper.pl "${mySIZE}".mapped.fa -p /RQusagers/johnsonr/miRNA_reference/hairpindna -c -m -s "${mySIZE}".mapped.fa.col -t "${mySIZE}".mapped.arf -n

                # Run quantification with most recent quantifier.pl script
                perl ~/mirdeep2_core-master/quantifier.pl -p "${myMIR}"/hairpindna.fa -m "${myMIR}"/mature.fa -r "${mySIZE}".mapped.fa.col -t mmu -y "${mySIZE}"v2

done < $FILENAME

