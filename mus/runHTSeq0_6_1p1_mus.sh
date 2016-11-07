
#!/bin/bash

# Add the modules needed for the analysis
module add HTSeq/0.6.1p1
module add SAMtools

# Variables
FILENAME=$1
myGTF="/RQexec/johnsonr/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"
count=0

while read mySAMPLE
        do
                let count++
                echo "$count $mySAMPLE"

                # Go to the STAR directory in the PARK_LAB folder
                cd $SCRATCH/STAR/"${mySAMPLE}"_STAR/

                python -m HTSeq.scripts.count -f sam "${mySAMPLE}"Aligned.out.sam "${myGTF}" >  "${mySAMPLE}".cnts


done < $FILENAME

echo -e "\nTotal $count lines read"
