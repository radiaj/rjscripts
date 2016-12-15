#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=47:00:00
#PBS -A bap-052-aa
#PBS -o out12152016
#PBS -e error121516
#PBS -V
#PBS -N vidjilRunBCR1

module add gcc/4.9.1
module add clang/3.5.0

FILENAME=$1
myDIR=$2

cd /gs/project/bap-052-aa/Park_data/RNAseq/vidjil-2016.03/
while read mySAMPLE
        do
                let count++
                echo "$count $mySAMPLE"

                                cd /gs/project/bap-052-aa/Park_data/RNAseq/vidjil-2016.03/

                                ./vidjil -g ../vidjil-2015.12/germline/germlines.data -i -w 30 -r 1 -U -y 1000 -z 1000 "${myDIR}"/"${mySAMPLE}".fa > "${mySAMPLE}".stdout

                done < $FILENAME

echo -e "\nTotal $count Lines read"


