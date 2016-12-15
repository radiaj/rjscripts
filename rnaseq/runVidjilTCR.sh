#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=47:00:00
#PBS -A bap-052-aa
#PBS -o outputfileB
#PBS -e errorfileB
#PBS -V
#PBS -N vidjilRunTCR002

module add gcc/4.9.1
module add clang/3.5.0

FILENAME=$1

cd /gs/project/bap-052-aa/Park_data/RNAseq/vidjil-2016.03/
while read mySAMPLE
        do
                let count++
                echo "$count $mySAMPLE"

                                cd /gs/project/bap-052-aa/Park_data/RNAseq/vidjil-2016.03/

                                ./vidjil -g ../vidjil-2015.12/germline/germlines-TR.data -i -w 30 -r 1 -U -y 1000 -z 1000 ../VIDJIL_SEGVDJ_FA/"${mySAMPLE}".segmented.vdj.fa > "${mySAMPLE}".stdout

                done < $FILENAME

echo -e "\nTotal $count Lines read"
