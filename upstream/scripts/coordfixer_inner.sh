#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=128G
set -xv
set -o errexit
set -o nounset

module use /usr/local/package/modulefiles/
module load samtools/1.19

CHR=$1
echo ${CHR}	
OUTPUT_FOLDER=$2
echo ${OUTPUT_FOLDER}	
python3.12 scripts/coordfixer.py ${CHR} ${OUTPUT_FOLDER}



