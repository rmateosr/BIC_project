#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l lmem,s_vmem=64G
#$ -pe def_slot 1
set -xv
set -o errexit
set -o nounset

CHR=$1
OUTPUT_FOLDER=$2
echo $CHR
Rscript scripts/linebylineextractingreadinfo_HG002.R $CHR ${OUTPUT_FOLDER}



