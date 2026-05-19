#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=32G
set -xv
set -o errexit
set -o nounset

module use /usr/local/package/modulefiles/
module load samtools/1.19
module load apptainer/
chr=$1
HP=$2
OUTPUT_FOLDER=$3
awk -F'\t' '$3 != "-1" && $12 == "m" {print}' ${OUTPUT_FOLDER}/modkit_splitbam_referenced_bychrandhap/modkit_extract_${chr}_${HP}.tsv > ${OUTPUT_FOLDER}/modkit_splitbam_referenced_bychrandhap_filtered/modkit_extract_${chr}_${HP}_filtered.tsv



