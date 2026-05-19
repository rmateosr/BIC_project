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
OUTPUT_FOLDER=$2
samtools view -h ${OUTPUT_FOLDER}/splitbam/phased_output.bam_${chr}.bam | awk '($0 ~ /^@/ || !($0 ~ /HP:i:1/ || $0 ~ /HP:i:2/))  {print}' | samtools view -b -o ${OUTPUT_FOLDER}/splitbam_bychrandhap/${chr}_noH.bam
samtools index ${OUTPUT_FOLDER}/splitbam_bychrandhap/${chr}_noH.bam


