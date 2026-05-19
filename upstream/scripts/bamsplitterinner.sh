#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=32G
set -xv
set -o errexit
set -o nounset
chr=$1
OUTPUT_FOLDER=$2

module use /usr/local/package/modulefiles/
module load samtools/1.19
module load apptainer/
mkdir -p ${OUTPUT_FOLDER}/splitbam/
samtools view -b ${OUTPUT_FOLDER}/phased_output.bam $chr > ${OUTPUT_FOLDER}/splitbam/phased_output.bam_$chr.bam
samtools index ${OUTPUT_FOLDER}/splitbam/phased_output.bam_$chr.bam



