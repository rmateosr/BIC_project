#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=4G
set -xv
set -o errexit
set -o nounset

module use /usr/local/package/modulefiles/
module load samtools/1.19
module load apptainer/
OUTPUT_FOLDER=$1
REF=$2
SAMPLE_ID=${3:-s}
SHARED_COORDS_DIR=${4:-}
mkdir -p ${OUTPUT_FOLDER}/splitbam_bychrandhap
CHROMS=${CHROMS:-$(seq -f "chr%g" 1 22) chrX chrY}
chromosomes=(${CHROMS})
for chr in "${chromosomes[@]}"
do
    qsub -N ${SAMPLE_ID}_bamsplit2 scripts/bamsplitter_part2_nowbyhap_inner_H1.sh $chr ${OUTPUT_FOLDER}
    qsub -N ${SAMPLE_ID}_bamsplit2 scripts/bamsplitter_part2_nowbyhap_inner_H2.sh $chr ${OUTPUT_FOLDER}
    qsub -N ${SAMPLE_ID}_bamsplit2 scripts/bamsplitter_part2_nowbyhap_inner_noH.sh $chr ${OUTPUT_FOLDER}
done
qsub -V -N ${SAMPLE_ID}_modkitref -hold_jid ${SAMPLE_ID}_bamsplit2 scripts/modkitouter_withreference_bychrandhap.sh ${OUTPUT_FOLDER} ${REF} ${SAMPLE_ID} "${SHARED_COORDS_DIR}"


