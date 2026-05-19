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
# NOTE: rm -rf of upstream modkit output disabled — intermediates kept for validation/resume.
#rm -rf ${OUTPUT_FOLDER}/modkit_splitbam_referenced_bychrandhap
mkdir -p ${OUTPUT_FOLDER}/modkit_referenced_splitbam_mergedbychr
CHROMS=${CHROMS:-$(seq -f "chr%g" 1 22) chrX chrY}
chromosomes=(${CHROMS})
for chr in "${chromosomes[@]}"
do
    qsub -N ${SAMPLE_ID}_fixer scripts/coordfixer_inner.sh $chr ${OUTPUT_FOLDER}
done

qsub -V -N ${SAMPLE_ID}_m2r -hold_jid ${SAMPLE_ID}_fixer scripts/merged2reads_outer.sh ${OUTPUT_FOLDER} ${REF} ${SAMPLE_ID} "${SHARED_COORDS_DIR}"

