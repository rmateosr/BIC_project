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
mkdir -p ${OUTPUT_FOLDER}/modkit_splitbam_referenced_bychrandhap_filtered
# NOTE: the next line targets a path the upstream never creates (stale no-op); disabled for clarity.
#rm -rf ${OUTPUT_FOLDER}/modkit_splitbam_bychrandhap
CHROMS=${CHROMS:-$(seq -f "chr%g" 1 22) chrX chrY}
chromosomes=(${CHROMS})
for chr in "${chromosomes[@]}"
do
    qsub -N ${SAMPLE_ID}_filter scripts/modkitfilterinner_bychrandhap.sh $chr H1 ${OUTPUT_FOLDER}
    qsub -N ${SAMPLE_ID}_filter scripts/modkitfilterinner_bychrandhap.sh $chr H2 ${OUTPUT_FOLDER}
    qsub -N ${SAMPLE_ID}_filter scripts/modkitfilterinner_bychrandhap.sh $chr noH ${OUTPUT_FOLDER}
done

qsub -V -N ${SAMPLE_ID}_fixer -hold_jid ${SAMPLE_ID}_filter scripts/coordfixer_outer.sh ${OUTPUT_FOLDER} ${REF} ${SAMPLE_ID} "${SHARED_COORDS_DIR}"


