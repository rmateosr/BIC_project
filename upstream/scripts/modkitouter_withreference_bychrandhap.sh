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
mkdir -p ${OUTPUT_FOLDER}/modkit_splitbam_referenced_bychrandhap
# NOTE: rm -rf of splitbam/ disabled — intermediates kept for validation/resume.
# Re-enable selectively after smoke test passes to save Lustre space.
#rm -rf ${OUTPUT_FOLDER}/splitbam
CHROMS=${CHROMS:-$(seq -f "chr%g" 1 22) chrX chrY}
chromosomes=(${CHROMS})
for chr in "${chromosomes[@]}"
do
    qsub -N ${SAMPLE_ID}_modkitref scripts/modkitinner_withreference_bychrandhap.sh $chr H1 ${OUTPUT_FOLDER} ${REF}
    qsub -N ${SAMPLE_ID}_modkitref scripts/modkitinner_withreference_bychrandhap.sh $chr H2 ${OUTPUT_FOLDER} ${REF}
    qsub -N ${SAMPLE_ID}_modkitref scripts/modkitinner_withreference_bychrandhap.sh $chr noH ${OUTPUT_FOLDER} ${REF}
done
qsub -V -N ${SAMPLE_ID}_filter -hold_jid ${SAMPLE_ID}_modkitref scripts/modkitfilterouter_bychrandhap.sh ${OUTPUT_FOLDER} ${REF} ${SAMPLE_ID} "${SHARED_COORDS_DIR}"
