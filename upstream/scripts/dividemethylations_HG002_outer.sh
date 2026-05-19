#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=4G
#$ -pe def_slot 1
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
# NOTE: rm -rf of merged modkit output disabled — intermediates kept for validation/resume.
#rm -rf ${OUTPUT_FOLDER}/modkit_referenced_splitbam_mergedbychr
CHROMS=${CHROMS:-$(seq -f "chr%g" 1 22) chrX chrY}
chromosomes=(${CHROMS})
for chr in "${chromosomes[@]}"
do
 qsub -N ${SAMPLE_ID}_div scripts/dividemethylations_HG002_inner.sh ${chr} ${OUTPUT_FOLDER} "${SHARED_COORDS_DIR}"
done

# Cross-sectional BIC (Step 9) auto-chain disabled — the time-based workflow runs
# its own EM driver separately after both samples complete Step 8.
#qsub -N BIC -hold_jid div scripts/qsub_chromosome_BIC_2024_HG002.sh ${OUTPUT_FOLDER} ${REF}