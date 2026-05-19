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
# NOTE: rm -rf of filtered modkit output disabled — intermediates kept for validation/resume.
#rm -rf ${OUTPUT_FOLDER}/modkit_splitbam_referenced_bychrandhap_filtered
mkdir -p ${OUTPUT_FOLDER}/read_format

CHROMS=${CHROMS:-$(seq -f "chr%g" 1 22) chrX chrY}
chromosomes=(${CHROMS})
for chr in "${chromosomes[@]}"
do
 qsub -N ${SAMPLE_ID}_m2r scripts/merged2reads_inner.sh ${chr} ${OUTPUT_FOLDER}
done

# Divide step chaining.
# Time-based mode: the runner must supply CHAIN_DIVIDE=0 so that divide waits for the
# shared-coords build across BOTH samples (scheduled separately by the runner).
# Cross-sectional / single-sample mode: CHAIN_DIVIDE=1 (default) keeps the original
# auto-chain behavior using per-sample CpG coords.
CHAIN_DIVIDE=${CHAIN_DIVIDE:-1}
if [ "${CHAIN_DIVIDE}" = "1" ]; then
    qsub -V -N ${SAMPLE_ID}_div -hold_jid ${SAMPLE_ID}_m2r scripts/dividemethylations_HG002_outer.sh ${OUTPUT_FOLDER} ${REF} ${SAMPLE_ID} "${SHARED_COORDS_DIR}"
else
    echo "CHAIN_DIVIDE=0 — runner will schedule divide after shared-coords build."
fi
