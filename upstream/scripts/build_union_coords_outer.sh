#!/bin/bash
# ABOUTME: Fan-out wrapper: launches per-chromosome union-coords builder jobs.
# ABOUTME: After completion, dividemethylations_HG002_outer.sh can use UNION_COORDS_DIR.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=4G
#$ -pe def_slot 1
set -xv
set -o errexit
set -o nounset

SAMPLE1_OUT=$1
SAMPLE2_OUT=$2
UNION_COORDS_DIR=$3
UNION_JOB_NAME=${4:-union_coords}

mkdir -p "${UNION_COORDS_DIR}"

CHROMS=${CHROMS:-$(seq -f "chr%g" 1 22) chrX chrY}
chromosomes=(${CHROMS})
for chr in "${chromosomes[@]}"; do
    qsub -N "${UNION_JOB_NAME}" scripts/build_union_coords_inner.sh \
        "${chr}" "${SAMPLE1_OUT}" "${SAMPLE2_OUT}" "${UNION_COORDS_DIR}"
done
