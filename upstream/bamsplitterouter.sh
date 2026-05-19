#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=4G
set -xv
set -o errexit
set -o nounset

# Usage: qsub bamsplitterouter.sh <OUTPUT_FOLDER> <REF> [SAMPLE_ID] [SHARED_COORDS_DIR]
# SAMPLE_ID prefixes SGE job names (avoids cross-sample hold_jid collisions when two samples run concurrently).
# SHARED_COORDS_DIR (optional) is forwarded to the divide step so regions align across time points.
# Override chromosomes via env var CHROMS (e.g. CHROMS="chr22" for smoke test).

module use /usr/local/package/modulefiles/
module load samtools/1.19
module load apptainer/
OUTPUT_FOLDER=$1
REF=$2
SAMPLE_ID=${3:-s}
SHARED_COORDS_DIR=${4:-}
CHROMS=${CHROMS:-$(seq -f "chr%g" 1 22) chrX chrY}
chromosomes=(${CHROMS})
for chr in "${chromosomes[@]}"
do
    qsub -N ${SAMPLE_ID}_bamsplit scripts/bamsplitterinner.sh $chr ${OUTPUT_FOLDER}
done
qsub -V -N ${SAMPLE_ID}_bamsplit2 -hold_jid ${SAMPLE_ID}_bamsplit scripts/bamsplitter_part2_nowbyhap_outer.sh ${OUTPUT_FOLDER} ${REF} ${SAMPLE_ID} "${SHARED_COORDS_DIR}"
