#!/bin/bash
# ABOUTME: Per-chromosome SGE job — builds union CpG coord set from two samples' read_format TSVs.
# ABOUTME: Output: ${UNION_COORDS_DIR}/${CHR}_union_coords.txt
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=8G
#$ -pe def_slot 1
set -xv
set -o errexit
set -o nounset

module use /usr/local/package/modulefiles/
module load python/3.12.0

CHR=$1
SAMPLE1_OUT=$2     # OUTPUT_FOLDER for sample 1
SAMPLE2_OUT=$3     # OUTPUT_FOLDER for sample 2
UNION_COORDS_DIR=$4

mkdir -p "${UNION_COORDS_DIR}"

python3 scripts/build_union_coords.py \
    --reads-file "${SAMPLE1_OUT}/read_format/${CHR}_reads.tsv" \
    --reads-file "${SAMPLE2_OUT}/read_format/${CHR}_reads.tsv" \
    --output "${UNION_COORDS_DIR}/${CHR}_union_coords.txt"
