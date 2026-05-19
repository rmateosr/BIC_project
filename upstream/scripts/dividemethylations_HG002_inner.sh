#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l lmem,s_vmem=64G
#$ -pe def_slot 1
set -xv
set -o errexit
set -o nounset

CHR=$1
echo $CHR
OUTPUT_FOLDER=$2
UNION_COORDS_DIR=${3:-}
# If UNION_COORDS_DIR is provided, pass the per-chromosome union coords file
# so region boundaries are identical across samples (required for time-based pairing).
if [ -n "${UNION_COORDS_DIR}" ] && [ -f "${UNION_COORDS_DIR}/${CHR}_union_coords.txt" ]; then
    Rscript scripts/dividemethylationintosmallerregions_updated_HG002.R $CHR ${OUTPUT_FOLDER} "${UNION_COORDS_DIR}/${CHR}_union_coords.txt"
else
    Rscript scripts/dividemethylationintosmallerregions_updated_HG002.R $CHR ${OUTPUT_FOLDER}
fi



