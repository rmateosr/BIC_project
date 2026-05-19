#!/bin/bash
# ABOUTME: SGE inner job — runs nanopolish_to_mergedbychr.py for one chromosome.
# ABOUTME: Produces modkit_referenced_splitbam_mergedbychr/modkit_extract_<chr>_merged_filtered_coordmod.tsv.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=16G
#$ -pe def_slot 1
set -euo pipefail

module use /usr/local/package/modulefiles/
module load samtools/1.19
module load python/3.12.0

CHR=$1
METH_TSV_GZ=$2
BAM=$3
OUTPUT_FOLDER=$4

# Under SGE, BASH_SOURCE[0] resolves to the spool path — use SGE_O_WORKDIR instead.
# The submitter (nanopolish_preprocess_outer.sh) runs from UPSTREAM_DIR with -cwd,
# so SGE_O_WORKDIR = UPSTREAM_DIR here.
UPSTREAM_DIR="${SGE_O_WORKDIR:-$(pwd)}"

echo "[$(date)] nanopolish_preprocess_inner: chr=${CHR} tsv=${METH_TSV_GZ} bam=${BAM} out=${OUTPUT_FOLDER}"
python3 "${UPSTREAM_DIR}/scripts/nanopolish_to_mergedbychr.py" "${CHR}" "${METH_TSV_GZ}" "${BAM}" "${OUTPUT_FOLDER}"
echo "[$(date)] done."
