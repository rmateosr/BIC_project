#!/bin/bash
# ABOUTME: SGE outer — submits per-chr nanopolish preprocessing then chains to merged2reads.
# ABOUTME: Replaces bamsplitterouter.sh (stages 1-5) for nanopolish-era samples without MM/ML tags.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=2G
set -euo pipefail

# Usage: qsub nanopolish_preprocess_outer.sh <OUTPUT_FOLDER> <REF> <SAMPLE_ID> <SHARED_COORDS_DIR> <METH_TSV_GZ> <BAM>
# CHROMS env var controls which chromosomes to process (default: chr1..22 chrX chrY).

OUTPUT_FOLDER=$1
REF=$2
SAMPLE_ID=${3:-s}
SHARED_COORDS_DIR=${4:-}
METH_TSV_GZ=$5
BAM=$6

CHROMS=${CHROMS:-$(seq -f "chr%g" 1 22) chrX chrY}
chromosomes=(${CHROMS})

echo "[$(date)] nanopolish_preprocess_outer: sample=${SAMPLE_ID} chroms=${CHROMS}"

mkdir -p "${OUTPUT_FOLDER}/modkit_referenced_splitbam_mergedbychr"
mkdir -p "${OUTPUT_FOLDER}/read_format"

PREPROC_JOB="${SAMPLE_ID}_npp"
for chr in "${chromosomes[@]}"; do
    qsub -N "${PREPROC_JOB}" \
        scripts/nanopolish_preprocess_inner.sh "${chr}" "${METH_TSV_GZ}" "${BAM}" "${OUTPUT_FOLDER}"
done

# Chain merged2reads_outer.sh after preprocessing completes.
# It submits ${SAMPLE_ID}_m2r per chr, which the coordinator polls for.
qsub -V -N "${SAMPLE_ID}_m2r_outer" -hold_jid "${PREPROC_JOB}" \
    scripts/merged2reads_outer.sh "${OUTPUT_FOLDER}" "${REF}" "${SAMPLE_ID}" "${SHARED_COORDS_DIR}"

echo "[$(date)] submitted ${#chromosomes[@]} preprocess + m2r chain for ${SAMPLE_ID}"
