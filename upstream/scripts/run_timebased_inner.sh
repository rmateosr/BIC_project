#!/bin/bash
# ABOUTME: Per-chromosome SGE job — runs time-based EM+BIC over matched regions from two samples.
# ABOUTME: Invokes the driver passed as 5th arg (em/ for production, em_unconstrained/ for ablation); 10 workers.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=8G
#$ -pe def_slot 10
set -xv
set -o errexit
set -o nounset

module use /usr/local/package/modulefiles/
module load python/3.12.0

# Pin each worker's BLAS backend to 1 thread.
# Rationale: with WORKERS=10 and BLAS defaults (one pool per worker, each ~= physical cores),
# the process would spawn >>10 threads and contend with itself. Empirical A/B on a chr22 subset
# (2026-04-20, jobs 121246279 vs 121246317) showed math-safe behaviour (BICsinglecomp bit-identical)
# and a ~65% peak-memory reduction (9.3 GB -> 3.2 GB). Wall-time impact at 2 workers was noise;
# the memory saving and reproducibility were the real wins.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

CHR=$1
SAMPLE1_OUT=$2       # OUTPUT_FOLDER for sample 1 (time point 1)
SAMPLE2_OUT=$3       # OUTPUT_FOLDER for sample 2 (time point 2)
TIMEBASED_OUT=$4     # output directory for EM TSVs
# 5th arg (driver) is required — caller (run_pipeline.sh) always passes it.
TIMEBASED_DRIVER=$5
WINDOWSIZE=${WINDOWSIZE:-10}
WORKERS=${WORKERS:-10}

mkdir -p "${TIMEBASED_OUT}"

SAMPLE1_DIR="${SAMPLE1_OUT}/read_format_split/${CHR}"
SAMPLE2_DIR="${SAMPLE2_OUT}/read_format_split/${CHR}"

if [ ! -d "${SAMPLE1_DIR}" ] || [ ! -d "${SAMPLE2_DIR}" ]; then
    echo "ERROR: missing per-chr region dirs: ${SAMPLE1_DIR} and/or ${SAMPLE2_DIR}" >&2
    exit 1
fi

python3 "${TIMEBASED_DRIVER}" \
    "${SAMPLE1_DIR}" \
    "${SAMPLE2_DIR}" \
    --windowsize "${WINDOWSIZE}" \
    --workers "${WORKERS}" \
    --output "${TIMEBASED_OUT}/${CHR}_timebased_ASM.tsv"
