#!/bin/bash
# ABOUTME: Time-based ASM pipeline starting from nanopolish methylation TSV files.
# ABOUTME: TSV + phased BAM -> per-read tables -> union coords -> paired regions -> EM+BIC.
set -euo pipefail

# ============================================================
# Usage
# ============================================================
# bash run_pipeline.sh \
#     <SAMPLE1_ID> <SAMPLE1_METH_TSV> <SAMPLE1_PHASED_BAM> <SAMPLE1_OUT> \
#     <SAMPLE2_ID> <SAMPLE2_METH_TSV> <SAMPLE2_PHASED_BAM> <SAMPLE2_OUT> \
#     <REF> <TIMEBASED_OUT>
#
# SAMPLE_METH_TSV: tabix-indexed nanopolish methylation calls (.tsv.gz + .tbi)
# SAMPLE_PHASED_BAM: Clair3+WhatsHap output with HP tags (for haplotype assignment)
# REF: reference genome FASTA (passed through to downstream scripts)
#
# Optional env overrides:
#   CHROMS      — space-separated chromosomes (default: chr1..22 chrX chrY)
#   WINDOWSIZE  — CpG sites per EM window (default: 10)
#   WORKERS     — parallel processes for EM driver (default: 10)
#   RUN_LABEL   — prefix for SGE job names (default: tsv)
#
# Example (chr22 smoke test):
#   CHROMS="chr22" RUN_LABEL=smoke \
#   bash run_pipeline.sh \
#       weekinit input_data/methylation_calls_weekinit_okada.tsv.gz \
#                /path/to/weekinit/phased_output.bam  output/weekinit \
#       week20   input_data/methylation_calls_okada_week20.tsv.gz \
#                /path/to/week20/phased_output.bam    output/week20 \
#       /path/to/ref.fasta output/timebased
# ============================================================

if [ "$#" -ne 10 ]; then
    cat <<EOF
Usage: $0 \\
    <SAMPLE1_ID> <SAMPLE1_METH_TSV> <SAMPLE1_PHASED_BAM> <SAMPLE1_OUT> \\
    <SAMPLE2_ID> <SAMPLE2_METH_TSV> <SAMPLE2_PHASED_BAM> <SAMPLE2_OUT> \\
    <REF> <TIMEBASED_OUT>
EOF
    exit 2
fi

SAMPLE1_ID=$1
SAMPLE1_METH_TSV=$2
SAMPLE1_BAM=$3
SAMPLE1_OUT=$4
SAMPLE2_ID=$5
SAMPLE2_METH_TSV=$6
SAMPLE2_BAM=$7
SAMPLE2_OUT=$8
REF=$9
TIMEBASED_OUT=${10}

export CHROMS="${CHROMS:-$(seq -f "chr%g" 1 22) chrX chrY}"
export WINDOWSIZE="${WINDOWSIZE:-10}"
export WORKERS="${WORKERS:-10}"
export UNION_COORDS_DIR="${TIMEBASED_OUT}/union_coords"
RUN_LABEL="${RUN_LABEL:-tsv}"

# Break auto-chain: divide waits on union-coords build (coordinator handles this).
export CHAIN_DIVIDE=0

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export UPSTREAM_DIR="${SCRIPT_DIR}"
TIMEBASED_DIR="${SCRIPT_DIR}/timebased_ASM"
# Default driver: the diploid-aware (constrained) model — production since 2026-04-21.
# Override with TIMEBASED_DRIVER=…/em_unconstrained/BIC_algorithm_timebased.py for the
# frozen unconstrained ablation (manuscript A/B comparison only).
TIMEBASED_DRIVER="${TIMEBASED_DRIVER:-${TIMEBASED_DIR}/BIC_algorithm_timebased_diploid_aware_model.py}"

# ------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------
for f in "${SAMPLE1_METH_TSV}" "${SAMPLE2_METH_TSV}" \
         "${SAMPLE1_BAM}" "${SAMPLE2_BAM}" "${REF}"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: file not found: ${f}" >&2
        exit 1
    fi
done
for f in "${SAMPLE1_METH_TSV}" "${SAMPLE2_METH_TSV}"; do
    if [ ! -f "${f}.tbi" ]; then
        echo "ERROR: tabix index not found: ${f}.tbi" >&2
        exit 1
    fi
done
if [ ! -f "${TIMEBASED_DRIVER}" ]; then
    echo "ERROR: EM driver not found: ${TIMEBASED_DRIVER}" >&2
    exit 1
fi

# ------------------------------------------------------------
# Set up per-sample output dirs with phased BAM symlinks
# ------------------------------------------------------------
setup_sample_dir () {
    local out=$1 bam=$2
    mkdir -p "${out}" "${out}/log"
    ln -sfn "$(readlink -f "${bam}")" "${out}/phased_output.bam"
    if [ -f "${bam}.bai" ]; then
        ln -sfn "$(readlink -f "${bam}.bai")" "${out}/phased_output.bam.bai"
    elif [ -f "${bam%.bam}.bai" ]; then
        ln -sfn "$(readlink -f "${bam%.bam}.bai")" "${out}/phased_output.bam.bai"
    else
        echo "WARN: no .bai index alongside ${bam}" >&2
    fi
}

setup_sample_dir "${SAMPLE1_OUT}" "${SAMPLE1_BAM}"
setup_sample_dir "${SAMPLE2_OUT}" "${SAMPLE2_BAM}"
mkdir -p "${TIMEBASED_OUT}" "${UNION_COORDS_DIR}" "${SCRIPT_DIR}/log"

# ------------------------------------------------------------
# Submit nanopolish preprocessing + merged2reads chain per sample
# ------------------------------------------------------------
pushd "${SCRIPT_DIR}" >/dev/null

S1_ROOT_JOB="${RUN_LABEL}_${SAMPLE1_ID}_npp_outer"
S2_ROOT_JOB="${RUN_LABEL}_${SAMPLE2_ID}_npp_outer"

qsub -V -N "${S1_ROOT_JOB}" scripts/nanopolish_preprocess_outer.sh \
    "${SAMPLE1_OUT}" "${REF}" "${RUN_LABEL}_${SAMPLE1_ID}" "${UNION_COORDS_DIR}" \
    "${SAMPLE1_METH_TSV}" "${SAMPLE1_BAM}"

qsub -V -N "${S2_ROOT_JOB}" scripts/nanopolish_preprocess_outer.sh \
    "${SAMPLE2_OUT}" "${REF}" "${RUN_LABEL}_${SAMPLE2_ID}" "${UNION_COORDS_DIR}" \
    "${SAMPLE2_METH_TSV}" "${SAMPLE2_BAM}"

popd >/dev/null

# ------------------------------------------------------------
# Submit coordinator job — polls for m2r outputs, then submits
# union_coords -> divide -> time-based EM.
# ------------------------------------------------------------
export RUN_LABEL SAMPLE1_ID SAMPLE2_ID SAMPLE1_OUT SAMPLE2_OUT
export REF TIMEBASED_OUT TIMEBASED_DRIVER

COORD_JOB="${RUN_LABEL}_coord"
qsub -V -N "${COORD_JOB}" \
    -hold_jid "${S1_ROOT_JOB},${S2_ROOT_JOB}" \
    -o "${SCRIPT_DIR}/log" -e "${SCRIPT_DIR}/log" \
    "${SCRIPT_DIR}/submit_downstream_after_m2r.sh"

cat <<EOF
Submitted time-based ASM pipeline (from TSV).
  sample1:        ${SAMPLE1_ID}
    meth TSV:     ${SAMPLE1_METH_TSV}
    phased BAM:   ${SAMPLE1_BAM}
    output:       ${SAMPLE1_OUT}
  sample2:        ${SAMPLE2_ID}
    meth TSV:     ${SAMPLE2_METH_TSV}
    phased BAM:   ${SAMPLE2_BAM}
    output:       ${SAMPLE2_OUT}
  reference:      ${REF}
  union coords:   ${UNION_COORDS_DIR}
  EM output:      ${TIMEBASED_OUT}
  chromosomes:    ${CHROMS}
  windowsize:     ${WINDOWSIZE}
  workers:        ${WORKERS}
  run label:      ${RUN_LABEL}

Monitor:  qstat -u \$USER | grep ${RUN_LABEL}
EOF
