#!/bin/bash
# ABOUTME: Time-based ASM pipeline starting from phased BAM files (modkit path).
# ABOUTME: BAM -> bamsplit -> modkit -> filter -> coordfix -> per-read tables -> union coords -> EM+BIC.
set -euo pipefail

# ============================================================
# Usage
# ============================================================
# bash run_pipeline.sh \
#     <SAMPLE1_ID> <SAMPLE1_PHASED_BAM> <SAMPLE1_OUT> \
#     <SAMPLE2_ID> <SAMPLE2_PHASED_BAM> <SAMPLE2_OUT> \
#     <REF> <TIMEBASED_OUT>
#
# SAMPLE_PHASED_BAM: Clair3+WhatsHap output with HP tags and MM/ML modification tags.
# REF: reference genome FASTA (required by modkit and downstream scripts).
#
# Optional env overrides:
#   CHROMS      — space-separated chromosomes (default: chr1..22 chrX chrY)
#   WINDOWSIZE  — CpG sites per EM window (default: 10)
#   WORKERS     — parallel processes for EM driver (default: 10)
#   RUN_LABEL   — prefix for SGE job names (default: bam)
#
# Example (chr22 smoke test):
#   BASE="$(cd "$(dirname "$0")/.." && pwd)"   # resolves to BIC_project/
#   REF=/path/to/Homo_sapiens_assembly38.fasta
#   CHROMS="chr22" RUN_LABEL=smoke \
#   bash run_pipeline.sh \
#       normal ${NORMAL_MODBAM} output/normal \
#       tumor  ${TUMOR_MODBAM}  output/tumor  \
#       ${REF} output/timebased
#
# NOTE: this pipeline requires modBAM input (Dorado-style MM/ML tags). The
# weekinit/ and week20/ BAMs shipped at BIC_project/ are Guppy/nanopolish-era
# and DO NOT carry MM/ML — they only work with timebased_from_tsv/. The
# startup guard below aborts if MM/ML are absent.
# ============================================================

if [ "$#" -ne 8 ]; then
    cat <<EOF
Usage: $0 \\
    <SAMPLE1_ID> <SAMPLE1_PHASED_BAM> <SAMPLE1_OUT> \\
    <SAMPLE2_ID> <SAMPLE2_PHASED_BAM> <SAMPLE2_OUT> \\
    <REF> <TIMEBASED_OUT>
EOF
    exit 2
fi

SAMPLE1_ID=$1
SAMPLE1_BAM=$2
SAMPLE1_OUT=$3
SAMPLE2_ID=$4
SAMPLE2_BAM=$5
SAMPLE2_OUT=$6
REF=$7
TIMEBASED_OUT=$8

export CHROMS="${CHROMS:-$(seq -f "chr%g" 1 22) chrX chrY}"
export WINDOWSIZE="${WINDOWSIZE:-10}"
export WORKERS="${WORKERS:-10}"
export UNION_COORDS_DIR="${TIMEBASED_OUT}/union_coords"
RUN_LABEL="${RUN_LABEL:-bam}"

# Break auto-chain: divide waits on union-coords build (coordinator handles this).
export CHAIN_DIVIDE=0

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export UPSTREAM_DIR="${SCRIPT_DIR}"
TIMEBASED_DIR="${SCRIPT_DIR}/timebased_ASM"
# Default driver: the diploid-aware (constrained) model — production since 2026-04-21.
# Override with TIMEBASED_DRIVER=…/em_unconstrained/BIC_algorithm_timebased.py for the
# frozen unconstrained ablation (manuscript A/B comparison only).
TIMEBASED_DRIVER="${TIMEBASED_DRIVER:-${TIMEBASED_DIR}/BIC_algorithm_timebased_diploid_aware_model.py}"

# modkit binary — exported so SGE inner scripts resolve correctly
# (under SGE, $0 / BASH_SOURCE[0] resolve to /var/spool/ge/.../job_scripts/).
export MODKIT_BIN="${SCRIPT_DIR}/bin/modkit"

# ------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------
for f in "${SAMPLE1_BAM}" "${SAMPLE2_BAM}" "${REF}"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: file not found: ${f}" >&2
        exit 1
    fi
done
if [ ! -x "${MODKIT_BIN}" ]; then
    echo "ERROR: modkit binary not executable: ${MODKIT_BIN}" >&2
    exit 1
fi
if [ ! -f "${TIMEBASED_DRIVER}" ]; then
    echo "ERROR: EM driver not found: ${TIMEBASED_DRIVER}" >&2
    exit 1
fi

# ------------------------------------------------------------
# modBAM guard: ensure inputs carry MM/ML methylation tags.
# Empty-tag BAMs (Guppy/nanopolish-era) silently produce empty modkit TSVs
# and crash coordfixer mid-chain. Sample the first 5000 mapped reads of each
# BAM; abort with a clear error if no MM tag is present.
# Skips the check if samtools is unavailable (warns and continues).
# Set SKIP_MODBAM_CHECK=1 to bypass entirely (eg if you know your BAM is fine
# but the first 5000 reads happen to lack MM).
# ------------------------------------------------------------
_check_modbam () {
    local bam="$1" label="$2"
    local found
    # Locally disable pipefail: head closes the pipe early once it's seen 5000
    # reads, which SIGPIPE-kills samtools. Under -o pipefail that propagates
    # as a failure even though we got the data we needed.
    set +o pipefail
    found=$(samtools view -F 260 "$bam" 2>/dev/null | head -5000 | awk '
        { for(i=12;i<=NF;i++) if($i ~ /^M[Mm]:/) { print 1; exit } }
    ')
    set -o pipefail
    if [ "${found:-0}" != "1" ]; then
        cat <<ERR >&2
ERROR: Input BAM has no MM/ML methylation tags in the first 5000 reads.
       label : ${label}
       file  : ${bam}
       timebased_from_bam REQUIRES modBAM input (Dorado MM/ML tags).
       For Guppy / nanopolish-era data, use timebased_from_tsv/ instead.
       To bypass this check, set SKIP_MODBAM_CHECK=1.
ERR
        exit 3
    fi
}

if [ "${SKIP_MODBAM_CHECK:-0}" != "1" ]; then
    if command -v samtools >/dev/null 2>&1; then
        _check_modbam "${SAMPLE1_BAM}" "${SAMPLE1_ID}"
        _check_modbam "${SAMPLE2_BAM}" "${SAMPLE2_ID}"
    else
        echo "WARN: samtools not on PATH — skipping modBAM tag check. Load samtools/1.19 before submission to enable." >&2
    fi
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
# Submit per-sample upstream chain: bamsplit -> modkit -> filter -> coordfix -> m2r
# ------------------------------------------------------------
pushd "${SCRIPT_DIR}" >/dev/null

S1_ROOT_JOB="${RUN_LABEL}_${SAMPLE1_ID}_bamsplit"
S2_ROOT_JOB="${RUN_LABEL}_${SAMPLE2_ID}_bamsplit"

qsub -V -N "${S1_ROOT_JOB}" bamsplitterouter.sh \
    "${SAMPLE1_OUT}" "${REF}" "${RUN_LABEL}_${SAMPLE1_ID}" "${UNION_COORDS_DIR}"
qsub -V -N "${S2_ROOT_JOB}" bamsplitterouter.sh \
    "${SAMPLE2_OUT}" "${REF}" "${RUN_LABEL}_${SAMPLE2_ID}" "${UNION_COORDS_DIR}"

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
Submitted time-based ASM pipeline (from BAM).
  sample1:        ${SAMPLE1_ID} -> ${SAMPLE1_OUT}
    BAM:          ${SAMPLE1_BAM}
  sample2:        ${SAMPLE2_ID} -> ${SAMPLE2_OUT}
    BAM:          ${SAMPLE2_BAM}
  reference:      ${REF}
  modkit:         ${MODKIT_BIN}
  union coords:   ${UNION_COORDS_DIR}
  EM output:      ${TIMEBASED_OUT}
  chromosomes:    ${CHROMS}
  windowsize:     ${WINDOWSIZE}
  workers:        ${WORKERS}
  run label:      ${RUN_LABEL}

Monitor:  qstat -u \$USER | grep ${RUN_LABEL}
EOF
