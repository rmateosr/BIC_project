#!/bin/bash
# ABOUTME: SGE coordinator — polls until both samples' _m2r jobs exist in queue, then submits union_coords/div/tbem.
# ABOUTME: Works around SGE's behavior where -hold_jid on non-existent job names is silently ignored.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=2G
#$ -pe def_slot 1
set -euo pipefail

# Required env (forwarded via qsub -V from run_pipeline.sh):
#   RUN_LABEL, SAMPLE1_ID, SAMPLE2_ID, SAMPLE1_OUT, SAMPLE2_OUT
#   REF, UNION_COORDS_DIR, TIMEBASED_OUT, TIMEBASED_DRIVER
#   CHROMS (required), WINDOWSIZE, WORKERS
# Under SGE, BASH_SOURCE[0] resolves to the spool dir (/var/spool/ge/...).
# Use the exported UPSTREAM_DIR from the runner (preferred), or fall back to
# SGE_O_WORKDIR/upstream.
if [ -z "${UPSTREAM_DIR:-}" ]; then
    SCRIPT_DIR="${SGE_O_WORKDIR:-$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )}"
    UPSTREAM_DIR="${SCRIPT_DIR}/upstream"
fi

for v in RUN_LABEL SAMPLE1_ID SAMPLE2_ID SAMPLE1_OUT SAMPLE2_OUT \
         REF UNION_COORDS_DIR TIMEBASED_OUT TIMEBASED_DRIVER CHROMS UPSTREAM_DIR; do
    if [ -z "${!v:-}" ]; then
        echo "ERROR: env var $v not set" >&2
        exit 1
    fi
done

S1_M2R="${RUN_LABEL}_${SAMPLE1_ID}_m2r"
S2_M2R="${RUN_LABEL}_${SAMPLE2_ID}_m2r"
UNION_JOB="${RUN_LABEL}_union_coords"
S1_DIV="${RUN_LABEL}_${SAMPLE1_ID}_div"
S2_DIV="${RUN_LABEL}_${SAMPLE2_ID}_div"
TB_JOB="${RUN_LABEL}_tbem"

# Polling (Bug B fix): check for terminal-state output files rather than
# queue state. `qstat -j <name>` returns "does not exist" for both
# not-yet-submitted AND already-completed jobs — if m2r finishes between
# poll iterations, the old logic would hang forever. File-existence is
# unambiguous: the files are either there or they're not.
#
# m2r writes ${SAMPLE_OUT}/read_format/${chr}_reads.tsv. Require all
# chromosomes' files present and non-empty for both samples before
# releasing downstream stages.
MAX_ITERS=360
iter=0
echo "[$(date)] coordinator: polling for read_format/*_reads.tsv in both sample dirs..."
while :; do
    all_present=true
    for chr in ${CHROMS}; do
        for out in "${SAMPLE1_OUT}" "${SAMPLE2_OUT}"; do
            f="${out}/read_format/${chr}_reads.tsv"
            if [ ! -s "${f}" ]; then
                all_present=false
                break 2
            fi
        done
    done
    if ${all_present}; then
        break
    fi
    iter=$((iter + 1))
    if [ "${iter}" -ge "${MAX_ITERS}" ]; then
        echo "ERROR: timed out after ${MAX_ITERS} minutes waiting for m2r output files" >&2
        echo "  expected: <sample_out>/read_format/<chr>_reads.tsv for chroms: ${CHROMS}" >&2
        exit 1
    fi
    sleep 60
done
echo "[$(date)] coordinator: m2r outputs present. Submitting downstream stages."

cd "${UPSTREAM_DIR}"

for chr in ${CHROMS}; do
    qsub -V -N "${UNION_JOB}" \
        -hold_jid "${S1_M2R},${S2_M2R}" \
        scripts/build_union_coords_inner.sh \
        "${chr}" "${SAMPLE1_OUT}" "${SAMPLE2_OUT}" "${UNION_COORDS_DIR}"
done

for chr in ${CHROMS}; do
    qsub -V -N "${S1_DIV}" -hold_jid "${UNION_JOB}" \
        scripts/dividemethylations_HG002_inner.sh "${chr}" "${SAMPLE1_OUT}" "${UNION_COORDS_DIR}"
    qsub -V -N "${S2_DIV}" -hold_jid "${UNION_JOB}" \
        scripts/dividemethylations_HG002_inner.sh "${chr}" "${SAMPLE2_OUT}" "${UNION_COORDS_DIR}"
done

for chr in ${CHROMS}; do
    qsub -V -N "${TB_JOB}" -hold_jid "${S1_DIV},${S2_DIV}" \
        scripts/run_timebased_inner.sh \
        "${chr}" "${SAMPLE1_OUT}" "${SAMPLE2_OUT}" "${TIMEBASED_OUT}" "${TIMEBASED_DRIVER}"
done

echo "[$(date)] coordinator: all downstream stages submitted."
