#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=64G
set -xv
set -o errexit
set -o nounset

module use /usr/local/package/modulefiles/
module load samtools/1.19
module load apptainer/
chr=$1
HP=$2
OUTPUT_FOLDER=$3
REF=$4
# Default to the modkit binary bundled in this repo (<UPSTREAM_DIR>/bin/modkit).
# Under SGE, $0 is the spooled job-script path (/var/spool/ge/.../job_scripts/<JID>),
# so $(readlink -f "$0") does NOT point back into the repo — it resolves to the spool
# dir and derives a path like /var/spool/ge/.../bin/modkit that doesn't exist.
# Use $SGE_O_WORKDIR instead: SGE auto-sets it to the submitter's cwd for every job
# (no -V required). The submitter (modkitouter) runs with -cwd from UPSTREAM_DIR, so
# $SGE_O_WORKDIR = UPSTREAM_DIR here. Fall back to $(pwd) for non-SGE invocations.
# Override with MODKIT_BIN env var if needed.
MODKIT_BIN=${MODKIT_BIN:-${SGE_O_WORKDIR:-$(pwd)}/bin/modkit}
# modkit v0.4.1 moved read-level extraction under the `calls` subcommand.
# Top-level `modkit extract` no longer accepts flags directly — it requires a
# subcommand (full | calls). `calls` preserves the --cpg / --mapped-only /
# --reference flags and emits a read-level TSV with the same column layout
# downstream filters expect ($3 = ref_position, $12 = call_code).
# --threads 1 because the SGE header allocates only 1 slot; modkit's default
# of 4 threads would oversubscribe and risk tripping queue policy.
${MODKIT_BIN} extract calls --threads 1 --cpg --mapped-only --reference ${REF} ${OUTPUT_FOLDER}/splitbam_bychrandhap/${chr}_${HP}.bam ${OUTPUT_FOLDER}/modkit_splitbam_referenced_bychrandhap/modkit_extract_${chr}_${HP}.tsv




