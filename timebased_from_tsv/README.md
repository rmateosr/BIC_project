# Time-Based ASM Pipeline (from nanopolish TSV)

Detects **Allele-Specific Methylation changes between two time points** using
nanopolish methylation calls as input (not modBAM/modkit).

This is a thin wrapper around the shared pipeline code in `../upstream/` and `../em/`.
The BAM-starting variant (stages 1-5: bamsplit -> modkit -> filter -> coordfix)
is replaced by `nanopolish_to_mergedbychr.py`, which reads nanopolish TSV +
phased BAM directly.

## Data flow

```
nanopolish TSV (.tsv.gz + .tbi)  ──┐
                                    ├──► nanopolish_to_mergedbychr.py ──► coordmod.tsv
phased BAM (Clair3+WhatsHap)    ──┘         (per chr, per sample)
                                                     │
                                    ┌────────────────┘
                                    ▼
                              merged2reads (R)         ──► read_format/<chr>_reads.tsv
                                    │
                   ┌────────────────┤  (both samples)
                   ▼                ▼
             build_union_coords.py     ──► union_coords/<chr>_union_coords.txt
                         │
                         ▼
              divide into regions (R)  ──► read_format_split/<chr>/methylationfraction_*_.tsv
                         │
                         ▼
               time-based EM+BIC       ──► <chr>_timebased_ASM.tsv
```

## Quick start (chr22 smoke test)

```bash
cd <repo>/timebased_from_tsv

BASE=<repo>   # path to a directory holding the per-timepoint BAMs + methylation TSVs
REF=/path/to/Homo_sapiens_assembly38.fasta
OUT=output_smoke

CHROMS="chr22" RUN_LABEL=tsmk \
bash run_pipeline.sh \
    weekinit input_data/methylation_calls_weekinit_okada.tsv.gz \
             ${BASE}/weekinit/phased_output.bam  ${OUT}/weekinit \
    week20   input_data/methylation_calls_okada_week20.tsv.gz \
             ${BASE}/week20/phased_output.bam    ${OUT}/week20 \
    ${REF}   ${OUT}/timebased
```

Monitor: `qstat -u $USER | grep tsmk`

## 4-way BIC (Shared / ASM / ASM-with-time / Shared-with-time)

The default driver runs the 3-way BIC (Shared / ASM / ASM-with-time-emergence).
To run the 4-way competition that also evaluates the Shared-with-time-emergence
model (no-h symmetric drift; LaTeX: `latex/EM_stepbystep_derivation_emergent_altered_model.tex`),
set `TIMEBASED_DRIVER` to the orchestrator before launching:

```bash
TIMEBASED_DRIVER=${BASE}/em/BIC_algorithm_timebased_4way.py \
CHROMS="chr22" RUN_LABEL=tsmk4 \
bash run_pipeline.sh \
    weekinit input_data/methylation_calls_weekinit_okada.tsv.gz \
             ${BASE}/weekinit/phased_output.bam  ${OUT}/weekinit \
    week20   input_data/methylation_calls_okada_week20.tsv.gz \
             ${BASE}/week20/phased_output.bam    ${OUT}/week20 \
    ${REF}   ${OUT}/timebased
```

Output adds four columns on top of the 3-way schema: `BIC_shared_time`,
`BIC_4way_winner` (string label from `{"shared", "asm", "asm_time", "shared_time"}`),
`em_iterations_shared_time`, `pi_altered_t2_shared`. The legacy 3-way columns
(`BICsinglecomp`, `BICmiddlecomp`, `BICpaircomp`, `BICresult`, `BIC_3way_winner`)
are emitted unchanged so existing downstream code keeps working.

## Directory structure

```
timebased_from_tsv/
├── run_pipeline.sh                        ← THE entry point (TSV-only)
├── submit_downstream_after_m2r.sh         → symlink (coordinator job)
├── input_data/                            → symlink (nanopolish TSV files)
├── timebased_ASM/                         → symlink (EM algorithm + tests)
├── log/                                   ← SGE logs
├── scripts/
│   ├── nanopolish_to_mergedbychr.py       → symlink (TSV preprocessor)
│   ├── nanopolish_preprocess_{inner,outer}.sh → symlinks
│   ├── merged2reads_{inner,outer}.sh      → symlinks (shared downstream)
│   ├── linebylineextractingreadinfo_HG002.R   → symlink
│   ├── build_union_coords{.py,_inner.sh}      → symlinks
│   ├── dividemethylations_HG002_inner.sh      → symlink
│   ├── dividemethylationintosmallerregions_updated_HG002.R → symlink
│   └── run_timebased_inner.sh             → symlink
└── README.md                              ← this file
```

All `→ symlink` entries point to `../upstream/scripts/`
or `../em/`. Changes to the shared code propagate automatically.

## Required inputs

| Input | Format | Example |
|-------|--------|---------|
| Nanopolish methylation TSV | `.tsv.gz` + `.tbi` (tabix-indexed) | `input_data/methylation_calls_weekinit_okada.tsv.gz` |
| Phased BAM | BAM with HP tags from Clair3+WhatsHap | `weekinit/phased_output.bam` |
| Reference genome | FASTA | `referencegenome/Homo_sapiens_assembly38.fasta` |

## Parameters (matches BIC ASM production)

| Parameter | Default | Override |
|-----------|---------|----------|
| Window size | 10 CpGs | `WINDOWSIZE=N` |
| Workers | 10 | `WORKERS=N` |
| Chromosomes | chr1-22 + chrX chrY | `CHROMS="chr22"` |
| Run label | `tsv` | `RUN_LABEL=name` |
| Pseudo count | 1e-10 | (in EMfunctions_timebased.py) |
| Max iterations | 1000 | (in EMfunctions_timebased.py) |
| Convergence tol | 1e-8 | (in EMfunctions_timebased.py) |

## Relationship to the rest of `BIC_project/`

This directory contains **no duplicated code**. The runner (`run_pipeline.sh`) is
the only new file; everything else is a symlink. The parent `em/` driver and
the `upstream/` script collection support both BAM-starting (modkit) and
TSV-starting (nanopolish) modes via conditional logic. This directory provides
a clean, TSV-only entry point.

| Component | Source |
|-----------|--------|
| EM algorithm | `em/` (reached as `timebased_ASM/`) |
| Upstream scripts | `upstream/scripts/` |
| Coordinator | `submit_downstream_after_m2r.sh` (at BIC_project root) |
| Input data | `input_data/` (at BIC_project root, symlinked to source) |
