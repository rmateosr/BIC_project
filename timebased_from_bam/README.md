# Time-Based ASM Pipeline (from phased BAM)

Detects **Allele-Specific Methylation changes between two time points** using
phased BAM files with modification tags (MM/ML) as input.

This is a thin wrapper around the shared pipeline code in `../upstream/` and `../em/`.
Runs the full 9-step chain: bamsplit -> modkit -> filter -> coordfix -> m2r ->
union_coords -> divide -> EM+BIC.

## Data flow

```
phased BAM (with HP + MM/ML tags)
    │
    ├──► bamsplit by chr
    ├──► bamsplit by haplotype (H1/H2/noH)
    ├──► modkit extract (per chr × hap)
    ├──► filter low-confidence calls
    ├──► coordfixer + merge haplotypes   ──► coordmod.tsv
    │
    ▼
merged2reads (R)                         ──► read_format/<chr>_reads.tsv
    │
    ├──  (both samples)
    ▼
build_union_coords.py                    ──► union_coords/<chr>_union_coords.txt
    │
    ▼
divide into regions (R)                  ──► read_format_split/<chr>/methylationfraction_*_.tsv
    │
    ▼
time-based EM+BIC                        ──► <chr>_timebased_ASM.tsv
```

## Quick start (chr22 smoke test)

Uses the HG008 PacBio HiFi Revio normal/tumor pair (already phased, with
MM/ML 5mC tags) as a real modBAM input.

```bash
cd <repo>/timebased_from_bam

TN=/path/to/TN_pancreas
NBAM=${TN}/Normal_HiFi_GRCh38/HG008-N-P_PacBio-HiFi-Revio_20240125_35x_GRCh38-GIABv3.bam
TBAM=${TN}/Tumor_HiFi_GRCh38/HG008-T_PacBio-HiFi-Revio_20240125_116x_GRCh38-GIABv3.bam
REF=${TN}/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
OUT=output_smoke

CHROMS="chr22" RUN_LABEL=bsmk \
bash run_pipeline.sh \
    normal ${NBAM}  ${OUT}/normal \
    tumor  ${TBAM}  ${OUT}/tumor  \
    ${REF} ${OUT}/timebased
```

Monitor: `qstat -u $USER | grep bsmk`

Note: the two "time points" here are normal vs tumor states, not literal time
points — the time-based EM treats any pair of samples symmetrically, so a
normal/tumor pair is a valid input. Use a true longitudinal pair (e.g.
treatment day 0 vs day N) when interpreting results as temporal change.

## Required inputs

| Input | Format | Notes |
|-------|--------|-------|
| Phased BAM | BAM with HP + MM/ML tags | Must be from modBAM-aware basecaller (Dorado, not old Guppy) |
| Reference genome | FASTA | GRCh38 for this project |

## Directory structure

```
timebased_from_bam/
├── run_pipeline.sh                         ← THE entry point (BAM-only)
├── submit_downstream_after_m2r.sh          → symlink (coordinator job)
├── bamsplitterouter.sh                     → symlink (BAM split entry)
├── bin/modkit                              → symlink (bundled v0.4.1)
├── timebased_ASM/                          → symlink (EM algorithm + tests)
├── log/                                    ← SGE logs
├── scripts/
│   ├── bamsplitter*.sh                     → symlinks (BAM-specific stages)
│   ├── modkit*.sh                          → symlinks (modkit stages)
│   ├── coordfixer*.{sh,py}                 → symlinks (coord fix stage)
│   ├── merged2reads_{inner,outer}.sh       → symlinks (shared downstream)
│   ├── linebylineextractingreadinfo_HG002.R    → symlink (shared)
│   ├── build_union_coords{.py,_inner.sh}       → symlink (shared)
│   ├── dividemethylations_HG002_inner.sh       → symlink (shared)
│   ├── dividemethylationintosmallerregions_updated_HG002.R → symlink (shared)
│   └── run_timebased_inner.sh              → symlink (shared)
└── README.md                               ← this file
```

All symlinks point to `../upstream/` (or `../em/` for `timebased_ASM/`).

## Data note: bundled `weekinit/`/`week20/` BAMs are NOT compatible with this pipeline

The BAMs symlinked at `BIC_project/weekinit/phased_output.bam` and
`BIC_project/week20/phased_output.bam` are Guppy/nanopolish-era — they carry
HP haplotype tags but **no MM/ML methylation tags**. `modkit extract calls`
will produce empty TSVs on them. The `run_pipeline.sh` startup guard will
abort immediately if you point this pipeline at them.

For the shTET3 weekinit/week20 longitudinal data, use **`timebased_from_tsv/`**
(it consumes the nanopolish `.tsv.gz` calls in `input_data/` directly and
bypasses modkit).

Use `timebased_from_bam/` only with modBAM input (Dorado, PacBio jasmine).

## Parameters (matches BIC ASM production)

| Parameter | Default | Override |
|-----------|---------|----------|
| Window size | 10 CpGs | `WINDOWSIZE=N` |
| Workers | 10 | `WORKERS=N` |
| Chromosomes | chr1-22 + chrX chrY | `CHROMS="chr22"` |
| Run label | `bam` | `RUN_LABEL=name` |
| Pseudo count | 1e-10 | (in EMfunctions_timebased.py) |
| Max iterations | 1000 | (in EMfunctions_timebased.py) |
| Convergence tol | 1e-8 | (in EMfunctions_timebased.py) |

## When to use this vs timebased_from_tsv/

| Use this (from BAM) | Use timebased_from_tsv/ |
|---------------------|------------------------|
| BAM has MM/ML tags (modBAM format) | Pre-modBAM data (nanopolish-era) |
| Basecalled with Dorado | Basecalled with old Guppy |
| Methylation calls embedded in BAM | Methylation calls in separate TSV |
