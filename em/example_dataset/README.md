# Example synthetic dataset — explore the generator

Generated 2026-04-27 with seed=42 to give Raul something concrete to poke at.

## Reproduce this exact dataset

```bash
cd em/
python -m synthetic.generate_synthetic \
    --n-regions 8 \
    --reads-per-t 50 \
    --cpgs-per-region 20 \
    --T 3 \
    --pi-final 0.7 \
    --class-mix 1,1,1,1 \
    --seed 42 \
    --chrom chrEXPLORE \
    --output-dir example_dataset
```

## What's in here

```
example_dataset/
├── manifest.tsv                                  # answer key (truth)
├── t1/methylationfraction_<start>_<end>_.tsv     # 8 regions, time point 1
├── t2/...                                        # 8 regions, time point 2
└── t3/...                                        # 8 regions, time point 3
```

- 8 regions, mixed across all 4 truth classes (3 NULL, 1 STATIC, 2 sym TE,
  2 ASYMMETRIC) → `manifest.tsv` tells you which is which.
- `T=3` time points (t1 = baseline; t2/t3 may carry altered signal).
- Each region has J=20 CpG sites, ~50 reads per time point.
- Per-region read files contain `H1` / `H2` / `noH` haplotype tags
  (mimics WhatsHap output).

## Per-region truth (this dataset)

| region | start | true_class | ramping_allele | θ_h0_k0 | θ_h1_k0 | π_h0_alt_t3 | π_h1_alt_t3 |
|---|---|---|---|---|---|---|---|
| 0 | 1000183 | STATIC_ASM | -1 | 0.852 | 0.063 | 0.000 | 0.000 |
| 1 | 1103055 | ASYMMETRIC | h0 | 0.818 | 0.127 | **0.350** | 0.000 |
| 2 | 1207968 | TE_sym | -1 | 0.824 | 0.114 | 0.350 | 0.350 |
| 3 | 1312645 | TE_sym | -1 | 0.945 | 0.226 | 0.350 | 0.350 |
| 4 | 1417329 | NULL | -1 | 0.310 | 0.310 | 0.000 | 0.000 |
| 5 | 1521419 | NULL | -1 | 0.531 | 0.531 | 0.000 | 0.000 |
| 6 | 1625250 | NULL | -1 | 0.346 | 0.346 | 0.000 | 0.000 |
| 7 | 1729168 | ASYMMETRIC | h0 | 0.908 | 0.198 | **0.350** | 0.000 |

> Note: `--pi-final 0.7` produces `π_alt = 0.350` in the manifest because the
> manifest reports the **joint π[h, k=1, t]** (with diploid marginal), which
> is half the **within-allele altered fraction**. So `pi_final` is the
> within-allele scale; the manifest column is `pi_final / 2`.

---

## Generator parameter cheatsheet

### Top-level knobs (`generate_synthetic_dataset` / CLI)

| CLI | Python kwarg | Default | What it controls |
|---|---|---|---|
| `--n-regions` | `n_regions` | 50 | Total regions in the dataset |
| `--reads-per-t` | `reads_per_t` | 80 | Reads per region per time point (coverage) |
| `--cpgs-per-region` | `J` | 50 | CpG sites per region |
| `--T` | `T` | 4 | Number of time points |
| `--pi-final` | `pi_final` | 0.6 | **Within-allele** altered fraction at t=T |
| `--pi-trajectory` | `pi_trajectory` | None | Custom per-t altered fraction (overrides linear ramp). e.g. `0,0.2,0.4,0.6` |
| `--pi-final-h0/h1` | `pi_final_h0/h1` | None | Per-allele final values (for asymmetric / non-uniform symmetric truth) |
| `--pi-trajectory-h0/h1` | `pi_trajectory_h0/h1` | None | Per-allele explicit trajectories |
| `--class-mix` | `class_mix` | equal | Weights for `NULL,STATIC,TIME_EMERGENT,ASYMMETRIC` (e.g. `1,0,0,4` = mostly asymmetric) |
| `--theta-noise` | `theta_noise` | 0.05 | Std-dev of per-CpG θ jitter around the per-region base value |
| `--windowsize` | `windowsize` | 10 | Sliding window size (only affects manifest's `coord_evaluated`) |
| `--chrom` | `chrom` | `chrSYN` | Chromosome name written into TSVs |
| `--seed` | `seed` | None | Reproducibility |
| `--output-dir` | `output_dir` | `synthetic_output` | Where to write `t*/` and `manifest.tsv` |

### What each truth class produces (θ structure)

| Class | θ pattern | π pattern |
|---|---|---|
| `NULL` | All four θ slots equal (random base ∈ [0.3, 0.7], jittered per CpG) | All π_alt = 0 |
| `STATIC_ASM` | θ[h=0,k=0] ∈ [0.75, 0.95]; θ[h=1,k=0] ∈ [0.05, 0.25]; θ[:,:,k=1] = random unused noise | All π_alt = 0 |
| `TIME_EMERGENT_ASM` | k=0 has the static-like allele difference; k=1 is a distinct altered pattern | Both alleles ramp π_alt 0 → `pi_final` |
| `ASYMMETRIC_TIME_EMERGENT_ASM` | Same θ structure as symmetric TE | Only the chosen `ramping_allele` ramps; other stays at 0 |

### Read-level mechanics (`generate_reads_for_region`)

Each read picks:

1. A **haplotype tag** from `(P(H1), P(H2), P(noH)) = (0.3, 0.3, 0.4)` —
   matching real WhatsHap output proportions.
2. A **true allele** `h ∈ {0, 1}` (forced by H1/H2; sampled from the joint
   `(h, k)` distribution for noH reads — important for asymmetric truth).
3. A **true component** `k ∈ {0=normal, 1=altered}` from the per-allele
   `pi_alt_t_h{h}` at this time point.
4. A **contiguous CpG span** (uniform random, ≥ `min_cpg_span` CpGs).
5. **Methylation calls** at each covered CpG, drawn `Bernoulli(θ[j, h, k])`.

### Manifest columns (answer key)

| Column | What it means |
|---|---|
| `region_idx` | 0-indexed region |
| `chrom`, `region_start`, `region_end` | Genomic coordinates |
| `coord_evaluated` | Per-window CpG position (one row per sliding window) |
| `true_class` | `NULL` / `STATIC_ASM` / `TIME_EMERGENT_ASM` / `ASYMMETRIC_TIME_EMERGENT_ASM` |
| `ramping_allele` | `0` or `1` for ASYMMETRIC; `-1` otherwise |
| `true_theta_h{0,1}_k{0,1}` | Region-level base θ (CpGs jitter around this by ±`theta_noise`) |
| `expected_BICresult` | What the legacy 2-way BIC should call (1=ASM, 0=not) |
| `pi_h{0,1}_alt_t{1..T}` | **Joint** per-allele π_alt at each time point. Half the within-allele altered fraction. Manifest reports under the diploid marginal constraint π[h,0,t] + π[h,1,t] = 0.5. |

---

## Things to try

```bash
# Inspect one region's read calls visually:
column -t -s$'\t' t2/methylationfraction_1103055_1107958_.tsv | less

# How many H1/H2/noH reads per file:
for f in t*/methylationfraction_*.tsv; do
    awk -F'\t' 'NR>1 {print $6}' "$f" | sort | uniq -c | tr '\n' ' '
    echo " $f"
done

# Run the 3-way BIC pipeline on a single region:
python -c "
import os, glob
from EMfunctions_timebased_diploid_aware_model import EMBIC_bin_path_timebased
paths = {t: f'example_dataset/t{t}/methylationfraction_1103055_1107958_.tsv' for t in (1,2,3)}
print(EMBIC_bin_path_timebased(paths, windowsize=10).to_string())
"

# Run parameter recovery on this dataset (compare M0/M1/M2 fits to truth):
python -c "
import pandas as pd
from parameter_recovery import recover_params
manifest = pd.read_csv('example_dataset/manifest.tsv', sep='\t', keep_default_na=False)
print(recover_params('example_dataset', manifest, T=3).to_string())
"

# Generate variants by tweaking just one knob:
python -m synthetic.generate_synthetic --seed 1 --T 3 --pi-final 0.9 \
    --class-mix 0,0,0,1 --n-regions 4 --output-dir asym_strong_only
python -m synthetic.generate_synthetic --seed 1 --T 5 --pi-trajectory 0,0,0.3,0.6,0.9 \
    --class-mix 0,0,1,0 --output-dir delayed_onset_TE
```
