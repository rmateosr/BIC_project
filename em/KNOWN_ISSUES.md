# Known issues — `em/` (constrained / diploid-aware EM+BIC driver)

Issues discovered against this driver that are **real** (not migration-cosmetic) and warrant a dedicated investigation+fix pass. Each entry: symptom, evidence, impact, location, suggested fix.

---

## 1. `total_reads` output column counts REGION rows, not WINDOW rows

**Status:** **resolved 2026-05-12** — fix log in [`em/FIX_LOG_window_n_count.md`](FIX_LOG_window_n_count.md). The original diagnosis below overstated the bug: investigation showed the BIC penalty was already using per-window `n`. Only the **reported `total_reads` column** was inflated. Kept here for traceability.

### Symptom

The `total_reads` column emitted by `EMBIC_bin_path_timebased` is roughly constant within a region (~5,277 in this dataset) regardless of which 10-CpG window inside the region is being evaluated. This is biologically impossible — a 10-CpG window spans ~500–2,000 bp on chr22, so for nanopolish-era shTET3 coverage the per-window depth should be on the order of tens of reads, not thousands.

### Root cause (confirmed)

`EMBIC_bin_path_timebased` loads the entire region file into a matrix `X_full[t]` of shape `(n_reads_in_region, n_CpGs_in_region)`. The sliding window only slices the CpG (column) axis:

```python
# em/EMfunctions_timebased_diploid_aware_model.py:117
X_w = X_full[t][:, cont:cont + windowsize]
```

So `X_window[t].shape[0] == X_full[t].shape[0]` for every window — it's always the **entire region's** read population, regardless of how many of those reads have a non-NaN call at the 10 CpGs of the current window.

Both `total_reads` (line 124) and `n_total` (line 237) inherit this inflated count:

```python
# em/EMfunctions_timebased_diploid_aware_model.py:124
total_reads += X_window[t].shape[0]

# em/EMfunctions_timebased_diploid_aware_model.py:237-239
n_total = sum(X_window[t].shape[0] for t in X_window)
p_alt = 4 * J + 2 * T - 2
bic_alt = p_alt * np.log(n_total) - 2 * ll_alt
```

The same pattern likely affects `compute_null_BIC` and `compute_M1_BIC` (need to verify — they take `X_window` and presumably compute their own `n` the same way).

### Evidence

On chr22 smoke run (job `121954972`, 60 regions, 188,640 reads total):

| Metric | Value |
|--------|-------|
| Reads on chr22 (weekinit) | 132,827 unique |
| Reads on chr22 (week20) | 55,813 unique |
| Regions after `divide` | 60 |
| Avg reads per region (one sample) | ~3,150 |
| Avg reads per region (both samples summed) | ~6,300 |
| Observed median `total_reads` across 45,879 windows | **5,277** |
| Observed p25 / p75 | 2,703 / 5,340 |
| Observed min / max | 2,703 / 6,516 |

The region file `methylationfraction_10513852_11923505_.tsv` (1.4 Mb span) contains 3,642 unique reads from weekinit alone — matching the per-region weekinit contribution. Adding week20 brings the sum to ~5,300, matching the observed median exactly.

In contrast, **realistic per-window depth** on shTET3 nanopolish data with ~14× weekinit + ~11× week20 coverage would be roughly **20–60 reads** for an average 10-CpG window — two orders of magnitude lower than what's reported.

### Impact

1. **`total_reads` is a misleading output column.** Users reading the TSV (or downstream filters that gate on read depth) interpret it as window depth and act on inflated numbers.

2. **BIC penalty term is too harsh.** Every model's BIC is `p * log(n) - 2 * ll`. With inflated `n ≈ 5,000` instead of true `n_window ≈ 30–60`, the penalty per parameter is `log(5000) ≈ 8.5` rather than `log(50) ≈ 3.9`. Each extra parameter in M1 vs M0 (or M2 vs M1) has to overcome ~2× more BIC penalty than it should.

3. **Bias toward M0 (no ASM).** Because M1 has `(J + 2T)` more parameters than M0 and M2 has another `2T-2` over M1, the inflated penalty preferentially crushes M1 and M2. The first 50k+ windows in the chr22 smoke output all showed `BIC_3way_winner = 0` — consistent with this bias.

4. **Cross-region comparisons are slightly off.** Because regions vary in length, `n_total` varies in proportion (2,703 → 6,516 observed) — the BIC penalty is harsher in larger regions than in smaller ones for the same per-window biological signal.

### Suggested fix

Define `n_effective_window` per window per timepoint as the number of rows of `X_window[t]` that have at least one non-NaN entry across the 10 CpGs:

```python
n_effective_t = int(np.sum(~np.isnan(X_w).all(axis=1)))
```

Use `sum(n_effective_t for t)` as both:
- the reported `total_reads` column, and
- the `n_total` argument to every BIC formula in this module.

Verify the same change is applied consistently in:
- `_compute_window_BIC` (line 237)
- `compute_null_BIC`
- `compute_M1_BIC`
- the single-window path `_process_single_window` (line 178)

### What to NOT do

- Don't replace `X_window[t].shape[0]` blindly. The EM itself (`EM`, `doE`, `doM_*`) needs the full matrix including all-NaN rows so its E-step probabilities are consistent — pre-filtering rows inside the EM would break the constraint propagation. Filter only when computing `n` for BIC and for the reported column.

- Don't change the unconstrained driver in `em_unconstrained/` to match — that's a frozen ablation and changing it would invalidate the manuscript A/B comparison.

### Repro

```bash
cd <repo>
# Inspect a region file (1.4 Mb of chr22 = 3,642 weekinit reads)
wc -l smokes/run_smoke4_full_pipeline/weekinit/read_format_split/chr22/methylationfraction_10513852_*.tsv
# Inspect the corresponding EM output windows — total_reads is ~5,277 across all of them
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="total_reads") c=i; next} $1 < 11923505 {print $1, $c}' \
    smokes/run_smoke4_full_pipeline/timebased/chr22_timebased_ASM.tsv.tmp | head
```
