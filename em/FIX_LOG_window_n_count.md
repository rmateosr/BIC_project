# Fix log — window `n` count in BIC / `total_reads`

Issue tracked in `em/KNOWN_ISSUES.md` §1.

## Phase 1 — Evaluation

### 1. Where does `n` come from in each model?

Traced through `em/EMfunctions_timebased_diploid_aware_model.py`:

| Code location | What it computes | Per-region or per-window? |
|---|---|---|
| `EMBIC_bin_path_timebased`, line 117–123 | `X_w = X_full[t][:, cont:cont+10]`; filters all-NaN-in-window rows → `X_window[t]` | **Per-window** (rows filtered) |
| `_compute_window_BIC`, line 237: `n_total = sum(X_window[t].shape[0] for t in X_window)` | n in BIC penalty `p_alt * log(n) - 2*ll_alt` | **Per-window** (X_window already filtered) |
| `compute_null_BIC`, line 644: `total_reads += X.shape[0]`; line 667 BIC | n in null BIC penalty | **Per-window** (receives X_window) |
| `compute_M1_BIC`, line 723: `n_total = sum(X_by_t[t].shape[0] for t in X_by_t)` | n in M1 BIC penalty | **Per-window** (receives X_window) |
| `_process_single_window`, line 178 / 213 | n for BIC and reporting | **Per-window** (region == window when J_full ≤ windowsize) |
| **`EMBIC_bin_path_timebased`, lines 164–165** | **`total_reads` output column** | **Per-region (inflated)** ← only bug |

### 2. Empirical verification

On region `chr22:10513852-11923505` (J=10,159 CpGs, 3,642 weekinit + 1,635 week20 = 5,277 reads):

- **What's reported in TSV `total_reads` column:** 5,277 (region total)
- **What's actually passed to BIC as `n`:** per-window real depth, verified by instrumenting `_compute_window_BIC`. First 5 windows had `n_total = 2`; median across 10 random windows: ~108 reads (range 17–164), ratio to region total ≈ 2%.

```
window call #0: n_total used in BIC = 2  per_t={1: 1, 2: 1}  J=10
window call #1: n_total used in BIC = 2  per_t={1: 1, 2: 1}  J=10
...
```

Empirical proof: the BIC penalty was **already using the correctly-filtered per-window n**, not the inflated 5,277.

### 3. Re-reading `KNOWN_ISSUES.md` §1

The diagnosis is correct that the reported `total_reads` column is inflated (uses `X_full` not `X_window`). But the diagnosis is **wrong** about the BIC penalty being affected:

- KNOWN_ISSUES.md claim: *"n_total inherits this inflated count [...] Each extra parameter in M1 vs M0 has to overcome ~2× more BIC penalty than it should."*
- Reality: `_compute_window_BIC` (and the null/M1 helpers it dispatches to) receive `X_window`, which is already row-filtered to per-window non-NaN reads at lines 121–123. The BIC penalty uses `log(n_per_window)`, not `log(n_per_region)`.

The "all M0 winners" pattern observed in the first 50k+ chr22 windows is therefore **not** caused by an inflated penalty. It comes from elsewhere — likely the LL deltas between M0/M1/M2 being genuinely small for low-depth windows, or a separate issue worth investigating.

### 4. Test baseline

`pytest em/tests/` → **15 passed in 95.66s**. Synthetic tests use full-coverage matrices so the row-filter is a no-op there; they will remain green under any fix that only changes the reported column.

### 5. Decision — chose (c)

Among the three options in the prompt:
- **(a) Fix n correctly in BIC** — there is nothing to fix; the BIC already uses per-window n.
- **(b) Add `n_effective_window` as a new column, leave old in place** — misleading; implies the BIC was using the wrong n when it wasn't.
- **(c) Fix only the misleading `total_reads` reporting column** — **chosen.**

The fix is therefore a small reporting-only patch:
1. Change the output `total_reads` column at line 165 to record the per-window count (already accumulated in the local `total_reads` variable at line 124, just needs to be captured per-window into an array).
2. `_process_single_window` (line 213) is already correct (region == window).
3. Correct `em/KNOWN_ISSUES.md` to reflect the actual scope.

---

## Phase 2 — Implementation

(see commit + test additions below)

## Phase 3 — Validation

### Setup

The pre-fix chr22 smoke (job `121954972`) was killed at ~63 min after I had completed Phase 1 + 2. By that time **20 of 60 regions** had been fully written to `chr22_timebased_ASM.tsv.tmp`. (The driver uses `Pool.imap_unordered`, so regions complete in arbitrary order; I confirmed completeness per region by intersecting expected `coord_evaluated` values from the union of `weekinit` + `week20` CpG coords against the .tmp contents.)

To get an apples-to-apples comparison without waiting ~2 more hours for the rest of chr22:

1. Extracted the 20 complete regions' rows (242,897 windows, 200,800 unique coords) into `smokes/run_smoke4_full_pipeline/timebased/chr22_timebased_ASM.PREFIX.tsv` as the pre-fix baseline.
2. Renamed the in-progress .tmp aside as `chr22_timebased_ASM.tsv.PREFIX_INCOMPLETE`.
3. Built `smokes/postfix_subset/{weekinit,week20}/` with symlinks to the same 20 region files.
4. Submitted `postfix_smk` (job `121959963`, SGE `def_slot 10`, `s_vmem=8G`) running the FIXED driver on those 20 regions. Wallclock 66 min.

### Comparison

Comparison is restricted to **158,703 windows whose `coord_evaluated` appears exactly once in each output** — this removes the ambiguity from ~10 kb region overlap zones where the same coord is evaluated under two different region CpG contexts.

#### `total_reads` column

| | min | median | max | distinct values |
|---|---|---|---|---|
| Pre-fix | 2,333 | 3,358 | 6,516 | **20** (one per region) |
| Post-fix | 1 | **78** | 292 | **272** (varies per window) |

The post-fix value is now the count of reads with at least one non-NaN call in the window (summed across both timepoints), as intended. Pre-fix value was per-region, repeated for every window in the region — confirmed by `distinct = 20` (one constant per region).

#### `BIC_3way_winner` flips: PRE vs POST

```
POST       0     1    2     All
PRE
0     150615   845   92  151552
1        852  5238  264    6354
2         84   292  421     797
All   151551  6375  777  158703
```

- 2,429 flips out of 158,703 (**1.53%**)
- Flips are roughly symmetric (e.g. 845 of `0→1` vs 852 of `1→0`)
- 3-way winner distribution barely shifts (M0: 151,552 → 151,551; M1: 6,354 → 6,375; M2: 797 → 777)

#### BIC value differences

| Column | model | max abs diff | median | p99 |
|---|---|---|---|---|
| `BICsinglecomp` | M0 (null, deterministic) | **0** | 0 | 0 |
| `BICmiddlecomp` | M1 (EM-based) | 753 | 4.7e-9 | 37.7 |
| `BICpaircomp` | M2/alt (EM-based) | 712 | 2.47 | 56.8 |

**`BICsinglecomp` is bit-identical across pre and post for every one of the 158,703 unique-coord windows.** Since the null BIC is the only fully deterministic of the three (compute_null_BIC has no random init), this is conclusive evidence that the fix did not change BIC computation logic — exactly as Phase 1 predicted.

M1 and M2 differ because the EM in `EMfunctions_timebased_diploid_aware_model.py` reseeds with `np.random.seed()` (no arg) at lines 534 and 699, i.e. fresh entropy on every run. The differences observed are pre-existing run-to-run EM stochasticity, not a consequence of the fix. (This is itself a finding worth recording — see "Side observation" below.)

### Conclusion

✅ **Fix does exactly what Phase 1 predicted: corrects the `total_reads` output column from region-level to per-window count, while leaving the BIC penalty unchanged (because it was already correct).**

The 1.5% winner flips and ~2-unit median M2 BIC differences are entirely attributable to EM random initialization, not to the fix. The `BICsinglecomp` zero-diff result is the smoking gun.

### Side observation — EM nondeterminism

`compute_M1_BIC` and `EM` both call `np.random.seed()` with no argument, reseeding from the OS each invocation. This means running the smoke twice on the same inputs gives different M1/M2 BIC values, and ~1.5% of windows flip 3-way winners between runs. This is independent of the `total_reads` fix and is worth a separate ticket if reproducibility of the per-window 3-way decision matters. Quick mitigations: pass a deterministic seed (e.g. derived from region coords) or run K=3 EM restarts and keep the best LL.

### Files touched

- `em/EMfunctions_timebased_diploid_aware_model.py` — fix in `EMBIC_bin_path_timebased`; the `total_reads` output column now uses `total_reads_per_window[cont]` (the count already accumulated locally) instead of the constant `sum(X_full[t].shape[0] ...)`.
- `em/tests/test_bic_n_correction.py` — 3 new tests:
  - `test_total_reads_and_bic_n_are_per_window`
  - `test_total_reads_column_varies_across_windows`
  - `test_bic_penalty_uses_per_window_n` (analytical null-BIC check)
- `em/KNOWN_ISSUES.md` — issue 1 marked resolved with link here.
- `em/FIX_LOG_window_n_count.md` — this file.

### Validation artifacts

| Path | Contents |
|---|---|
| `smokes/run_smoke4_full_pipeline/timebased/chr22_timebased_ASM.PREFIX.tsv` | Pre-fix baseline, 20 complete regions, 242,897 rows |
| `smokes/run_smoke4_full_pipeline/timebased/chr22_timebased_ASM.tsv.PREFIX_INCOMPLETE` | Original .tmp at kill time (kept for archaeology) |
| `smokes/postfix_subset/chr22_postfix_ASM.tsv` | Post-fix output on the same 20 regions, 242,897 rows |
| `smokes/postfix_subset/weekinit/`, `smokes/postfix_subset/week20/` | Symlinks to the 20 region inputs |
| `smokes/postfix_subset/log/postfix_smk.{o,e}121959963` | SGE logs |

### Tests

`pytest em/tests/` — **18 passed** (15 prior + 3 new). No regressions.
