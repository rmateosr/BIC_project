# Future-session prompt — fix the inflated `n` in BIC penalty

Paste this prompt into a fresh Claude Code session running under the repo root. It is self-contained — no prior context needed beyond `em/KNOWN_ISSUES.md`.

---

## Prompt

> I need you to evaluate and fix issue #1 in `em/KNOWN_ISSUES.md` ("`total_reads` and BIC `n` count REGION rows, not WINDOW rows"). Read that file in full before starting — it has the diagnosis, evidence, and a suggested fix sketch.
>
> ### Phase 1 — Evaluate (no code changes yet)
>
> 1. **Confirm the bug applies to all three models.** Open `em/EMfunctions_timebased_diploid_aware_model.py` and trace `n` through:
>    - `_compute_window_BIC` (line ~237)
>    - `compute_null_BIC` — find it, check whether it derives `n` from `X_window.shape[0]` or already filters
>    - `compute_M1_BIC` — same check
>    - `_process_single_window` (line ~178)
>
>    For each, report: where `n` comes from, whether it's per-region or per-window, what BIC formula uses it.
>
> 2. **Quantify the bias on real data.** Use the chr22 smoke output at `smokes/run_smoke4_full_pipeline/timebased/chr22_timebased_ASM.tsv` (or `.tmp` if the run isn't finished). Compute:
>    - Distribution of `total_reads` (already known: median ~5277)
>    - Per-window real depth = number of rows in `X_window[t]` with any non-NaN call. You'll have to load the region files and re-derive this. Sample 10 random regions × 10 random windows each. Report the median ratio `real_depth / total_reads`.
>    - The BIC penalty difference between models at the inflated vs corrected `n`. Specifically:
>      - `Δ_penalty_M1_vs_M0 = (p_M1 − p_M0) × log(n)` at the two values of `n`
>      - `Δ_penalty_M2_vs_M1 = (p_M2 − p_M1) × log(n)` at the two values of `n`
>
>    Conclude: how many of the current `BIC_3way_winner=0` windows would FLIP to M1 or M2 if `n` were corrected? Estimate from log-likelihood deltas that are already in the output columns. If the answer is "a meaningful fraction", the fix is high-value; if it's "almost none because BIC differences are huge anyway", flag that too.
>
> 3. **Check tests.** Will any existing test in `em/tests/` break under the fix? Specifically `test_3way_BIC.py` and `test_synthetic_diploid_aware.py` — they use synthetic data where region == window, so the bug might not manifest there at all. Run them now to establish the baseline.
>
> 4. **Decide.** Three possible outcomes:
>    - (a) Fix `n` correctly and accept that BIC values shift across the existing output.
>    - (b) Add the correct `n` as a new column (`n_effective_window`) and leave the old `n` in place to avoid breaking downstream consumers; flag the old column as deprecated.
>    - (c) Decide the bias is harmless because `n` is constant within a region and only relative BIC matters → just fix the misleading `total_reads` reporting, leave the penalty alone.
>
>    Recommend one. Explain why.
>
> Report findings as a short markdown summary (under 400 words). DO NOT change any code yet.
>
> ### Phase 2 — Implement (only after I approve the Phase 1 plan)
>
> Apply the chosen fix. Constraints:
>
> - Do NOT pre-filter rows inside the EM functions (`EM`, `doE`, `doM_pi`, etc.). The EM needs the full matrix.
> - Add unit tests in `em/tests/test_bic_n_correction.py` that:
>   - Construct a synthetic region with 1,000 reads but only 30 reads having calls at a specific 10-CpG window.
>   - Assert that the BIC penalty uses 30, not 1,000.
>   - Assert that `total_reads` in the output column is 30, not 1,000.
> - Do NOT touch `em_unconstrained/` — it's a frozen manuscript ablation.
> - Run `pytest em/tests/` and confirm zero failures before declaring done.
>
> ### Phase 3 — Validate
>
> First, **confirm the pre-fix smoke (job 121954972) has finished** by checking `qstat -u $USER | grep dip_smk` returns empty AND `smokes/run_smoke4_full_pipeline/timebased/chr22_timebased_ASM.tsv` exists (the `.tmp` rename signals clean exit). Do NOT proceed if the job is still running — Phase 3 would clobber the in-flight output.
>
> Then preserve the pre-fix baseline:
>
> ```bash
> cd <repo>
> mv smokes/run_smoke4_full_pipeline smokes/run_smoke4_full_pipeline_PREFIX
> ```
>
> Re-run the chr22 smoke (`bash smokes/run_diploid_full_pipeline_smoke_chr22.sh`) and compare the new output against `smokes/run_smoke4_full_pipeline_PREFIX/timebased/chr22_timebased_ASM.tsv`. Report:
> - How many windows changed `BIC_3way_winner` (0/1/2 flips)
> - Distribution of `total_reads` before vs after (should drop by ~50–200×)
> - Whether the 3-way winner distribution shifted toward M1/M2 (expected) or stayed the same (interesting — would indicate the BIC LL term dominates the penalty regardless)
>
> Save findings to `em/FIX_LOG_window_n_count.md` and update `em/KNOWN_ISSUES.md` to mark the issue resolved (with a link to the fix log).

---

## What this prompt assumes

- The chr22 smoke output (or `.tmp`) is still on disk at `smokes/run_smoke4_full_pipeline/timebased/`.
- The region files are still at `smokes/run_smoke4_full_pipeline/{weekinit,week20}/read_format_split/chr22/`.
- Python 3.12 + numpy/pandas/scipy/pytest are available (on an HPC with environment modules, e.g. `module load python/3.12.0`).
- The `em_unconstrained/` exclusion principle from the original migration is still in force.

## What this prompt does NOT do

- It doesn't try to redesign the BIC criterion itself (the use of BIC vs AIC vs cross-validation is a separate methodological question, out of scope).
- It doesn't try to fix the `compute_M1_BIC` and other helpers if Phase 1 reveals they have unrelated bugs — only the `n` inflation.
- It doesn't auto-rerun the full chr22 smoke without operator approval (Phase 3 is gated on Phase 1+2 sign-off).
