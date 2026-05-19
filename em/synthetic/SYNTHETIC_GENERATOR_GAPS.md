# Where `em/synthetic/generate_synthetic.py` falls short of the diploid-aware model

> **TL;DR.** The generator is a faithful but **conservative** simulator: every read it
> writes is consistent with the production model, but it only ever produces the easy
> *symmetric* subcase. The hardest and most diagnostic regimes — the ones the
> diploid-aware constraint was designed for — never appear in the synthetic data.
> This doc explains each gap in plain language so we know what to fix.

---

## 0. Quick glossary (read first)

The model has a handful of one-letter symbols. They're easy once you've seen them:

| Symbol | Stands for | What it actually means |
|---|---|---|
| `h` | **h**aplotype / **a**llele index | `h ∈ {0, 1}`. `h=0` is "the maternal copy of the chromosome", `h=1` is "the paternal copy" (or vice versa — we don't actually know which is which, just that they're the two parental copies). When a read is tagged `H1` by WhatsHap, we set `h=0`; tagged `H2` → `h=1`; untagged `noH` → unknown, EM has to guess. |
| `k` | "**k**omponent" of the mixture | `k ∈ {0, 1}`. `k=0` is the **normal** methylation program ("baseline"). `k=1` is the **altered** program ("the new state that emerges over time"). At `t=1` we *force* `k=0` for everyone — there's no altered state at the start. |
| `t` | **t**ime point | `t ∈ {1, 2, …, T}`. We sample at multiple time points (e.g. `weekinit`, `week20`). |
| `j` | CpG site index within a window | `j ∈ {0, …, J-1}`. A window has `J=10` CpGs by default. |
| `i` | read index | one read = one Oxford-Nanopore molecule that crossed the window. |
| `θ[j, h, k]` | per-site methylation prob | "if a read comes from allele `h` and is in component `k`, what's the chance the CpG at position `j` reads as methylated?" |
| `π[h, k, t]` | mixture weights | "at time `t`, what fraction of all reads come from (allele `h`, component `k`)?" |
| diploid marginal constraint | the new thing in this model | for every allele `h` and every time `t`: `π[h,0,t] + π[h,1,t] = 1/2`. In words: each allele contributes exactly half of all reads, no matter how it splits between normal and altered programs. |
| symmetric vs asymmetric π | whether both alleles change in lockstep | **symmetric:** `π[0,1,t] == π[1,1,t]` (both alleles drift to "altered" at the same rate). **Asymmetric:** `π[0,1,t] != π[1,1,t]` (only one allele changes — *this is the genuinely interesting ASM case*). |

If you remember nothing else: **`h` = which parental chromosome copy, `k` = normal vs altered, `t` = time point, `θ` = methylation probability, `π` = how reads are distributed over (allele, component).**

---

## 1. The seven gaps

### Gap 1 — The generator only produces **symmetric** π trajectories

**The code.** In `generate_reads_for_region`:
```python
if t == 1 or rng.rand() >= pi_alt_t:
    true_k = 0
else:
    true_k = 1
```
The probability `pi_alt_t` of being in the altered component (`k=1`) is **the same number** whether the read is on allele `h=0` or allele `h=1`.

**What that means in plain English.** When we generate a TIME_EMERGENT region, *both copies of the chromosome* drift toward the altered state at the **exact same rate**. If 30% of allele-0 reads are altered at week 20, then 30% of allele-1 reads are altered at week 20 too.

**Why that's wrong / weak.** The whole point of calling this *Allele-Specific* Methylation
is that the two alleles behave **differently**. A real ASM event looks like:
- Allele 0: stays methylated (normal) the whole time. `π[0,1,t]` ≈ 0 always.
- Allele 1: progressively un-methylates (drifts to altered). `π[1,1,t]` ramps up.

The diploid-aware model was *specifically* designed to capture this asymmetry while
still pinning the per-allele marginal at 1/2 (because we still expect 50/50 reads
between the two chromosomes regardless of methylation). The synthetic data **never
exercises this asymmetry**, so:
- The diploid-aware EM and the unconstrained EM look equally good on this data.
- We can't actually demonstrate that the diploid constraint *helps detection*.
- We can't tell if `pi_altered_t2` (one of our key reported metrics) is meaningful
  for asymmetric cases.

**ADHD-friendly mental picture.** Imagine you're testing a stereo system that has
separate left and right volume knobs. The current synthetic data is like only ever
testing it with both knobs turned to the same position. You never find out whether
the knobs actually work independently. That's what we're doing to the diploid model.

**Fix sketch.** Let the trajectory be per-allele:
```python
pi_alt_t_h0, pi_alt_t_h1 = pi_alt_traj_h0[t-1], pi_alt_traj_h1[t-1]
prob_k1 = pi_alt_t_h0 if true_h == 0 else pi_alt_t_h1
true_k = 0 if (t == 1 or rng.rand() >= prob_k1) else 1
```
…and add at least one new region class, e.g. `ASYMMETRIC_TIME_EMERGENT_ASM`,
where one allele ramps and the other doesn't.

---

### Gap 2 — `noH` reads draw `h` first, then `k` independently — this only works because π is symmetric

**The code.**
```python
hap = rng.choice(["H1", "H2", "noH"], p=hap_probs)
true_h = 0 if hap == "H1" else (1 if hap == "H2" else rng.choice([0, 1]))
# ... then later, independently:
true_k = 0 if (t == 1 or rng.rand() >= pi_alt_t) else 1
```
For a `noH` read, the generator draws `h` (50/50), then draws `k` *separately* from
`pi_alt_t` — as if `h` and `k` were independent.

**What that means in plain English.** "First flip a coin to pick which chromosome
this read came from. Then, independently, flip another coin to decide if it's normal
or altered."

**Why that's a hidden landmine.** This factorized two-coin draw only matches the
true joint distribution `π[h, k, t]` *if π is symmetric across h*. In other words:
**Gap 2 is invisible right now only because of Gap 1.** The moment we fix Gap 1 and
start generating asymmetric data, this code becomes incorrect for `noH` reads.

**ADHD-friendly mental picture.** Imagine a deck where red cards (allele 0) are
50% aces (altered) and black cards (allele 1) are 0% aces. The probability of drawing
an ace depends on the colour — they're *not independent*. If you draw colour first,
then "is it an ace" with a single shared probability, you'll mis-simulate the deck.

**Fix sketch.** For `noH` reads, do a single 4-way draw over the joint
`π[·, ·, t]` (which sums to 1). For tagged reads, condition on the fixed `h` and
draw `k` from `π[h, ·, t] / 0.5`.
```python
if hap == "noH":
    flat = pi_t.flatten()  # length 4: (h0,k0), (h0,k1), (h1,k0), (h1,k1)
    idx = rng.choice(4, p=flat)
    true_h, true_k = idx // 2, idx % 2
else:
    true_h = 0 if hap == "H1" else 1
    p_k1 = pi_t[true_h, 1] / 0.5
    true_k = 0 if (t == 1 or rng.rand() >= p_k1) else 1
```

---

### Gap 3 — Phasing tags are perfect

**The code.** A read tagged `H1` is *always* generated from allele `h=0`. Tagged
`H2` is *always* from `h=1`. There's no probability of mis-tagging.

**What that means in plain English.** WhatsHap is treated as an oracle that never
makes mistakes.

**Why that's wrong / weak.** Real WhatsHap output **does** have errors — typically
a small fraction of reads (≈1–5%) get tagged to the wrong haplotype, especially in
homozygous regions or low-coverage spans. Our EM uses the tag as a *hard constraint*
in the E-step (`gamma[h1_mask, 1, :]` is forced to 0). Mis-tagged reads therefore
get plugged into the wrong allele's likelihood and bias the θ estimates.

We don't currently know how robust the diploid model is to this. The synthetic
generator can't tell us.

**ADHD-friendly mental picture.** It's like training a self-driving car only on
sunny days. The car looks perfect in the lab. Then you put it on a real road in
fog (tag errors) and have no idea what happens.

**Fix sketch.** Add a `tag_error_rate` parameter (default 0 to preserve current
behaviour). With probability `tag_error_rate`, flip the tag of a fraction of reads
**after** their true `h` is sampled. Document this clearly in the manifest so
evaluation knows which reads are mis-tagged.

---

### Gap 4 — `STATIC_ASM` regions are "doubly null" on the altered component

**The code.** For STATIC_ASM:
```python
theta[:, 0, :] = vals_h0[:, None]   # h=0, both k get the SAME values
theta[:, 1, :] = vals_h1[:, None]   # h=1, both k get the SAME values
# ... and:
region_traj = np.zeros(T)            # pi_alt_t is 0 at every time point
```
Two things are simultaneously degenerate:
1. `θ[:, h, 0] == θ[:, h, 1]` — the normal and altered components emit identical methylation.
2. `π[:, 1, t] = 0` for all `t` — the altered component is *never used* anyway.

**What that means in plain English.** For STATIC ASM, "altered" is set up to be
both invisible (no reads come from it) AND indistinguishable (even if reads did
come from it, they'd look the same as normal). Two safety nets when one would do.

**Why that's wrong / weak.** This is **over-engineered to make the test easy**.
The generator could be more discriminating: keep the allele difference at `k=0`,
and put the altered component somewhere different (with `pi_alt_t > 0`), but in a
way that doesn't violate "static" — e.g. *both* alleles drift toward the altered
component at the same rate, and the altered component's θ is also allele-split but
the same shape across time. Right now STATIC is essentially "no time dynamics,
period," which conflates two distinct things our model is trying to separate.

This matters because Pending work item *"3-way BIC comparison (M0/M1/M2)"* in
`CLAUDE.md` wants to distinguish *stable* ASM from *time-emergent diseased* ASM.
The current STATIC generator produces a degenerate version of stable ASM where
the EM has nothing at all to fit on the altered side.

**ADHD-friendly mental picture.** It's like testing whether a metal detector can
find coins by burying *no* coin and also hiding the spot under concrete. Sure, the
detector says "no coin" — but you didn't really test it.

**Fix sketch.** Make STATIC mean "real allele difference, present at all `t`,
captured by `k=0` only, with `π[h,1,t]=0`." Drop the `θ[:,h,0]=θ[:,h,1]` redundancy
— let `θ[:,h,1]` be drawn freely so the EM has to learn that `π[h,1,t]` is what
suppresses it, not the θ similarity.

---

### Gap 5 — Every TIME_EMERGENT region uses the **same** trajectory shape

**The code.**
```python
traj = make_pi_trajectory(T, pi_final, pi_trajectory)
# ... for every region:
if rc == TIME_EMERGENT_ASM:
    region_traj = traj.copy()
```
Either a linear ramp `0 → pi_final` or a single `pi_trajectory` provided on the
command line. **All TIME_EMERGENT regions in the dataset get this one shape.**

**What that means in plain English.** If you simulate 100 time-emergent regions,
all 100 of them ramp from 0% altered to (say) 60% altered, perfectly linearly,
with the same end-point.

**Why that's wrong / weak.** Real biology has *different* regions changing on
*different* timescales. Some flip fast, some slowly, some plateau, some are
non-monotonic. The model fits a per-region `π[h, k, t]` so it can in principle
handle all of these — but our test data only ever shows it the linear ramp case.
We have no idea whether the EM handles a step function, an S-curve, or a partial
reversal correctly.

**ADHD-friendly mental picture.** It's like a dating app that's been tested only
on people who reply within 5 minutes. You don't know how it handles slow repliers,
silent treatment, or "yes then no then yes again" until you see real users.

**Fix sketch.** Per-region random trajectory shape. E.g., draw `pi_final` per
region from a uniform `[0.2, 0.8]`, and choose between linear / step /
saturating-exponential shapes. Or accept a list of named trajectory types and
sample one per region.

---

### Gap 6 — Within a region, all CpGs share one base θ value (no spatial structure)

**The code.**
```python
base = rng.uniform(0.3, 0.7)
vals = np.clip(base + rng.normal(0, theta_noise, J), 0.01, 0.99)
theta[:, :, :] = vals[:, None, None]
```
Each CpG gets the base value plus i.i.d. Gaussian noise (σ=0.05). Position along
the genome is irrelevant — CpG `j=0` and CpG `j=49` are statistically independent.

**What that means in plain English.** The 50 CpGs in a region all hover around a
single methylation level, with random per-CpG noise. There's no notion of "this
half of the region is methylated, that half isn't" or "methylation correlates with
distance to the nearest CpG island."

**Why that's wrong / weak.** Real CpG methylation is **spatially correlated** —
neighbouring CpGs are more similar than distant ones (this is literally what defines
a CpG island vs. an open-sea CpG). The sliding-window approach we use (10 CpGs at
a time) implicitly assumes signal is *contiguous*: an ASM event spans a stretch
of nearby CpGs. The current generator can't produce regions where, say, CpGs 0–20
show ASM and CpGs 21–49 don't.

This means our windowed BIC is being tested only on regions where the signal is
uniform across all windows — we never test the boundary case where some windows
should call ASM and others shouldn't *within the same region*.

**ADHD-friendly mental picture.** It's like testing a smoke detector by filling
the entire house with smoke or none of it. You never test the realistic case where
smoke is in one room and the detector down the hall has to decide.

**Fix sketch.** Add an option for spatially-structured θ: e.g. an autocorrelated
Gaussian process along the CpG positions, or piecewise-constant θ with a few change
points, so some sliding windows within a region carry ASM signal and others don't.

---

### Gap 7 — Coverage is constant: same `reads_per_t` per region per time point

**The code.**
```python
for t in range(1, T + 1):
    reads = generate_reads_for_region(
        coords, theta, region_traj[t - 1], reads_per_t, t, chrom, rng,
    )
```
`reads_per_t = 80` by default, applied to every region and every time point.

**What that means in plain English.** Every region at every time point gets exactly
the same number of reads (80, default). No variation.

**Why that's wrong / weak.** Real long-read coverage varies enormously:
- Across regions: telomeric, centromeric, and repeat-rich regions have lower
  effective coverage. Some regions have 5×, others 50×.
- Across time points: different sequencing runs have different depth.
- Across haplotypes: heterozygous SVs cause uneven coverage between H1 and H2.

The pipeline already has machinery for coping with low-coverage windows
(`DEPTH = 1` cutoff, all-NaN-row removal, etc.), but the synthetic data never
stresses it. We don't know how the EM behaves at, say, 3 reads/region or with
80% NaNs in a window.

**ADHD-friendly mental picture.** It's like a restaurant testing its menu only
during weekdays at 2pm with three regulars. You never see what happens at the
Saturday-night rush or at Tuesday 10am with one walk-in.

**Fix sketch.** Sample `reads_per_t` per (region, time) from a distribution
(e.g. `Poisson(λ=80)` or `NegativeBinomial` for over-dispersion). Optionally also
add a `coverage_dropout_rate` that occasionally yields tiny-coverage regions to
test the depth-cutoff path.

---

## 2. Priority of fixes (proposed)

| # | Gap | Severity | Effort | Reasoning |
|---|---|---|---|---|
| 1 | Symmetric-only π | 🔴 high | medium | This is the **headline gap** — without it we can't validate the diploid model's reason for existing. |
| 2 | `noH` factorized draw | 🔴 high | small | Becomes **wrong** as soon as Gap 1 is fixed. Fix in the same PR. |
| 4 | STATIC doubly null | 🟡 med | small | Cleaner test design; needed for the planned 3-way BIC. |
| 5 | One trajectory shape | 🟡 med | medium | Stress-tests EM convergence on diverse dynamics. |
| 7 | Constant coverage | 🟡 med | small | Stress-tests low-coverage code paths. |
| 6 | No spatial θ structure | 🟢 low | medium | Real but a separate, larger generator-design conversation. |
| 3 | Perfect phasing tags | 🟢 low | small | Useful but orthogonal — robustness probe, not correctness. |

---

### Gap 8 — Class assignment is a **random multinomial draw**, not stratified

**The code.** In `generate_synthetic_dataset`:
```python
weights = np.array([class_mix.get(c, 0) for c in CLASSES], dtype=float)
weights /= weights.sum()
region_classes = rng.choice(CLASSES, size=n_regions, p=weights)
```
With the typical "equal mix" config (`class_mix={c: 1 for c in CLASSES}`), this
is `n_regions` independent uniform-over-4 draws. The realised per-class counts
follow a multinomial distribution, **not** a fixed `n_regions / 4` per class.

**What that means in plain English.** If you ask for 32 regions with equal
class mix, you don't get exactly 8 of each. You get 32 i.i.d. draws — typically
something like 7/9/8/8 but occasionally 5/11/8/8 or 4/12/9/7.

**Why that's wrong / weak.** Two real-world consequences:

1. **Some classes can in principle be missing entirely from a dataset.**
   `P(class absent | 32 draws, p=1/4) = (3/4)^32 ≈ 0.0001`. Across a 600-cell
   sweep × 4 classes that's an expected ~0.24 missing-class events, so usually
   nothing happens — but the floor isn't zero. The
   `predictive_power_depth_2026-05/` sweep got lucky (all 2 400 cell×class
   slots populated, min count 4); a smaller `n_per_class` would not.
2. **Per-seed accuracy resolution varies by cell.** A seed that drew 4
   ASYMMETRIC regions can only report accuracy on the grid {0.00, 0.25, 0.50,
   0.75, 1.00}; a seed that drew 14 has 15 possible values. The bootstrap CI
   picks this up as legitimate variance, but it's design noise, not signal.

**Miss-rate scaling.** For `n_per_class = K` (so `n_regions = 4K` regions with
equiprobable class mix):

| K | n_regions | P(any one class absent) | Expected missing per 600 cells × 4 classes |
|---|---|---|---|
| 4 | 16 | 1.0 % | ~24 |
| 8 | 32 | 0.01 % | ~0.24 |
| 16 | 64 | 1.6e-8 | effectively 0 |

So the current default (`n_per_class=8`) is safe; smaller values are not.

**ADHD-friendly mental picture.** Imagine you're packing a lunchbox to test
a kid's food preferences. You want them to try one of each: apple, banana,
sandwich, juice. The current code is "reach into a bag containing piles of
each item and pull 32 things out blindly." Usually you get a balanced mix
but sometimes the kid gets 15 apples and zero juices. Stratified sampling
is "put exactly 8 of each in the box on purpose."

**Fix sketch.** Add a `stratified` flag (or detect "all weights equal") and
pre-build the class assignment without sampling:
```python
if all(w == weights[0] for w in weights) and n_regions % len(CLASSES) == 0:
    # Stratified: exactly n_regions / len(CLASSES) of each class, shuffled
    per_class = n_regions // len(CLASSES)
    region_classes = np.repeat(CLASSES, per_class)
    rng.shuffle(region_classes)
else:
    region_classes = rng.choice(CLASSES, size=n_regions, p=weights)
```
This preserves backward compatibility (any non-uniform `class_mix` keeps the
random draw) while making the common "equal mix" case deterministic in its
per-class count.

**Priority.** 🟢 low — current `n_per_class=8` default keeps the miss rate at
~0.01 % per cell, so the gap is mostly cosmetic for runs of that size. Worth
fixing before anyone shrinks `n_per_class` to 4 or below, or runs sweeps
where exact per-class counts matter (e.g. paired-region comparisons across
seeds).

---

## 2b. Downstream footgun — `"NULL"` is in pandas' default NA list

**Not a gap in the generator's math; a sharp edge in the *output format* that
bit us in `predictive_power_depth_2026-05/` (2026-05-13).**

**What happened.** `generate_synthetic.py:33` defines `NULL = "NULL"` and the
manifest stores it as the literal string `"NULL"` in the `true_class` column.
Any downstream consumer that loads a CSV/TSV produced from this manifest with
`pandas.read_csv(...)` *without* extra arguments silently converts those cells
to `NaN`, because `"NULL"` is in pandas' default `na_values` list (alongside
`""`, `"NA"`, `"N/A"`, `"NaN"`, `"nan"`, `"null"`, etc.).

**Symptom seen in production.** `evaluate.py` did
`df.groupby("true_class")`, and pandas drops NaN groups by default, so the
entire NULL class (4 740 of 19 200 rows) silently disappeared from
`accuracy_per_class.csv` and from `fig1_accuracy_vs_depth.png`. The bug was
caught only because the resulting figure had 3 panels instead of 4 — there
was no exception, no warning, no error log.

**The fix that consumers must apply.** Read with `keep_default_na=False`:
```python
pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[""])
```
(`na_values=[""]` keeps the empty string as the only NA marker, which is the
behaviour you almost always want for these TSVs.)

**Why not rename `NULL` to `NO_ASM` in the generator?** Cleaner, but
non-local — would force the prior `parameter_recovery_sweep.py`, the
2026-04-27 report's harness, every test in `em/tests/`, and every external
reader to update at the same time. The read-side workaround is one line per
consumer and breaks nothing.

**If you're writing a new aggregator against these CSVs:** use the
`keep_default_na=False` pattern from the start. Don't assume the column is
clean just because there are no genuinely missing cells — pandas isn't
testing for missingness, it's pattern-matching strings.

---

## 3. What this doc is NOT

- Not a critique of the generator's *correctness*. Every read it produces is a
  legitimate sample from a special case of the diploid-aware model. The issue
  is that the special cases it samples are too benign.
- Not a list of bugs. Nothing here is broken. Each gap is a conscious-or-accidental
  simplification.
- Not a plan. See the next planning step for what we actually do about it.
