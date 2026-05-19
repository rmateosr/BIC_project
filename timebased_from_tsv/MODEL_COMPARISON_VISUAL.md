# Visual map: time-based vs. no-time-based ASM models

> Goal: lay out, side by side, what a sample is expected to look like when it
> truly behaves under the **null** (single-component) model vs. the
> **alternative** (4-component, time-resolved) model — so we can sanity-check
> our intuition before staring at BIC outputs.

---

## 1. Model architecture (what each one assumes about a window)

### Null model (no time, no allele, no mixture)

```
            ┌──────────────────────────────────────────┐
            │  ONE methylation profile θ_j (J probs)   │
            │                                          │
            │     all reads, all alleles, all t ──►    │
            │     drawn iid Bernoulli(θ_j) per CpG     │
            └──────────────────────────────────────────┘

  parameters:    p_null  =  J
```

There is **no t dimension and no h dimension**. Time-point label and
haplotag are ignored entirely. Tagged H1 / H2 / noH reads are pooled.

---

### Alternative model (time-based, 2 alleles × 2 models)

```
                   ┌─────────────── allele h=1 ────────────────┐
                   │                                           │
                   │   k=1 normal  : θ_{j,1,1}                 │
                   │   k=2 altered : θ_{j,1,2}                 │
                   │                                           │
                   └───────────────────────────────────────────┘
                   ┌─────────────── allele h=2 ────────────────┐
                   │                                           │
                   │   k=1 normal  : θ_{j,2,1}                 │
                   │   k=2 altered : θ_{j,2,2}                 │
                   │                                           │
                   └───────────────────────────────────────────┘

   four components share J×2×2 = 4J methylation rates θ_{j,h,k}
   that are FIXED across time
                          │
                          ▼
   what changes across time is the MIXING WEIGHT π_{h,k,t}
   (one 2×2 simplex per time point t = 1..T)
                          │
   ┌──────────────────────┼─────────────────────────────────────┐
   │   t=1 (baseline)     │   π_{h,2,1} = 0  by construction    │
   │                      │   only the (h, k=1) cells exist     │
   │   t=2..T             │   altered model can switch on       │
   └──────────────────────┴─────────────────────────────────────┘

  parameters:   p_alt  =  4J  +  (3T - 2)
                            ▲          ▲
                            │          └── π has 3 free / time at t≥2,
                            │              1 free at t=1  →  3T - 2 total
                            └── fixed θ pool across time
```

Read-type constraints (E-step):

```
   read tag       allowed (h,k) cells          normalised over
   ─────────      ────────────────────────     ───────────────────
   H1             (1,1), (1,2)                 k only
   H2             (2,1), (2,2)                 k only
   noH            (1,1),(1,2),(2,1),(2,2)      all 4
   any read,t=1   (h, k=1) only                reduced subset
```

---

## 2. The 2×2×T grid that the alternative model is fitting

For T = 4 (illustrated), each cell is π_{h,k,t} and the column must sum to 1:

```
                       t=1        t=2        t=3        t=4
                    ┌────────┬──────────┬──────────┬──────────┐
   h=1,  k=1 normal │  π11,1 │  π11,2   │  π11,3   │  π11,4   │
                    ├────────┼──────────┼──────────┼──────────┤
   h=1,  k=2 alter. │   0    │  π12,2   │  π12,3   │  π12,4   │
                    ├────────┼──────────┼──────────┼──────────┤
   h=2,  k=1 normal │  π21,1 │  π21,2   │  π21,3   │  π21,4   │
                    ├────────┼──────────┼──────────┼──────────┤
   h=2,  k=2 alter. │   0    │  π22,2   │  π22,3   │  π22,4   │
                    └────────┴──────────┴──────────┴──────────┘
                      sum=1     sum=1      sum=1      sum=1
```

A whole sample's "behavior" is essentially **the trajectory of this grid
across t** (plus how distinct the four θ profiles end up being).

---

## 3. The four prototypical sample behaviors

I'm picking the four canonical regimes that span what the pipeline can see.
Each box shows the expected π trajectory and the expected θ profiles.
"●●●" thick bars are large weights; "·" near zero.

### Class A — pure null (boring window)

> No ASM, no temporal drift. Both alleles share the same methylation profile
> at every time point. Nothing emerges.

```
   π over time                          θ profiles (across CpGs in window)
   ─────────────────────                ──────────────────────────────────
   (1,1) ●●●●●  ●●●●●  ●●●●●  ●●●●●     θ_{j,1,1}  ≈  θ_{j,2,1}
   (1,2)   ·      ·      ·      ·       θ_{j,1,2}  ≈  θ_{j,1,1}    (degenerate)
   (2,1) ●●●●●  ●●●●●  ●●●●●  ●●●●●     θ_{j,2,2}  ≈  θ_{j,2,1}    (degenerate)
   (2,2)   ·      ·      ·      ·
            t=1    t=2    t=3    t=4

   → all four θ collapse to one shared curve
   → BIC verdict:  null wins  (BICresult = 0)
   → readable proxies: pi_altered_t2 ≈ 0,  mean_theta_diff ≈ 0
```

---

### Class B — static ASM, no temporal change

> The alleles really do have different methylation, but nothing changes over
> time. This is what the *original* (non-time-based) model was built to find;
> the time-based model still recovers it because the (h=1) and (h=2) normal
> cells just have different θ.

```
   π over time                          θ profiles
   ─────────────────────                ──────────────────────────────────
   (1,1) ●●●●●  ●●●●●  ●●●●●  ●●●●●     θ_{j,1,1}  =/=  θ_{j,2,1}   ← ASM
   (1,2)   ·      ·      ·      ·       θ_{*,*,2}   degenerate
   (2,1) ●●●●●  ●●●●●  ●●●●●  ●●●●●
   (2,2)   ·      ·      ·      ·
            t=1    t=2    t=3    t=4

   → altered cells stay empty,  normal cells stay near 0.5 each
   → BIC verdict:  alt wins on the strength of θ_{j,1,1} ≠ θ_{j,2,1} alone
   → readable proxies: pi_altered_t2 ≈ 0,  mean_theta_diff ≈ 0  (k=2 ≈ k=1)

   ⚠ this is the case where time-based ≈ classical ASM —
     temporal information adds nothing.
```

---

### Class C — bilateral temporal change (no ASM)

> Both alleles drift toward the altered state together (e.g., global LOI,
> bulk hyper/hypo-methylation event affecting both copies). No allele
> asymmetry; methylation just changes with time.

```
   π over time                          θ profiles
   ─────────────────────                ──────────────────────────────────
   (1,1) ●●●●●  ●●●●   ●●●    ●●        θ_{j,1,1}  ≈  θ_{j,2,1}    (normal,
   (1,2)   ·    ●●     ●●●    ●●●●        e.g. unmethylated)
   (2,1) ●●●●●  ●●●●   ●●●    ●●        θ_{j,1,2}  ≈  θ_{j,2,2}    (altered,
   (2,2)   ·    ●●     ●●●    ●●●●        e.g. methylated)
            t=1    t=2    t=3    t=4

   → π_{*,2,t} grows symmetrically on both alleles
   → BIC verdict:  alt wins because mixture explains the bimodal
     methylation seen at later t much better than a single θ
   → readable proxies: pi_altered_t2 > 0  on BOTH h,
                       mean_theta_diff large (the two k profiles really differ)

   note: an unconstrained π_{h,k,t} cannot tell C from D unless we look at
   the SPLIT between h=1 and h=2 within the k=2 row — see Class D.
```

---

### Class D — allele-specific emergent ASM   ← *the regime the pipeline is designed to detect*

> One allele stays normal across time; the other slides into the altered
> state. This is the canonical "ASM appears over time" story (e.g., one
> allele acquires aberrant methylation in a tumor population).

```
   π over time                          θ profiles
   ─────────────────────                ──────────────────────────────────
   (1,1) ●●●●●  ●●●●●  ●●●●●  ●●●●●     θ_{j,1,1}  (normal allele 1)
   (1,2)   ·      ·      ·      ·       θ_{j,2,1}  ≈ θ_{j,1,1}  (h=2 starts normal)
   (2,1) ●●●●●  ●●●●   ●●●    ●●        θ_{j,2,2}  =/=  θ_{j,2,1}   ← emerges
   (2,2)   ·    ●●     ●●●    ●●●●      θ_{j,1,2}  unconstrained (small mass)
            t=1    t=2    t=3    t=4

   → mass leaks from (2,1) into (2,2) over time;  (1,*) row stays put
   → BIC verdict:  alt wins, and *needs* both allele and model dimensions
   → readable proxies:
        pi_altered_t2 ≈ π_{2,2,2}  (most weight on h=2 only)
        mean_theta_diff large
        H1 reads remain consistent across t,
        H2 reads bifurcate across t
```

Mirror image (h=1 drifts, h=2 stays) is also Class D — same structure.

---

## 4. Decision matrix — which model wins, and why

| Behavior | π pattern | θ pattern | BIC pick | What `pi_altered_t2` reads | What `mean_theta_diff` reads |
|---|---|---|---|---|---|
| A. Null | flat 0.5/0.5 normal cells, altered ≈ 0 | one shared profile | **null** (BICresult = 0) | ≈ 0 | ≈ 0 |
| B. Static ASM | altered ≈ 0 across t | θ differs across h, normal only | **alt** | ≈ 0 | small (k=1 vs k=2 collapsed) |
| C. Bilateral drift | altered grows on both h | k=1 vs k=2 distinct, h-symmetric | **alt** | > 0 on both h | large |
| D. Allele-specific emergent ASM | altered grows on one h only | θ_{j,h,2} differs from θ_{j,h,1} on the drifting allele | **alt** (strongest signal) | > 0, asymmetric | large |

---

## 5. What the time-based extension genuinely buys you

Looking at the matrix, **only Class C and Class D need the time axis**:

- A and B are detectable by the original (no-time) model;
  the time-based model does not lose them, but it doesn't *need* the t-axis to find them.
- C is invisible to the cross-sectional model when run on a *single* time
  point of pooled data — the bimodality at later t looks like noise unless
  you know it grew over t.
- D is the case the time-based formalism is uniquely good at: it ties the
  emerging altered component to a specific allele AND to a specific time
  slice via the structured π_{h,k,t} grid.

So when we look at a sample's per-window outputs we should expect:

- **Most windows → Class A** (BICresult = 0).
- **Some windows → Class B** (constitutive ASM regions like imprinted loci,
  XCI on chrX in females). Here `pi_altered_t2` should stay near 0 and
  `mean_theta_diff` should also be small — the win comes from h, not from k.
- **A few windows → Class C** (e.g., bulk methylation drift around
  treatment). `pi_altered_t2` non-zero on both alleles, `mean_theta_diff`
  large.
- **The interesting hits → Class D** (asymmetric emergence). `pi_altered_t2`
  non-zero, dominated by one of the two alleles, `mean_theta_diff` large.

If empirical results don't follow this taxonomy — e.g. high
`pi_altered_t2` on windows that visually look static, or no Class-D
windows anywhere in regions we expect them — that's the diagnostic
signal that something in the model fit (or the data construction) is off.

---

## 6. One-line check against the code

The visual above is consistent with `EMfunctions_timebased.py`:

- `doE` zeroes `log_gamma_prime[:, :, 1]` at t=1   → enforces Class-A/B/D's empty altered row at baseline.
- `doM_theta` pools across `t` for each (h, k)     → θ is fixed across t (the four θ profiles in §1).
- `doM_pi` recomputes per-t                        → π grids in §3 evolve freely across t.
- `pi_altered_t2 = π[0,1,1] + π[1,1,1]`            → matches the "altered mass at t=2" intuition used in §4.
- `mean_theta_diff = mean |θ[:,:,1] − θ[:,:,0]|`   → matches the "k=1 vs k=2 spread" axis used in §4.
