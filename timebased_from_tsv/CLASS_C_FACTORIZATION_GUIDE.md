# Why Class C "factors" — an ADHD-friendly deep dive

> Companion to `MODEL_COMPARISON_VISUAL.md`. You spotted that in Class C the
> mixing weights split as 0.5/0.5 *within* each model and only the
> normal-vs-altered share moves over time. This guide unpacks why that
> matters — biologically, mathematically, and for BIC behavior.

---

## TL;DR (read this first)

> **The 0.5/0.5 within each model is the fingerprint of a perturbation that
> cannot tell the two alleles apart.** It implies the mixture is really
> happening at the *cell* level (some cells normal, some cells altered),
> not at the *allele* level. Class D, by contrast, requires a mechanism
> that *can* read allele identity (SNPs, imprints, XCI). Your derivations
> already encoded this distinction: the old `_ASM.tex` assumed Class C by
> construction; the clean `_ASM_clean.tex` is the strict generalisation
> that handles both. There is one BIC caveat: when the truth is Class C,
> the alt model is overspending parameters and may lose to null in
> low-coverage windows.

---

## 1. What you actually noticed (the factorisation)

> **TL;DR:** π_{h,k,t} can be written as `(allele share) × (state share at time t)`.
> The fact that it factors at all is the whole story.

Inside Class C the four π values at each t look like this:

```
                 t=2 example with altered share = 0.30

                       k=1 normal     k=2 altered
                     ┌──────────────┬──────────────┐
              h=1    │     0.35     │     0.15     │   row sum 0.50
                     ├──────────────┼──────────────┤
              h=2    │     0.35     │     0.15     │   row sum 0.50
                     └──────────────┴──────────────┘
                       col sum 0.70   col sum 0.30
```

Two things to notice:

1. **Row sums are 0.50** at every t. (Allele-1 reads always make up half of all reads, allele-2 reads the other half.)
2. **Column ratios stay 0.5/0.5 within each column.** Inside the normal column, h=1 and h=2 split evenly; same inside the altered column.

> **Key idea:** Both observations are the same observation. They say
> π_{h,k,t} = ψ_h × φ_{k,t} with ψ_h = 0.5. The 2×2 grid carries no
> *interaction* between h and k — only marginals.

▷ **Analogy:** imagine a bag of marbles where each marble has two
independent labels stamped on it: a *colour* (red/blue) and a *finish*
(matte/glossy). If colour and finish were assigned by independent coin
flips, then 50% of matte marbles are red and 50% are blue, and the same
for glossy. Knowing one label tells you nothing about the other. That
independence is what "factorises" looks like.

---

## 2. Class D is the version that does **not** factor

> **TL;DR:** In Class D the marble's colour and finish are correlated.
> Knowing the colour predicts the finish. That's the interaction term.

The Class D grid at the same t = 2 example looks like:

```
                 t=2 example with h=2 altered share = 0.20

                       k=1 normal     k=2 altered
                     ┌──────────────┬──────────────┐
              h=1    │     0.50     │     0.00     │   row sum 0.50
                     ├──────────────┼──────────────┤
              h=2    │     0.30     │     0.20     │   row sum 0.50
                     └──────────────┴──────────────┘
                       col sum 0.80   col sum 0.20
```

Row sums still 0.5/0.5 (diploidy is preserved — both alleles still
contribute equal numbers of reads). But **the column ratios differ**:
inside the altered column, 0% comes from h=1 and 100% comes from h=2.

> **Key idea:** "Knowing a read is altered tells me which allele it came
> from" is the operational definition of allele-specific dynamics. That
> is exactly the conditional information the factorised model cannot carry.

---

## 3. The cell-level mental model (this is the load-bearing intuition)

> **TL;DR:** The mixture isn't really over reads, it's over **cells**.
> Class C = mixture of cell *types*. Class D = mixture of allele *states*
> within cells.

A long-read sample is reads pooled from many cells. Each cell is diploid
(two alleles, h=1 and h=2). What can vary is **what state each cell is
in**.

### Class C: cells in different *states*, both alleles drift together inside each state

```
                ╭────────────────────────────╮     ╭────────────────────────────╮
   normal cell  │  h=1  ──────────  normal   │     │ proportion at t=1: 100%    │
                │  h=2  ──────────  normal   │     │ proportion at t=4:  25%    │
                ╰────────────────────────────╯     ╰────────────────────────────╯

                ╭────────────────────────────╮     ╭────────────────────────────╮
   altered cell │  h=1  ──────────  altered  │     │ proportion at t=1:   0%    │
                │  h=2  ──────────  altered  │     │ proportion at t=4:  75%    │
                ╰────────────────────────────╯     ╰────────────────────────────╯
```

If you sequence reads from a 1:1 mix of these two cell types, you'll get:
50% from normal cells (split 50/50 across alleles inside) and 50% from
altered cells (split 50/50 across alleles inside). That's exactly the
factorisation in §1. Over time the cell-type proportion shifts, the
allele split inside each cell type does not.

### Class D: cells all of one *type*, but one allele is in a different state

```
                ╭────────────────────────────╮
       cell     │  h=1  ──────────  normal   │     ← every cell looks like this
                │  h=2  ──────────  altered  │
                ╰────────────────────────────╯
```

Now there's no cell-type heterogeneity — every cell has the same two
states on the two alleles. The "altered" reads are *always* h=2. Knowing
state determines allele.

▷ **Analogy:** Class C is like a gym where some lockers are open and
some closed (random which); within each open or closed locker, the left
and right shelves are equally likely to be in the same state. Class D is
a gym where the **left shelf** is always open and the **right shelf**
always closed — the door state is glued to the shelf side.

---

## 4. What kind of biology produces each pattern?

> **TL;DR:** Class C arises when the perturbation cannot read DNA. Class D
> requires the perturbation to discriminate the alleles somehow.

| Pattern | Mechanism must be… | Real-world examples |
|---|---|---|
| **Class C — factorises** | sequence-agnostic, cell-state-driven | aging-associated drift, replicative methylation noise, polyclonal tumour with mixed epigenetic clones, global hypomethylation in cancer that affects both copies, chromatin-context-driven changes (e.g. heterochromatin spreading) that don't read the sequence |
| **Class D — does not factor** | sequence-aware (must distinguish alleles) | imprinting (parent-of-origin tag), random monoallelic expression that has stabilised, X-inactivation choice, allele-specific silencing of a tumour-suppressor where one allele carries a regulatory variant, loss-of-imprinting events |

> **Key idea:** the question "is this Class C or Class D?" is biologically
> the question **"does the mechanism causing the change have any way of
> telling the two parental copies apart?"** Most generic stress, ageing,
> and tumour-progression epigenetic noise cannot. Imprinting and XCI
> can, by definition.

---

## 5. How this maps onto the two derivations you wrote

> **TL;DR:** Your old `_ASM.tex` *assumed* the factorisation. Your clean
> `_ASM_clean.tex` *generalises* it. You upgraded from Class-C-only to
> Class-C-and-D capable.

### The old derivation (`_ASM.tex`)

Used two separate latent variables:

- `w_{i,t,h}`: which allele does read i come from?
- `z_{i,t,k}`: which model is read i drawn from?

with parameters `ψ_{h,t}` and `π_{k,t}` *separately*. The factorisation
`π_{h,k,t} = ψ_{h,t} × π_{k,t}` was baked in. **That model can only see
Class C.** It literally cannot represent Class D, because the joint
"allele × model" probability is forced to be the product of marginals.

### The clean derivation (`_ASM_clean.tex`)

Used a unified latent variable `z_{i,t,h,k}` with a joint mixing weight
`π_{h,k,t}`. The factorisation is no longer assumed — it's allowed but
not required. So the clean model can:

- collapse onto Class C when the data look like the factorised pattern, **or**
- pick up the interaction term and fit Class D when one allele is
  driving the change.

> **Key idea:** Your insight that "C factors nicely" is exactly the
> insight that made the old parameterisation tempting in the first place.
> Moving to the unified `z` was the principled way to keep that
> capability while not closing the door on Class D — which is the regime
> the time-based extension is most uniquely useful for.

▷ **Analogy:** the old model is a coordinate system that can only
describe rectangles aligned with the axes. The new model is a coordinate
system that can also describe tilted rectangles. Anything you could
describe before, you can still describe. New things are now possible.

---

## 6. The BIC caveat — Class C "wastes" parameters

> **TL;DR:** When reality is Class C, the alt model is paying for
> ~`2T − 1` extra parameters it didn't actually need. In low-coverage
> windows this can flip the BIC verdict to null even when there's real
> signal.

### Counting

Under the full unified model, π contributes `3T − 2` free parameters
(3 free per t at t ≥ 2, 1 free at t = 1). Under the factorised
Class-C-only model, the truly free π parameters are:

- ψ_h: 0 free (pinned at 0.5 by diploidy)
- φ_{k,t}: T − 1 free (one altered share per t ≥ 2; pinned to 0 at t=1)

So `T − 1` parameters under independence vs `3T − 2` under the full
joint. For T = 4 that's **3 vs 10**.

### Consequence for BIC

```
BIC = p × log(n) − 2 × log L
```

The `p × log(n)` penalty is paid in full whether or not the data make
the extra parameters useful. So a window that genuinely is Class C will:

1. Drive the extra π parameters toward the factorised values during EM
   (ψ_h ≈ 0.5 falls out automatically since reads are 50/50 by tag), but
2. Still pay the full `3T − 2` penalty in `p_alt = 4J + 3T − 2`.

In moderate-coverage windows this is fine — `log L` improvement
dominates. In low-coverage windows or weak Class-C signal it can *flip*
the BIC verdict to null.

> **Key idea:** if we systematically see windows that look Class-C-shaped
> in a manual π-plot but report `BICresult = 0`, that's diagnostic of
> the penalty eating the win, not of the EM failing. The fix would be a
> third BIC contender — a constrained "Class-C-only" model with the
> factorised π — and pick the lowest of {null, factored-alt, full-alt}.
> Worth flagging now, but only worth implementing if the empirical
> outputs show this happening.

---

## 7. Diagnostic checklist — how to tell C from D in real output

> **TL;DR:** Look at the H1-only and H2-only π trajectories separately.
> If both bifurcate over time, it's C. If only one does, it's D.

For a window that won under the alt model, do these checks:

| Check | Class C looks like… | Class D looks like… |
|---|---|---|
| Sum of altered π across alleles, `π_{1,2,t} + π_{2,2,t}` | grows monotonically; this is `pi_altered_t2` for t=2 | grows monotonically (same!) |
| Ratio of allele-1 altered to allele-2 altered, `π_{1,2,t} / π_{2,2,t}` | ≈ 1 at every t (50/50 inside altered) | very far from 1 (one is ~0) |
| Within-tag bimodality: do H1 reads alone show two methylation modes at late t? | yes — H1 reads bifurcate too | no — H1 reads stay clean |
| Same check on H2 reads alone | yes — H2 reads bifurcate too | yes — only H2 reads bifurcate |
| Reads tagged H1 vs H2: do their per-read posterior γ over (h,k) look symmetric? | yes — mirror-image distributions | no — only one tag's reads have non-trivial altered γ |

> **Key idea:** `pi_altered_t2` alone cannot distinguish C from D. The
> ratio `π_{1,2,t} / π_{2,2,t}` (or equivalently the per-tag bifurcation
> pattern) is what does it. If we want a single summary column for this,
> something like `allele_specificity = |π_{1,2,t} − π_{2,2,t}| / (π_{1,2,t} + π_{2,2,t})`
> reads ~0 for C and ~1 for D and would be a useful addition to the output TSV.

---

## 8. When the factorisation assumption breaks (edge cases worth knowing)

> **TL;DR:** Three situations make the simple Class-C picture less clean.

1. **Tag bias.** If WhatsHap tags H1 reads more reliably than H2 reads
   (or vice versa) in a region, the row sums won't be exactly 0.5/0.5.
   The interaction we see in the joint π would then partly reflect
   tagging bias, not biology. Worth checking total H1 vs H2 read counts
   per window.

2. **Aneuploidy / allele dropout.** If one allele is deleted in a subset
   of cells (common in cancer), the row sums become time-dependent —
   you would see ψ_h drift with t, breaking the "row sums are 0.5"
   assumption used in the cell-level mental model.

3. **Mixed Class C + Class D.** Real biology can layer them: a region
   that is constitutively imprinted (Class B baseline ASM) and *also*
   drifts toward complete loss of imprinting in a fraction of cells over
   time. The π grid in that case has both an interaction term *and* a
   shifting marginal. The unified model can fit it; interpretation
   requires looking at both `mean_theta_diff` and the allele-specificity
   ratio above.

> **Key idea:** the four-class taxonomy is a *coordinate system*, not a
> mutually exclusive set of buckets. Real windows can sit between
> classes. The diagnostics in §7 are the right way to localise where a
> window sits in that 2-D space (allele-asymmetry × temporal-drift).

---

## 9. One-paragraph executive summary

The within-model 50/50 you noticed in Class C is the signature of
independence between allele identity and methylation state. It comes
from the mixture being at the cell level (some cells normal, some
altered) rather than at the allele level inside cells. Biologically it
means the perturbation is sequence-agnostic — it cannot tell the two
parental copies apart. Class D requires the opposite: a mechanism that
can. Mathematically, the factorisation `π_{h,k,t} = 0.5 × φ_{k,t}` is
exactly what your old `_ASM.tex` derivation hard-coded; your clean
unified-`z` derivation generalises that to also support Class D, at the
cost of paying for `2T − 1` extra π parameters that are wasted when the
truth is Class C. If small windows ever look visually Class-C but report
`BICresult = 0`, that's the parameter penalty, not a bug.
