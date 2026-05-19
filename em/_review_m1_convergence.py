# ABOUTME: Review-only sanity check for item 5 — M1 EM convergence with pinned pi.
# ABOUTME: Verifies gamma[:,:,1] stays 0, ll is monotone, and converges fast.
"""One-shot diagnostic for the parameter-recovery review (2026-04-27 report).

Runs the constrained M1 EM loop on a synthetic STATIC region, instrumenting
each iteration to record:
  - max(gamma[:, :, 1])   (must be 0 because pi[:, 1, :] is pinned)
  - observed-data log-likelihood  (must be non-decreasing)
  - iteration count to convergence (vs maxIter=1000)

This is read-only — does not modify any pipeline file.
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from EMfunctions_timebased_diploid_aware_model import (
    _observed_loglikelihood,
    doE,
    doM_theta,
)
from parameter_recovery import load_region
from synthetic.generate_synthetic import STATIC_ASM, generate_synthetic_dataset


def m1_em_instrumented(X_by_t, tags_by_t, T, maxIter=1000, tol=1e-8):
    """Run M1 EM with per-iteration diagnostics."""
    J = X_by_t[next(iter(X_by_t))].shape[1]
    np.random.seed(2026)
    theta = np.random.uniform(0.01, 0.99, size=(J, 2, 2))
    pi = np.zeros((2, 2, T))
    pi[:, 0, :] = 0.5

    history = []
    previous_ll = -np.inf
    n_iters = 0
    for iteration in range(maxIter):
        gamma_by_t = {
            t: doE(X_by_t[t], theta, pi[:, :, t - 1], tags_by_t[t], t)
            for t in range(1, T + 1)
        }
        # Instrument: max gamma at k=1 across all time points.
        max_g_k1 = max(
            float(np.max(gamma_by_t[t][:, :, 1])) for t in range(1, T + 1)
        )
        theta = doM_theta(X_by_t, gamma_by_t)
        ll = _observed_loglikelihood(X_by_t, tags_by_t, theta, pi, T)
        history.append((iteration, ll, max_g_k1))
        n_iters = iteration + 1
        if abs(ll - previous_ll) < tol:
            break
        previous_ll = ll
    return history, n_iters


def main():
    with tempfile.TemporaryDirectory() as out:
        manifest = generate_synthetic_dataset(
            n_regions=4,
            reads_per_t=100,
            J=20,
            T=2,
            class_mix={STATIC_ASM: 1},
            pi_final=0.6,
            output_dir=out,
            seed=2026,
        )
        import glob
        paths = {
            t: sorted(glob.glob(os.path.join(out, f"t{t}", "methylationfraction_*_.tsv")))[0]
            for t in (1, 2)
        }
        X_by_t, tags_by_t, _ = load_region(paths)

    history, n_iters = m1_em_instrumented(X_by_t, tags_by_t, T=2)

    print(f"n_iters_to_convergence = {n_iters} / 1000")
    print()
    print("iter |       ll       | max gamma[:,:,1]")
    print("-----+----------------+-----------------")
    for it, ll, gk1 in history[:12]:
        print(f"{it:4d} | {ll:14.6f} | {gk1:.3e}")
    if len(history) > 12:
        print(f" ... ({len(history) - 12} more iterations not shown)")
        it, ll, gk1 = history[-1]
        print(f"{it:4d} | {ll:14.6f} | {gk1:.3e}  (last)")
    print()

    # Monotonicity check.
    lls = [h[1] for h in history]
    diffs = np.diff(lls)
    n_decreases = int((diffs < -1e-9).sum())
    max_decrease = float(diffs.min()) if len(diffs) else 0.0
    print(f"ll monotone non-decreasing: n_decreases (Δ<-1e-9) = {n_decreases}")
    print(f"  worst (most negative) Δll = {max_decrease:.3e}")

    # gamma[:,:,1] must be exactly 0.
    max_gk1_overall = max(h[2] for h in history)
    print(f"max gamma[:,:,1] across all iterations = {max_gk1_overall:.3e}")
    print()
    print(
        "PASS criteria:\n"
        "  - gamma[:,:,1] == 0 every iteration\n"
        "  - ll monotone non-decreasing within EM tolerance\n"
        "  - n_iters << maxIter=1000"
    )


if __name__ == "__main__":
    main()
