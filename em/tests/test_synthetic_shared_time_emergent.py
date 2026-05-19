# ABOUTME: Tests the no-h model on pure-Class-C synthetic data (shared baseline + symmetric drift).
# ABOUTME: The shared generator's TIME_EMERGENT_ASM bakes in allele difference at k=0, so we synthesise inline.

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from EMfunctions_timebased_shared_time_emergent_model import (
    EM,
    compute_shared_BIC,
    compute_shared_time_BIC,
)
import EMfunctions_timebased_diploid_aware_model as em_dip


def _generate_class_c_pure(I, J, T, pi_alt_final, seed):
    """Pure Class C: alleles share theta at baseline; a subpopulation drifts symmetrically.

    Both alleles are indistinguishable — theta has no h-axis dependency at
    either k. The fraction of reads in the altered program ramps linearly
    from 0 at t=1 to pi_alt_final at t=T. Returns X_by_t and tags_by_t.

    Generator matches the no-h LaTeX model exactly:
      theta_{j, k=0} fixed across alleles, theta_{j, k=1} distinct profile.
      pi_{k=1, t} = pi_alt_final * (t-1)/(T-1)
    """
    rng = np.random.default_rng(seed)
    theta_k0 = rng.uniform(0.6, 0.85, size=J)
    theta_k1 = rng.uniform(0.1, 0.35, size=J)

    pi_alt_by_t = np.zeros(T)
    if T > 1:
        pi_alt_by_t = np.linspace(0.0, pi_alt_final, T)

    X_by_t = {}
    tags_by_t = {}
    for t in range(1, T + 1):
        pi_alt = pi_alt_by_t[t - 1]
        k_assignments = (rng.random(size=I) < pi_alt).astype(int)
        X = np.zeros((I, J), dtype=float)
        for i in range(I):
            theta = theta_k1 if k_assignments[i] == 1 else theta_k0
            X[i] = (rng.random(size=J) < theta).astype(float)
        X_by_t[t] = X
        # Tags don't matter to the no-h EM; pass dummies for the diploid EM.
        tags_by_t[t] = np.array(["noH"] * I)
    return X_by_t, tags_by_t


def test_recovered_pi_t2_tracks_truth():
    """The recovered pi[1, T-1] should land within sampling noise of the truth."""
    T, I, J = 2, 200, 10
    pi_alt_final = 0.5
    X_by_t, _ = _generate_class_c_pure(I=I, J=J, T=T, pi_alt_final=pi_alt_final, seed=0)
    _, pi, _, _ = EM(X_by_t, T=T)
    # Within ±0.10 of the truth at I=200 reads, J=10 sites — well above the
    # 2*sqrt(p*(1-p)/I) ≈ 0.07 binomial std at p=0.5.
    assert abs(pi[1, T - 1] - pi_alt_final) < 0.10, (
        f"Recovered pi[1, T-1]={pi[1, T-1]:.3f} far from truth {pi_alt_final:.3f}"
    )
    # And the t=1 anchor still holds.
    assert pi[1, 0] == 0.0
    assert pi[0, 0] == 1.0


def test_shared_time_beats_shared_on_real_drift():
    """When the truth has a real second component, shared_time should beat shared.

    Class C is nested in the shared null model. With a strong real drift
    (pi_alt_final=0.6, I=200) the extra (J + T - 1) parameters pay back.
    """
    T, I, J = 2, 200, 10
    X_by_t, _ = _generate_class_c_pure(I=I, J=J, T=T, pi_alt_final=0.6, seed=1)
    bic_shared = compute_shared_BIC(X_by_t, T)
    bic_shared_time, _ = compute_shared_time_BIC(X_by_t, T)
    assert bic_shared_time < bic_shared, (
        f"Expected shared_time to beat shared on Class-C-pure drift: "
        f"shared={bic_shared:.2f}, shared_time={bic_shared_time:.2f}"
    )


def test_shared_time_beats_diploid_M2_on_pure_class_c():
    """On Class-C-pure data the diploid M2 overspends parameters and loses.

    Both EMs can fit the bimodality, but M2 has 4J + 2T - 2 params vs
    shared_time's 2J + T - 1. Extra penalty = (2J + T - 1) * log(n).
    At J=10, T=2, n=400 (=2 * I): penalty diff ≈ 21 * log(400) ≈ 126 in
    M2's disfavour. Unless M2 fits ~63 LL units better — which it can't
    on truly shared baseline data — shared_time wins.
    """
    T, I, J = 2, 200, 10
    X_by_t, tags_by_t = _generate_class_c_pure(
        I=I, J=J, T=T, pi_alt_final=0.6, seed=2,
    )
    n_total = sum(X_by_t[t].shape[0] for t in X_by_t)

    # Diploid M2 BIC.
    theta_dip, pi_dip, _, _ = em_dip.EM(X_by_t, tags_by_t, T)
    ll_dip = em_dip.compute_alt_loglikelihood(X_by_t, tags_by_t, theta_dip, pi_dip, T)
    p_asm_time = 4 * J + 2 * T - 2
    bic_asm_time = p_asm_time * np.log(n_total) - 2 * ll_dip

    bic_shared_time, _ = compute_shared_time_BIC(X_by_t, T)

    assert bic_shared_time < bic_asm_time, (
        f"Expected shared_time to beat asm_time on pure Class C: "
        f"shared_time={bic_shared_time:.2f}, asm_time={bic_asm_time:.2f}"
    )
