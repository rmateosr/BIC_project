# ABOUTME: Param-count + nesting sanity for the shared-time BIC (LaTeX §11).
# ABOUTME: p = 2J + (T-1); on K=1-collapse data the gap to compute_shared_BIC equals (J + T-1)*log(n).

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from EMfunctions_timebased_shared_time_emergent_model import (
    EM,
    _observed_loglikelihood,
    compute_shared_BIC,
    compute_shared_time_BIC,
)


def _make_X_by_t(T=2, I=80, J=10, p=0.3, seed=0):
    rng = np.random.default_rng(seed)
    return {t: rng.binomial(1, p, size=(I, J)).astype(float) for t in range(1, T + 1)}


def test_param_count_in_BIC_formula():
    """p_shared_time = 2J + (T-1) (LaTeX §11). Verify by reconstructing BIC from LL.

    BIC = p*log(n) - 2*LL. Refit, recompute LL, and check the relation holds
    with p = 2J + (T-1) on a known input.
    """
    T, J = 2, 10
    X_by_t = _make_X_by_t(T=T, J=J)
    n = sum(X_by_t[t].shape[0] for t in X_by_t)
    bic, _ = compute_shared_time_BIC(X_by_t, T)

    theta, pi, _, _ = EM(X_by_t, T)
    ll = _observed_loglikelihood(X_by_t, theta, pi, T)
    p_expected = 2 * J + (T - 1)
    bic_recomputed = p_expected * np.log(n) - 2 * ll
    np.testing.assert_allclose(bic, bic_recomputed, rtol=1e-10)


def test_shared_time_BIC_finite():
    """Sanity on random input — no NaN, no inf, no divide-by-zero blowups."""
    X_by_t = _make_X_by_t(T=2, J=10, seed=42)
    bic, n_iters = compute_shared_time_BIC(X_by_t, T=2)
    assert np.isfinite(bic)
    assert n_iters >= 1


def test_shared_time_BIC_loses_to_shared_on_null_data():
    """On NULL-like data (no real second component), Shared should win.

    The K=1 collapse (LaTeX §9) makes the shared-time model nested in the
    shared model. With p_shared = J and p_shared_time = 2J + (T-1), the
    extra (J + T - 1) parameters cost (J + T - 1) * log(n) on the BIC
    side. The LL gain from a spurious second component on noise is bounded,
    so Shared should win at moderate n.
    """
    T, J, I = 2, 10, 100
    X_by_t = _make_X_by_t(T=T, J=J, I=I, p=0.3, seed=7)
    bic_shared = compute_shared_BIC(X_by_t, T)
    bic_shared_time, _ = compute_shared_time_BIC(X_by_t, T)
    assert bic_shared < bic_shared_time, (
        f"Expected Shared to win on null-like data: "
        f"shared={bic_shared:.2f}, shared_time={bic_shared_time:.2f}"
    )


def test_reproducible_BIC_on_same_input():
    """Deterministic seed-from-X => identical BIC on re-runs."""
    X_by_t = _make_X_by_t(seed=11)
    a, na = compute_shared_time_BIC(X_by_t, T=2)
    b, nb = compute_shared_time_BIC(X_by_t, T=2)
    assert a == b
    assert na == nb
