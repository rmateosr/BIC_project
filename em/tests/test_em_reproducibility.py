# ABOUTME: Verify EM() and compute_M1_BIC() produce identical results when
# ABOUTME: re-run on the same input (deterministic seed from X_by_t content).

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from em.EMfunctions_timebased_diploid_aware_model import EM, compute_M1_BIC


def _make_X_by_t(T=2, I=40, J=10, p=0.3, seed=123):
    rng = np.random.default_rng(seed)
    X_by_t = {}
    tags_by_t = {}
    for t in range(1, T + 1):
        X_by_t[t] = rng.binomial(1, p, size=(I, J)).astype(float)
        tags_by_t[t] = np.array(["noH"] * I)
    return X_by_t, tags_by_t


def test_EM_is_reproducible():
    X_by_t, tags_by_t = _make_X_by_t()
    theta_a, pi_a, _, _ = EM(X_by_t, tags_by_t, T=2)
    theta_b, pi_b, _, _ = EM(X_by_t, tags_by_t, T=2)
    np.testing.assert_array_equal(theta_a, theta_b)
    np.testing.assert_array_equal(pi_a, pi_b)


def test_compute_M1_BIC_is_reproducible():
    X_by_t, tags_by_t = _make_X_by_t()
    bic_a, n_a = compute_M1_BIC(X_by_t, tags_by_t, T=2)
    bic_b, n_b = compute_M1_BIC(X_by_t, tags_by_t, T=2)
    assert bic_a == bic_b
    assert n_a == n_b


def test_EM_differs_on_different_input():
    """Sanity: the seed actually changes when the input changes."""
    X1, tags = _make_X_by_t(seed=1)
    X2, _ = _make_X_by_t(seed=2)
    theta_1, _, _, _ = EM(X1, tags, T=2)
    theta_2, _, _, _ = EM(X2, tags, T=2)
    assert not np.array_equal(theta_1, theta_2)
