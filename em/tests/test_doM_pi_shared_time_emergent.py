# ABOUTME: Pure-math tests for the no-h M-step on pi (LaTeX Eq. 16).
# ABOUTME: Single simplex per t, no factor of 1/2, t=1 boundary emerges from soft counts.

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from EMfunctions_timebased_shared_time_emergent_model import doM_pi


def test_doM_pi_satisfies_simplex_per_t():
    """sum_k pi[k, t] = 1 to machine precision at every t (LaTeX §10.5 sanity 5)."""
    rng = np.random.default_rng(0)
    T = 3
    gamma_by_t = {}
    for t in range(1, T + 1):
        # Random soft assignments, then renormalize so each row sums to 1.
        raw = rng.random(size=(25, 2))
        gamma_by_t[t] = raw / raw.sum(axis=1, keepdims=True)
    pi = doM_pi(gamma_by_t, T=T)
    assert pi.shape == (2, T)
    np.testing.assert_allclose(pi.sum(axis=0), 1.0, atol=1e-12)


def test_doM_pi_matches_N_over_I_formula():
    """pi[k, t] = N_{k,t} / I_t exactly on a hand-built gamma matrix (LaTeX Eq. 16)."""
    # I_2 = 4 reads at t=2, evenly split.
    gamma_2 = np.array([
        [1.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [0.0, 1.0],
    ])
    # t=1 is all-baseline (k=1 forbidden via gamma).
    gamma_1 = np.array([[1.0, 0.0]] * 3)
    gamma_by_t = {1: gamma_1, 2: gamma_2}
    pi = doM_pi(gamma_by_t, T=2)
    # t=1: 3/3 in k=0, 0/3 in k=1
    np.testing.assert_allclose(pi[:, 0], [1.0, 0.0])
    # t=2: 2/4 each
    np.testing.assert_allclose(pi[:, 1], [0.5, 0.5])


def test_doM_pi_boundary_at_t1_from_soft_counts_alone():
    """t=1 baseline anchor is preserved by the formula with no special case (LaTeX §10.6).

    If the E-step has zeroed gamma[:, 1] at t=1 (which it does — see the
    PSEUDO short-circuit in doE), then N_{2,1} = 0 and the formula gives
    pi[1, 0] = 0 / I_1 = 0 without any branching.
    """
    rng = np.random.default_rng(1)
    I_1 = 10
    # E-step output at t=1 with the boundary applied: all mass on k=0.
    gamma_1 = np.zeros((I_1, 2))
    gamma_1[:, 0] = 1.0
    # Arbitrary t=2.
    raw_2 = rng.random(size=(8, 2))
    gamma_2 = raw_2 / raw_2.sum(axis=1, keepdims=True)
    pi = doM_pi({1: gamma_1, 2: gamma_2}, T=2)
    assert pi[1, 0] == 0.0
    assert pi[0, 0] == 1.0
    # And t=2 is unaffected.
    np.testing.assert_allclose(pi[:, 1].sum(), 1.0, atol=1e-12)


def test_doM_pi_zero_reads_falls_back_to_baseline():
    """Degenerate I_t=0 case (LaTeX §10.7): keep the constraint-feasible default."""
    gamma_by_t = {1: np.zeros((0, 2)), 2: np.array([[0.5, 0.5]] * 4)}
    pi = doM_pi(gamma_by_t, T=2)
    assert pi[0, 0] == 1.0 and pi[1, 0] == 0.0
    np.testing.assert_allclose(pi[:, 1], [0.5, 0.5])
