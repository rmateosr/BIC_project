# ABOUTME: Regression test — BIC alternative-model parameter count is 4J + 2T - 2.
# ABOUTME: Verifies via the BIC expression; the constant is exercised in _compute_window_BIC.
from unittest.mock import patch

import numpy as np

import EMfunctions_timebased_diploid_aware_model as em


def test_bic_alternative_parameter_count_is_4J_plus_2T_minus_2():
    """Mock EM to fix theta/pi, then recompute expected BIC and compare.

    _compute_window_BIC internally runs EM with random init; mocking it removes
    the stochasticity so we can assert an exact match against the 4J+2T-2 formula.
    """
    rng = np.random.RandomState(0)
    T, J = 3, 5
    X_by_t = {t: (rng.rand(30, J) > 0.5).astype(float) for t in range(1, T + 1)}
    tags_by_t = {
        t: np.array(["H1"] * 10 + ["H2"] * 10 + ["noH"] * 10)
        for t in range(1, T + 1)
    }

    theta_fixed = np.full((J, 2, 2), 0.5)
    pi_fixed = np.ones((2, 2, T)) / 4
    pi_fixed[:, 1, 0] = 0.0
    pi_fixed[:, 0, 0] = 0.5

    with patch.object(em, "EM", return_value=(theta_fixed, pi_fixed, {}, 1)):
        _, _, bic_alt, _, _, _, _, _ = em._compute_window_BIC(X_by_t, tags_by_t, T)

    ll_alt = em.compute_alt_loglikelihood(X_by_t, tags_by_t, theta_fixed, pi_fixed, T)
    n_total = sum(X_by_t[t].shape[0] for t in range(1, T + 1))
    expected_p_alt = 4 * J + 2 * T - 2
    expected_bic_alt = expected_p_alt * np.log(n_total) - 2 * ll_alt

    np.testing.assert_allclose(bic_alt, expected_bic_alt, atol=1e-8)


def test_bic_formula_string_is_present_in_source():
    """Cheap regression check: the source of _compute_window_BIC uses 4J+2T-2, not 4J+3T-2."""
    import inspect

    src = inspect.getsource(em._compute_window_BIC)
    assert "4 * J + 2 * T - 2" in src
    assert "4 * J + 3 * T - 2" not in src
