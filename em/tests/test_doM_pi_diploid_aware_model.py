# ABOUTME: Unit tests for doM_pi under the diploid per-allele marginal constraint.
# ABOUTME: Covers: marginal invariant, asymmetry preservation, degenerate-allele fallback.
import numpy as np

from EMfunctions_timebased_diploid_aware_model import doM_pi


class TestMarginalInvariant:
    """After any doM_pi call: sum_k pi[h, k, t] == 0.5 for each allele h and time t."""

    def test_random_gamma(self):
        rng = np.random.RandomState(0)
        T, I = 4, 25
        gamma_by_t = {t: rng.rand(I, 2, 2) for t in range(1, T + 1)}
        gamma_by_t[1][:, :, 1] = 0.0  # t=1 boundary: k=1 absent

        pi = doM_pi(gamma_by_t, T)

        for t in range(T):
            assert np.isclose(pi[0, :, t].sum(), 0.5), f"allele 0 marginal != 0.5 at t={t+1}"
            assert np.isclose(pi[1, :, t].sum(), 0.5), f"allele 1 marginal != 0.5 at t={t+1}"


class TestAsymmetryPreserved:
    """Toy dataset from the educational guide — the signature asymmetry is preserved."""

    def test_toy_from_educational_guide(self):
        # t=2: 10 allele-1 normal, 0 altered; 6 allele-2 normal, 4 altered.
        # Expected: pi[:, :, 1] = [[0.5, 0.0], [0.3, 0.2]].
        gamma_t2 = np.zeros((20, 2, 2))
        gamma_t2[:10, 0, 0] = 1.0
        gamma_t2[10:16, 1, 0] = 1.0
        gamma_t2[16:20, 1, 1] = 1.0

        # Trivial t=1 consistent with the boundary (all reads on their tagged allele, k=0).
        gamma_t1 = np.zeros((20, 2, 2))
        gamma_t1[:10, 0, 0] = 1.0
        gamma_t1[10:, 1, 0] = 1.0

        pi = doM_pi({1: gamma_t1, 2: gamma_t2}, T=2)

        np.testing.assert_allclose(pi[:, :, 1], [[0.5, 0.0], [0.3, 0.2]], atol=1e-12)
        # Allele-specific altered fraction is recovered, not erased.
        assert pi[0, 1, 1] != pi[1, 1, 1]


class TestDegenerateAllele:
    """If S_h == 0 at some (h, t), fall back to (pi[h,0,t], pi[h,1,t]) = (0.5, 0) — no NaN/Inf."""

    def test_zero_mass_on_allele_1(self):
        # All reads at t=2 assign zero posterior to allele 1.
        gamma_t2 = np.zeros((10, 2, 2))
        gamma_t2[:, 0, 0] = 0.8
        gamma_t2[:, 0, 1] = 0.2
        gamma_t1 = np.zeros((10, 2, 2))
        gamma_t1[:, 0, 0] = 1.0

        pi = doM_pi({1: gamma_t1, 2: gamma_t2}, T=2)

        assert np.isfinite(pi).all()
        np.testing.assert_allclose(pi[1, 0, 1], 0.5, atol=1e-12)
        np.testing.assert_allclose(pi[1, 1, 1], 0.0, atol=1e-12)
        # And allele 0 still respects the marginal.
        np.testing.assert_allclose(pi[0, :, 1].sum(), 0.5, atol=1e-12)
