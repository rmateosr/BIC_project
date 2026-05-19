# ABOUTME: Asserts the --theta-iid generator branch draws each θ[j, h, k] independently.
# ABOUTME: Within-region variance must be substantially higher than the noise-mode baseline.

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from synthetic.generate_synthetic import (  # noqa: E402
    NULL,
    STATIC_ASM,
    TIME_EMERGENT_ASM,
    ASYMMETRIC_TIME_EMERGENT_ASM,
    SHARED_TIME_EMERGENT,
    generate_region_theta,
)


# Expected per-(h, k) Uniform ranges baked into generate_region_theta.
# Keep in lockstep with the generator if the bases are ever changed.
RANGES = {
    STATIC_ASM: {(0, 0): (0.75, 1.00), (1, 0): (0.00, 0.25)},
    TIME_EMERGENT_ASM: {
        (0, 0): (0.75, 1.00), (1, 0): (0.00, 0.25),
        (0, 1): (0.25, 0.55), (1, 1): (0.45, 0.75),
    },
    ASYMMETRIC_TIME_EMERGENT_ASM: {
        (0, 0): (0.75, 1.00), (1, 0): (0.00, 0.25),
        (0, 1): (0.25, 0.55), (1, 1): (0.45, 0.75),
    },
    SHARED_TIME_EMERGENT: {
        (0, 0): (0.6, 0.85), (1, 0): (0.6, 0.85),
        (0, 1): (0.1, 0.35), (1, 1): (0.1, 0.35),
    },
}


def _gen(region_class, J, seed, theta_iid):
    rng = np.random.RandomState(seed)
    return generate_region_theta(region_class, J=J, rng=rng, theta_noise=0.05,
                                 theta_iid=theta_iid)


def test_iid_mode_per_region_means_concentrate_vs_noise_mode():
    # The key behavioural difference: in iid mode each region's mean θ is an
    # average of J independent draws from Uniform(lo, hi), so the per-region
    # mean concentrates around the midpoint as 1/sqrt(J). In noise mode each
    # region's mean is "one base + tiny noise" so the per-region mean inherits
    # the full Uniform(lo, hi) spread of the base. For J=50 the ratio of
    # between-region std-of-means is ~sqrt(50) ≈ 7x.
    J = 50
    n_regions = 80
    for rc, cells in RANGES.items():
        for (h, k), (lo, hi) in cells.items():
            iid_means = np.array([
                _gen(rc, J, seed=1000 + s, theta_iid=True)[:, h, k].mean()
                for s in range(n_regions)
            ])
            noise_means = np.array([
                _gen(rc, J, seed=1000 + s, theta_iid=False)[:, h, k].mean()
                for s in range(n_regions)
            ])
            iid_std = float(np.std(iid_means))
            noise_std = float(np.std(noise_means))
            assert iid_std < 0.5 * noise_std, (
                f"{rc} (h={h}, k={k}): per-region mean std iid={iid_std:.3f} "
                f"vs noise={noise_std:.3f} — iid mode is not concentrating "
                f"as expected"
            )


def test_iid_mode_means_land_in_specified_ranges():
    # Averaging many regions, the mean θ for each (h, k) cell should be close
    # to (lo + hi) / 2.
    J = 20
    n_regions = 200
    for rc, cells in RANGES.items():
        thetas = np.stack([_gen(rc, J, seed=s, theta_iid=True)
                           for s in range(n_regions)])  # (R, J, 2, 2)
        for (h, k), (lo, hi) in cells.items():
            mean = float(thetas[:, :, h, k].mean())
            mid = (lo + hi) / 2
            # Tolerance: 5% of full [0, 1] range.
            assert abs(mean - mid) < 0.05, (
                f"{rc} (h={h}, k={k}): mean θ {mean:.3f} not near "
                f"midpoint {mid:.3f}"
            )


def test_shared_time_emergent_still_shared_across_alleles_in_iid_mode():
    # Even with iid per-CpG draws, both alleles must share the same θ at every
    # CpG and every k for SHARED_TIME_EMERGENT. The "shared" property is about
    # cross-allele equality, not per-CpG noise.
    theta = _gen(SHARED_TIME_EMERGENT, J=30, seed=7, theta_iid=True)
    assert np.array_equal(theta[:, 0, 0], theta[:, 1, 0])
    assert np.array_equal(theta[:, 0, 1], theta[:, 1, 1])
    # k=0 and k=1 still distinct programs
    assert not np.allclose(theta[:, 0, 0], theta[:, 0, 1])


def test_iid_mode_values_in_unit_interval():
    # No clipping is applied in iid mode (uniform draws can't go out of range),
    # but the values must still be valid probabilities.
    for rc in RANGES:
        theta = _gen(rc, J=30, seed=11, theta_iid=True)
        assert theta.min() > 0.0
        assert theta.max() < 1.0


def test_null_iid_keeps_alleles_and_components_identical():
    # NULL must have theta identical across h and k even in iid mode (one
    # per-CpG draw shared across all 4 cells).
    theta = _gen(NULL, J=40, seed=99, theta_iid=True)
    assert np.array_equal(theta[:, 0, 0], theta[:, 1, 0])
    assert np.array_equal(theta[:, 0, 0], theta[:, 0, 1])
    assert np.array_equal(theta[:, 0, 0], theta[:, 1, 1])
