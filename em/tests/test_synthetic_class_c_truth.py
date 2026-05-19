# ABOUTME: Asserts the SHARED_TIME_EMERGENT generator branch produces pure Class C truth.
# ABOUTME: Both alleles must share theta at every k and ramp pi[h,1,t] together.

import os
import shutil
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from synthetic.generate_synthetic import (  # noqa: E402
    SHARED_TIME_EMERGENT,
    CLASSES,
    generate_region_theta,
    generate_synthetic_dataset,
)


def test_shared_time_emergent_is_registered():
    assert SHARED_TIME_EMERGENT in CLASSES


def test_theta_shared_across_alleles_at_both_components():
    rng = np.random.RandomState(0)
    theta = generate_region_theta(SHARED_TIME_EMERGENT, J=20, rng=rng, theta_noise=0.05)
    assert theta.shape == (20, 2, 2)
    for k in (0, 1):
        assert np.allclose(theta[:, 0, k], theta[:, 1, k]), (
            f"theta[:, 0, {k}] and theta[:, 1, {k}] differ — not pure Class C"
        )
    # Sanity: the two components are not collapsed to the same profile.
    assert not np.allclose(theta[:, 0, 0], theta[:, 0, 1]), (
        "k=0 and k=1 should be distinct programs"
    )


def test_pi_truth_symmetric_in_manifest():
    """Per-allele pi truth columns must match for SHARED_TIME_EMERGENT regions."""
    out = tempfile.mkdtemp(prefix="class_c_truth_")
    try:
        T = 4
        manifest = generate_synthetic_dataset(
            n_regions=5, reads_per_t=40, J=10, T=T,
            class_mix={SHARED_TIME_EMERGENT: 1},
            pi_final=0.6,
            output_dir=out, seed=123,
        )
        rows = manifest[manifest["true_class"] == SHARED_TIME_EMERGENT]
        assert len(rows) > 0, "no SHARED_TIME_EMERGENT rows generated"
        for ti in range(1, T + 1):
            h0 = rows[f"pi_h0_alt_t{ti}"].to_numpy()
            h1 = rows[f"pi_h1_alt_t{ti}"].to_numpy()
            assert np.allclose(h0, h1), (
                f"pi_h0_alt_t{ti} != pi_h1_alt_t{ti}: not bilateral symmetric"
            )
        # t=1 anchor: both alleles must be zero.
        assert (rows["pi_h0_alt_t1"] == 0).all()
        assert (rows["pi_h1_alt_t1"] == 0).all()
        # Ramping: t=T must be nonzero somewhere.
        assert (rows[f"pi_h0_alt_t{T}"] > 0).any()
    finally:
        shutil.rmtree(out, ignore_errors=True)


def test_pi_truth_is_per_region_random_not_linear():
    """Trajectories must vary across regions (per-region random draws, not a
    fixed linear ramp). Take one row per region; the t=2 values across regions
    should not all coincide.
    """
    out = tempfile.mkdtemp(prefix="class_c_random_")
    try:
        T = 4
        manifest = generate_synthetic_dataset(
            n_regions=20, reads_per_t=20, J=10, T=T,
            class_mix={SHARED_TIME_EMERGENT: 1},
            pi_final=0.6,
            output_dir=out, seed=42,
        )
        # One row per region (manifest has one row per sliding window).
        per_region = manifest.drop_duplicates("region_idx")
        t2 = per_region["pi_h0_alt_t2"].to_numpy()
        tT = per_region[f"pi_h0_alt_t{T}"].to_numpy()
        # If trajectories were linear with the same pi_final, t=2 would be
        # constant across regions (= pi_final/(T-1)/2).
        assert t2.std() > 1e-6, "pi_h0_alt_t2 is constant across regions — looks linear"
        # And t=T should also vary (random endpoint, not a fixed pi_final/2).
        assert tT.std() > 1e-6, f"pi_h0_alt_t{T} is constant — fixed endpoint, not random"
        # Bounded: every draw must lie in [0, pi_final/2] (the /2 is the
        # diploid marginal: pi[h,1,t] = 0.5 * region_traj[t-1, h]).
        for ti in range(2, T + 1):
            col = per_region[f"pi_h0_alt_t{ti}"].to_numpy()
            assert (col >= 0).all() and (col <= 0.6 / 2 + 1e-9).all()
    finally:
        shutil.rmtree(out, ignore_errors=True)


def test_class_mix_dispatch_does_not_leak_other_classes():
    out = tempfile.mkdtemp(prefix="class_c_dispatch_")
    try:
        manifest = generate_synthetic_dataset(
            n_regions=10, reads_per_t=20, J=10, T=2,
            class_mix={SHARED_TIME_EMERGENT: 1},
            pi_final=0.5,
            output_dir=out, seed=7,
        )
        assert (manifest["true_class"] == SHARED_TIME_EMERGENT).all()
    finally:
        shutil.rmtree(out, ignore_errors=True)
