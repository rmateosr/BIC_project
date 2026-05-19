# ABOUTME: Sanity tests for the diploid-aware time-based ASM pipeline + synthetic generator.
# ABOUTME: Validates NULL/STATIC/symmetric/asymmetric detection, per-allele pi recovery, joint (h,k) draw.
"""Diploid-aware sanity tests.

Generates synthetic data, runs EMBIC_bin_path_timebased on it, and checks:
  - detection rate per region class (NULL/STATIC/TIME_EMERGENT/ASYMMETRIC)
  - per-allele pi_h{0,1}_alt_t2 recovery (the new EM output columns)
  - manifest format (incl. the new ramping_allele column)
  - the noH read sampler does a true joint draw over (h, k)

Detection-rate thresholds are loose because EM(...) reseeds via np.random.seed()
in EMfunctions_timebased_diploid_aware_model.py — generator seed alone does not
fully determinize results.
"""

import glob
import os
import tempfile

import numpy as np
import pandas as pd

from EMfunctions_timebased_diploid_aware_model import EMBIC_bin_path_timebased
from synthetic.generate_synthetic import (
    ASYMMETRIC_TIME_EMERGENT_ASM,
    NULL,
    STATIC_ASM,
    TIME_EMERGENT_ASM,
    generate_reads_for_region,
    generate_synthetic_dataset,
)


def _run_pipeline(output_dir, T, windowsize=10):
    """Run EMBIC over every region in a synthetic dataset and concatenate results."""
    files_per_t = {}
    for t in range(1, T + 1):
        files_per_t[t] = sorted(
            glob.glob(os.path.join(output_dir, f"t{t}", "methylationfraction_*_.tsv"))
        )
    n_regions = len(files_per_t[1])
    results = []
    for ridx in range(n_regions):
        paths = {t: files_per_t[t][ridx] for t in range(1, T + 1)}
        result = EMBIC_bin_path_timebased(paths, windowsize=windowsize)
        result["region_idx"] = ridx
        results.append(result)
    return pd.concat(results, ignore_index=True)


# Test sizing — small enough for fast pytest, big enough for stable detection rates.
TEST_N_REGIONS = 8
TEST_READS = 40
TEST_J = 20
TEST_T = 4
TEST_WS = 10
TEST_SEED = 42


class TestNullRegions:
    def test_low_detection_rate(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            generate_synthetic_dataset(
                n_regions=TEST_N_REGIONS, reads_per_t=TEST_READS, J=TEST_J, T=TEST_T,
                class_mix={NULL: 1, STATIC_ASM: 0, TIME_EMERGENT_ASM: 0,
                           ASYMMETRIC_TIME_EMERGENT_ASM: 0},
                windowsize=TEST_WS, output_dir=tmpdir, seed=TEST_SEED,
            )
            results = _run_pipeline(tmpdir, TEST_T, TEST_WS)
        rate = results["BICresult"].mean()
        assert rate < 0.15, f"NULL detection rate {rate:.2%} exceeds 15%"


class TestStaticASM:
    """STATIC_ASM: detection should fire (allele Δθ at k=0); per-allele altered π near 0."""

    def test_high_detection_and_per_allele_pi_low(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            generate_synthetic_dataset(
                n_regions=TEST_N_REGIONS, reads_per_t=TEST_READS, J=TEST_J, T=TEST_T,
                class_mix={NULL: 0, STATIC_ASM: 1, TIME_EMERGENT_ASM: 0,
                           ASYMMETRIC_TIME_EMERGENT_ASM: 0},
                windowsize=TEST_WS, output_dir=tmpdir, seed=TEST_SEED,
            )
            results = _run_pipeline(tmpdir, TEST_T, TEST_WS)

        rate = results["BICresult"].mean()
        assert rate > 0.70, f"STATIC_ASM detection rate {rate:.2%} below 70%"

        detected = results[results["BICresult"] == 1]
        if len(detected) > 0:
            # Both per-allele altered fractions should be low (truth is 0 on each).
            # Use median to absorb the few EM runs where init lands in a poor basin.
            assert detected["pi_h0_alt_t2"].median() < 0.25, (
                f"STATIC_ASM pi_h0_alt_t2 median "
                f"{detected['pi_h0_alt_t2'].median():.3f} exceeds 0.25"
            )
            assert detected["pi_h1_alt_t2"].median() < 0.25, (
                f"STATIC_ASM pi_h1_alt_t2 median "
                f"{detected['pi_h1_alt_t2'].median():.3f} exceeds 0.25"
            )


class TestSymmetricTimeEmergent:
    """TIME_EMERGENT (symmetric): both alleles ramp; per-allele π should be roughly equal."""

    def test_detection_and_balanced_pi(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            generate_synthetic_dataset(
                n_regions=TEST_N_REGIONS, reads_per_t=TEST_READS, J=TEST_J, T=TEST_T,
                class_mix={NULL: 0, STATIC_ASM: 0, TIME_EMERGENT_ASM: 1,
                           ASYMMETRIC_TIME_EMERGENT_ASM: 0},
                pi_final=0.6,
                windowsize=TEST_WS, output_dir=tmpdir, seed=TEST_SEED,
            )
            results = _run_pipeline(tmpdir, TEST_T, TEST_WS)

        rate = results["BICresult"].mean()
        assert rate > 0.50, f"Symmetric TIME_EMERGENT detection rate {rate:.2%} below 50%"

        detected = results[results["BICresult"] == 1]
        if len(detected) > 0:
            assert 0.05 < detected["pi_altered_t2"].mean() < 0.90
            # Symmetric truth: |π_h0 - π_h1| should be small on average.
            mean_imbalance = (
                detected["pi_h0_alt_t2"] - detected["pi_h1_alt_t2"]
            ).abs().mean()
            assert mean_imbalance < 0.25, (
                f"Symmetric mean |π_h0 - π_h1| = {mean_imbalance:.3f} exceeds 0.25"
            )


class TestAsymmetricTimeEmergent:
    """ASYMMETRIC: one allele ramps, the other stays at 0. EM should recover the asymmetry.

    Uses T=2 so the full pi_final magnitude lands at t=2 (the time point exposed
    in the per-window output). With T=4 + linear ramp, only ~0.15 within-allele
    altered fraction reaches t=2, which is too weak for reliable allele
    identification under random EM init.
    """

    def test_detection_and_asymmetric_pi_recovered(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest = generate_synthetic_dataset(
                n_regions=12, reads_per_t=TEST_READS, J=TEST_J, T=2,
                class_mix={NULL: 0, STATIC_ASM: 0, TIME_EMERGENT_ASM: 0,
                           ASYMMETRIC_TIME_EMERGENT_ASM: 1},
                pi_final=0.6,
                windowsize=TEST_WS, output_dir=tmpdir, seed=TEST_SEED,
            )
            results = _run_pipeline(tmpdir, T=2, windowsize=TEST_WS)

        rate = results["BICresult"].mean()
        assert rate > 0.50, f"ASYMMETRIC detection rate {rate:.2%} below 50%"

        # Join EM output to manifest truth on (region_idx, coord_evaluated).
        truth = (
            manifest[["region_idx", "coord_evaluated", "ramping_allele"]]
            .drop_duplicates(subset=["region_idx", "coord_evaluated"])
        )
        merged = results.merge(truth, on=["region_idx", "coord_evaluated"], how="left")
        detected = merged[merged["BICresult"] == 1]

        if len(detected) > 0:
            # Asymmetry magnitude: |π_h0 - π_h1| should be substantial.
            mean_imbalance = (
                detected["pi_h0_alt_t2"] - detected["pi_h1_alt_t2"]
            ).abs().mean()
            assert mean_imbalance > 0.10, (
                f"ASYMMETRIC mean |π_h0 - π_h1| = {mean_imbalance:.3f} below 0.10"
            )

            # Direction: EM's elevated π should land on the truly-ramping allele
            # in the majority of detected windows.
            ramps_h1 = detected[detected["ramping_allele"] == 1]
            ramps_h0 = detected[detected["ramping_allele"] == 0]
            correct = 0
            total = 0
            if len(ramps_h1) > 0:
                correct += int((ramps_h1["pi_h1_alt_t2"] > ramps_h1["pi_h0_alt_t2"]).sum())
                total += len(ramps_h1)
            if len(ramps_h0) > 0:
                correct += int((ramps_h0["pi_h0_alt_t2"] > ramps_h0["pi_h1_alt_t2"]).sum())
                total += len(ramps_h0)
            if total > 0:
                accuracy = correct / total
                # Loose threshold: random EM init causes occasional misattribution
                # even with strong t=2 signal.
                assert accuracy > 0.55, (
                    f"ASYMMETRIC ramping-allele identification accuracy "
                    f"{accuracy:.2%} below 55%"
                )


class TestManifestFormat:
    """Manifest must carry per-allele truth and the new ramping_allele column."""

    def test_columns_and_asymmetric_truth(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest = generate_synthetic_dataset(
                n_regions=12, reads_per_t=10, J=15, T=4,
                class_mix={NULL: 1, STATIC_ASM: 1, TIME_EMERGENT_ASM: 1,
                           ASYMMETRIC_TIME_EMERGENT_ASM: 2},
                pi_final=0.6,
                output_dir=tmpdir, seed=TEST_SEED,
            )

        for col in ["ramping_allele", "pi_h0_alt_t1", "pi_h1_alt_t4"]:
            assert col in manifest.columns, f"Missing column: {col}"

        # Boundary: pi at t=1 always 0 for every class.
        assert (manifest["pi_h0_alt_t1"] == 0).all()
        assert (manifest["pi_h1_alt_t1"] == 0).all()

        # ramping_allele: -1 for NULL/STATIC/TIME_EMERGENT, 0 or 1 for ASYMMETRIC.
        non_asym = manifest[manifest["true_class"] != ASYMMETRIC_TIME_EMERGENT_ASM]
        assert (non_asym["ramping_allele"] == -1).all()
        asym = manifest[manifest["true_class"] == ASYMMETRIC_TIME_EMERGENT_ASM]
        assert asym["ramping_allele"].isin([0, 1]).all()

        # Asymmetric truth: exactly one of pi_h0_alt_t4, pi_h1_alt_t4 is > 0.
        for _, row in asym.iterrows():
            assert (row["pi_h0_alt_t4"] > 0) ^ (row["pi_h1_alt_t4"] > 0), (
                "ASYMMETRIC region must have exactly one allele ramping at t=T"
            )

        # Symmetric TIME_EMERGENT truth: both alleles equal at t=T.
        sym = manifest[manifest["true_class"] == TIME_EMERGENT_ASM]
        np.testing.assert_allclose(sym["pi_h0_alt_t4"], sym["pi_h1_alt_t4"])


class TestNoHReadsJointDraw:
    """Catch regression of Gap 2: the noH sampler must do a joint (h, k) draw."""

    def test_marginal_consistency_under_asymmetry(self):
        """With pi_alt_h0=0, pi_alt_h1=0.6, no noH read should be (h=0, k=1).

        Truth probabilities for noH reads:
            P(h=0, k=0) = 0.5,    P(h=0, k=1) = 0.0
            P(h=1, k=0) = 0.2,    P(h=1, k=1) = 0.3

        Engineered θ makes each (h, k) identifiable from one read's mean
        methylation:
            (h=0, k=0): θ=0  → meth = 0 exactly
            (h=1, k=0): θ=1  → meth = 1 exactly
            (h=0, k=1): θ=0.7 → meth ≈ 0.7 (noisy, but ∈ [0.5, 1) w.h.p.)
            (h=1, k=1): θ=0.3 → meth ≈ 0.3 (noisy, but ∈ (0, 0.5) w.h.p.)
        Under a broken factorized (h,k) draw, P(h=0, k=1) would be nonzero;
        under the correct joint draw it is exactly 0.
        """
        rng = np.random.RandomState(TEST_SEED)
        J = 30
        coords = np.arange(1000, 1000 + 10 * J, 10)
        theta = np.zeros((J, 2, 2))
        theta[:, 0, 0] = 0.0
        theta[:, 0, 1] = 0.7
        theta[:, 1, 0] = 1.0
        theta[:, 1, 1] = 0.3

        n = 5000
        rows = generate_reads_for_region(
            coords, theta,
            pi_alt_t_h0=0.0, pi_alt_t_h1=0.6,
            n_reads=n, chrom="chrTEST", rng=rng,
            hap_probs=(0.0, 0.0, 1.0),  # only noH reads
            min_cpg_span=J,
        )
        df = pd.DataFrame(rows)
        meth_means = np.array(
            [np.mean([int(s) for s in r.split(",")]) for r in df["status"]]
        )

        n_h0_k0 = int(np.sum(meth_means == 0.0))
        n_h1_k0 = int(np.sum(meth_means == 1.0))
        n_h0_k1 = int(np.sum((meth_means >= 0.5) & (meth_means < 1.0)))
        n_h1_k1 = int(np.sum((meth_means > 0.0) & (meth_means < 0.5)))

        # Expected: 0.5, 0.2, 0.0, 0.3.
        assert abs(n_h0_k0 / n - 0.5) < 0.03, f"P(h=0,k=0) ≈ {n_h0_k0/n:.3f}"
        assert abs(n_h1_k0 / n - 0.2) < 0.03, f"P(h=1,k=0) ≈ {n_h1_k0/n:.3f}"
        # KEY assertion — Gap-2 regression test:
        # Correct joint draw → 0. Broken factorized draw with averaged
        # pi_alt_t = 0.3 would give 0.5 * 0.3 = 0.15 here.
        assert n_h0_k1 / n < 0.02, (
            f"(h=0, k=1) noH reads = {n_h0_k1/n:.3f}, expected ~0 — "
            f"Gap-2 regression: noH (h,k) draw is no longer joint."
        )
        assert abs(n_h1_k1 / n - 0.3) < 0.04, f"P(h=1,k=1) ≈ {n_h1_k1/n:.3f}"
