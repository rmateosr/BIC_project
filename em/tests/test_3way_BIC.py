# ABOUTME: 3-way BIC (M0/M1/M2) tests on synthetic data — NULL→M0, STATIC→M1, ASYMMETRIC→M2.
# ABOUTME: Loose thresholds; EM(...) reseeds via np.random.seed() so single runs are stochastic.
"""3-way BIC sanity tests.

Generates synthetic data for each of the four region classes and checks that
the new BIC_3way_winner column lands in the expected model column (0=M0,
1=M1, 2=M2) at a reasonable rate.

Thresholds are loose for the same reason as test_synthetic_diploid_aware.py:
EMfunctions_timebased_diploid_aware_model.EM() reseeds via np.random.seed()
on every call, so identical generator seeds still yield stochastic EM
outputs. Run pytest 3× before assuming a green is real.

Symmetric TIME_EMERGENT_ASM is intentionally not tested here: bulk drift on
both alleles is biologically ambiguous (not really ASM) and the prompt doc
flagged its expected M1-vs-M2 outcome as an open question.
"""

import glob
import os
import tempfile

import numpy as np
import pandas as pd

from EMfunctions_timebased_diploid_aware_model import (
    EMBIC_bin_path_timebased,
    compute_M1_BIC,
)
from synthetic.generate_synthetic import (
    ASYMMETRIC_TIME_EMERGENT_ASM,
    NULL,
    STATIC_ASM,
    TIME_EMERGENT_ASM,
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


# Test sizing — small enough for fast pytest, big enough for stable rates.
TEST_N_REGIONS = 8
TEST_READS = 50
TEST_J = 20
TEST_T = 2
TEST_WS = 10
TEST_SEED = 2026


class Test3wayBIC_OnNull:
    """NULL truth: BIC_3way_winner should land on column 0 (M0) most of the time."""

    def test_M0_dominates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            generate_synthetic_dataset(
                n_regions=TEST_N_REGIONS, reads_per_t=TEST_READS, J=TEST_J, T=TEST_T,
                class_mix={NULL: 1, STATIC_ASM: 0, TIME_EMERGENT_ASM: 0,
                           ASYMMETRIC_TIME_EMERGENT_ASM: 0},
                windowsize=TEST_WS, output_dir=tmpdir, seed=TEST_SEED,
            )
            results = _run_pipeline(tmpdir, TEST_T, TEST_WS)

        m0_rate = (results["BIC_3way_winner"] == 0).mean()
        assert m0_rate > 0.80, (
            f"NULL: M0-winner rate {m0_rate:.2%} below 80% — "
            f"distribution {results['BIC_3way_winner'].value_counts().to_dict()}"
        )


class Test3wayBIC_OnStatic:
    """STATIC truth: BIC_3way_winner should land on column 1 (M1) most of the time.

    This is the headline empirical question — does the BIC penalty correctly
    punish M2's spurious altered-component fit on stable allele-specific
    methylation? See docs/NEXT_SESSION_PROMPT_3way_BIC.md.
    """

    def test_M1_dominates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            generate_synthetic_dataset(
                n_regions=TEST_N_REGIONS, reads_per_t=TEST_READS, J=TEST_J, T=TEST_T,
                class_mix={NULL: 0, STATIC_ASM: 1, TIME_EMERGENT_ASM: 0,
                           ASYMMETRIC_TIME_EMERGENT_ASM: 0},
                windowsize=TEST_WS, output_dir=tmpdir, seed=TEST_SEED,
            )
            results = _run_pipeline(tmpdir, TEST_T, TEST_WS)

        m1_rate = (results["BIC_3way_winner"] == 1).mean()
        assert m1_rate > 0.70, (
            f"STATIC: M1-winner rate {m1_rate:.2%} below 70% — "
            f"distribution {results['BIC_3way_winner'].value_counts().to_dict()}"
        )


class Test3wayBIC_OnAsymmetric:
    """ASYMMETRIC truth: BIC_3way_winner should fire on M2 in a meaningful
    fraction of windows.

    Empirical finding (32-region mixed smoke, 2026-04-27 follow-up): the BIC
    penalty for M2 over M1 is `(2J + 2T - 2) * log(n)` ≈ 100 units at J=10,
    T=2, n≈100. The LL gain from a single ramping allele at pi_final=0.9 is
    in the same range, so M2 fires only ~25% of the time on ASYMMETRIC truth
    even with strong signal. M1 absorbs most of the asymmetric signal because
    one-allele drift looks similar to allele-specific θ when averaged over
    time. This is a property of BIC, not a bug — captured here as a
    permissive lower bound until/unless a mitigation (smaller M2, prior on
    π, or AIC for M1-vs-M2) is added.
    """

    def test_M2_fires(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            generate_synthetic_dataset(
                n_regions=12, reads_per_t=TEST_READS, J=TEST_J, T=TEST_T,
                class_mix={NULL: 0, STATIC_ASM: 0, TIME_EMERGENT_ASM: 0,
                           ASYMMETRIC_TIME_EMERGENT_ASM: 1},
                pi_final=0.9,
                windowsize=TEST_WS, output_dir=tmpdir, seed=TEST_SEED,
            )
            results = _run_pipeline(tmpdir, TEST_T, TEST_WS)

        m2_rate = (results["BIC_3way_winner"] == 2).mean()
        assert m2_rate > 0.15, (
            f"ASYMMETRIC: M2-winner rate {m2_rate:.2%} below 15% — "
            f"distribution {results['BIC_3way_winner'].value_counts().to_dict()}"
        )


class TestBIC_M1_finite:
    """Sanity: compute_M1_BIC returns a finite float on small fixed input."""

    def test_M1_BIC_finite_on_random_input(self):
        rng = np.random.RandomState(0)
        T, J = 2, 10
        X_by_t = {t: (rng.rand(40, J) > 0.5).astype(float) for t in (1, 2)}
        tags_by_t = {
            t: np.array(["H1"] * 15 + ["H2"] * 15 + ["noH"] * 10) for t in (1, 2)
        }

        bic_M1, n_iters = compute_M1_BIC(X_by_t, tags_by_t, T)
        assert np.isfinite(bic_M1), f"BIC_M1 is non-finite: {bic_M1}"
        assert n_iters >= 1
