# ABOUTME: End-to-end 4-way BIC pipeline tests through BIC_algorithm_timebased_4way.
# ABOUTME: Confirms shared/asm/shared_time/asm_time each win on their intended truth class.

import glob
import os
import sys
import tempfile

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from BIC_algorithm_timebased_4way import EMBIC_bin_path_timebased_4way
from synthetic.generate_synthetic import (
    ASYMMETRIC_TIME_EMERGENT_ASM,
    NULL,
    STATIC_ASM,
    TIME_EMERGENT_ASM,
    generate_synthetic_dataset,
)


TEST_N_REGIONS = 8
TEST_READS = 60
TEST_J = 20
TEST_T = 2
TEST_WS = 10
TEST_SEED = 2026


def _run_pipeline(output_dir, T, windowsize=TEST_WS):
    files_per_t = {}
    for t in range(1, T + 1):
        files_per_t[t] = sorted(
            glob.glob(os.path.join(output_dir, f"t{t}", "methylationfraction_*_.tsv"))
        )
    n_regions = len(files_per_t[1])
    out = []
    for ridx in range(n_regions):
        paths = {t: files_per_t[t][ridx] for t in range(1, T + 1)}
        df = EMBIC_bin_path_timebased_4way(paths, windowsize=windowsize)
        df["region_idx"] = ridx
        out.append(df)
    return pd.concat(out, ignore_index=True)


def test_shared_wins_on_NULL():
    with tempfile.TemporaryDirectory() as tmp:
        generate_synthetic_dataset(
            n_regions=TEST_N_REGIONS, reads_per_t=TEST_READS, J=TEST_J, T=TEST_T,
            class_mix={NULL: 1, STATIC_ASM: 0, TIME_EMERGENT_ASM: 0,
                       ASYMMETRIC_TIME_EMERGENT_ASM: 0},
            windowsize=TEST_WS, output_dir=tmp, seed=TEST_SEED,
        )
        results = _run_pipeline(tmp, TEST_T)
    rate = (results["BIC_4way_winner"] == "shared").mean()
    assert rate > 0.80, (
        f"NULL: shared-winner rate {rate:.2%} below 80% — "
        f"distribution {results['BIC_4way_winner'].value_counts().to_dict()}"
    )


def test_asm_wins_on_STATIC():
    with tempfile.TemporaryDirectory() as tmp:
        generate_synthetic_dataset(
            n_regions=TEST_N_REGIONS, reads_per_t=TEST_READS, J=TEST_J, T=TEST_T,
            class_mix={NULL: 0, STATIC_ASM: 1, TIME_EMERGENT_ASM: 0,
                       ASYMMETRIC_TIME_EMERGENT_ASM: 0},
            windowsize=TEST_WS, output_dir=tmp, seed=TEST_SEED,
        )
        results = _run_pipeline(tmp, TEST_T)
    rate = (results["BIC_4way_winner"] == "asm").mean()
    assert rate > 0.60, (
        f"STATIC: asm-winner rate {rate:.2%} below 60% — "
        f"distribution {results['BIC_4way_winner'].value_counts().to_dict()}"
    )


def _build_pure_class_c_tsvs(tmp_dir, n_regions, reads_per_t, J, T, pi_alt_final, seed):
    """Write per-time TSVs for pure-Class-C regions (shared baseline + symmetric drift).

    The shared generator's TIME_EMERGENT_ASM bakes in allele difference at
    k=0 — see em/synthetic/generate_synthetic.py:117-128. Pure Class C
    requires both alleles to share theta at every k. We synthesise inline
    in the format the 4-way driver consumes.
    """
    import uuid

    rng = np.random.default_rng(seed)
    for t in range(1, T + 1):
        os.makedirs(os.path.join(tmp_dir, f"t{t}"), exist_ok=True)

    pi_alt_by_t = np.linspace(0.0, pi_alt_final, T) if T > 1 else np.array([0.0])
    offset = 1_000_000
    for ridx in range(n_regions):
        coords = np.cumsum(rng.integers(80, 250, size=J)) + offset
        offset = int(coords[-1]) + 50_000
        theta_k0 = rng.uniform(0.65, 0.85, size=J)
        theta_k1 = rng.uniform(0.10, 0.30, size=J)

        fname = f"methylationfraction_{int(coords[0])}_{int(coords[-1])}_.tsv"
        for t in range(1, T + 1):
            pi_alt = pi_alt_by_t[t - 1]
            rows = []
            for _ in range(reads_per_t):
                k = 1 if rng.random() < pi_alt else 0
                hap = rng.choice(["H1", "H2", "noH"], p=[0.4, 0.4, 0.2])
                theta = theta_k1 if k == 1 else theta_k0
                # Full-coverage reads keep the test fast and the signal clean.
                status = (rng.random(J) < theta).astype(int)
                rows.append({
                    "readlabel": uuid.uuid4().hex[:16],
                    "chrom": "chrSYN",
                    "strand": "+",
                    "startcoord": ",".join(str(c) for c in coords),
                    "status": ",".join(str(s) for s in status),
                    "haplotype": hap,
                })
            pd.DataFrame(rows).to_csv(
                os.path.join(tmp_dir, f"t{t}", fname), sep="\t", index=False,
            )


def test_shared_time_wins_on_pure_class_c():
    """The headline test: when truth is pure Class C, shared_time wins BIC."""
    with tempfile.TemporaryDirectory() as tmp:
        _build_pure_class_c_tsvs(
            tmp, n_regions=6, reads_per_t=120, J=20, T=2,
            pi_alt_final=0.6, seed=TEST_SEED,
        )
        results = _run_pipeline(tmp, T=2)
    rate = (results["BIC_4way_winner"] == "shared_time").mean()
    assert rate > 0.60, (
        f"Pure Class C: shared_time-winner rate {rate:.2%} below 60% — "
        f"distribution {results['BIC_4way_winner'].value_counts().to_dict()}"
    )


def test_asm_time_still_fires_on_ASYMMETRIC():
    """ASYMMETRIC (allele-specific drift) should still pull asm_time wins.

    Loose lower bound mirrors the existing test_3way_BIC test: BIC penalty
    for asm_time is large, so it wins ~15-25% even at strong signal. The
    point of this test is to confirm shared_time doesn't *steal* every
    asymmetric region from asm_time.
    """
    with tempfile.TemporaryDirectory() as tmp:
        generate_synthetic_dataset(
            n_regions=12, reads_per_t=TEST_READS, J=TEST_J, T=TEST_T,
            class_mix={NULL: 0, STATIC_ASM: 0, TIME_EMERGENT_ASM: 0,
                       ASYMMETRIC_TIME_EMERGENT_ASM: 1},
            pi_final=0.9,
            windowsize=TEST_WS, output_dir=tmp, seed=TEST_SEED,
        )
        results = _run_pipeline(tmp, TEST_T)
    rate = (results["BIC_4way_winner"] == "asm_time").mean()
    assert rate > 0.10, (
        f"ASYMMETRIC: asm_time-winner rate {rate:.2%} below 10% — "
        f"distribution {results['BIC_4way_winner'].value_counts().to_dict()}"
    )


def test_BICsinglecomp_identical_to_shared_null():
    """The diploid null and the no-h null compute the same number.

    The driver emits the diploid value as `BICsinglecomp` for backward
    compatibility; this asserts that's not silently lossy.
    """
    from EMfunctions_timebased_diploid_aware_model import compute_null_BIC as dip_null
    from EMfunctions_timebased_shared_time_emergent_model import compute_shared_BIC

    rng = np.random.default_rng(0)
    X_by_t = {t: rng.binomial(1, 0.4, size=(50, 10)).astype(float) for t in (1, 2)}
    bic_dip = dip_null(X_by_t, T=2)
    bic_shared = compute_shared_BIC(X_by_t, T=2)
    assert bic_dip == pytest.approx(bic_shared, rel=1e-12)
