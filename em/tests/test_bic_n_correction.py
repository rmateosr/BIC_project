# ABOUTME: Verify BIC penalty and reported total_reads use the per-window
# ABOUTME: non-NaN read count, not the inflated per-region read count.

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from em.EMfunctions_timebased_diploid_aware_model import EMBIC_bin_path_timebased


WINDOWSIZE = 10


def _make_region_file(path, n_reads, region_cpgs, window_cpgs, p_meth, rng):
    """Write a synthetic region where only `len(window_cpgs)` of the reads
    have calls at `window_cpgs`. The remaining reads cover only CpGs outside
    that window. All reads cover at least one CpG in the region (so they
    survive the region-level all-NaN filter).

    Used to assert that the BIC `n` and reported `total_reads` reflect the
    per-window coverage, not the per-region coverage.
    """
    rows = []
    n_covered = len(window_cpgs)
    # First n_covered reads cover only the window CpGs.
    for _ in range(n_covered):
        statuses = rng.binomial(1, p_meth, size=len(window_cpgs)).astype(int)
        rows.append({
            "readlabel": f"r{len(rows)}",
            "haplotype": "noH",
            "chrom": "chrTEST",
            "startcoord": ",".join(str(c) for c in window_cpgs),
            "status": ",".join(str(s) for s in statuses),
        })
    # Remaining reads cover only the non-window CpGs.
    non_window = [c for c in region_cpgs if c not in window_cpgs]
    for _ in range(n_reads - n_covered):
        statuses = rng.binomial(1, p_meth, size=len(non_window)).astype(int)
        rows.append({
            "readlabel": f"r{len(rows)}",
            "haplotype": "noH",
            "chrom": "chrTEST",
            "startcoord": ",".join(str(c) for c in non_window),
            "status": ",".join(str(s) for s in statuses),
        })
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False)


def test_total_reads_and_bic_n_are_per_window(tmp_path):
    """1000 reads in the region but only 30 reads cover the target window.
    The output column `total_reads` and the BIC penalty must use 30, not 1000."""
    rng = np.random.default_rng(42)
    region_cpgs = list(range(1000, 1000 + 200))  # 200 CpGs in the region
    target_window = region_cpgs[50:60]            # 10 contiguous CpGs

    paths = {}
    for t in (1, 2):
        p = tmp_path / f"region_t{t}.tsv"
        n_covered_at_window = 15  # per timepoint; sum = 30
        _make_region_file(
            p,
            n_reads=500,
            region_cpgs=region_cpgs,
            window_cpgs=target_window,
            p_meth=0.3 if t == 1 else 0.7,
            rng=rng,
        )
        # Force the per-window-covered count: rewrite manually to control it.
        # (The helper above gives n_covered = len(window_cpgs) = 10; override.)
        rows = []
        for _ in range(n_covered_at_window):
            statuses = rng.binomial(1, 0.5, size=len(target_window)).astype(int)
            rows.append({
                "readlabel": f"r{len(rows)}",
                "haplotype": "noH",
                "chrom": "chrTEST",
                "startcoord": ",".join(str(c) for c in target_window),
                "status": ",".join(str(s) for s in statuses),
            })
        non_window = [c for c in region_cpgs if c not in target_window]
        for _ in range(500 - n_covered_at_window):
            statuses = rng.binomial(1, 0.5, size=len(non_window)).astype(int)
            rows.append({
                "readlabel": f"r{len(rows)}",
                "haplotype": "noH",
                "chrom": "chrTEST",
                "startcoord": ",".join(str(c) for c in non_window),
                "status": ",".join(str(s) for s in statuses),
            })
        pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
        paths[t] = str(p)

    df = EMBIC_bin_path_timebased(paths, windowsize=WINDOWSIZE)

    # Locate the row for the target window's leading CpG.
    target_start = target_window[0]
    target_row = df[df["coord_evaluated"] == target_start]
    assert len(target_row) == 1, (
        f"Expected exactly one row at coord {target_start}; got {len(target_row)}"
    )
    reported = int(target_row["total_reads"].iloc[0])

    expected_per_window = 15 + 15  # 15 per timepoint
    region_total = 500 + 500

    assert reported == expected_per_window, (
        f"total_reads at target window should be per-window count "
        f"({expected_per_window}), got {reported}. Region-level inflated total "
        f"would have been {region_total}."
    )
    # Sanity: must not be the inflated region count.
    assert reported != region_total


def test_total_reads_column_varies_across_windows(tmp_path):
    """The `total_reads` column must vary across windows when window-level
    coverage genuinely varies — i.e. it is not a per-region constant."""
    rng = np.random.default_rng(0)
    region_cpgs = list(range(1, 1 + 30))  # 30 CpGs → 21 sliding windows of size 10

    # Reads: some cover all 30, some cover only the first 10, some only the last 10.
    rows_by_t = {1: [], 2: []}
    for t in (1, 2):
        for _ in range(20):  # full-coverage reads
            statuses = rng.binomial(1, 0.4, size=30).astype(int)
            rows_by_t[t].append({
                "readlabel": f"full_{len(rows_by_t[t])}",
                "haplotype": "noH",
                "chrom": "chrTEST",
                "startcoord": ",".join(str(c) for c in region_cpgs),
                "status": ",".join(str(s) for s in statuses),
            })
        for _ in range(10):  # only first-10 CpGs
            statuses = rng.binomial(1, 0.4, size=10).astype(int)
            rows_by_t[t].append({
                "readlabel": f"first_{len(rows_by_t[t])}",
                "haplotype": "noH",
                "chrom": "chrTEST",
                "startcoord": ",".join(str(c) for c in region_cpgs[:10]),
                "status": ",".join(str(s) for s in statuses),
            })

    paths = {}
    for t in (1, 2):
        p = tmp_path / f"region_t{t}.tsv"
        pd.DataFrame(rows_by_t[t]).to_csv(p, sep="\t", index=False)
        paths[t] = str(p)

    df = EMBIC_bin_path_timebased(paths, windowsize=WINDOWSIZE)

    first_window_n = int(df.iloc[0]["total_reads"])    # covers CpGs[0:10]
    last_window_n = int(df.iloc[-1]["total_reads"])    # covers CpGs[20:30]

    # First window has full-coverage reads (20 per t = 40) + first-only reads
    # (10 per t = 20). Last window has only full-coverage reads (40).
    assert first_window_n == 60
    assert last_window_n == 40
    assert first_window_n != last_window_n, (
        "total_reads must vary across windows when coverage varies"
    )


def test_bic_penalty_uses_per_window_n(tmp_path):
    """Direct check: with 30 reads at the window and 1000 in the region, the
    BIC values produced must be consistent with `log(30)`-scaled penalties,
    not `log(1000)`-scaled penalties.

    We assert this by comparing the BIC against an analytically reconstructed
    value computed under the assumption n=30."""
    from em.EMfunctions_timebased_diploid_aware_model import (
        compute_null_BIC,
        Reshape2Matrix,
    )

    rng = np.random.default_rng(7)
    region_cpgs = list(range(1, 1 + 100))
    target_window = region_cpgs[40:50]

    rows = []
    n_covered = 30
    for _ in range(n_covered):
        statuses = rng.binomial(1, 0.5, size=10).astype(int)
        rows.append({
            "readlabel": f"r{len(rows)}",
            "haplotype": "noH",
            "chrom": "chrTEST",
            "startcoord": ",".join(str(c) for c in target_window),
            "status": ",".join(str(s) for s in statuses),
        })
    non_window = [c for c in region_cpgs if c not in target_window]
    for _ in range(1000 - n_covered):
        statuses = rng.binomial(1, 0.5, size=len(non_window)).astype(int)
        rows.append({
            "readlabel": f"r{len(rows)}",
            "haplotype": "noH",
            "chrom": "chrTEST",
            "startcoord": ",".join(str(c) for c in non_window),
            "status": ",".join(str(s) for s in statuses),
        })

    df = pd.DataFrame(rows)
    p = tmp_path / "region_t1.tsv"
    df.to_csv(p, sep="\t", index=False)
    # Run a 1-timepoint pipeline-equivalent by hand: slice to the target window
    coords_array = np.array(sorted(region_cpgs))
    X_full = Reshape2Matrix(df, coords_array)
    # Drop reads that are all-NaN across the whole region (none here, but mirror pipeline)
    X_full = X_full[~np.isnan(X_full).all(axis=1)]

    # Find window starting at target_window[0]
    j0 = list(coords_array).index(target_window[0])
    X_w = X_full[:, j0:j0 + WINDOWSIZE]
    valid = ~np.isnan(X_w).all(axis=1)
    X_window = X_w[valid]

    assert X_window.shape[0] == n_covered, (
        f"per-window filter should leave {n_covered} reads, got {X_window.shape[0]}"
    )

    bic = compute_null_BIC({1: X_window}, T=1)
    # The BIC penalty for the null model is J * log(n). With n=30 and J=10:
    # penalty = 10 * log(30) ≈ 34.0; with the inflated n=1000 it would be
    # ≈ 69.1. Reconstruct the penalty by re-computing the LL portion.
    observed = ~np.isnan(X_window)
    X_clean = np.where(np.isnan(X_window), 0.0, X_window)
    total_meth = X_clean.sum(axis=0)
    total_obs = observed.astype(float).sum(axis=0)
    PSEUDO = 1e-10
    theta_null = total_meth / (total_obs + PSEUDO)
    log_t = np.log(theta_null + PSEUDO)
    log_1mt = np.log(1 - theta_null + PSEUDO)
    ll = float((X_clean * log_t).sum() + ((observed.astype(float) - X_clean) * log_1mt).sum())

    expected_with_window_n = 10 * np.log(n_covered) - 2 * ll
    expected_with_region_n = 10 * np.log(1000) - 2 * ll

    assert abs(bic - expected_with_window_n) < 1e-6, (
        f"BIC ({bic}) does not match the per-window-n expectation "
        f"({expected_with_window_n})"
    )
    assert abs(bic - expected_with_region_n) > 1.0, (
        f"BIC ({bic}) is suspiciously close to the inflated region-n "
        f"expectation ({expected_with_region_n})"
    )
