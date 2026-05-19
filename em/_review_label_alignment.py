# ABOUTME: Review-only test for item 4 — does per-CpG alignment change recovery numbers?
# ABOUTME: Compares global-swap (current) vs per-CpG argmin swap on per-(j,h,k) MAE.
"""Test for partial-label-swap risk in `_align_M2_components`.

The prompt's concern: EM might converge with the k=0/k=1 labels flipped on
some CpGs but not others within the same region. The current alignment picks
ONE swap decision for the whole region; per-CpG argmin would pick the optimal
swap per CpG.

All three MAEs below use abs-first, mean-second (over j × h × k), so they
are directly comparable:
  - identity_mae = mean |theta_M2[j,h,k] - truth[h,k]|             (no swap)
  - global_mae   = min over {identity, full-swap} of identity_mae  (current)
  - per_cpg_mae  = mean_j of min over {identity, swap} of mean_(h,k) |..|
                                                        per CpG j (best)

If per_cpg_mae << global_mae, partial swap is real and the report's recovery
numbers are biased high. If per_cpg_mae == global_mae, partial swap is
structurally impossible — consistent with the EM update sharing one γ[i,h,k]
across all CpGs in a window.
"""
import glob
import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from parameter_recovery import fit_M2, load_region
from synthetic.generate_synthetic import (
    ASYMMETRIC_TIME_EMERGENT_ASM,
    NULL,
    STATIC_ASM,
    TIME_EMERGENT_ASM,
    generate_synthetic_dataset,
)


def per_jhk_mae(theta_M2, truth_arr):
    """Both alignments computed in the same per-(j,h,k) abs-first MAE units."""
    err_identity_jhk = np.abs(theta_M2 - truth_arr[None, :, :])           # (J,2,2)
    err_swap_jhk     = np.abs(theta_M2[:, :, ::-1] - truth_arr[None, :, :])  # (J,2,2)

    identity_mae = float(err_identity_jhk.mean())
    full_swap_mae = float(err_swap_jhk.mean())
    global_mae = min(identity_mae, full_swap_mae)

    # Per-CpG argmin: for each j, pick the swap orientation that minimizes
    # the (h,k)-mean error at that CpG; then average that minimum across j.
    cpg_id   = err_identity_jhk.mean(axis=(1, 2))   # (J,)
    cpg_swap = err_swap_jhk.mean(axis=(1, 2))       # (J,)
    per_cpg_mae = float(np.minimum(cpg_id, cpg_swap).mean())

    n_cpgs_prefer_swap = int((cpg_swap < cpg_id).sum())
    return identity_mae, full_swap_mae, global_mae, per_cpg_mae, n_cpgs_prefer_swap


def main():
    rows = []
    for seed in (2026, 2027, 2028):
        for cls in (NULL, STATIC_ASM, TIME_EMERGENT_ASM, ASYMMETRIC_TIME_EMERGENT_ASM):
            with tempfile.TemporaryDirectory() as out:
                manifest = generate_synthetic_dataset(
                    n_regions=8,
                    reads_per_t=100,
                    J=20,
                    T=2,
                    class_mix={cls: 1},
                    pi_final=0.9,
                    output_dir=out,
                    seed=seed,
                )
                files_per_t = {
                    t: sorted(glob.glob(os.path.join(out, f"t{t}", "methylationfraction_*_.tsv")))
                    for t in (1, 2)
                }
                for ridx in range(len(files_per_t[1])):
                    paths = {t: files_per_t[t][ridx] for t in (1, 2)}
                    X_by_t, tags_by_t, _ = load_region(paths)
                    truth = manifest[manifest["region_idx"] == ridx].iloc[0]
                    truth_arr = np.array([
                        [truth["true_theta_h0_k0"], truth["true_theta_h0_k1"]],
                        [truth["true_theta_h1_k0"], truth["true_theta_h1_k1"]],
                    ])

                    theta_M2_raw, _ = fit_M2(X_by_t, tags_by_t, T=2)

                    id_mae, swap_mae, g_mae, p_mae, n_swap = per_jhk_mae(theta_M2_raw, truth_arr)
                    rows.append(dict(
                        seed=seed, true_class=cls, region_idx=ridx,
                        identity_mae=id_mae, full_swap_mae=swap_mae,
                        global_mae=g_mae, per_cpg_mae=p_mae,
                        gap=g_mae - p_mae,
                        n_cpgs_prefer_swap=n_swap,
                    ))

    df = pd.DataFrame(rows)

    print("=== Per-class summary (3 seeds × 8 regions = 24 regions per class) ===")
    summary = (
        df.groupby("true_class")[["global_mae", "per_cpg_mae", "gap", "n_cpgs_prefer_swap"]]
        .agg(["mean", "max"])
        .round(4)
    )
    print(summary.to_string())
    print()

    # Partial-swap regions: per-CpG swap-counts away from {0, J}.
    J = 20
    df["swap_pattern"] = df["n_cpgs_prefer_swap"].apply(
        lambda n: "all_id" if n == 0 else ("all_swap" if n == J else "partial")
    )
    pat_counts = df["swap_pattern"].value_counts()
    print("=== Per-CpG swap-pattern counts across 96 regions ===")
    print(pat_counts.to_string())
    print()

    print("=== Largest per-CpG vs global gaps ===")
    print(df.nlargest(10, "gap")[
        ["seed", "true_class", "region_idx", "global_mae", "per_cpg_mae", "gap",
         "n_cpgs_prefer_swap"]
    ].to_string(index=False))
    print()

    max_gap = df["gap"].max()
    mean_gap = df["gap"].mean()
    print(f"max gap (global − per_cpg) across 96 regions: {max_gap:.5f}")
    print(f"mean gap                                   : {mean_gap:.5f}")
    print(
        "\nINTERPRETATION:\n"
        "  - per_cpg_mae ≤ global_mae always (optimal per-CpG ≤ best of two global).\n"
        "  - If max gap is small (e.g. < 0.005), partial swap is essentially absent\n"
        "    and the existing alignment is correct.\n"
        "  - If max gap > ~0.01, partial swap is real; report's MAE is biased high.\n"
        "  - Many 'partial' patterns with tiny gap = ties due to noise, not real swap."
    )


if __name__ == "__main__":
    main()
