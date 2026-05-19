# ABOUTME: Fits M0/M1/M2 to each synthetic region and reports per-class θ and π recovery error.
# ABOUTME: Run from em/. CLI flags: --n-per-class, --reads-per-t, --J, --T, --pi-final, --seed.
"""Parameter-recovery experiment for the 3-way BIC models.

For each region in a synthetic dataset, fits all three models and compares the
recovered parameters to the truth carried in the manifest:
    M0 — shared theta_j (closed-form mean methylation across all reads × time)
    M1 — per-allele theta_j at k=0 (constrained M2 with π[:, 1, :] = 0)
    M2 — full diploid-aware time-based mixture (production EM)

Per-class summary reports mean |θ_fit - θ_true| and signed π error at t=T.
A model that "catches" the inner parameters used in generation should land
near zero on the columns relevant to its own structure.
"""

import argparse
import glob
import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from EMfunctions_timebased_diploid_aware_model import (
    EM,
    PSEUDO,
    Reshape2Matrix,
    _observed_loglikelihood,
    doE,
    doM_theta,
)
from synthetic.generate_synthetic import (
    ASYMMETRIC_TIME_EMERGENT_ASM,
    NULL,
    STATIC_ASM,
    TIME_EMERGENT_ASM,
    generate_synthetic_dataset,
)


def load_region(paths_by_t):
    """Load per-t TSVs for one region; return X_by_t, tags_by_t, coords_array.

    Mirrors the matrix-construction step at the top of EMBIC_bin_path_timebased.
    """
    T = len(paths_by_t)
    data_by_t = {t: pd.read_csv(p, sep="\t") for t, p in paths_by_t.items()}

    all_coords = set()
    for df in data_by_t.values():
        coords_str = ",".join(df["startcoord"].astype(str))
        all_coords.update(int(float(c)) for c in coords_str.split(","))
    coords_array = np.array(sorted(all_coords))

    X_by_t = {}
    tags_by_t = {}
    for t, df in data_by_t.items():
        X = Reshape2Matrix(df, coords_array)
        tags = np.array(df["haplotype"])
        valid = ~np.isnan(X).all(axis=1)
        X_by_t[t] = X[valid]
        tags_by_t[t] = tags[valid]
    return X_by_t, tags_by_t, coords_array


def fit_M0(X_by_t):
    """Fit M0: shared theta_j across alleles and times. Returns theta of shape (J,)."""
    any_t = next(iter(X_by_t))
    J = X_by_t[any_t].shape[1]
    total_meth = np.zeros(J)
    total_obs = np.zeros(J)
    for X in X_by_t.values():
        observed = ~np.isnan(X)
        X_clean = np.where(np.isnan(X), 0, X)
        total_meth += X_clean.sum(axis=0)
        total_obs += observed.astype(float).sum(axis=0)
    return total_meth / (total_obs + PSEUDO)


def fit_M1(X_by_t, tags_by_t, T, maxIter=1000, tol=1e-8):
    """Fit M1 — constrained M2 with π[:, 1, :] = 0. Returns (theta, pi).

    Same compute as compute_M1_BIC but returns parameters. theta[:, :, 1] is
    unfit (stays at the doM_theta default of 0.5) and unused for the LL.
    """
    J = X_by_t[next(iter(X_by_t))].shape[1]
    np.random.seed()
    theta = np.random.uniform(0.01, 0.99, size=(J, 2, 2))
    pi = np.zeros((2, 2, T))
    pi[:, 0, :] = 0.5

    previous_ll = -np.inf
    for _ in range(maxIter):
        gamma_by_t = {
            t: doE(X_by_t[t], theta, pi[:, :, t - 1], tags_by_t[t], t)
            for t in range(1, T + 1)
        }
        theta = doM_theta(X_by_t, gamma_by_t)
        ll = _observed_loglikelihood(X_by_t, tags_by_t, theta, pi, T)
        if abs(ll - previous_ll) < tol:
            break
        previous_ll = ll
    return theta, pi


def fit_M2(X_by_t, tags_by_t, T):
    """Fit M2 — production EM(...). Returns (theta, pi)."""
    theta, pi, _, _ = EM(X_by_t, tags_by_t, T)
    return theta, pi


def _align_M2_components(theta_M2, pi_M2, truth):
    """M2's k=0/k=1 labels are arbitrary — swap if the alternative ordering is closer to truth.

    Compares both label assignments (identity vs swap) by total |θ - truth| over
    all (h, k) and picks the smaller. Without this, a region where EM happens
    to land with k=0 ↔ k=1 swapped would report large recovery error from
    label permutation, not real misfit.
    """
    truth_arr = np.array([
        [truth["true_theta_h0_k0"], truth["true_theta_h0_k1"]],
        [truth["true_theta_h1_k0"], truth["true_theta_h1_k1"]],
    ])  # shape (2, 2): [h, k]

    err_identity = np.mean(np.abs(theta_M2.mean(axis=0) - truth_arr))
    swapped_theta = theta_M2[:, :, ::-1]
    err_swap = np.mean(np.abs(swapped_theta.mean(axis=0) - truth_arr))
    if err_swap < err_identity:
        return swapped_theta, pi_M2[:, ::-1, :]
    return theta_M2, pi_M2


def recover_params(output_dir, manifest, T):
    """Per-region: load data, fit each model, return per-region recovery metrics."""
    files_per_t = {
        t: sorted(glob.glob(os.path.join(output_dir, f"t{t}", "methylationfraction_*_.tsv")))
        for t in range(1, T + 1)
    }
    n_regions = len(files_per_t[1])

    rows = []
    for ridx in range(n_regions):
        paths = {t: files_per_t[t][ridx] for t in range(1, T + 1)}
        X_by_t, tags_by_t, _ = load_region(paths)

        truth = manifest[manifest["region_idx"] == ridx].iloc[0]

        theta_M0 = fit_M0(X_by_t)
        theta_M1, _ = fit_M1(X_by_t, tags_by_t, T)
        theta_M2_raw, pi_M2_raw = fit_M2(X_by_t, tags_by_t, T)
        theta_M2, pi_M2 = _align_M2_components(theta_M2_raw, pi_M2_raw, truth)

        # M0's "best" target: the allele-and-component-average θ. For NULL it
        # equals θ_h0_k0; for STATIC it's the per-allele mean at k=0.
        m0_target = (truth["true_theta_h0_k0"] + truth["true_theta_h1_k0"]) / 2

        rows.append({
            "region_idx": ridx,
            "true_class": truth["true_class"],
            "M0_theta_mae": float(np.mean(np.abs(theta_M0 - m0_target))),
            "M1_theta_h0_mae": float(np.mean(np.abs(theta_M1[:, 0, 0] - truth["true_theta_h0_k0"]))),
            "M1_theta_h1_mae": float(np.mean(np.abs(theta_M1[:, 1, 0] - truth["true_theta_h1_k0"]))),
            "M2_theta_h0_k0_mae": float(np.mean(np.abs(theta_M2[:, 0, 0] - truth["true_theta_h0_k0"]))),
            "M2_theta_h1_k0_mae": float(np.mean(np.abs(theta_M2[:, 1, 0] - truth["true_theta_h1_k0"]))),
            "M2_theta_h0_k1_mae": float(np.mean(np.abs(theta_M2[:, 0, 1] - truth["true_theta_h0_k1"]))),
            "M2_theta_h1_k1_mae": float(np.mean(np.abs(theta_M2[:, 1, 1] - truth["true_theta_h1_k1"]))),
            "M2_pi_h0_alt_t2_abs_err": abs(float(pi_M2[0, 1, T - 1] - truth["pi_h0_alt_t2"])),
            "M2_pi_h1_alt_t2_abs_err": abs(float(pi_M2[1, 1, T - 1] - truth["pi_h1_alt_t2"])),
            "M2_pi_max_abs_err": max(
                abs(float(pi_M2[0, 1, T - 1] - truth["pi_h0_alt_t2"])),
                abs(float(pi_M2[1, 1, T - 1] - truth["pi_h1_alt_t2"])),
            ),
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-per-class", type=int, default=8)
    ap.add_argument("--reads-per-t", type=int, default=50)
    ap.add_argument("--J", type=int, default=20)
    ap.add_argument("--T", type=int, default=2)
    ap.add_argument("--pi-final", type=float, default=0.6)
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    with tempfile.TemporaryDirectory(prefix="param_recovery_") as out:
        manifest = generate_synthetic_dataset(
            n_regions=args.n_per_class * 4,
            reads_per_t=args.reads_per_t,
            J=args.J,
            T=args.T,
            class_mix={NULL: 1, STATIC_ASM: 1, TIME_EMERGENT_ASM: 1, ASYMMETRIC_TIME_EMERGENT_ASM: 1},
            pi_final=args.pi_final,
            output_dir=out,
            seed=args.seed,
        )
        recovery = recover_params(out, manifest, args.T)

    print(
        f"\n=== Parameter recovery "
        f"(n_per_class={args.n_per_class}, reads_per_t={args.reads_per_t}, "
        f"J={args.J}, T={args.T}, pi_final={args.pi_final}, seed={args.seed}) ===\n"
    )
    summary = (
        recovery.drop(columns=["region_idx"])
        .groupby("true_class")
        .mean(numeric_only=True)
        .round(3)
    )
    print(summary.to_string())
    print(
        "\nLegend: *_mae = mean |θ_fit - θ_true| across CpGs (lower = better).  "
        "*_abs_err = |π_fit - π_true| at t=T per allele.  "
        "*_max_abs_err = worst of the two alleles per region (use this for ASYMMETRIC)."
    )
    print(
        "\nClass truth structure (constant over CpGs within a region):\n"
        "  NULL: θ_h0_k0 = θ_h1_k0 = θ_h0_k1 = θ_h1_k1; π_alt = 0 always\n"
        "  STATIC_ASM: θ_h0_k0 ≠ θ_h1_k0 (per-allele difference); π_alt = 0 always\n"
        "  TIME_EMERGENT_ASM (sym): θ_h_k differ; π_h0_alt and π_h1_alt both ramp to pi_final\n"
        "  ASYMMETRIC_TIME_EMERGENT_ASM: only one allele's π ramps (the other stays at 0)"
    )


if __name__ == "__main__":
    main()
