# ABOUTME: Sweep parameter recovery + 3-way BIC classification across pi_final × reads_per_t.
# ABOUTME: Writes long-form CSVs and PNG plots to em/sweep_artifacts/ for the recovery report.
"""Sweep harness for the parameter-recovery + 3-way BIC experiments.

For every (pi_final, reads_per_t, seed) cell in the configured grid:
  1. Generate a fresh synthetic dataset (8 regions per class, 4 classes).
  2. Fit M0 / M1 / M2 to each region; compute parameter recovery error
     against the manifest truth (delegates to parameter_recovery.recover_params).
  3. Run the production EMBIC_bin_path_timebased pipeline per region;
     tabulate BIC_3way_winner % per class.

Outputs (default em/sweep_artifacts/):
  - sweep_recovery.csv  long-form per-cell, per-class θ + π recovery error
  - sweep_bic.csv       long-form per-cell, per-class BIC 3-way winner %
  - figures/*.png       summary plots (one per panel)

Runtime: ~60–90 s per cell on this cluster, depending on reads_per_t.
"""

import argparse
import glob
import os
import sys
import tempfile
import time
from itertools import product

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from EMfunctions_timebased_diploid_aware_model import EMBIC_bin_path_timebased
from parameter_recovery import recover_params
from synthetic.generate_synthetic import (
    ASYMMETRIC_TIME_EMERGENT_ASM,
    NULL,
    STATIC_ASM,
    TIME_EMERGENT_ASM,
    generate_synthetic_dataset,
)


CLASS_ORDER = [NULL, STATIC_ASM, TIME_EMERGENT_ASM, ASYMMETRIC_TIME_EMERGENT_ASM]
CLASS_SHORT = {
    NULL: "NULL",
    STATIC_ASM: "STATIC",
    TIME_EMERGENT_ASM: "TE_sym",
    ASYMMETRIC_TIME_EMERGENT_ASM: "ASYM",
}


def run_one_cell(pi_final, reads_per_t, seed, n_per_class=8, J=20, T=2, windowsize=10):
    """Generate data, run BIC pipeline + recovery; return (recovery_df, bic_df)."""
    with tempfile.TemporaryDirectory(
        prefix=f"sweep_pi{pi_final}_n{reads_per_t}_s{seed}_"
    ) as out:
        manifest = generate_synthetic_dataset(
            n_regions=n_per_class * 4,
            reads_per_t=reads_per_t,
            J=J,
            T=T,
            class_mix={
                NULL: 1, STATIC_ASM: 1, TIME_EMERGENT_ASM: 1,
                ASYMMETRIC_TIME_EMERGENT_ASM: 1,
            },
            pi_final=pi_final,
            output_dir=out,
            seed=seed,
        )
        recovery = recover_params(out, manifest, T)

        files = {
            t: sorted(glob.glob(os.path.join(out, f"t{t}", "methylationfraction_*_.tsv")))
            for t in range(1, T + 1)
        }
        bic_results = []
        for ridx in range(len(files[1])):
            paths = {t: files[t][ridx] for t in range(1, T + 1)}
            r = EMBIC_bin_path_timebased(paths, windowsize=windowsize)
            r["region_idx"] = ridx
            bic_results.append(r)
        bic_results = pd.concat(bic_results, ignore_index=True)
        truth_per_region = manifest[["region_idx", "true_class"]].drop_duplicates()
        bic_results = bic_results.merge(truth_per_region, on="region_idx", how="left")

    rec_summary = (
        recovery.drop(columns=["region_idx"], errors="ignore")
        .groupby("true_class")
        .mean(numeric_only=True)
        .round(4)
        .reset_index()
    )
    rec_summary["pi_final"] = pi_final
    rec_summary["reads_per_t"] = reads_per_t
    rec_summary["seed"] = seed

    bic_summary = (
        bic_results.groupby("true_class")["BIC_3way_winner"]
        .value_counts(normalize=True)
        .unstack(fill_value=0.0)
        .reset_index()
    )
    for col_int, col_name in [(0, "pct_M0"), (1, "pct_M1"), (2, "pct_M2")]:
        if col_int in bic_summary.columns:
            bic_summary = bic_summary.rename(columns={col_int: col_name})
        if col_name not in bic_summary.columns:
            bic_summary[col_name] = 0.0
    bic_summary = bic_summary[["true_class", "pct_M0", "pct_M1", "pct_M2"]]
    bic_summary["pi_final"] = pi_final
    bic_summary["reads_per_t"] = reads_per_t
    bic_summary["seed"] = seed

    return rec_summary, bic_summary


def make_plots(rec_long, bic_long, fig_dir):
    """Generate the report's headline plots."""
    os.makedirs(fig_dir, exist_ok=True)

    # Plot 1: M2 max-abs π recovery error vs reads_per_t, by truth class,
    # at the highest pi_final. Lines are TRUTH CLASSES, not models — only
    # M2 is being fit and evaluated here.
    pi_max = rec_long["pi_final"].max()
    rec_top = rec_long[rec_long["pi_final"] == pi_max].copy()
    rec_top["class_short"] = rec_top["true_class"].map(CLASS_SHORT)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    class_label = {
        NULL: "NULL truth (no ASM)",
        STATIC_ASM: "STATIC truth (imprinted-like)",
        TIME_EMERGENT_ASM: "TE_sym truth (both alleles ramp)",
        ASYMMETRIC_TIME_EMERGENT_ASM: "ASYM truth (one allele ramps)",
    }
    for cls in CLASS_ORDER:
        sub = rec_top[rec_top["true_class"] == cls].sort_values("reads_per_t")
        ax.plot(sub["reads_per_t"], sub["M2_pi_max_abs_err"],
                marker="o", label=class_label[cls])
    ax.set_xscale("log")
    ax.set_xlabel("reads_per_t (log scale)")
    ax.set_ylabel("M2 π max |fit − truth| at t=T   (lower = better)")
    ax.set_title(
        f"How well does M2 alone recover π?  (only M2 evaluated; pi_final={pi_max})\n"
        "Lines = synthetic truth classes the data was generated under, NOT models."
    )
    ax.axhline(0.10, color="gray", linestyle=":", linewidth=0.8,
               label="0.10 reference")
    ax.legend(loc="upper right", fontsize=8.5, title="Truth class")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "fig1_pi_recovery_vs_reads.png"), dpi=120)
    plt.close(fig)

    # Plot 2: BIC M2-winner rate on ASYMMETRIC, heatmap pi_final × reads_per_t
    asym = bic_long[bic_long["true_class"] == ASYMMETRIC_TIME_EMERGENT_ASM]
    pivot_M2 = asym.pivot_table(
        index="pi_final", columns="reads_per_t", values="pct_M2", aggfunc="mean"
    )
    fig, ax = plt.subplots(figsize=(6.5, 4))
    im = ax.imshow(pivot_M2.values, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(pivot_M2.columns)))
    ax.set_xticklabels(pivot_M2.columns)
    ax.set_yticks(range(len(pivot_M2.index)))
    ax.set_yticklabels([f"{v:.2f}" for v in pivot_M2.index])
    ax.set_xlabel("reads_per_t")
    ax.set_ylabel("pi_final (signal strength)")
    ax.set_title("ASYMMETRIC truth — fraction of windows where 3-way BIC picks M2\n"
                 "(higher = M2 fires more often on the model that matches truth)")
    for i in range(pivot_M2.shape[0]):
        for j in range(pivot_M2.shape[1]):
            ax.text(j, i, f"{pivot_M2.values[i, j]:.2f}",
                    ha="center", va="center",
                    color="white" if pivot_M2.values[i, j] < 0.5 else "black",
                    fontsize=10)
    fig.colorbar(im, ax=ax, label="fraction M2-winner")
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "fig2_bic_M2_rate_heatmap.png"), dpi=120)
    plt.close(fig)

    # Plot 3: per-class BIC 3-way winner % stacked bars at pi_final=max, varying reads
    fig, axes = plt.subplots(1, len(CLASS_ORDER), figsize=(13, 4), sharey=True)
    for ax, cls in zip(axes, CLASS_ORDER):
        sub = (
            bic_long[(bic_long["true_class"] == cls) & (bic_long["pi_final"] == pi_max)]
            .groupby("reads_per_t")[["pct_M0", "pct_M1", "pct_M2"]]
            .mean()
            .sort_index()
        )
        bottom = np.zeros(len(sub))
        for col, color in zip(["pct_M0", "pct_M1", "pct_M2"], ["#888", "#1f77b4", "#d62728"]):
            ax.bar(range(len(sub)), sub[col].values, bottom=bottom,
                   label=col.replace("pct_", ""), color=color)
            bottom += sub[col].values
        ax.set_xticks(range(len(sub)))
        ax.set_xticklabels(sub.index, rotation=0)
        ax.set_title(f"{CLASS_SHORT[cls]}\ntruth")
        ax.set_xlabel("reads_per_t")
    axes[0].set_ylabel(f"BIC 3-way winner share (pi_final={pi_max})")
    axes[-1].legend(loc="lower right", fontsize=9)
    fig.suptitle("3-way BIC winner share by class — does the right model win?")
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "fig3_bic_winner_stacks.png"), dpi=120)
    plt.close(fig)

    # Plot 4: θ recovery — M0/M1/M2 (k=0) vs reads, on STATIC truth.
    # IMPORTANT: each line is a different model fit, judged against the
    # truth target that model can structurally recover. M0 vs allele
    # average; M1 vs per-allele; M2 vs per-(allele, component) at k=0.
    # Lines being similar means each model successfully fits its own
    # target — NOT that they are equally appropriate models.
    fig, ax = plt.subplots(figsize=(7.5, 5))
    static = rec_long[rec_long["true_class"] == STATIC_ASM]
    static = static[static["pi_final"] == pi_max].sort_values("reads_per_t")
    ax.plot(static["reads_per_t"], static["M0_theta_mae"],
            marker="o", label="M0  vs allele-avg truth")
    ax.plot(static["reads_per_t"],
            (static["M1_theta_h0_mae"] + static["M1_theta_h1_mae"]) / 2,
            marker="s", label="M1  vs per-allele truth")
    ax.plot(static["reads_per_t"],
            (static["M2_theta_h0_k0_mae"] + static["M2_theta_h1_k0_mae"]) / 2,
            marker="^", label="M2 k=0  vs per-allele truth")
    ax.set_xscale("log")
    ax.set_xlabel("reads_per_t (log scale)")
    ax.set_ylabel("θ MAE — each model vs its own structural target")
    ax.set_title(
        f"θ recovery on STATIC truth — every model hits its own target\n"
        f"(pi_final={pi_max}; lower = better; lines NOT directly comparable — see caption)"
    )
    ax.legend(title="model  vs  truth target it can recover", fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "fig4_theta_recovery_static.png"), dpi=120)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pi-final-vals", nargs="+", type=float, default=[0.3, 0.6, 0.9])
    ap.add_argument("--reads-per-t-vals", nargs="+", type=int, default=[25, 50, 100, 200])
    ap.add_argument("--seeds", nargs="+", type=int, default=[2026])
    ap.add_argument("--n-per-class", type=int, default=8)
    ap.add_argument("--out-dir", type=str, default="sweep_artifacts")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    fig_dir = os.path.join(args.out_dir, "figures")

    cells = list(product(args.pi_final_vals, args.reads_per_t_vals, args.seeds))
    print(f"Running {len(cells)} cells "
          f"({len(args.pi_final_vals)} pi_final × {len(args.reads_per_t_vals)} reads_per_t "
          f"× {len(args.seeds)} seed(s))...")

    rec_rows = []
    bic_rows = []
    for i, (pi_f, reads, seed) in enumerate(cells, 1):
        t0 = time.time()
        rec, bic = run_one_cell(pi_f, reads, seed, n_per_class=args.n_per_class)
        rec_rows.append(rec)
        bic_rows.append(bic)
        elapsed = time.time() - t0
        print(f"  [{i}/{len(cells)}] pi_final={pi_f} reads_per_t={reads} "
              f"seed={seed} took {elapsed:.1f}s", flush=True)

    rec_long = pd.concat(rec_rows, ignore_index=True)
    bic_long = pd.concat(bic_rows, ignore_index=True)

    rec_path = os.path.join(args.out_dir, "sweep_recovery.csv")
    bic_path = os.path.join(args.out_dir, "sweep_bic.csv")
    rec_long.to_csv(rec_path, index=False)
    bic_long.to_csv(bic_path, index=False)
    print(f"\nWrote {rec_path}, {bic_path}.")

    print("Generating plots...")
    make_plots(rec_long, bic_long, fig_dir)
    print(f"Wrote plots to {fig_dir}/")


if __name__ == "__main__":
    main()
