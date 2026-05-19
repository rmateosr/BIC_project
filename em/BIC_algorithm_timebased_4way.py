# ABOUTME: 4-way BIC orchestrator — Shared / ASM / ASM-with-time / Shared-with-time.
# ABOUTME: Imports both EM modules, runs all four BICs per window, emits unified TSV.

import argparse
import glob
import os
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd

import EMfunctions_timebased_diploid_aware_model as em_dip
import EMfunctions_timebased_shared_time_emergent_model as em_shared_time


DEPTH = em_dip.DEPTH

# Verdict labels for BIC_4way_winner, in argmin order over
# [BICsinglecomp, BICmiddlecomp, BICpaircomp, BIC_shared_time].
WINNER_LABELS = ("shared", "asm", "asm_time", "shared_time")


def main():
    warnings.filterwarnings("ignore", category=FutureWarning, module="numpy")

    parser = argparse.ArgumentParser(
        description=(
            "4-way BIC for time-based methylation: Shared, ASM, "
            "ASM-with-time-emergence, Shared-with-time-emergence."
        )
    )
    parser.add_argument(
        'input_paths', nargs='+',
        help='Input directories, one per time period (in chronological order)',
    )
    parser.add_argument('--output', '-o', type=str, required=True, help='Output TSV path')
    parser.add_argument('--windowsize', type=int, default=10, help='CpG sites per window')
    parser.add_argument(
        '--workers', type=int, default=10,
        help='Parallel worker processes (default 10 matches BIC ASM production).',
    )
    parser.add_argument(
        '--glob-pattern', type=str, default='methylationfraction_*_.tsv',
        help='Glob pattern for per-region TSVs.',
    )
    args = parser.parse_args()

    T = len(args.input_paths)

    files_by_t = {}
    for t, dir_path in enumerate(args.input_paths, 1):
        files_by_t[t] = sorted(glob.glob(os.path.join(dir_path, args.glob_pattern)))
        print(f"Time period {t}: {len(files_by_t[t])} region files in {dir_path}")

    n_regions = len(files_by_t[1])
    for t in range(2, T + 1):
        if len(files_by_t[t]) != n_regions:
            raise ValueError(
                f"Time period {t} has {len(files_by_t[t])} files, but time "
                f"period 1 has {n_regions}. Region counts must match."
            )

    tasks = []
    for region_idx in range(n_regions):
        paths_for_region = {t: files_by_t[t][region_idx] for t in range(1, T + 1)}
        tasks.append((paths_for_region, args.windowsize))

    tmp_output = args.output + ".tmp"
    header_written = False
    with Pool(processes=args.workers) as pool, open(tmp_output, 'w') as f:
        for df in pool.imap_unordered(_process_region, tasks):
            df.to_csv(f, sep='\t', index=False, header=not header_written)
            header_written = True
            f.flush()
    os.replace(tmp_output, args.output)
    print(f"Output written to {args.output}")


def _process_region(task):
    paths_by_t, windowsize = task
    return EMBIC_bin_path_timebased_4way(paths_by_t, windowsize=windowsize)


def EMBIC_bin_path_timebased_4way(paths_by_t, windowsize=10):
    """4-way BIC over sliding CpG windows for one region.

    Reuses the diploid module's data-loading layer (Reshape2Matrix and the
    glob/load logic). The per-window BIC computation is local to this driver
    so that we can run the diploid EM and the no-h EM side by side.
    """
    T = len(paths_by_t)

    data_by_t = {}
    for t, path in paths_by_t.items():
        data_by_t[t] = pd.read_csv(path, sep='\t')

    all_coords = set()
    for t, df in data_by_t.items():
        coords_str = ','.join(df['startcoord'].astype(str))
        all_coords.update(int(float(c)) for c in coords_str.split(','))
    coords_array = np.array(sorted(all_coords))

    X_full = {}
    tags_full = {}
    for t, df in data_by_t.items():
        X_full[t] = em_dip.Reshape2Matrix(df, coords_array)
        tags_full[t] = np.array(df['haplotype'])
        valid = ~np.isnan(X_full[t]).all(axis=1)
        X_full[t] = X_full[t][valid]
        tags_full[t] = tags_full[t][valid]

    J_full = len(coords_array)

    if J_full <= windowsize:
        return _process_single_window(X_full, tags_full, T, coords_array, windowsize)

    n_windows = J_full - windowsize + 1
    BICsinglecomp = np.full(n_windows, np.nan)
    BICmiddlecomp = np.full(n_windows, np.nan)
    BICpaircomp = np.full(n_windows, np.nan)
    BIC_shared_time = np.full(n_windows, np.nan)
    em_iters = np.full(n_windows, 0, dtype=int)
    em_iters_shared = np.full(n_windows, 0, dtype=int)
    pi_alt_t2 = np.full(n_windows, np.nan)
    pi_h0_alt_t2 = np.full(n_windows, np.nan)
    pi_h1_alt_t2 = np.full(n_windows, np.nan)
    pi_alt_t2_shared = np.full(n_windows, np.nan)
    theta_diff = np.full(n_windows, np.nan)
    total_reads_per_window = np.zeros(n_windows, dtype=int)

    for cont in range(n_windows):
        X_window = {}
        tags_window = {}
        total_reads = 0
        for t in range(1, T + 1):
            X_w = X_full[t][:, cont:cont + windowsize]
            tags_w = tags_full[t]
            valid = ~np.isnan(X_w).all(axis=1)
            X_window[t] = X_w[valid]
            tags_window[t] = tags_w[valid]
            total_reads += X_window[t].shape[0]

        total_reads_per_window[cont] = total_reads

        if total_reads < DEPTH or any(X_window[t].shape[0] == 0 for t in range(1, T + 1)):
            # Sentinel: argmin([0, 2, 1, 3]) = 0 -> "shared". Matches the
            # diploid driver's degenerate-window encoding for slots 0-2 and
            # places shared_time at a value strictly greater than all others.
            BICsinglecomp[cont] = 0
            BICmiddlecomp[cont] = 2
            BICpaircomp[cont] = 1
            BIC_shared_time[cont] = 3
            continue

        result = _compute_window_BIC_4way(X_window, tags_window, T)
        (BICsinglecomp[cont], BICmiddlecomp[cont], BICpaircomp[cont],
         BIC_shared_time[cont], em_iters[cont], em_iters_shared[cont],
         pi_alt_t2[cont], pi_h0_alt_t2[cont], pi_h1_alt_t2[cont],
         pi_alt_t2_shared[cont], theta_diff[cont]) = result

    bic_stack_3 = np.stack([BICsinglecomp, BICmiddlecomp, BICpaircomp], axis=1)
    BIC_3way_winner = np.argmin(bic_stack_3, axis=1)

    bic_stack_4 = np.stack(
        [BICsinglecomp, BICmiddlecomp, BICpaircomp, BIC_shared_time], axis=1,
    )
    winner_idx = np.argmin(bic_stack_4, axis=1)
    BIC_4way_winner = np.array([WINNER_LABELS[i] for i in winner_idx])

    return pd.DataFrame({
        "coord_evaluated": coords_array[:n_windows],
        "BICsinglecomp": BICsinglecomp,
        "BICmiddlecomp": BICmiddlecomp,
        "BICpaircomp": BICpaircomp,
        "BIC_shared_time": BIC_shared_time,
        "BICresult": (BICpaircomp < BICsinglecomp).astype(int),
        "BIC_3way_winner": BIC_3way_winner,
        "BIC_4way_winner": BIC_4way_winner,
        "windows_size": windowsize,
        "T": T,
        "total_reads": total_reads_per_window,
        "em_iterations": em_iters,
        "em_iterations_shared_time": em_iters_shared,
        "pi_altered_t2": pi_alt_t2,
        "pi_h0_alt_t2": pi_h0_alt_t2,
        "pi_h1_alt_t2": pi_h1_alt_t2,
        "pi_altered_t2_shared": pi_alt_t2_shared,
        "mean_theta_diff": theta_diff,
    })


def _process_single_window(X_full, tags_full, T, coords_array, windowsize):
    """Single-window path (region <= windowsize CpGs)."""
    total_reads = sum(X_full[t].shape[0] for t in X_full)

    if total_reads < DEPTH or any(X_full[t].shape[0] == 0 for t in X_full):
        BICsinglecomp = 0
        BICmiddlecomp = 2
        BICpaircomp = 1
        BIC_shared_time_val = 3
        n_iters = 0
        n_iters_shared = 0
        pi_altered_t2 = np.nan
        pi_h0_alt_t2 = np.nan
        pi_h1_alt_t2 = np.nan
        pi_altered_t2_shared = np.nan
        mean_theta_diff = np.nan
    else:
        (BICsinglecomp, BICmiddlecomp, BICpaircomp, BIC_shared_time_val,
         n_iters, n_iters_shared, pi_altered_t2, pi_h0_alt_t2, pi_h1_alt_t2,
         pi_altered_t2_shared, mean_theta_diff) = _compute_window_BIC_4way(
            X_full, tags_full, T,
        )

    bic_stack_3 = np.array([BICsinglecomp, BICmiddlecomp, BICpaircomp])
    BIC_3way_winner = int(np.argmin(bic_stack_3))
    bic_stack_4 = np.array(
        [BICsinglecomp, BICmiddlecomp, BICpaircomp, BIC_shared_time_val],
    )
    BIC_4way_winner = WINNER_LABELS[int(np.argmin(bic_stack_4))]

    return pd.DataFrame({
        "coord_evaluated": coords_array[:1],
        "BICsinglecomp": [BICsinglecomp],
        "BICmiddlecomp": [BICmiddlecomp],
        "BICpaircomp": [BICpaircomp],
        "BIC_shared_time": [BIC_shared_time_val],
        "BICresult": [int(BICpaircomp < BICsinglecomp)],
        "BIC_3way_winner": [BIC_3way_winner],
        "BIC_4way_winner": [BIC_4way_winner],
        "windows_size": [windowsize],
        "T": [T],
        "total_reads": [total_reads],
        "em_iterations": [n_iters],
        "em_iterations_shared_time": [n_iters_shared],
        "pi_altered_t2": [pi_altered_t2],
        "pi_h0_alt_t2": [pi_h0_alt_t2],
        "pi_h1_alt_t2": [pi_h1_alt_t2],
        "pi_altered_t2_shared": [pi_altered_t2_shared],
        "mean_theta_diff": [mean_theta_diff],
    })


def _compute_window_BIC_4way(X_window, tags_window, T):
    """Run both EMs and compute all four BICs for one window.

    The diploid null and the no-h null are the same pooled-theta MLE — we
    use the diploid value as the single `BICsinglecomp`. The shared-time
    BIC reports the no-h M3.
    """
    theta_dip, pi_dip, _, n_iters_dip = em_dip.EM(X_window, tags_window, T)

    bic_shared = em_dip.compute_null_BIC(X_window, T)
    bic_asm, _ = em_dip.compute_M1_BIC(X_window, tags_window, T)

    ll_alt_dip = em_dip.compute_alt_loglikelihood(
        X_window, tags_window, theta_dip, pi_dip, T,
    )
    J = X_window[next(iter(X_window))].shape[1]
    n_total = sum(X_window[t].shape[0] for t in X_window)
    p_asm_time = 4 * J + 2 * T - 2
    bic_asm_time = p_asm_time * np.log(n_total) - 2 * ll_alt_dip

    theta_shared, pi_shared, _, n_iters_shared = em_shared_time.EM(X_window, T)
    ll_shared = em_shared_time.compute_alt_loglikelihood(
        X_window, theta_shared, pi_shared, T,
    )
    p_shared_time = 2 * J + (T - 1)
    bic_shared_time_val = p_shared_time * np.log(n_total) - 2 * ll_shared

    # Diagnostic columns. pi_dip shape: (2, 2, T) -> pi[h, k, t-1].
    if T >= 2:
        pi_h0_alt_t2 = float(pi_dip[0, 1, 1])
        pi_h1_alt_t2 = float(pi_dip[1, 1, 1])
        pi_altered_t2 = pi_h0_alt_t2 + pi_h1_alt_t2
        pi_altered_t2_shared = float(pi_shared[1, T - 1])
    else:
        pi_h0_alt_t2 = 0.0
        pi_h1_alt_t2 = 0.0
        pi_altered_t2 = 0.0
        pi_altered_t2_shared = 0.0
    mean_theta_diff = float(np.mean(np.abs(theta_dip[:, :, 1] - theta_dip[:, :, 0])))

    return (
        float(bic_shared), float(bic_asm), float(bic_asm_time),
        float(bic_shared_time_val),
        n_iters_dip, n_iters_shared,
        pi_altered_t2, pi_h0_alt_t2, pi_h1_alt_t2,
        pi_altered_t2_shared, mean_theta_diff,
    )


if __name__ == "__main__":
    main()
