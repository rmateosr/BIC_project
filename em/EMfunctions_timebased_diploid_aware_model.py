# ABOUTME: EM algorithm for time-based ASM with diploid per-allele marginal constraint.
# ABOUTME: pi[h,:,t] pinned at 0.5 for each allele; see latex/EM_algorithm_derivation_timebased_ASM_optionA.tex.

import hashlib

import numpy as np
import pandas as pd


def _seed_from_X(X_by_t):
    """Deterministic 32-bit seed derived from the per-timepoint read matrices.

    Same input matrices -> same seed -> same random theta init -> same EM
    trajectory. Used by EM() and compute_M1_BIC() to make per-window BIC
    values reproducible across runs.
    """
    h = hashlib.md5()
    for t in sorted(X_by_t):
        h.update(np.ascontiguousarray(X_by_t[t]).tobytes())
    return int(h.hexdigest()[:8], 16)

np.seterr(all='raise')

PSEUDO = 1e-10
DEPTH = 1


def Reshape2Matrix(binfullset, coords_array):
    """Convert read-level data to a binary matrix (n_reads x J).

    Each row is a read; each column is a CpG site in coords_array.
    Entries are 0 or 1 (methylation status), NaN if the read doesn't cover that site.

    Parameters
    ----------
    binfullset : DataFrame
        Must have columns 'startcoord' (comma-separated positions) and 'status' (comma-separated 0/1).
    coords_array : ndarray of int
        Sorted CpG positions defining the columns.

    Returns
    -------
    matrix : ndarray, shape (n_reads, len(coords_array))
    """
    nrows = binfullset.shape[0]
    ncols = len(coords_array)
    matrix = np.full((nrows, ncols), np.nan)

    for i in range(nrows):
        row = binfullset.iloc[i]
        read_coords = np.array([int(float(c)) for c in str(row['startcoord']).split(',')])
        read_values = np.array([float(v) for v in str(row['status']).split(',')])

        # Remove duplicate coordinates (keep first)
        _, unique_idx = np.unique(read_coords, return_index=True)
        read_coords = read_coords[unique_idx]
        read_values = read_values[unique_idx]

        # Map read coordinates to matrix columns
        coord_in_set = np.isin(coords_array, read_coords)
        read_in_set = np.isin(read_coords, coords_array)

        matrix[i, np.where(coord_in_set)[0]] = read_values[read_in_set]

    return matrix


def EMBIC_bin_path_timebased(paths_by_t, windowsize=10):
    """Run the time-based ASM pipeline on input files for multiple time periods.

    Parameters
    ----------
    paths_by_t : dict {t: str}
        File paths for each time period (1-indexed keys).
    windowsize : int
        Number of CpG sites per sliding window.

    Returns
    -------
    BIC_summary : DataFrame
    """
    T = len(paths_by_t)

    # Load all time periods
    data_by_t = {}
    for t, path in paths_by_t.items():
        data_by_t[t] = pd.read_csv(path, sep='\t')

    # Find union of all CpG coordinates across time periods
    all_coords = set()
    for t, df in data_by_t.items():
        coords_str = ','.join(df['startcoord'].astype(str))
        all_coords.update(int(float(c)) for c in coords_str.split(','))
    coords_array = np.array(sorted(all_coords))

    # Build matrices per time period
    X_full = {}
    tags_full = {}
    for t, df in data_by_t.items():
        X_full[t] = Reshape2Matrix(df, coords_array)
        tags_full[t] = np.array(df['haplotype'])

        # Remove all-NaN rows (reads that don't overlap this region)
        valid = ~np.isnan(X_full[t]).all(axis=1)
        X_full[t] = X_full[t][valid]
        tags_full[t] = tags_full[t][valid]

    J_full = len(coords_array)

    # If the region is smaller than or equal to the window size, process as one window
    if J_full <= windowsize:
        return _process_single_window(X_full, tags_full, T, coords_array, windowsize)

    # Slide window across CpG sites
    n_windows = J_full - windowsize + 1
    BICsinglecomp = np.full(n_windows, np.nan)
    BICmiddlecomp = np.full(n_windows, np.nan)
    BICpaircomp = np.full(n_windows, np.nan)
    em_iters = np.full(n_windows, 0, dtype=int)
    pi_alt_t2 = np.full(n_windows, np.nan)
    pi_h0_alt_t2 = np.full(n_windows, np.nan)
    pi_h1_alt_t2 = np.full(n_windows, np.nan)
    theta_diff = np.full(n_windows, np.nan)
    total_reads_per_window = np.zeros(n_windows, dtype=int)

    for cont in range(n_windows):
        # Slice each time period to the window
        X_window = {}
        tags_window = {}
        total_reads = 0

        for t in range(1, T + 1):
            X_w = X_full[t][:, cont:cont + windowsize]
            tags_w = tags_full[t]

            # Remove reads that are all-NaN within this window
            valid = ~np.isnan(X_w).all(axis=1)
            X_window[t] = X_w[valid]
            tags_window[t] = tags_w[valid]
            total_reads += X_window[t].shape[0]

        total_reads_per_window[cont] = total_reads

        if total_reads < DEPTH:
            # Sentinels: argmin of (0, 2, 1) = 0 → BIC_3way_winner = M0.
            # Order also preserves the legacy 2-way `BICresult` (=0).
            BICsinglecomp[cont] = 0
            BICmiddlecomp[cont] = 2
            BICpaircomp[cont] = 1
            continue

        # Skip windows where any time point has 0 reads — EM can't estimate
        # mixing weights (pi) for a time point with no observations.
        if any(X_window[t].shape[0] == 0 for t in range(1, T + 1)):
            BICsinglecomp[cont] = 0
            BICmiddlecomp[cont] = 2
            BICpaircomp[cont] = 1
            continue

        result = _compute_window_BIC(X_window, tags_window, T)
        BICsinglecomp[cont] = result[0]
        BICmiddlecomp[cont] = result[1]
        BICpaircomp[cont] = result[2]
        em_iters[cont] = result[3]
        pi_alt_t2[cont] = result[4]
        pi_h0_alt_t2[cont] = result[5]
        pi_h1_alt_t2[cont] = result[6]
        theta_diff[cont] = result[7]

    bic_stack = np.stack([BICsinglecomp, BICmiddlecomp, BICpaircomp], axis=1)
    BIC_3way_winner = np.argmin(bic_stack, axis=1)

    BIC_summary = pd.DataFrame({
        "coord_evaluated": coords_array[:n_windows],
        "BICsinglecomp": BICsinglecomp,
        "BICmiddlecomp": BICmiddlecomp,
        "BICpaircomp": BICpaircomp,
        "BICresult": (BICpaircomp < BICsinglecomp).astype(int),
        "BIC_3way_winner": BIC_3way_winner,
        "windows_size": windowsize,
        "T": T,
        "total_reads": total_reads_per_window,
        "em_iterations": em_iters,
        "pi_altered_t2": pi_alt_t2,
        "pi_h0_alt_t2": pi_h0_alt_t2,
        "pi_h1_alt_t2": pi_h1_alt_t2,
        "mean_theta_diff": theta_diff,
    })
    return BIC_summary


def _process_single_window(X_full, tags_full, T, coords_array, windowsize):
    """Handle the case where the region fits in a single window."""
    total_reads = sum(X_full[t].shape[0] for t in X_full)

    if total_reads < DEPTH or any(X_full[t].shape[0] == 0 for t in X_full):
        # Same sentinel ordering as the sliding-window path: argmin → M0.
        BICsinglecomp = 0
        BICmiddlecomp = 2
        BICpaircomp = 1
        n_iters = 0
        pi_altered_t2 = np.nan
        pi_h0_alt_t2 = np.nan
        pi_h1_alt_t2 = np.nan
        mean_theta_diff = np.nan
    else:
        result = _compute_window_BIC(X_full, tags_full, T)
        BICsinglecomp = result[0]
        BICmiddlecomp = result[1]
        BICpaircomp = result[2]
        n_iters = result[3]
        pi_altered_t2 = result[4]
        pi_h0_alt_t2 = result[5]
        pi_h1_alt_t2 = result[6]
        mean_theta_diff = result[7]

    bic_stack = np.array([BICsinglecomp, BICmiddlecomp, BICpaircomp])
    BIC_3way_winner = int(np.argmin(bic_stack))

    BIC_summary = pd.DataFrame({
        "coord_evaluated": coords_array[:1],
        "BICsinglecomp": [BICsinglecomp],
        "BICmiddlecomp": [BICmiddlecomp],
        "BICpaircomp": [BICpaircomp],
        "BICresult": [int(BICpaircomp < BICsinglecomp)],
        "BIC_3way_winner": [BIC_3way_winner],
        "windows_size": [windowsize],
        "T": [T],
        "total_reads": [total_reads],
        "em_iterations": [n_iters],
        "pi_altered_t2": [pi_altered_t2],
        "pi_h0_alt_t2": [pi_h0_alt_t2],
        "pi_h1_alt_t2": [pi_h1_alt_t2],
        "mean_theta_diff": [mean_theta_diff],
    })
    return BIC_summary


def _compute_window_BIC(X_window, tags_window, T):
    """Run EM and compute BIC for a single window.

    Returns (bic_null, bic_M1, bic_alt, n_iters, pi_altered_t2, pi_h0_alt_t2,
    pi_h1_alt_t2, mean_theta_diff). bic_M1 is the constrained allele-specific,
    no-time-emergence model — see compute_M1_BIC.
    """
    theta, pi, gamma_by_t, n_iters = EM(X_window, tags_window, T)

    bic_null = compute_null_BIC(X_window, T)
    bic_M1, _ = compute_M1_BIC(X_window, tags_window, T)

    ll_alt = compute_alt_loglikelihood(X_window, tags_window, theta, pi, T)
    J = X_window[next(iter(X_window))].shape[1]
    n_total = sum(X_window[t].shape[0] for t in X_window)
    p_alt = 4 * J + 2 * T - 2
    bic_alt = p_alt * np.log(n_total) - 2 * ll_alt

    # Per-allele altered weights at t=2 expose asymmetric ASM events: an event
    # where only one parental copy drifts to the altered methylation program
    # has pi_h0_alt_t2 ≈ 0 and pi_h1_alt_t2 elevated (or vice versa), while
    # their sum (pi_altered_t2) cannot distinguish that from a balanced drift.
    # pi shape: (2, 2, T) — pi[h, k, t-1]
    if T >= 2:
        pi_h0_alt_t2 = float(pi[0, 1, 1])
        pi_h1_alt_t2 = float(pi[1, 1, 1])
        pi_altered_t2 = pi_h0_alt_t2 + pi_h1_alt_t2
    else:
        pi_h0_alt_t2 = 0.0
        pi_h1_alt_t2 = 0.0
        pi_altered_t2 = 0.0

    # mean_theta_diff: average |theta[:,:,k=1] - theta[:,:,k=0]| across sites
    # and alleles. Measures how different the altered methylation pattern is
    # from the normal one. High = the two states are very different.
    # theta shape: (J, 2, 2) — theta[j, h, k]
    mean_theta_diff = float(np.mean(np.abs(theta[:, :, 1] - theta[:, :, 0])))

    return (
        float(bic_null), float(bic_M1), float(bic_alt), n_iters,
        pi_altered_t2, pi_h0_alt_t2, pi_h1_alt_t2,
        mean_theta_diff,
    )


def doE(X_t, theta, pi_t, tags_t, t):
    """E-step for a single time period.

    Parameters
    ----------
    X_t : ndarray, shape (I_t, J)
        Binary methylation matrix for time period t. NaN for missing.
    theta : ndarray, shape (J, 2, 2)
        Current theta estimates. theta[j, h, k] = P(methylated at site j | allele h, model k).
    pi_t : ndarray, shape (2, 2)
        Current pi estimates for time period t. pi_t[h, k] = P(allele h, model k | time t).
    tags_t : ndarray of str, shape (I_t,)
        Haplotype tags: "H1", "H2", or "noH".
    t : int
        Time period index (1-indexed). Used for boundary condition at t=1.

    Returns
    -------
    gamma : ndarray, shape (I_t, 2, 2)
        Posterior probabilities. gamma[i, h, k] = P(z_{i,t,h,k}=1 | X, theta, pi).
    """
    I_t, J = X_t.shape
    gamma = np.zeros((I_t, 2, 2))

    # Handle NaN: split into methylated/unmethylated contribution matrices
    # Where NaN, both are 0 so the site contributes 0 to the log-likelihood
    xformeth = X_t.copy()
    xformeth[np.isnan(xformeth)] = 0

    xforunmeth = X_t.copy()
    xforunmeth[np.isnan(xforunmeth)] = 1
    xforunmeth = 1 - xforunmeth

    # Compute unnormalized log-posteriors for each (h, k) combination
    # log_gamma_prime[i, h, k] = log(pi_t[h,k]) + sum_j [x*log(theta) + (1-x)*log(1-theta)]
    log_gamma_prime = np.full((I_t, 2, 2), -np.inf)

    for h in range(2):
        for k in range(2):
            if pi_t[h, k] < PSEUDO:
                # Component has zero prior — leave at -inf
                continue

            log_theta = np.log(theta[:, h, k] + PSEUDO)        # (J,)
            log_1_minus_theta = np.log(1 - theta[:, h, k] + PSEUDO)  # (J,)

            # x*log(theta) + (1-x)*log(1-theta), with NaN handled via xformeth/xforunmeth
            ll_per_read = xformeth @ log_theta + xforunmeth @ log_1_minus_theta  # (I_t,)

            log_gamma_prime[:, h, k] = np.log(pi_t[h, k] + PSEUDO) + ll_per_read

    # Boundary condition at t=1: model k=1 (altered, 0-indexed) does not exist
    if t == 1:
        log_gamma_prime[:, :, 1] = -np.inf

    # Vectorized normalization by tag group. Replaces per-read Python loop.
    # Groups: H1 -> normalize over k for h=0; H2 -> over k for h=1; noH -> over all 4.
    h1_mask = tags_t == "H1"
    h2_mask = tags_t == "H2"
    noh_mask = ~(h1_mask | h2_mask)

    if np.any(h1_mask):
        gamma[h1_mask, 0, :] = _normalize_log_batch(
            log_gamma_prime[h1_mask, 0, :], axis=-1
        )
        # gamma[h1_mask, 1, :] remains 0 (h=1 forbidden for H1 reads)

    if np.any(h2_mask):
        gamma[h2_mask, 1, :] = _normalize_log_batch(
            log_gamma_prime[h2_mask, 1, :], axis=-1
        )
        # gamma[h2_mask, 0, :] remains 0 (h=0 forbidden for H2 reads)

    if np.any(noh_mask):
        n_noh = int(noh_mask.sum())
        flat = log_gamma_prime[noh_mask].reshape(n_noh, 4)
        gamma[noh_mask] = _normalize_log_batch(flat, axis=-1).reshape(n_noh, 2, 2)

    return gamma


def _normalize_log(log_values):
    """Convert unnormalized log-probabilities to normalized probabilities.

    Uses log-sum-exp with clipping for numerical stability.
    """
    finite_mask = np.isfinite(log_values)
    if not np.any(finite_mask):
        # All components are -inf: return uniform over the positions
        # (degenerate case, shouldn't happen in practice)
        return np.ones_like(log_values) / len(log_values)

    L_star = np.max(log_values[finite_mask])
    shifted = np.clip(log_values - L_star, -700, None)
    exp_shifted = np.exp(shifted)
    # Components that were -inf must be exactly 0, not exp(-700)
    exp_shifted[~finite_mask] = 0.0
    total = exp_shifted.sum()

    if total == 0:
        return np.ones_like(log_values) / len(log_values)

    return exp_shifted / total


def _normalize_log_batch(log_values, axis=-1):
    """Vectorized softmax over `axis` with -inf masking.

    Per-slice semantics match `_normalize_log`: slices that are all -inf
    (or whose exp-shifted sum is 0) return uniform; otherwise returns
    normalized probabilities summing to 1 along `axis`.
    """
    finite_mask = np.isfinite(log_values)
    any_finite = np.any(finite_mask, axis=axis, keepdims=True)

    # Max over finite entries; for all-(-inf) slices, force L_star to 0 so
    # the subtraction below is well-defined (those slices get uniform at the end).
    safe_for_max = np.where(finite_mask, log_values, -np.inf)
    L_star = np.max(safe_for_max, axis=axis, keepdims=True)
    L_star = np.where(any_finite, L_star, 0.0)

    shifted = np.clip(log_values - L_star, -700, None)
    exp_shifted = np.where(finite_mask, np.exp(shifted), 0.0)
    total = exp_shifted.sum(axis=axis, keepdims=True)

    n_along = log_values.shape[axis]
    uniform = np.full_like(log_values, 1.0 / n_along)
    safe_total = np.where(total > 0, total, 1.0)
    normalized = exp_shifted / safe_total
    degenerate = (total == 0) | (~any_finite)
    return np.where(degenerate, uniform, normalized)


def _logsumexp_batch(log_values, axis=-1):
    """Vectorized log-sum-exp over `axis` with -inf masking.

    Matches the per-read log-sum-exp used in `_observed_loglikelihood`:
    `L_star + log(sum(exp(shifted)) + PSEUDO)`. For slices that are
    entirely -inf the contribution is 0 (caller summed them and skipped).
    """
    finite_mask = np.isfinite(log_values)
    any_finite = np.any(finite_mask, axis=axis, keepdims=True)

    safe_for_max = np.where(finite_mask, log_values, -np.inf)
    L_star = np.max(safe_for_max, axis=axis, keepdims=True)
    L_star = np.where(any_finite, L_star, 0.0)

    shifted = np.clip(log_values - L_star, -700, None)
    exp_shifted = np.where(finite_mask, np.exp(shifted), 0.0)
    total = exp_shifted.sum(axis=axis, keepdims=True)

    result = np.squeeze(L_star + np.log(total + PSEUDO), axis=axis)
    any_finite_squeezed = np.squeeze(any_finite, axis=axis)
    return np.where(any_finite_squeezed, result, 0.0)


def doM_theta(X_by_t, gamma_by_t):
    """M-step for theta: pools data across ALL time periods.

    theta[j, h, k] = sum_{t,i} gamma[i,h,k] * x_{i,j} / sum_{t,i} gamma[i,h,k] * observed_{i,j}

    Parameters
    ----------
    X_by_t : dict {t: ndarray (I_t, J)}
        Binary methylation matrices per time period.
    gamma_by_t : dict {t: ndarray (I_t, 2, 2)}
        Posterior probabilities per time period.

    Returns
    -------
    theta : ndarray, shape (J, 2, 2)
    """
    # Determine J from any time period
    any_t = next(iter(X_by_t))
    J = X_by_t[any_t].shape[1]
    theta = np.full((J, 2, 2), 0.5)

    for h in range(2):
        for k in range(2):
            numerator = np.zeros(J)
            denominator = np.zeros(J)

            for t in X_by_t:
                X_t = X_by_t[t]
                g = gamma_by_t[t][:, h, k]  # (I_t,)

                observed = ~np.isnan(X_t)  # (I_t, J)
                X_clean = X_t.copy()
                X_clean[np.isnan(X_clean)] = 0

                # Weighted sum of methylation values (NaN sites contribute 0)
                numerator += (g[:, np.newaxis] * X_clean).sum(axis=0)
                # Weighted count of observed sites
                denominator += (g[:, np.newaxis] * observed.astype(float)).sum(axis=0)

            # Where denominator > 0, compute weighted mean; else keep default 0.5
            has_data = denominator > PSEUDO
            theta[has_data, h, k] = numerator[has_data] / denominator[has_data]

    return theta


def doM_pi(gamma_by_t, T):
    """M-step for pi under the diploid per-allele marginal constraint.

    pi[h, k, t] = 0.5 * sum_i gamma[i, h, k] / sum_{k'} sum_i gamma[i, h, k']

    Enforces sum_k pi[h, k, t] = 0.5 for each allele h at every time t, derived
    by a two-Lagrange-multiplier M-step (one per allele). See §3.4 of
    latex/EM_algorithm_derivation_timebased_ASM_optionA.tex.

    Parameters
    ----------
    gamma_by_t : dict {t: ndarray (I_t, 2, 2)}
    T : int
        Number of time periods.

    Returns
    -------
    pi : ndarray, shape (2, 2, T)
    """
    pi = np.zeros((2, 2, T))

    for t in range(1, T + 1):
        gamma_t = gamma_by_t[t]            # (I_t, 2, 2)
        numer = gamma_t.sum(axis=0)         # (2, 2)  sum over reads
        S_h = numer.sum(axis=1)             # (2,)    per-allele posterior mass

        # Guard: if all reads at time t put zero posterior mass on allele h,
        # the update is undefined. Fall back to (pi[h,0,t], pi[h,1,t]) = (0.5, 0).
        S_h_safe = np.where(S_h > 0, S_h, 1.0)
        pi[:, :, t - 1] = 0.5 * numer / S_h_safe[:, None]
        for h in range(2):
            if S_h[h] == 0:
                pi[h, 0, t - 1] = 0.5
                pi[h, 1, t - 1] = 0.0

    return pi


def EM(X_by_t, tags_by_t, T, maxIter=1000, tol=1e-8):
    """Full EM algorithm for the time-based ASM model.

    Parameters
    ----------
    X_by_t : dict {t: ndarray (I_t, J)}
        Binary methylation matrices per time period (1-indexed keys).
    tags_by_t : dict {t: ndarray of str (I_t,)}
        Haplotype tags per time period.
    T : int
        Number of time periods.
    maxIter : int
        Maximum number of EM iterations.
    tol : float
        Convergence tolerance on log-likelihood change.

    Returns
    -------
    theta : ndarray, shape (J, 2, 2)
    pi : ndarray, shape (2, 2, T)
    gamma_by_t : dict {t: ndarray (I_t, 2, 2)}
    n_iters : int
        Number of EM iterations actually run before convergence or maxIter.
    """
    J = X_by_t[next(iter(X_by_t))].shape[1]

    # Initialize theta randomly, seeded from the input so re-runs reproduce.
    np.random.seed(_seed_from_X(X_by_t))
    theta = np.random.uniform(0.01, 0.99, size=(J, 2, 2))

    # Initialize pi: uniform, with t=1 boundary enforced
    pi = np.ones((2, 2, T)) / 4
    pi[:, 1, 0] = 0.0      # k=1 (altered) does not exist at t=1
    pi[:, 0, 0] = 0.5       # at t=1, only k=0, split equally between alleles

    previous_ll = -np.inf
    gamma_by_t = {}

    for iteration in range(maxIter):
        # E-step: compute gamma for each time period
        gamma_by_t = {}
        for t in range(1, T + 1):
            gamma_by_t[t] = doE(X_by_t[t], theta, pi[:, :, t - 1], tags_by_t[t], t)

        # M-step: update theta (pooled across t) and pi (per t)
        theta = doM_theta(X_by_t, gamma_by_t)
        pi = doM_pi(gamma_by_t, T)

        # Enforce t=1 boundary on pi (should already hold via gamma, but be explicit)
        pi[:, 1, 0] = 0.0

        # Convergence check: observed-data log-likelihood
        ll = _observed_loglikelihood(X_by_t, tags_by_t, theta, pi, T)
        if abs(ll - previous_ll) < tol:
            break
        previous_ll = ll

    n_iters = iteration + 1
    return theta, pi, gamma_by_t, n_iters


def _observed_loglikelihood(X_by_t, tags_by_t, theta, pi, T):
    """Compute observed-data log-likelihood for convergence monitoring.

    For each read, marginalizes over the relevant (h, k) components.
    Uses log-sum-exp for numerical stability.
    """
    ll = 0.0

    for t in range(1, T + 1):
        X_t = X_by_t[t]
        tags_t = tags_by_t[t]
        pi_t = pi[:, :, t - 1]
        I_t, J = X_t.shape

        # Precompute NaN-safe matrices
        xformeth = X_t.copy()
        xformeth[np.isnan(xformeth)] = 0
        xforunmeth = X_t.copy()
        xforunmeth[np.isnan(xforunmeth)] = 1
        xforunmeth = 1 - xforunmeth

        # Compute log P(x_i | h, k) + log pi_{h,k,t} for each (h, k)
        # Shape: (I_t, 2, 2)
        log_terms = np.full((I_t, 2, 2), -np.inf)

        for h in range(2):
            for k in range(2):
                if pi_t[h, k] < PSEUDO:
                    continue
                log_theta = np.log(theta[:, h, k] + PSEUDO)
                log_1_minus_theta = np.log(1 - theta[:, h, k] + PSEUDO)
                ll_per_read = xformeth @ log_theta + xforunmeth @ log_1_minus_theta
                log_terms[:, h, k] = np.log(pi_t[h, k] + PSEUDO) + ll_per_read

        # Vectorized per-read log-sum-exp by tag group.
        h1_mask = tags_t == "H1"
        h2_mask = tags_t == "H2"
        noh_mask = ~(h1_mask | h2_mask)

        if np.any(h1_mask):
            ll += float(_logsumexp_batch(log_terms[h1_mask, 0, :], axis=-1).sum())
        if np.any(h2_mask):
            ll += float(_logsumexp_batch(log_terms[h2_mask, 1, :], axis=-1).sum())
        if np.any(noh_mask):
            n_noh = int(noh_mask.sum())
            flat = log_terms[noh_mask].reshape(n_noh, 4)
            ll += float(_logsumexp_batch(flat, axis=-1).sum())

    return ll


def compute_null_BIC(X_by_t, T):
    """BIC for the null model: single theta_j shared across all reads and time periods.

    Parameters
    ----------
    X_by_t : dict {t: ndarray (I_t, J)}
    T : int

    Returns
    -------
    BIC_null : float
    """
    any_t = next(iter(X_by_t))
    J = X_by_t[any_t].shape[1]

    total_meth = np.zeros(J)
    total_obs = np.zeros(J)
    total_reads = 0

    for t in X_by_t:
        X = X_by_t[t]
        observed = ~np.isnan(X)
        X_clean = X.copy()
        X_clean[np.isnan(X_clean)] = 0
        total_meth += X_clean.sum(axis=0)
        total_obs += observed.astype(float).sum(axis=0)
        total_reads += X.shape[0]

    theta_null = total_meth / (total_obs + PSEUDO)

    # Log-likelihood under null model
    ll_null = 0.0
    for t in X_by_t:
        X = X_by_t[t]
        observed = ~np.isnan(X)
        X_clean = X.copy()
        X_clean[np.isnan(X_clean)] = 0
        X_unmeth = observed.astype(float) - X_clean

        log_theta = np.log(theta_null + PSEUDO)
        log_1_minus_theta = np.log(1 - theta_null + PSEUDO)

        # Handle 0*log(0): where x=0 and theta=0, or x=0 and theta=1
        meth_contrib = X_clean * log_theta
        unmeth_contrib = X_unmeth * log_1_minus_theta

        ll_null += np.sum(meth_contrib) + np.sum(unmeth_contrib)

    p_null = J
    BIC_null = p_null * np.log(total_reads) - 2 * ll_null
    return float(BIC_null)


def compute_M1_BIC(X_by_t, tags_by_t, T, maxIter=1000, tol=1e-8):
    """BIC for the allele-specific, no-time-emergence model (M1).

    Constrained version of the M2 EM with pi[h, 1, t] = 0 pinned for all
    (h, t) so the altered component is disallowed. Components collapse to a
    single per-allele state at k=0; theta[:, :, 1] is unused. Reuses doE,
    doM_theta, and _observed_loglikelihood from the M2 path so M0 / M1 / M2
    BICs are computed by the same code with identical numerical safeguards.

    Parameter count is 2J (per-allele theta at k=0), matching the legacy
    single-time-period M1 in scripts/EMfunctions_HG002.py.

    Parameters
    ----------
    X_by_t : dict {t: ndarray (I_t, J)}
    tags_by_t : dict {t: ndarray of str (I_t,)}
    T : int
    maxIter : int
    tol : float

    Returns
    -------
    BIC_M1 : float
    n_iters : int
    """
    J = X_by_t[next(iter(X_by_t))].shape[1]

    np.random.seed(_seed_from_X(X_by_t))
    theta = np.random.uniform(0.01, 0.99, size=(J, 2, 2))

    # pi pinned: only k=0 has mass, split equally between alleles, every t.
    # doE short-circuits on `pi_t[h, k] < PSEUDO` so gamma[:, :, 1] = 0
    # automatically; doM_theta then leaves theta[:, :, 1] at its default 0.5
    # (unused, never enters the LL).
    pi = np.zeros((2, 2, T))
    pi[:, 0, :] = 0.5

    previous_ll = -np.inf
    ll = -np.inf
    for iteration in range(maxIter):
        gamma_by_t = {}
        for t in range(1, T + 1):
            gamma_by_t[t] = doE(X_by_t[t], theta, pi[:, :, t - 1], tags_by_t[t], t)

        # M-step for theta only; pi stays pinned across iterations.
        theta = doM_theta(X_by_t, gamma_by_t)

        ll = _observed_loglikelihood(X_by_t, tags_by_t, theta, pi, T)
        if abs(ll - previous_ll) < tol:
            break
        previous_ll = ll

    n_total = sum(X_by_t[t].shape[0] for t in X_by_t)
    p_M1 = 2 * J
    BIC_M1 = p_M1 * np.log(n_total) - 2 * ll
    return float(BIC_M1), iteration + 1


def compute_alt_loglikelihood(X_by_t, tags_by_t, theta, pi, T):
    """Observed-data log-likelihood for the alternative (mixture) model.

    Same computation as _observed_loglikelihood but exposed as public API
    for BIC computation.

    Parameters
    ----------
    X_by_t : dict {t: ndarray (I_t, J)}
    tags_by_t : dict {t: ndarray of str (I_t,)}
    theta : ndarray, shape (J, 2, 2)
    pi : ndarray, shape (2, 2, T)
    T : int

    Returns
    -------
    ll : float
    """
    return _observed_loglikelihood(X_by_t, tags_by_t, theta, pi, T)
