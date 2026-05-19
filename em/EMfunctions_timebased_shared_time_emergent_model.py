# ABOUTME: EM for the "Shared with time emergence" model — no haplotype axis, single simplex per t.
# ABOUTME: Derivation in latex/EM_stepbystep_derivation_emergent_altered_model.tex. p = 2J + (T-1).

import hashlib

import numpy as np


def _seed_from_X(X_by_t):
    """Deterministic 32-bit seed derived from the per-timepoint read matrices.

    Same input matrices -> same seed -> same random theta init -> same EM
    trajectory. Mirrors the diploid module's helper so reproducibility tests
    can share the pattern.
    """
    h = hashlib.md5()
    for t in sorted(X_by_t):
        h.update(np.ascontiguousarray(X_by_t[t]).tobytes())
    return int(h.hexdigest()[:8], 16)


np.seterr(all='raise')

PSEUDO = 1e-10


def doE(X_t, theta, pi_t, t):
    """E-step for one time period (LaTeX Eqs. 2–4).

    Parameters
    ----------
    X_t : ndarray, shape (I_t, J)
        Binary methylation matrix; NaN for sites the read does not cover.
    theta : ndarray, shape (J, 2)
        theta[j, k] = P(methylated at site j | program k). k=0 normal, k=1 altered.
    pi_t : ndarray, shape (2,)
        Mixing weights at time t. pi_t[1] = 0 enforces the baseline anchor at t=1.
    t : int
        Time index, 1-indexed. Used only for the t=1 boundary safety net.

    Returns
    -------
    gamma : ndarray, shape (I_t, 2)
        Posterior probabilities. gamma[i, k] = P(z_{i,t,k}=1 | X, theta, pi).
    """
    I_t, J = X_t.shape

    # NaN-safe Bernoulli factor: methylated sites contribute via xformeth,
    # unmethylated via xforunmeth; NaN sites contribute 0 in both.
    xformeth = X_t.copy()
    xformeth[np.isnan(xformeth)] = 0
    xforunmeth = X_t.copy()
    xforunmeth[np.isnan(xforunmeth)] = 1
    xforunmeth = 1 - xforunmeth

    log_gamma_prime = np.full((I_t, 2), -np.inf)
    for k in range(2):
        if pi_t[k] < PSEUDO:
            continue
        log_theta = np.log(theta[:, k] + PSEUDO)
        log_1_minus_theta = np.log(1 - theta[:, k] + PSEUDO)
        ll_per_read = xformeth @ log_theta + xforunmeth @ log_1_minus_theta
        log_gamma_prime[:, k] = np.log(pi_t[k] + PSEUDO) + ll_per_read

    # Boundary at t=1: altered program does not exist. pi_t[1] = 0 already
    # forces gamma[:, 1] = 0 via the PSEUDO short-circuit above, but be
    # explicit so the invariant is visible.
    if t == 1:
        log_gamma_prime[:, 1] = -np.inf

    return _normalize_log_batch(log_gamma_prime, axis=-1)


def _normalize_log_batch(log_values, axis=-1):
    """Vectorized softmax over `axis` with -inf masking.

    Slices that are all -inf (or whose exp-shifted sum is 0) return uniform;
    otherwise returns normalized probabilities summing to 1 along `axis`.
    """
    finite_mask = np.isfinite(log_values)
    any_finite = np.any(finite_mask, axis=axis, keepdims=True)

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

    Slices entirely -inf contribute 0 (caller is summing across reads).
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
    """M-step for theta (LaTeX Eq. 9): weighted Bernoulli MLE pooled across t.

    theta[j, k] = sum_{t,i} gamma[i,k] * x_{i,j} / sum_{t,i} gamma[i,k] * observed_{i,j}

    Parameters
    ----------
    X_by_t : dict {t: ndarray (I_t, J)}
    gamma_by_t : dict {t: ndarray (I_t, 2)}

    Returns
    -------
    theta : ndarray, shape (J, 2)
    """
    any_t = next(iter(X_by_t))
    J = X_by_t[any_t].shape[1]
    theta = np.full((J, 2), 0.5)

    for k in range(2):
        numerator = np.zeros(J)
        denominator = np.zeros(J)
        for t in X_by_t:
            X_t = X_by_t[t]
            g = gamma_by_t[t][:, k]                # (I_t,)
            observed = ~np.isnan(X_t)              # (I_t, J)
            X_clean = X_t.copy()
            X_clean[np.isnan(X_clean)] = 0
            numerator += (g[:, np.newaxis] * X_clean).sum(axis=0)
            denominator += (g[:, np.newaxis] * observed.astype(float)).sum(axis=0)

        has_data = denominator > PSEUDO
        theta[has_data, k] = numerator[has_data] / denominator[has_data]

    return theta


def doM_pi(gamma_by_t, T):
    """M-step for pi (LaTeX Eq. 16): pi[k, t] = (sum_i gamma[i, k]) / I_t.

    Single simplex per t (sum_k pi[k, t] = 1), no factor of 1/2. Baseline
    anchor pi[1, 0] = 0 emerges automatically: at t=1 the E-step forces
    gamma[:, 1] = 0, so N_{2,1} = 0 and pi[1, 0] = 0 / I_1 = 0 (LaTeX §10.6).

    Parameters
    ----------
    gamma_by_t : dict {t: ndarray (I_t, 2)}
    T : int

    Returns
    -------
    pi : ndarray, shape (2, T)
    """
    pi = np.zeros((2, T))

    for t in range(1, T + 1):
        gamma_t = gamma_by_t[t]                   # (I_t, 2)
        I_t = gamma_t.shape[0]
        if I_t == 0:
            # Degenerate: no reads at this time. Fall back to (1, 0) — the
            # constraint-feasible default — and let the caller deal with it.
            pi[0, t - 1] = 1.0
            pi[1, t - 1] = 0.0
            continue
        N_kt = gamma_t.sum(axis=0)                # (2,)
        pi[:, t - 1] = N_kt / I_t

    return pi


def EM(X_by_t, T, maxIter=1000, tol=1e-8):
    """Full EM loop for the Shared-with-time-emergence model.

    Parameters
    ----------
    X_by_t : dict {t: ndarray (I_t, J)}
    T : int
    maxIter, tol : convergence controls

    Returns
    -------
    theta : ndarray, shape (J, 2)
    pi : ndarray, shape (2, T)
    gamma_by_t : dict {t: ndarray (I_t, 2)}
    n_iters : int
    """
    J = X_by_t[next(iter(X_by_t))].shape[1]

    np.random.seed(_seed_from_X(X_by_t))
    theta = np.random.uniform(0.01, 0.99, size=(J, 2))

    # pi init: uniform at t>=2, anchored at t=1.
    pi = np.full((2, T), 0.5)
    pi[0, 0] = 1.0
    pi[1, 0] = 0.0

    previous_ll = -np.inf
    gamma_by_t = {}
    for iteration in range(maxIter):
        gamma_by_t = {}
        for t in range(1, T + 1):
            gamma_by_t[t] = doE(X_by_t[t], theta, pi[:, t - 1], t)

        theta = doM_theta(X_by_t, gamma_by_t)
        pi = doM_pi(gamma_by_t, T)
        # Re-assert the baseline anchor — doM_pi already produces it for free
        # (LaTeX §10.6) but explicit beats implicit.
        pi[1, 0] = 0.0
        pi[0, 0] = 1.0

        ll = _observed_loglikelihood(X_by_t, theta, pi, T)
        if abs(ll - previous_ll) < tol:
            break
        previous_ll = ll

    n_iters = iteration + 1
    return theta, pi, gamma_by_t, n_iters


def _observed_loglikelihood(X_by_t, theta, pi, T):
    """Observed-data log-likelihood (LaTeX Eq. 19)."""
    ll = 0.0
    for t in range(1, T + 1):
        X_t = X_by_t[t]
        pi_t = pi[:, t - 1]
        I_t, J = X_t.shape

        xformeth = X_t.copy()
        xformeth[np.isnan(xformeth)] = 0
        xforunmeth = X_t.copy()
        xforunmeth[np.isnan(xforunmeth)] = 1
        xforunmeth = 1 - xforunmeth

        log_terms = np.full((I_t, 2), -np.inf)
        for k in range(2):
            if pi_t[k] < PSEUDO:
                continue
            log_theta = np.log(theta[:, k] + PSEUDO)
            log_1_minus_theta = np.log(1 - theta[:, k] + PSEUDO)
            ll_per_read = xformeth @ log_theta + xforunmeth @ log_1_minus_theta
            log_terms[:, k] = np.log(pi_t[k] + PSEUDO) + ll_per_read

        ll += float(_logsumexp_batch(log_terms, axis=-1).sum())

    return ll


def compute_shared_BIC(X_by_t, T):
    """BIC for the Shared (null) model: one theta_j pooled across reads and t.

    Mathematically identical to compute_null_BIC in the diploid module — the
    pooled-theta MLE has no h-axis dependency. Re-implemented here so this
    module is self-contained.

    p = J. Returns float.
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

    ll_null = 0.0
    for t in X_by_t:
        X = X_by_t[t]
        observed = ~np.isnan(X)
        X_clean = X.copy()
        X_clean[np.isnan(X_clean)] = 0
        X_unmeth = observed.astype(float) - X_clean

        log_theta = np.log(theta_null + PSEUDO)
        log_1_minus_theta = np.log(1 - theta_null + PSEUDO)
        ll_null += np.sum(X_clean * log_theta) + np.sum(X_unmeth * log_1_minus_theta)

    p_null = J
    BIC_null = p_null * np.log(total_reads) - 2 * ll_null
    return float(BIC_null)


def compute_shared_time_BIC(X_by_t, T, maxIter=1000, tol=1e-8):
    """BIC for the Shared-with-time-emergence model.

    Runs the no-h EM, evaluates the observed log-likelihood, and computes
    BIC with p = 2J + (T - 1) (LaTeX §11): 2J for theta_{j,k} with
    k in {0,1}, T-1 for pi (one free param per t>=2 on the simplex; t=1 is
    anchored).

    Parameters
    ----------
    X_by_t : dict {t: ndarray (I_t, J)}
    T : int

    Returns
    -------
    BIC : float
    n_iters : int
    """
    theta, pi, _, n_iters = EM(X_by_t, T, maxIter=maxIter, tol=tol)
    ll = _observed_loglikelihood(X_by_t, theta, pi, T)
    J = X_by_t[next(iter(X_by_t))].shape[1]
    n_total = sum(X_by_t[t].shape[0] for t in X_by_t)
    p = 2 * J + (T - 1)
    BIC = p * np.log(n_total) - 2 * ll
    return float(BIC), n_iters


def compute_alt_loglikelihood(X_by_t, theta, pi, T):
    """Public alias for the observed log-likelihood under (theta, pi)."""
    return _observed_loglikelihood(X_by_t, theta, pi, T)
