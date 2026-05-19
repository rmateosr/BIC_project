# ABOUTME: Generates synthetic methylation data for evaluating the time-based ASM pipeline.
# ABOUTME: Produces per-region TSVs in t1/.../tT/ directories plus an answer-key manifest.
"""Synthetic data generator for time-based ASM pipeline evaluation.

Generates region files in the format consumed by BIC_algorithm_timebased.py:
  - Tab-separated with header: readlabel, chrom, strand, startcoord, status, haplotype
  - One file per region: methylationfraction_<start>_<end>_.tsv
  - One directory per time period: t1/, t2/, ..., tT/

Four region classes:
  NULL                        -- identical methylation on both alleles (no ASM)
  STATIC_ASM                  -- alleles differ at k=0; π[h,1,t]=0; theta[:,:,1] is
                                 free (EM must learn k=1 unused via π, not via θ)
  TIME_EMERGENT_ASM           -- altered pattern emerges symmetrically: both alleles
                                 drift π[h,1,t] together
  ASYMMETRIC_TIME_EMERGENT_ASM -- only one allele's π[h,1,t] drifts; the other stays
                                 at zero. The canonical ASM time-emergence event.
                                 Per-region random choice of which allele drifts.
  SHARED_TIME_EMERGENT        -- alleles share theta at both k=0 and k=1, and both
                                 drift π[h,1,t] together. Pure Class C: no allele
                                 difference at any layer (theta or pi).

All three emergence classes use a NON-LINEAR per-region random trajectory:
t=1 anchored at 0, t>=2 drawn iid from Uniform(0, pi_final) — so the truth
spans the full parametric freedom of the models (free π[h,1,t] per t) instead
of baking in a linear monotonic ramp. Explicit --pi-trajectory* CLI overrides
still produce a deterministic shape (escape hatch for tests).

Also writes manifest.tsv (answer key) for downstream evaluation. The manifest's
pi_h0_alt_t{ti} and pi_h1_alt_t{ti} columns now carry per-allele truth — for
ASYMMETRIC regions one is 0 and the other is nonzero; for SHARED_TIME_EMERGENT
both columns are equal.
"""

import argparse
import os
import uuid

import numpy as np
import pandas as pd

# Region classes
NULL = "NULL"
STATIC_ASM = "STATIC_ASM"
TIME_EMERGENT_ASM = "TIME_EMERGENT_ASM"
ASYMMETRIC_TIME_EMERGENT_ASM = "ASYMMETRIC_TIME_EMERGENT_ASM"
SHARED_TIME_EMERGENT = "SHARED_TIME_EMERGENT"
CLASSES = [
    NULL, STATIC_ASM, TIME_EMERGENT_ASM, ASYMMETRIC_TIME_EMERGENT_ASM,
    SHARED_TIME_EMERGENT,
]


def _random_emergent_trajectory(T, pi_final, rng):
    """Per-region random emergence trajectory: t=1 anchored at 0; t>=2 iid Uniform(0, pi_final).

    The Class B / Class C time-emergent models leave π[h,1,t] (or its no-h
    sibling π[1,t]) free for every t>=2. Sampling each entry independently
    from a uniform — rather than a linear ramp 0 → pi_final — means the truth
    spans the full parametric freedom of the model: monotonic, non-monotonic,
    early-burst, late-burst, partial reversal, etc. all appear at equal rates.
    """
    traj = np.zeros(T)
    if T > 1:
        traj[1:] = rng.uniform(0.0, pi_final, size=T - 1)
    return traj


def _build_single_trajectory(T, pi_final, trajectory, label):
    """Return shape (T,) per-time-point pi_altered for one allele.

    Default: linear ramp 0 → pi_final. Custom: explicit list of length T with
    trajectory[0]==0 (the t=1 boundary condition).
    """
    if trajectory is not None:
        traj = np.array(trajectory, dtype=float)
        if len(traj) != T:
            raise ValueError(f"trajectory_{label} length {len(traj)} != T={T}")
        if traj[0] != 0:
            raise ValueError(
                f"trajectory_{label}[0] must be 0 (boundary condition at t=1)"
            )
        return traj
    if T == 1:
        return np.array([0.0])
    return np.linspace(0, pi_final, T)


def make_pi_trajectory_per_allele(
    T, pi_final_h0, pi_final_h1, trajectory_h0=None, trajectory_h1=None,
):
    """Build per-allele, per-time-point pi_altered values.

    Returns
    -------
    traj : ndarray, shape (T, 2)
        traj[t-1, h] = within-allele probability that a read on allele h is
        from the altered component (k=1) at time t. Equivalent to
        2*π[h, 1, t] under the diploid marginal constraint.
        Boundary: traj[0, :] == 0 (no altered component at t=1).
    """
    col_h0 = _build_single_trajectory(T, pi_final_h0, trajectory_h0, "h0")
    col_h1 = _build_single_trajectory(T, pi_final_h1, trajectory_h1, "h1")
    return np.stack([col_h0, col_h1], axis=1)


def generate_coords(J, chrom_offset, rng):
    """Generate J sorted CpG coordinates with jittered spacing.

    Spacing ~ Exponential(mean=200bp), mimicking CpG-island-like density.
    """
    spacings = rng.exponential(200, size=J).astype(int).clip(min=10)
    coords = np.cumsum(spacings) + chrom_offset
    return coords.astype(int)


def generate_region_theta(region_class, J, rng, theta_noise=0.05, theta_iid=False):
    """Generate true theta parameters for a region.

    Parameters
    ----------
    region_class : str
    J : int
    rng : np.random.RandomState
    theta_noise : float
        Per-CpG Gaussian noise added on top of the per-(h, k) region base value.
        Ignored when theta_iid=True.
    theta_iid : bool
        If False (default), every CpG in a region shares one base value per
        (h, k) cell with iid Gaussian noise of width theta_noise. If True, each
        θ[j, h, k] is drawn independently from the (h, k)-specific uniform range
        and theta_noise is ignored. The iid mode produces much higher within-
        region variance in θ — useful for testing the EM under realistic per-CpG
        heterogeneity rather than the "flat-region + noise" simplification.

    Returns
    -------
    theta : ndarray, shape (J, 2, 2)
        theta[j, h, k] = P(methylated at CpG j | allele h, component k).
    """
    theta = np.full((J, 2, 2), 0.5)

    if region_class == NULL:
        # Both alleles identical, both components identical.
        if theta_iid:
            vals = rng.uniform(0.3, 0.7, J)
        else:
            base = rng.uniform(0.3, 0.7)
            vals = np.clip(base + rng.normal(0, theta_noise, J), 0.01, 0.99)
        theta[:, :, :] = vals[:, None, None]

    elif region_class == STATIC_ASM:
        # k=0: strong allele difference (the static ASM signal).
        # Ranges extend to the endpoints (0.0 and 1.0) so the truth can be
        # fully unmethylated / fully methylated.
        if theta_iid:
            theta[:, 0, 0] = rng.uniform(0.75, 1.00, J)
            theta[:, 1, 0] = rng.uniform(0.00, 0.25, J)
        else:
            base_h0_k0 = rng.uniform(0.75, 1.00)
            base_h1_k0 = rng.uniform(0.00, 0.25)
            theta[:, 0, 0] = np.clip(base_h0_k0 + rng.normal(0, theta_noise, J), 0.01, 0.99)
            theta[:, 1, 0] = np.clip(base_h1_k0 + rng.normal(0, theta_noise, J), 0.01, 0.99)
        # k=1: free per-CpG values. Truth has π[:,1,t]=0 everywhere, so this
        # component is unused by the data generator. We let θ[:,:,1] vary freely
        # so the EM has to learn "k=1 is unused" from π alone, not from the
        # degenerate θ similarity that the previous version baked in.
        theta[:, 0, 1] = rng.uniform(0.01, 0.99, J)
        theta[:, 1, 1] = rng.uniform(0.01, 0.99, J)

    elif region_class in (TIME_EMERGENT_ASM, ASYMMETRIC_TIME_EMERGENT_ASM):
        # k=0 (normal): strong allele difference; k=1 (altered): distinct pattern.
        # ASYMMETRIC differs from TIME_EMERGENT only in the trajectory (one allele
        # ramps, the other stays at 0) — θ structure is identical.
        # k=0 ranges extend to 0.0/1.0 to allow fully unmethylated/methylated
        # baselines; k=1 ranges remain in the intermediate band so the altered
        # program is identifiably distinct from the baseline.
        if theta_iid:
            theta[:, 0, 0] = rng.uniform(0.75, 1.00, J)
            theta[:, 1, 0] = rng.uniform(0.00, 0.25, J)
            theta[:, 0, 1] = rng.uniform(0.25, 0.55, J)
            theta[:, 1, 1] = rng.uniform(0.45, 0.75, J)
        else:
            base_h0_k0 = rng.uniform(0.75, 1.00)
            base_h1_k0 = rng.uniform(0.00, 0.25)
            base_h0_k1 = rng.uniform(0.25, 0.55)
            base_h1_k1 = rng.uniform(0.45, 0.75)
            theta[:, 0, 0] = np.clip(base_h0_k0 + rng.normal(0, theta_noise, J), 0.01, 0.99)
            theta[:, 1, 0] = np.clip(base_h1_k0 + rng.normal(0, theta_noise, J), 0.01, 0.99)
            theta[:, 0, 1] = np.clip(base_h0_k1 + rng.normal(0, theta_noise, J), 0.01, 0.99)
            theta[:, 1, 1] = np.clip(base_h1_k1 + rng.normal(0, theta_noise, J), 0.01, 0.99)

    elif region_class == SHARED_TIME_EMERGENT:
        # Pure Class C: both alleles share theta at k=0 AND k=1. The k=1 program
        # is a genuinely distinct profile from k=0 (so the mixture is identifiable),
        # but neither allele can be told apart at any single read. Combined with a
        # symmetric pi trajectory in generate_synthetic_dataset, this is the truth
        # the Shared-with-time-emergence model is designed to recover.
        if theta_iid:
            shared_k0 = rng.uniform(0.6, 0.85, J)
            shared_k1 = rng.uniform(0.1, 0.35, J)
        else:
            base_k0 = rng.uniform(0.6, 0.85)
            base_k1 = rng.uniform(0.1, 0.35)
            shared_k0 = np.clip(base_k0 + rng.normal(0, theta_noise, J), 0.01, 0.99)
            shared_k1 = np.clip(base_k1 + rng.normal(0, theta_noise, J), 0.01, 0.99)
        theta[:, 0, 0] = shared_k0
        theta[:, 1, 0] = shared_k0
        theta[:, 0, 1] = shared_k1
        theta[:, 1, 1] = shared_k1

    return theta


def generate_reads_for_region(
    coords, theta, pi_alt_t_h0, pi_alt_t_h1, n_reads, chrom, rng,
    hap_probs=(0.3, 0.3, 0.4), min_cpg_span=3,
):
    """Generate reads for one region at one time point.

    Parameters
    ----------
    coords : ndarray of int, shape (J,)
        CpG positions.
    theta : ndarray, shape (J, 2, 2)
        True methylation probabilities.
    pi_alt_t_h0, pi_alt_t_h1 : float
        Within-allele probability that a read on allele h is from the altered
        component (k=1) at this time. Equivalent to 2*π[h, 1, t] under the
        diploid marginal constraint. The t=1 boundary (k=1 forbidden) is
        enforced by the caller passing zeros at t=1.
    n_reads : int
    chrom : str
    rng : np.random.RandomState
    hap_probs : tuple of 3 floats
        (P(H1), P(H2), P(noH)).
    min_cpg_span : int
        Minimum CpGs a read must cover.

    Returns
    -------
    rows : list of dict
    """
    J = len(coords)
    rows = []

    # Joint distribution over (h, k) for noH reads — single 4-way draw, not a
    # factorized (h then k) draw. Required for asymmetric π: when one allele
    # ramps and the other doesn't, P(k=1 | h) depends on h, so the two
    # marginals do not factorize.
    # Order: (h=0,k=0), (h=0,k=1), (h=1,k=0), (h=1,k=1). Sums to 1.
    joint_pi = np.array([
        0.5 * (1 - pi_alt_t_h0),
        0.5 * pi_alt_t_h0,
        0.5 * (1 - pi_alt_t_h1),
        0.5 * pi_alt_t_h1,
    ])

    for _ in range(n_reads):
        hap = rng.choice(["H1", "H2", "noH"], p=hap_probs)
        if hap == "H1":
            true_h = 0
            true_k = 1 if rng.rand() < pi_alt_t_h0 else 0
        elif hap == "H2":
            true_h = 1
            true_k = 1 if rng.rand() < pi_alt_t_h1 else 0
        else:
            idx = rng.choice(4, p=joint_pi)
            true_h, true_k = int(idx // 2), int(idx % 2)

        # Contiguous CpG span (uniform coverage model)
        span = rng.randint(min_cpg_span, J + 1)
        start = rng.randint(0, J - span + 1)
        sel = slice(start, start + span)
        sel_coords = coords[sel]

        # Sample methylation from Bernoulli(theta)
        probs = theta[sel, true_h, true_k]
        status = (rng.rand(span) < probs).astype(int)

        rows.append({
            "readlabel": uuid.uuid4().hex[:16],
            "chrom": chrom,
            "strand": "+",
            "startcoord": ",".join(str(c) for c in sel_coords),
            "status": ",".join(str(s) for s in status),
            "haplotype": hap,
        })

    return rows


def generate_synthetic_dataset(
    n_regions=50,
    reads_per_t=80,
    J=50,
    T=4,
    class_mix=None,
    pi_final=0.6,
    pi_trajectory=None,
    pi_final_h0=None,
    pi_final_h1=None,
    pi_trajectory_h0=None,
    pi_trajectory_h1=None,
    theta_noise=0.05,
    theta_iid=False,
    windowsize=10,
    chrom="chrSYN",
    output_dir=None,
    seed=None,
    hap_probs=None,
):
    """Generate a complete synthetic dataset for pipeline evaluation.

    Writes per-region TSV files into t1/...tT/ directories and a manifest.tsv
    answer key.

    Parameters
    ----------
    n_regions : int
    reads_per_t : int
    J : int
        CpG sites per region.
    T : int
        Number of time periods.
    class_mix : dict, optional
        {class_name: weight} over CLASSES. Default: equal weights for all four.
    pi_final : float
        Default within-allele altered fraction at t=T. Used for TIME_EMERGENT
        regions when neither --pi-final-h0/h1 nor --pi-trajectory-h0/h1 is set,
        and for ASYMMETRIC regions on whichever allele is selected to ramp.
    pi_trajectory : list of float, optional
        Explicit per-t altered fraction (length T, first must be 0). Overrides
        pi_final for the ASYMMETRIC ramping allele.
    pi_final_h0, pi_final_h1 : float, optional
        Per-allele overrides for TIME_EMERGENT regions. If unset, fall back to
        pi_final. ASYMMETRIC regions ignore these (they always pin one allele
        to zero).
    pi_trajectory_h0, pi_trajectory_h1 : list of float, optional
        Explicit per-allele trajectories for TIME_EMERGENT regions. Override
        the corresponding pi_final_h*.
    theta_noise : float
        Per-CpG Gaussian noise around each region-level (h, k) base value.
        Ignored when theta_iid=True.
    theta_iid : bool
        If True, each θ[j, h, k] is drawn independently from the (h, k)-specific
        uniform range with no Gaussian noise added. See generate_region_theta
        for per-class behavior. Default False (legacy "flat base + noise").
    windowsize : int
        Sliding window size (for manifest coord_evaluated computation).
    chrom : str
    output_dir : str
    seed : int, optional
    hap_probs : tuple of 3 floats, optional
        (P(H1), P(H2), P(noH)) used for every read in every region. Must sum
        to 1. Default (0.3, 0.3, 0.4) preserves prior behavior.

    Returns
    -------
    manifest : DataFrame
        Answer key with per-allele pi truth columns. For ASYMMETRIC regions one
        of pi_h{0,1}_alt_t{ti} is 0 and the other is nonzero.
    """
    if class_mix is None:
        class_mix = {c: 1 for c in CLASSES}
    if output_dir is None:
        output_dir = "synthetic_output"
    if hap_probs is None:
        hap_probs = (0.3, 0.3, 0.4)
    hap_probs = tuple(float(x) for x in hap_probs)
    if len(hap_probs) != 3 or abs(sum(hap_probs) - 1.0) > 1e-6:
        raise ValueError(
            f"hap_probs must be a length-3 tuple summing to 1; got {hap_probs}"
        )

    rng = np.random.RandomState(seed)

    # Explicit-trajectory escape hatches. These produce deterministic ramps
    # only when the user passes --pi-trajectory* on the CLI (or kwargs to this
    # function). The default emergence path below uses per-region random draws
    # — see _random_emergent_trajectory and the project rule against assuming
    # linearity of π emergence in BIC_project/CLAUDE.md.
    explicit_base = pi_trajectory is not None
    explicit_te = (
        pi_trajectory_h0 is not None or pi_trajectory_h1 is not None
        or pi_final_h0 is not None or pi_final_h1 is not None
    )
    base_traj_explicit = (
        _build_single_trajectory(T, pi_final, pi_trajectory, "base")
        if explicit_base else None
    )
    te_traj_explicit = (
        make_pi_trajectory_per_allele(
            T,
            pi_final_h0 if pi_final_h0 is not None else pi_final,
            pi_final_h1 if pi_final_h1 is not None else pi_final,
            pi_trajectory_h0 if pi_trajectory_h0 is not None else pi_trajectory,
            pi_trajectory_h1 if pi_trajectory_h1 is not None else pi_trajectory,
        )
        if explicit_te else None
    )

    # Assign classes to regions
    weights = np.array([class_mix.get(c, 0) for c in CLASSES], dtype=float)
    weights /= weights.sum()
    region_classes = rng.choice(CLASSES, size=n_regions, p=weights)

    # Create time-period directories
    t_dirs = {}
    for t in range(1, T + 1):
        d = os.path.join(output_dir, f"t{t}")
        os.makedirs(d, exist_ok=True)
        t_dirs[t] = d

    manifest_rows = []
    offset = 1_000_000

    for ridx in range(n_regions):
        rc = region_classes[ridx]
        coords = generate_coords(J, offset, rng)
        rstart = int(coords[0])
        rend = int(coords[-1])
        offset = rend + 100_000

        theta = generate_region_theta(rc, J, rng, theta_noise, theta_iid)

        # Per-class trajectory. Shape (T, 2): per-time, per-allele within-allele
        # altered fractions. Recall pi[h, 1, t] = 0.5 * region_traj[t-1, h] under
        # the diploid marginal.
        if rc == TIME_EMERGENT_ASM:
            if te_traj_explicit is not None:
                region_traj = te_traj_explicit.copy()
            else:
                # Symmetric in π, asymmetric in θ: both alleles drift on the
                # SAME per-region random trajectory. The allele difference is
                # carried entirely by θ (see generate_region_theta). Drawing
                # one trajectory and sharing it across alleles keeps the class
                # semantically distinct from SHARED_TIME_EMERGENT only in θ.
                shared_traj = _random_emergent_trajectory(T, pi_final, rng)
                region_traj = np.stack([shared_traj, shared_traj], axis=1)
            ramping_allele = -1
        elif rc == ASYMMETRIC_TIME_EMERGENT_ASM:
            region_traj = np.zeros((T, 2))
            ramping_allele = int(rng.choice([0, 1]))
            if base_traj_explicit is not None:
                region_traj[:, ramping_allele] = base_traj_explicit
            else:
                region_traj[:, ramping_allele] = _random_emergent_trajectory(
                    T, pi_final, rng,
                )
        elif rc == SHARED_TIME_EMERGENT:
            # Pure Class C: both alleles drift on the same per-region random
            # trajectory. Sharing the draw across alleles is what makes it
            # bilaterally symmetric (vs. TIME_EMERGENT_ASM, where each allele
            # draws independently).
            shared_traj = _random_emergent_trajectory(T, pi_final, rng)
            region_traj = np.stack([shared_traj, shared_traj], axis=1)
            ramping_allele = -1
        else:  # NULL, STATIC_ASM
            region_traj = np.zeros((T, 2))
            ramping_allele = -1

        # Write reads for each time period
        fname = f"methylationfraction_{rstart}_{rend}_.tsv"
        for t in range(1, T + 1):
            reads = generate_reads_for_region(
                coords, theta,
                region_traj[t - 1, 0], region_traj[t - 1, 1],
                reads_per_t, chrom, rng,
                hap_probs=hap_probs,
            )
            df = pd.DataFrame(reads)
            df.to_csv(os.path.join(t_dirs[t], fname), sep="\t", index=False)

        # Manifest: one row per sliding window
        n_windows = max(1, J - windowsize + 1)
        for w in range(n_windows):
            coord_eval = int(coords[w]) if J > windowsize else int(coords[0])

            row = {
                "region_idx": ridx,
                "chrom": chrom,
                "region_start": rstart,
                "region_end": rend,
                "coord_evaluated": coord_eval,
                "true_class": rc,
                "ramping_allele": ramping_allele,
                "true_theta_h0_k0": float(np.mean(theta[:, 0, 0])),
                "true_theta_h1_k0": float(np.mean(theta[:, 1, 0])),
                "true_theta_h0_k1": float(np.mean(theta[:, 0, 1])),
                "true_theta_h1_k1": float(np.mean(theta[:, 1, 1])),
                "expected_BICresult": 0 if rc == NULL else 1,
            }
            # Per-allele pi truth: pi[h, 1, t] = 0.5 * region_traj[t-1, h].
            for ti in range(1, T + 1):
                row[f"pi_h0_alt_t{ti}"] = region_traj[ti - 1, 0] / 2
                row[f"pi_h1_alt_t{ti}"] = region_traj[ti - 1, 1] / 2

            manifest_rows.append(row)

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(os.path.join(output_dir, "manifest.tsv"), sep="\t", index=False)

    return manifest


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic methylation data for time-based ASM evaluation."
    )
    parser.add_argument("--n-regions", type=int, default=50)
    parser.add_argument("--reads-per-t", type=int, default=80)
    parser.add_argument("--cpgs-per-region", type=int, default=50, dest="J")
    parser.add_argument("--T", type=int, default=4)
    parser.add_argument("--pi-final", type=float, default=0.6)
    parser.add_argument(
        "--pi-trajectory", type=str, default=None,
        help="Comma-separated per-t pi values (e.g. 0,0.2,0.4,0.6). Overrides --pi-final.",
    )
    parser.add_argument(
        "--pi-final-h0", type=float, default=None,
        help="Per-allele pi_final for h=0 (TIME_EMERGENT only). Default: --pi-final.",
    )
    parser.add_argument(
        "--pi-final-h1", type=float, default=None,
        help="Per-allele pi_final for h=1 (TIME_EMERGENT only). Default: --pi-final.",
    )
    parser.add_argument(
        "--pi-trajectory-h0", type=str, default=None,
        help="Per-allele explicit trajectory for h=0 (TIME_EMERGENT only). Overrides --pi-final-h0.",
    )
    parser.add_argument(
        "--pi-trajectory-h1", type=str, default=None,
        help="Per-allele explicit trajectory for h=1 (TIME_EMERGENT only). Overrides --pi-final-h1.",
    )
    parser.add_argument("--theta-noise", type=float, default=0.05)
    parser.add_argument(
        "--theta-iid", action="store_true",
        help=(
            "Draw each θ[j, h, k] independently from the (h, k) uniform range; "
            "skip the legacy 'one base + Gaussian noise' scheme. --theta-noise "
            "is ignored when this flag is set."
        ),
    )
    parser.add_argument("--windowsize", type=int, default=10)
    parser.add_argument("--chrom", type=str, default="chrSYN")
    parser.add_argument("--output-dir", type=str, default="synthetic_output")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--class-mix", type=str, default=None,
        help=(
            "Comma-separated weights for "
            "NULL,STATIC_ASM,TIME_EMERGENT_ASM,ASYMMETRIC_TIME_EMERGENT_ASM,"
            "SHARED_TIME_EMERGENT (e.g. 1,1,1,1,1). Trailing classes default to "
            "0 if fewer values given."
        ),
    )
    parser.add_argument(
        "--hap-probs", type=str, default=None,
        help=(
            "Comma-separated read-haplotype probabilities P(H1),P(H2),P(noH); "
            "must sum to 1. Default: 0.3,0.3,0.4. Mutually exclusive with --noh-frac."
        ),
    )
    parser.add_argument(
        "--noh-frac", type=float, default=None,
        help=(
            "Fraction of reads tagged noH (untagged); H1 and H2 split the remainder "
            "evenly. Mutually exclusive with --hap-probs."
        ),
    )
    args = parser.parse_args()

    pi_traj = None
    if args.pi_trajectory:
        pi_traj = [float(x) for x in args.pi_trajectory.split(",")]
    pi_traj_h0 = None
    if args.pi_trajectory_h0:
        pi_traj_h0 = [float(x) for x in args.pi_trajectory_h0.split(",")]
    pi_traj_h1 = None
    if args.pi_trajectory_h1:
        pi_traj_h1 = [float(x) for x in args.pi_trajectory_h1.split(",")]

    class_mix = None
    if args.class_mix:
        vals = [float(x) for x in args.class_mix.split(",")]
        class_mix = dict(zip(CLASSES, vals))

    if args.hap_probs is not None and args.noh_frac is not None:
        parser.error("--hap-probs and --noh-frac are mutually exclusive")
    hap_probs = None
    if args.hap_probs is not None:
        hap_probs = tuple(float(x) for x in args.hap_probs.split(","))
    elif args.noh_frac is not None:
        f = args.noh_frac
        hap_probs = ((1 - f) / 2, (1 - f) / 2, f)

    manifest = generate_synthetic_dataset(
        n_regions=args.n_regions,
        reads_per_t=args.reads_per_t,
        J=args.J,
        T=args.T,
        class_mix=class_mix,
        pi_final=args.pi_final,
        pi_trajectory=pi_traj,
        pi_final_h0=args.pi_final_h0,
        pi_final_h1=args.pi_final_h1,
        pi_trajectory_h0=pi_traj_h0,
        pi_trajectory_h1=pi_traj_h1,
        theta_noise=args.theta_noise,
        theta_iid=args.theta_iid,
        windowsize=args.windowsize,
        chrom=args.chrom,
        output_dir=args.output_dir,
        seed=args.seed,
        hap_probs=hap_probs,
    )

    print(f"Generated {len(manifest)} manifest rows in {args.output_dir}/")
    print("Class distribution:")
    for cls in CLASSES:
        n_regions = (manifest.groupby("region_idx")["true_class"].first() == cls).sum()
        n_windows = (manifest["true_class"] == cls).sum()
        print(f"  {cls}: {n_regions} regions, {n_windows} windows")


if __name__ == "__main__":
    main()
