# ABOUTME: Plot the expected pi_{h,k,t} dynamics for the four canonical regimes
# ABOUTME: (null, static ASM, bilateral drift, allele-specific emergent ASM).

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

T = 4
t = np.arange(1, T + 1)

# pi[h, k, t-1]; columns must sum to 1 at each t.
# h in {0,1} = allele 1,2;  k in {0,1} = normal, altered;  k=1 forced 0 at t=1.

def make_pi(scenario):
    pi = np.zeros((2, 2, T))
    if scenario == "A_null":
        pi[0, 0, :] = 0.5
        pi[1, 0, :] = 0.5
    elif scenario == "B_static_asm":
        pi[0, 0, :] = 0.5
        pi[1, 0, :] = 0.5
    elif scenario == "C_bilateral_drift":
        # both alleles drift symmetrically into altered
        altered_share = np.array([0.0, 0.30, 0.55, 0.75])
        pi[0, 1, :] = altered_share / 2
        pi[1, 1, :] = altered_share / 2
        pi[0, 0, :] = (1 - altered_share) / 2
        pi[1, 0, :] = (1 - altered_share) / 2
    elif scenario == "D_allele_specific":
        # only h=2 (allele 2) drifts into altered
        h2_altered = np.array([0.0, 0.20, 0.40, 0.60])
        pi[0, 0, :] = 0.5                  # h=1 stays normal
        pi[1, 0, :] = 0.5 - h2_altered     # h=2 normal share shrinks
        pi[1, 1, :] = h2_altered           # h=2 altered share grows
        pi[0, 1, :] = 0.0
    return pi

def make_theta(scenario, J=10):
    rng = np.random.default_rng(0)
    base_low = rng.uniform(0.05, 0.20, J)
    base_high = rng.uniform(0.80, 0.95, J)
    theta = np.zeros((J, 2, 2))
    if scenario == "A_null":
        shared = rng.uniform(0.4, 0.6, J)
        for h in range(2):
            for k in range(2):
                theta[:, h, k] = shared
    elif scenario == "B_static_asm":
        # alleles differ; altered profile collapses onto normal (no real k=2)
        theta[:, 0, 0] = base_low
        theta[:, 1, 0] = base_high
        theta[:, 0, 1] = base_low
        theta[:, 1, 1] = base_high
    elif scenario == "C_bilateral_drift":
        # both alleles share normal=low and altered=high
        theta[:, 0, 0] = base_low
        theta[:, 1, 0] = base_low
        theta[:, 0, 1] = base_high
        theta[:, 1, 1] = base_high
    elif scenario == "D_allele_specific":
        # h=1 normal=low; h=2 normal=low, altered=high
        theta[:, 0, 0] = base_low
        theta[:, 1, 0] = base_low
        theta[:, 0, 1] = base_low      # cell carries no real mass
        theta[:, 1, 1] = base_high
    return theta

scenarios = [
    ("A_null",            "Class A · pure null (no ASM, no drift)"),
    ("B_static_asm",      "Class B · static ASM (no temporal change)"),
    ("C_bilateral_drift", "Class C · bilateral temporal drift (no ASM)"),
    ("D_allele_specific", "Class D · allele-specific emergent ASM"),
]

colors = {
    (0, 0): "#1f77b4",  # h=1 normal  - blue
    (0, 1): "#aec7e8",  # h=1 altered - light blue
    (1, 0): "#d62728",  # h=2 normal  - red
    (1, 1): "#ff9896",  # h=2 altered - light red
}
labels = {
    (0, 0): r"$\pi_{1,1,t}$  allele 1, normal",
    (0, 1): r"$\pi_{1,2,t}$  allele 1, altered",
    (1, 0): r"$\pi_{2,1,t}$  allele 2, normal",
    (1, 1): r"$\pi_{2,2,t}$  allele 2, altered",
}

# ----------------------------------------------------------------------
# Figure 1: pi trajectories
# ----------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
for ax, (key, title) in zip(axes.flat, scenarios):
    pi = make_pi(key)
    for h in range(2):
        for k in range(2):
            linestyle = "-" if k == 0 else "--"
            ax.plot(
                t, pi[h, k, :],
                marker="o", linewidth=2.2, linestyle=linestyle,
                color=colors[(h, k)], label=labels[(h, k)],
            )
    ax.set_title(title, fontsize=11)
    ax.set_xticks(t)
    ax.set_xlabel("time point t")
    ax.set_ylabel(r"$\pi_{h,k,t}$")
    ax.set_ylim(-0.03, 1.03)
    ax.axhline(0, color="grey", linewidth=0.5)
    ax.grid(alpha=0.3)

handles, lab = axes[0, 0].get_legend_handles_labels()
fig.legend(
    handles, lab,
    loc="lower center", ncol=4, frameon=False,
    bbox_to_anchor=(0.5, -0.02), fontsize=10,
)
fig.suptitle(r"Expected mixing-weight dynamics $\pi_{h,k,t}$ per regime", fontsize=13)
fig.tight_layout(rect=(0, 0.03, 1, 0.96))
fig.savefig("model_dynamics_pi.png", dpi=150, bbox_inches="tight")

# ----------------------------------------------------------------------
# Figure 2: theta profiles per regime (CpG site x methylation prob)
# ----------------------------------------------------------------------
fig2, axes2 = plt.subplots(2, 2, figsize=(12, 7), sharex=True, sharey=True)
J = 10
sites = np.arange(1, J + 1)
for ax, (key, title) in zip(axes2.flat, scenarios):
    theta = make_theta(key, J=J)
    for h in range(2):
        for k in range(2):
            linestyle = "-" if k == 0 else "--"
            ax.plot(
                sites, theta[:, h, k],
                marker="s", linewidth=1.8, linestyle=linestyle,
                color=colors[(h, k)], label=labels[(h, k)],
            )
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("CpG site index j (within window)")
    ax.set_ylabel(r"$\theta_{j,h,k}$")
    ax.set_ylim(-0.03, 1.03)
    ax.set_xticks(sites)
    ax.grid(alpha=0.3)

handles2, lab2 = axes2[0, 0].get_legend_handles_labels()
fig2.legend(
    handles2, lab2,
    loc="lower center", ncol=4, frameon=False,
    bbox_to_anchor=(0.5, -0.02), fontsize=10,
)
fig2.suptitle(r"Expected methylation profiles $\theta_{j,h,k}$ per regime", fontsize=13)
fig2.tight_layout(rect=(0, 0.03, 1, 0.96))
fig2.savefig("model_dynamics_theta.png", dpi=150, bbox_inches="tight")

# ----------------------------------------------------------------------
# Figure 3: stacked-bar view of pi at each t (easier to read the simplex)
# ----------------------------------------------------------------------
fig3, axes3 = plt.subplots(1, 4, figsize=(15, 4.2), sharey=True)
order = [(0, 0), (1, 0), (0, 1), (1, 1)]   # normals first, altered on top
for ax, (key, title) in zip(axes3.flat, scenarios):
    pi = make_pi(key)
    bottom = np.zeros(T)
    for (h, k) in order:
        ax.bar(
            t, pi[h, k, :], bottom=bottom,
            color=colors[(h, k)], edgecolor="white", linewidth=0.7,
            label=labels[(h, k)],
        )
        bottom = bottom + pi[h, k, :]
    ax.set_title(title.split("·")[0].strip() + "\n" + title.split("·")[1].strip(), fontsize=10)
    ax.set_xticks(t)
    ax.set_xlabel("t")
    ax.set_ylim(0, 1.0)
    ax.grid(axis="y", alpha=0.3)
axes3[0].set_ylabel("mixing weight (stacked, sums to 1)")
handles3, lab3 = axes3[0].get_legend_handles_labels()
fig3.legend(
    handles3, lab3,
    loc="lower center", ncol=4, frameon=False,
    bbox_to_anchor=(0.5, -0.05), fontsize=10,
)
fig3.suptitle(r"$\pi_{h,k,t}$ as a stacked simplex per time point", fontsize=12)
fig3.tight_layout(rect=(0, 0.05, 1, 0.95))
fig3.savefig("model_dynamics_stacked.png", dpi=150, bbox_inches="tight")

print("wrote: model_dynamics_pi.png, model_dynamics_theta.png, model_dynamics_stacked.png")
