# ABOUTME: Renders sweep CSVs into Markdown tables for the parameter-recovery report appendix.
# ABOUTME: Run after parameter_recovery_sweep.py completes; pipes table fragments into the report.
"""Helper: turn sweep CSVs into 3 Markdown tables for the report appendix.

Outputs to stdout. Patch into `docs/PARAMETER_RECOVERY_REPORT_2026-04-27.md`
by replacing the `<!-- AUTO-FILL: ... -->` markers.

Usage:
    cd em
    python _render_report_tables.py \
        --recovery sweep_artifacts/sweep_recovery.csv \
        --bic       sweep_artifacts/sweep_bic.csv
"""

import argparse

import pandas as pd

CLASS_SHORT = {
    "NULL": "NULL",
    "STATIC_ASM": "STATIC",
    "TIME_EMERGENT_ASM": "TE_sym",
    "ASYMMETRIC_TIME_EMERGENT_ASM": "ASYM",
}


def df_to_md(df, fmt="{:.3f}"):
    """Render a DataFrame as a GitHub-flavored Markdown table."""
    cols = list(df.columns)
    header = "| " + " | ".join(str(c) for c in cols) + " |"
    sep = "|" + "|".join("---" for _ in cols) + "|"
    rows = []
    for _, row in df.iterrows():
        cells = []
        for c in cols:
            v = row[c]
            if isinstance(v, float):
                cells.append(fmt.format(v))
            else:
                cells.append(str(v))
        rows.append("| " + " | ".join(cells) + " |")
    return "\n".join([header, sep, *rows])


def render_theta_table(rec):
    rec = rec.copy()
    rec["class"] = rec["true_class"].map(CLASS_SHORT)
    rec["M1_theta_avg_mae"] = (rec["M1_theta_h0_mae"] + rec["M1_theta_h1_mae"]) / 2
    rec["M2_theta_k0_avg_mae"] = (rec["M2_theta_h0_k0_mae"] + rec["M2_theta_h1_k0_mae"]) / 2
    out = (
        rec[["class", "pi_final", "reads_per_t",
             "M0_theta_mae", "M1_theta_avg_mae", "M2_theta_k0_avg_mae"]]
        .sort_values(["class", "pi_final", "reads_per_t"])
        .reset_index(drop=True)
    )
    return df_to_md(out)


def render_pi_table(rec):
    rec = rec.copy()
    rec["class"] = rec["true_class"].map(CLASS_SHORT)
    out = (
        rec[["class", "pi_final", "reads_per_t",
             "M2_pi_h0_alt_t2_abs_err", "M2_pi_h1_alt_t2_abs_err",
             "M2_pi_max_abs_err"]]
        .sort_values(["class", "pi_final", "reads_per_t"])
        .reset_index(drop=True)
    )
    return df_to_md(out)


def render_bic_table(bic):
    bic = bic.copy()
    bic["class"] = bic["true_class"].map(CLASS_SHORT)
    out = (
        bic[["class", "pi_final", "reads_per_t",
             "pct_M0", "pct_M1", "pct_M2"]]
        .sort_values(["class", "pi_final", "reads_per_t"])
        .reset_index(drop=True)
    )
    return df_to_md(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recovery", required=True)
    ap.add_argument("--bic", required=True)
    args = ap.parse_args()

    # keep_default_na=False so the literal class name "NULL" isn't parsed as NaN
    rec = pd.read_csv(args.recovery, keep_default_na=False)
    bic = pd.read_csv(args.bic, keep_default_na=False)
    for df in (rec, bic):
        for c in df.select_dtypes(include="object").columns:
            if c == "true_class":
                continue
            df[c] = pd.to_numeric(df[c])

    print("\n=== AUTO-FILL: theta_recovery_table ===\n")
    print(render_theta_table(rec))
    print("\n=== AUTO-FILL: pi_recovery_table ===\n")
    print(render_pi_table(rec))
    print("\n=== AUTO-FILL: bic_winner_table ===\n")
    print(render_bic_table(bic))


if __name__ == "__main__":
    main()
