#!/usr/bin/env python3
"""
ZERO_DOWNLOAD Phase 3 -- cell-level AC020656.1 vs LYZ independence assessment,
using the ALREADY-EXISTING monocyte h5ad (no new download, no full-matrix load
-- column-sliced to just the 2 genes + metadata needed).

Outputs small CSV/text summaries only (no large files).
"""
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
import anndata as ad

H5AD = "data/processed/GSE268210_monocytes_final.h5ad"
OUTDIR = "results/tables/ac0206561_lyz"

GENES = ["AC020656.1", "LYZ"]


def main():
    print("Loading monocyte h5ad in BACKED mode (no full matrix load)...", file=sys.stderr)
    adata = ad.read_h5ad(H5AD, backed="r")
    print(f"  {adata.shape}", file=sys.stderr)

    # column-slice only the 2 genes we need, then bring just that slice into memory
    sub = adata[:, GENES].to_memory()
    counts = sub.layers["counts"]
    counts = counts.toarray() if sp.issparse(counts) else np.asarray(counts)
    lognorm = sub.layers["log_norm"]
    lognorm = lognorm.toarray() if sp.issparse(lognorm) else np.asarray(lognorm)

    ac_counts, lyz_counts = counts[:, 0], counts[:, 1]
    ac_ln, lyz_ln = lognorm[:, 0], lognorm[:, 1]
    subtype = adata.obs["mono_subtype"].astype(str).values
    donor = adata.obs["sample_id"].astype(str).values
    n_cells = len(ac_counts)

    ## ---- A. Detection ----
    rows = []
    for name, c in [("AC020656.1", ac_counts), ("LYZ", lyz_counts)]:
        detected = c > 0
        rows.append({
            "gene": name, "n_cells_total": n_cells,
            "n_detected": int(detected.sum()), "fraction_detected": detected.mean(),
            "mean_count_all_cells": c.mean(), "mean_count_detected_cells": c[detected].mean() if detected.sum() else np.nan,
            "median_count_detected_cells": np.median(c[detected]) if detected.sum() else np.nan,
        })
    detection_df = pd.DataFrame(rows)
    detection_df.to_csv(f"{OUTDIR}/03A_detection_summary.csv", index=False)
    print("=== A. Detection ==="); print(detection_df.to_string(index=False))

    ## ---- B. Cell-type distribution ----
    rows = []
    for st in sorted(set(subtype)):
        mask = subtype == st
        rows.append({
            "subtype": st, "n_cells": int(mask.sum()),
            "AC020656.1_detection_rate": (ac_counts[mask] > 0).mean(),
            "AC020656.1_mean_count": ac_counts[mask].mean(),
            "LYZ_detection_rate": (lyz_counts[mask] > 0).mean(),
            "LYZ_mean_count": lyz_counts[mask].mean(),
        })
    subtype_df = pd.DataFrame(rows)
    subtype_df.to_csv(f"{OUTDIR}/03B_subtype_distribution.csv", index=False)
    print("\n=== B. Cell-type distribution ==="); print(subtype_df.to_string(index=False))

    ## ---- C. Co-detection 2x2 table ----
    ac_pos, lyz_pos = ac_counts > 0, lyz_counts > 0
    table = pd.DataFrame({
        "category": ["AC+/LYZ+", "AC+/LYZ-", "AC-/LYZ+", "AC-/LYZ-"],
        "n_cells": [
            int((ac_pos & lyz_pos).sum()), int((ac_pos & ~lyz_pos).sum()),
            int((~ac_pos & lyz_pos).sum()), int((~ac_pos & ~lyz_pos).sum()),
        ],
    })
    table["fraction"] = table["n_cells"] / n_cells
    table.to_csv(f"{OUTDIR}/03C_codetection_table.csv", index=False)
    print("\n=== C. Co-detection ==="); print(table.to_string(index=False))
    # chi-square test of independence (descriptive only -- cells are not independent
    # biological replicates, so this is NOT used for any disease-level inference)
    contingency = np.array([[table.n_cells[0], table.n_cells[1]], [table.n_cells[2], table.n_cells[3]]])
    chi2, chi2_p, _, _ = stats.chi2_contingency(contingency)
    print(f"Chi-square test of AC+/LYZ+ association (descriptive, cell-level, NOT a biological-replicate test): "
          f"chi2={chi2:.1f}, p={chi2_p:.3e}")

    ## ---- D. Within-cell relationship ----
    rows = []
    # overall
    rho, p = stats.spearmanr(ac_ln, lyz_ln)
    rows.append({"stratum": "ALL_CELLS", "n": n_cells, "spearman_rho": rho, "spearman_p": p})
    # within subtype
    for st in sorted(set(subtype)):
        mask = subtype == st
        rho, p = stats.spearmanr(ac_ln[mask], lyz_ln[mask])
        rows.append({"stratum": f"subtype={st}", "n": int(mask.sum()), "spearman_rho": rho, "spearman_p": p})
    # among cells where at least one gene is detected
    either = ac_pos | lyz_pos
    rho, p = stats.spearmanr(ac_ln[either], lyz_ln[either])
    rows.append({"stratum": "AC+_OR_LYZ+_cells", "n": int(either.sum()), "spearman_rho": rho, "spearman_p": p})
    # among double-positive cells only
    both = ac_pos & lyz_pos
    if both.sum() > 10:
        rho, p = stats.spearmanr(ac_ln[both], lyz_ln[both])
        rows.append({"stratum": "AC+_AND_LYZ+_cells", "n": int(both.sum()), "spearman_rho": rho, "spearman_p": p})
    within_df = pd.DataFrame(rows)
    within_df.to_csv(f"{OUTDIR}/03D_within_cell_correlation.csv", index=False)
    print("\n=== D. Within-cell correlation ==="); print(within_df.to_string(index=False))

    ## ---- E. Conditional structure (descriptive residual variance, cell-level) ----
    # does AC020656.1 retain variation after regressing out LYZ + subtype + donor
    # (dummy-coded)? Purely descriptive -- explicitly NOT a causal/independence proof.
    df = pd.DataFrame({"ac": ac_ln, "lyz": lyz_ln, "subtype": subtype, "donor": donor})
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    m0 = smf.ols("ac ~ lyz", data=df).fit()
    m1 = smf.ols("ac ~ lyz + C(subtype)", data=df).fit()
    m2 = smf.ols("ac ~ lyz + C(subtype) + C(donor)", data=df).fit()
    resid_var_summary = pd.DataFrame({
        "model": ["ac ~ lyz", "ac ~ lyz + subtype", "ac ~ lyz + subtype + donor"],
        "r_squared": [m0.rsquared, m1.rsquared, m2.rsquared],
        "residual_variance": [m0.resid.var(), m1.resid.var(), m2.resid.var()],
        "residual_variance_pct_of_total": [m0.resid.var() / df["ac"].var() * 100,
                                            m1.resid.var() / df["ac"].var() * 100,
                                            m2.resid.var() / df["ac"].var() * 100],
    })
    resid_var_summary.to_csv(f"{OUTDIR}/03E_conditional_residual_variance.csv", index=False)
    print("\n=== E. Conditional structure (cell-level, descriptive only) ===")
    print(resid_var_summary.to_string(index=False))
    print(f"\nTotal AC020656.1 variance explained by LYZ+subtype+donor: {100 - resid_var_summary.residual_variance_pct_of_total.iloc[-1]:.1f}%")
    print(f"Remaining (unexplained) variance: {resid_var_summary.residual_variance_pct_of_total.iloc[-1]:.1f}%")
    print("NOTE: remaining variance is NOT proof of independent transcription -- it includes technical")
    print("noise, dropout, and any unmodeled biological or technical structure. This is descriptive only.")

    print(f"\nDone. Outputs in {OUTDIR}/ (03A-03E files)")


if __name__ == "__main__":
    main()
