#!/usr/bin/env python3
"""
06_ac020656_validation.py
Phase 6: Dedicated robustness analysis for AC020656.1, the leading
donor-aware trajectory candidate -- per-donor breakdown, leave-one-donor-out
sensitivity, and monocyte-subtype specificity.

Usage:
  python scripts/analysis/06_ac020656_validation.py \
      --mono_h5ad data/processed/GSE268210_monocytes_final.h5ad \
      --outdir results/tables
"""
import argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
import anndata as ad

GENE = "AC020656.1"


def fisher_z_meta_hksj(rhos, ns):
    rhos = np.clip(np.asarray(rhos, dtype=float), -0.999, 0.999)
    ns = np.asarray(ns, dtype=float)
    ok = ns > 3
    if ok.sum() < 2:
        return np.nan, np.nan, int(ok.sum())
    z = np.arctanh(rhos[ok]); n = ns[ok]
    v = 1.0 / (n - 3); w_fe = 1.0 / v
    z_fe = np.sum(w_fe * z) / np.sum(w_fe)
    Q = np.sum(w_fe * (z - z_fe) ** 2)
    df = ok.sum() - 1
    C = np.sum(w_fe) - np.sum(w_fe ** 2) / np.sum(w_fe)
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    w_re = 1.0 / (v + tau2)
    z_bar = np.sum(w_re * z) / np.sum(w_re)
    k = ok.sum()
    if k > 2:
        resid_ss = np.sum(w_re * (z - z_bar) ** 2) / (k - 1)
        se_hksj = np.sqrt(resid_ss / np.sum(w_re))
        t_stat = z_bar / se_hksj if se_hksj > 0 else np.nan
        p = 2 * stats.t.sf(abs(t_stat), df=k - 1) if not np.isnan(t_stat) else np.nan
    else:
        se = 1.0 / np.sqrt(np.sum(w_re)); t_stat = z_bar / se
        p = 2 * stats.norm.sf(abs(t_stat))
    return np.tanh(z_bar), p, int(k)


def main():
    parser = argparse.ArgumentParser(description=f"{GENE} donor-level robustness analysis")
    parser.add_argument("--mono_h5ad", required=True)
    parser.add_argument("--outdir", default="results/tables")
    args = parser.parse_args()

    adata = ad.read_h5ad(args.mono_h5ad)
    expr = adata[:, GENE].layers["log_norm"]
    expr = np.asarray(expr.todense()).flatten() if sp.issparse(expr) else np.asarray(expr).flatten()
    counts = adata[:, GENE].layers["counts"]
    counts = np.asarray(counts.todense()).flatten() if sp.issparse(counts) else np.asarray(counts).flatten()

    pt = adata.obs["dpt_pseudotime"].values.astype(float)
    donors = adata.obs["sample_id"].astype(str).values
    subtype = adata.obs["mono_subtype"].astype(str).values
    uniq_donors = sorted(set(donors))

    rows, per_donor_rho, per_donor_n = [], [], []
    for d in uniq_donors:
        idx = donors == d
        e, t, c = expr[idx], pt[idx], counts[idx]
        rho, pval = stats.spearmanr(t, e)
        per_donor_rho.append(rho); per_donor_n.append(idx.sum())
        rows.append({"donor": d, "n_cells": int(idx.sum()), "detection_rate": float((c > 0).mean()),
                      "mean_log_norm_expr": float(e.mean()), "spearman_rho_vs_pseudotime": rho, "p_value": pval})
    per_donor_df = pd.DataFrame(rows)
    per_donor_df.to_csv(f"{args.outdir}/ac020656_per_donor_summary.csv", index=False)
    print(f"=== {GENE} per-donor summary ===\n{per_donor_df.to_string(index=False)}")

    overall_rho, overall_p, k = fisher_z_meta_hksj(per_donor_rho, per_donor_n)
    print(f"\nDonor-aware combined rho={overall_rho:.4f}  p={overall_p:.3e}  n_donors={k}")

    loo_rows = []
    for excl in uniq_donors:
        keep_rho = [r for d, r in zip(uniq_donors, per_donor_rho) if d != excl]
        keep_n = [n for d, n in zip(uniq_donors, per_donor_n) if d != excl]
        rho_loo, p_loo, k_loo = fisher_z_meta_hksj(keep_rho, keep_n)
        loo_rows.append({"excluded_donor": excl, "donor_aware_rho_without": rho_loo, "p_value_without": p_loo})
    loo_df = pd.DataFrame(loo_rows)
    loo_df.to_csv(f"{args.outdir}/ac020656_leave_one_donor_out.csv", index=False)
    rho_range = loo_df["donor_aware_rho_without"].max() - loo_df["donor_aware_rho_without"].min()
    print(f"\nLeave-one-donor-out: rho range={rho_range:.4f}  worst-case p={loo_df['p_value_without'].max():.3e}")

    sub_rows = []
    for st in ["CD14_Classical", "Intermediate", "CD16_NonClassical"]:
        idx = subtype == st
        sub_rows.append({"subtype": st, "n_cells": int(idx.sum()),
                          "detection_rate": float((counts[idx] > 0).mean()),
                          "mean_log_norm_expr": float(expr[idx].mean())})
    sub_df = pd.DataFrame(sub_rows)
    sub_df.to_csv(f"{args.outdir}/ac020656_subtype_breakdown.csv", index=False)
    print(f"\n=== {GENE} expression by monocyte subtype ===\n{sub_df.to_string(index=False)}")


if __name__ == "__main__":
    main()
