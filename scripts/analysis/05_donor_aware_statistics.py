#!/usr/bin/env python3
"""
04_donor_aware_statistics.py
Phase 4: Donor-aware trajectory-lncRNA association testing.

Cells from the same donor are not independent observations: they share
genetic background, sample handling, and technical batch. For each candidate
lncRNA (GENCODE-confirmed, >=2% detection in monocytes), Spearman correlation
with diffusion pseudotime is computed SEPARATELY within each donor, then the
per-donor estimates are combined via a DerSimonian-Laird random-effects
meta-analysis on the Fisher z scale, with a Hartung-Knapp-Sidik-Jonkman (HKSJ)
small-sample correction applied to the final test (t-distribution, n_donors-1
degrees of freedom). This correctly reflects an effective sample size of
n_donors, not n_cells, and is the statistic used for every trajectory-lncRNA
significance claim in this project.

A locus is called significant at combined |rho| >= 0.15, Benjamini-Hochberg
FDR < 0.05 (across all tested loci), and direction agreement in at least
6 of 9 donors.

Usage:
  python scripts/analysis/04_donor_aware_statistics.py \
      --mono_h5ad data/processed/GSE268210_monocytes_final.h5ad \
      --lncrna_table results/tables/lncrna_annotation_gencode_v32.csv \
      --outdir results/tables
"""
import argparse
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import anndata as ad

MIN_DETECTION = 0.02
MIN_CELLS_PER_DONOR = 20
SEED = 42
np.random.seed(SEED)


def fisher_z_meta(rhos, ns):
    """DerSimonian-Laird random-effects Fisher z meta-analysis with HKSJ correction."""
    rhos = np.clip(np.asarray(rhos, dtype=float), -0.999, 0.999)
    ns = np.asarray(ns, dtype=float)
    ok = ns > 3
    if ok.sum() < 2:
        return dict(combined_rho=np.nan, z_stat=np.nan, p_value=np.nan,
                     n_donors=int(ok.sum()), Q=np.nan, Q_p=np.nan, I2=np.nan, tau2=np.nan)

    z = np.arctanh(rhos[ok])
    n = ns[ok]
    v = 1.0 / (n - 3)
    w_fe = 1.0 / v

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
        z_stat = t_stat
    else:
        se = 1.0 / np.sqrt(np.sum(w_re))
        z_stat = z_bar / se
        p = 2 * stats.norm.sf(abs(z_stat))

    Q_p = stats.chi2.sf(Q, df) if df > 0 else np.nan
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    return dict(combined_rho=np.tanh(z_bar), z_stat=z_stat, p_value=p,
                n_donors=int(k), Q=Q, Q_p=Q_p, I2=I2, tau2=tau2)


def main():
    parser = argparse.ArgumentParser(description="Donor-aware trajectory-lncRNA statistics")
    parser.add_argument("--mono_h5ad", required=True)
    parser.add_argument("--lncrna_table", required=True, help="lncrna_annotation_gencode_v32.csv")
    parser.add_argument("--outdir", default="results/tables")
    args = parser.parse_args()

    print("Loading GENCODE-confirmed lncRNA list...", file=sys.stderr)
    gencode_lnc = pd.read_csv(args.lncrna_table)
    gencode_lnc_genes = set(gencode_lnc["gene_symbol"])

    print("Loading monocyte h5ad (log_norm layer)...", file=sys.stderr)
    adata = ad.read_h5ad(args.mono_h5ad)
    print(f"  {adata.shape[0]:,} cells x {adata.shape[1]:,} genes", file=sys.stderr)

    candidate_genes = [g for g in gencode_lnc_genes if g in adata.var_names]
    sub = adata[:, candidate_genes]
    X = sub.layers["log_norm"]
    det_frac = np.asarray((X > 0).mean(axis=0)).flatten() if sp.issparse(X) else (X > 0).mean(axis=0)
    usable_mask = det_frac >= MIN_DETECTION
    usable_genes = [g for g, u in zip(candidate_genes, usable_mask) if u]
    print(f"  Usable (>{MIN_DETECTION*100:.0f}% detection) lncRNAs: {len(usable_genes):,}", file=sys.stderr)

    pt = adata.obs["dpt_pseudotime"].values.astype(float)
    donors = adata.obs["sample_id"].astype(str).values
    uniq_donors = sorted(set(donors))

    Xu = adata[:, usable_genes].layers["log_norm"]
    Xu = Xu.toarray() if sp.issparse(Xu) else np.asarray(Xu)

    donor_idx = {d: np.where(donors == d)[0] for d in uniq_donors}

    results = []
    for j, gene in enumerate(usable_genes):
        expr = Xu[:, j]
        per_donor_rho, per_donor_n = [], []
        for d in uniq_donors:
            idx = donor_idx[d]
            if len(idx) < MIN_CELLS_PER_DONOR:
                continue
            e, t = expr[idx], pt[idx]
            if e.std() < 1e-9 or t.std() < 1e-9:
                continue
            r, _ = stats.spearmanr(t, e)
            if np.isnan(r):
                continue
            per_donor_rho.append(r)
            per_donor_n.append(len(idx))

        meta = fisher_z_meta(per_donor_rho, per_donor_n)
        dirs = np.sign(per_donor_rho) if len(per_donor_rho) else np.array([])
        combined_sign = np.sign(meta["combined_rho"]) if not np.isnan(meta.get("combined_rho", np.nan)) else 0
        dir_consistency = float(np.mean(dirs == combined_sign)) if len(dirs) and combined_sign != 0 else np.nan

        results.append({
            "gene": gene, "donor_aware_rho": meta["combined_rho"], "donor_aware_pvalue": meta["p_value"],
            "n_donors_contributing": meta["n_donors"], "direction_consistency": dir_consistency,
            "heterogeneity_I2_pct": meta["I2"], "detection_frac": det_frac[candidate_genes.index(gene)],
        })
        if (j + 1) % 100 == 0:
            print(f"  ... {j+1}/{len(usable_genes)} genes processed", file=sys.stderr)

    res = pd.DataFrame(results).dropna(subset=["donor_aware_pvalue"])
    _, donor_fdr, _, _ = multipletests(res["donor_aware_pvalue"], method="fdr_bh")
    res["donor_aware_fdr"] = donor_fdr
    res["significant"] = (
        (res["donor_aware_fdr"] < 0.05) & (res["donor_aware_rho"].abs() >= 0.15) & (res["direction_consistency"] >= 0.66)
    )
    res = res.sort_values("donor_aware_rho", key=lambda s: s.abs(), ascending=False)
    res.to_csv(f"{args.outdir}/trajectory_donor_aware_results.csv", index=False)

    n_sig = int(res["significant"].sum())
    print(f"\n{n_sig} donor-aware-significant loci (of {len(res)} tested)")
    print(res[res["significant"]][["gene", "donor_aware_rho", "donor_aware_fdr",
                                     "direction_consistency"]].to_string(index=False))
    print(f"\nOutput: {args.outdir}/trajectory_donor_aware_results.csv")


if __name__ == "__main__":
    main()
