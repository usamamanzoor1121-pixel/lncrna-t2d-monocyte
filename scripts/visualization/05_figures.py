#!/usr/bin/env python3
"""
05_figures.py
Regenerate the FINAL publication figures from already-computed results
(no recomputation of pseudotime, trajectory statistics, or DE).

REPAIR NOTE (2026-09-22): this file previously contained a byte-identical
copy of the Phase-2 QC/annotation script (mismatched name and a CLI that did
not match this file's own README documentation). That content has been moved
to `scripts/preprocessing/02_cell_annotation.py`, where it belongs. This file
now actually does what the README's Quick Start Step 4 says it does. See
`REPAIR/REPORTS/01_GAP_RESOLUTION_MATRIX.md` (finding H1).

Regenerates:
  Fig9_pseudotime_final.png            — monocyte subtype + pseudotime UMAPs
  Fig10_trajectory_heatmap_final.png   — trajectory lncRNA expression heatmap
  Fig10b_trajectory_lollipop_final.png — lncRNA-pseudotime correlation lollipop
  Fig13_lncrna_celltype_specificity.png— lncRNA atlas across PBMC subsets
  Fig14_bulk_validation.png            — scRNA vs. bulk concordance (panels A/B only —
                                          panel C per-sample boxplots need the raw bulk
                                          count matrix, which is not part of this script's
                                          4 documented inputs; see note printed at runtime)

Usage:
  python scripts/visualization/05_figures.py \\
      --full_h5ad data/processed/GSE268210_phase1_full.h5ad \\
      --mono_h5ad results/phase3c/GSE268210_monocytes_final.h5ad \\
      --traj_csv  results/phase3c/trajectory_lncrna_final.csv \\
      --bulk_csv  results/phase4/bulk_concordance.csv \\
      --outdir    results/figures
"""
import argparse
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import to_rgba
from pathlib import Path

warnings.filterwarnings('ignore')
sc.settings.verbosity = 0

SUBTYPE_ORDER = ['CD14_Classical', 'Intermediate', 'CD16_NonClassical']
SUBTYPE_LABELS = ['CD14\nClassical', 'Intermediate', 'CD16\nNon-Classical']
SUBTYPE_COLS = {'CD14_Classical': '#E53E3E', 'Intermediate': '#D69E2E', 'CD16_NonClassical': '#3182CE'}


def savefig(fig, path, dpi=150):
    fig.savefig(path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: {Path(path).name}")


def fig9_pseudotime(adata_mono, outdir):
    umap = adata_mono.obsm['X_umap']
    pt = adata_mono.obs['dpt_pseudotime'].values
    counts = {st: (adata_mono.obs['mono_subtype'] == st).sum() for st in SUBTYPE_ORDER}
    total = sum(counts.values())

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'Monocyte Pseudotime Trajectory — T2D PBMCs\nn={total:,} cells', fontsize=11, fontweight='bold')

    ax = axes[0]
    for st in SUBTYPE_ORDER:
        mask = adata_mono.obs['mono_subtype'] == st
        ax.scatter(umap[mask, 0], umap[mask, 1], c=SUBTYPE_COLS[st], s=0.5, alpha=0.5, rasterized=True,
                   label=f"{st.replace('_', ' ')} (n={mask.sum():,})")
    ax.set_title('Monocyte Subtypes', fontsize=11, fontweight='bold')
    ax.legend(markerscale=6, fontsize=8, loc='best', framealpha=0.8)
    ax.spines[['top', 'right']].set_visible(False)

    ax = axes[1]
    sc_plot = ax.scatter(umap[:, 0], umap[:, 1], c=pt, cmap='viridis', s=0.5, alpha=0.6, rasterized=True)
    plt.colorbar(sc_plot, ax=ax, shrink=0.8, label='Pseudotime')
    ax.set_title('Diffusion Pseudotime', fontsize=11, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)

    ax = axes[2]
    for i, st in enumerate(SUBTYPE_ORDER):
        mask = adata_mono.obs['mono_subtype'] == st
        if mask.sum() > 0:
            vp = ax.violinplot([pt[mask.values]], positions=[i], showmedians=True, showextrema=True)
            for pc in vp['bodies']:
                pc.set_facecolor(SUBTYPE_COLS[st]); pc.set_alpha(0.75)
            med = np.median(pt[mask.values])
            ax.text(i, med + 0.015, f'{med:.3f}', ha='center', fontsize=9, fontweight='bold')
    ax.set_xticks(range(3)); ax.set_xticklabels(SUBTYPE_LABELS, fontsize=9)
    ax.set_ylabel('Pseudotime', fontsize=10)
    ax.set_title('Pseudotime per Subtype', fontsize=11, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)

    plt.tight_layout()
    savefig(fig, outdir / 'Fig9_pseudotime_final.png', dpi=200)


def fig10_heatmap_and_lollipop(adata_mono, traj_df, outdir):
    sig = traj_df[(traj_df['padj'] < 0.05) & (traj_df['abs_rho'] >= 0.15)].copy()
    genes = [g for g in sig.sort_values('abs_rho', ascending=False).head(30)['gene'] if g in adata_mono.var_names]
    if len(genes) < 3:
        print("  Too few significant genes for heatmap — skipping Fig10")
        return

    pseudotime = adata_mono.obs['dpt_pseudotime'].values
    pt_order = np.argsort(pseudotime)
    n_show = min(5000, len(pt_order))
    step = max(1, len(pt_order) // n_show)
    show_idx = pt_order[::step][:n_show]

    X = adata_mono[show_idx, :][:, genes].layers['log_norm']
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    X = X.T.astype(np.float32)
    win = max(5, n_show // 50)
    X_sm = pd.DataFrame(X).T.rolling(win, center=True, min_periods=1).mean().T.values
    rmean = X_sm.mean(axis=1, keepdims=True); rstd = X_sm.std(axis=1, keepdims=True)
    rstd[rstd < 1e-6] = 1.0
    X_z = np.clip((X_sm - rmean) / rstd, -2.5, 2.5)

    n_genes = len(genes)
    fig = plt.figure(figsize=(14, max(7, n_genes * 0.32 + 3)))
    gs = gridspec.GridSpec(3, 2, height_ratios=[0.04, 0.04, 1], width_ratios=[3, 1], hspace=0.05, wspace=0.08)

    ax_st = fig.add_subplot(gs[0, 0])
    st_vals = adata_mono.obs['mono_subtype'].values[show_idx]
    rgba = np.array([[to_rgba(SUBTYPE_COLS.get(str(s), '#9E9E9E')) for s in st_vals]])
    ax_st.imshow(rgba, aspect='auto', interpolation='nearest')
    ax_st.set_xticks([]); ax_st.set_yticks([])

    ax_pt = fig.add_subplot(gs[1, 0])
    ax_pt.imshow(pseudotime[show_idx].reshape(1, -1), aspect='auto', cmap='viridis', interpolation='nearest')
    ax_pt.set_xticks([]); ax_pt.set_yticks([])

    ax_hm = fig.add_subplot(gs[2, 0])
    im = ax_hm.imshow(X_z, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2, interpolation='nearest')
    ax_hm.set_yticks(range(n_genes)); ax_hm.set_yticklabels(genes, fontsize=8)
    ax_hm.set_xlabel('Cells ordered by pseudotime', fontsize=9); ax_hm.set_xticks([])
    plt.colorbar(im, ax=ax_hm, shrink=0.4, label='Z-score', pad=0.01)

    ax_rho = fig.add_subplot(gs[2, 1])
    rl = sig.set_index('gene')['rho'].to_dict()
    rvals = [rl.get(g, 0) for g in genes]
    ax_rho.barh(range(n_genes), rvals, color=['#E53E3E' if r > 0 else '#3182CE' for r in rvals], alpha=0.85)
    ax_rho.set_yticks(range(n_genes)); ax_rho.set_yticklabels([])
    ax_rho.axvline(x=0, color='black', linewidth=1.0)
    ax_rho.set_xlabel("Spearman rho", fontsize=8)
    ax_rho.spines[['top', 'right']].set_visible(False)

    fig.suptitle(f'Trajectory-Associated lncRNAs in T2D Monocytes (n={n_genes})', fontsize=11, fontweight='bold', y=1.01)
    savefig(fig, outdir / 'Fig10_trajectory_heatmap_final.png', dpi=180)

    df = sig.copy().sort_values('rho')
    fig, ax = plt.subplots(figsize=(8, max(6, len(df) * 0.35 + 2)))
    for i, (_, row) in enumerate(df.iterrows()):
        col = '#3182CE' if row['rho'] < 0 else '#E53E3E'
        ax.plot([0, row['rho']], [i, i], color=col, lw=2, alpha=0.75)
        ax.scatter(row['rho'], i, color=col, s=60, zorder=3)
        ha = 'left' if row['rho'] > 0 else 'right'
        ax.text(row['rho'] + (0.005 if row['rho'] > 0 else -0.005), i, row['gene'], va='center', ha=ha, fontsize=8.5)
    ax.axvline(x=0, color='black', lw=1.0)
    ax.set_yticks([])
    ax.set_xlabel('Spearman rho with pseudotime', fontsize=11)
    ax.set_title(f'Trajectory-Associated lncRNAs — T2D Monocytes (n={len(df)})', fontsize=11, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    savefig(fig, outdir / 'Fig10b_trajectory_lollipop_final.png')


def fig13_celltype_specificity(adata_full, outdir):
    lnc_mask = adata_full.var.get('is_lncrna', pd.Series(False, index=adata_full.var_names)).astype(bool)
    lncrna_genes = adata_full.var_names[lnc_mask].tolist()
    X_lnc = adata_full[:, lncrna_genes].layers.get('log_norm', adata_full[:, lncrna_genes].X)
    det = np.asarray((X_lnc > 0).mean(axis=0)).flatten() if sp.issparse(X_lnc) else (X_lnc > 0).mean(axis=0)
    usable = [g for g, d in zip(lncrna_genes, det) if d >= 0.01]

    celltypes = ['T_cell', 'NK', 'Monocyte', 'B_cell', 'DC', 'Platelet']
    mean_expr = {}
    for ct in celltypes:
        mask = adata_full.obs['broad_celltype'] == ct
        if mask.sum() < 50:
            continue
        X_ct = adata_full[mask, :][:, usable].layers.get('log_norm', adata_full[mask, :][:, usable].X)
        mean_expr[ct] = np.asarray(X_ct.mean(axis=0)).flatten() if sp.issparse(X_ct) else X_ct.mean(axis=0)

    mean_df = pd.DataFrame(mean_expr, index=usable)
    var_genes = mean_df.var(axis=1).nlargest(40).index.tolist()
    plot_data = mean_df.loc[var_genes]
    row_mean = plot_data.mean(axis=1); row_std = plot_data.std(axis=1).replace(0, 1)
    plot_z = plot_data.sub(row_mean, axis=0).div(row_std, axis=0)

    fig, ax = plt.subplots(figsize=(10, 13))
    im = ax.imshow(plot_z.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2, interpolation='nearest')
    ax.set_xticks(range(len(mean_df.columns))); ax.set_xticklabels(mean_df.columns, fontsize=11, fontweight='bold')
    ax.set_yticks(range(len(var_genes))); ax.set_yticklabels(var_genes, fontsize=8)
    for x in np.arange(-0.5, len(mean_df.columns), 1):
        ax.axvline(x, color='white', lw=1.5)
    plt.colorbar(im, ax=ax, shrink=0.35, pad=0.02, label='Z-score (across cell types)')
    ax.set_title(f'Cell-Type-Specific lncRNA Expression in T2D PBMCs\n(top 40 of {len(usable)} lncRNAs)',
                 fontsize=12, fontweight='bold', pad=12)
    plt.tight_layout()
    savefig(fig, outdir / 'Fig13_lncrna_celltype_specificity.png')


def fig14_bulk_validation(conc_df, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    conc = conc_df[conc_df['concordant']]
    disc = conc_df[~conc_df['concordant']]

    ax = axes[0]
    ax.scatter(conc['scrna_rho'], conc['bulk_lfc'], c='#2ECC71', s=60, alpha=0.8, label='Concordant')
    ax.scatter(disc['scrna_rho'], disc['bulk_lfc'], c='#E74C3C', s=60, alpha=0.8, label='Discordant')
    ax.axhline(0, color='gray', lw=0.8, ls='--', alpha=0.5); ax.axvline(0, color='gray', lw=0.8, ls='--', alpha=0.5)
    r, p = stats.pearsonr(conc_df['scrna_rho'], conc_df['bulk_lfc'])
    ax.text(0.05, 0.95, f'r = {r:.3f}\np = {p:.3f}', transform=ax.transAxes, va='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    ax.set_xlabel('scRNA-seq Spearman rho (trajectory)', fontsize=10)
    ax.set_ylabel('Bulk log2FC (T2D vs Control)', fontsize=10)
    ax.set_title('scRNA-seq Trajectory vs Bulk DE', fontsize=10, fontweight='bold')
    ax.legend(fontsize=8); ax.spines[['top', 'right']].set_visible(False)

    ax = axes[1]
    n_conc, n_disc = conc_df['concordant'].sum(), (~conc_df['concordant']).sum()
    from scipy.stats import binomtest
    bt = binomtest(int(n_conc), n_conc + n_disc, 0.5, alternative='greater')
    ax.bar(['Concordant', 'Discordant'], [n_conc, n_disc], color=['#2ECC71', '#E74C3C'], alpha=0.85, width=0.5)
    ax.axhline(y=(n_conc + n_disc) / 2, color='gray', ls='--', lw=1.5, alpha=0.7, label='50% chance')
    pct = n_conc / (n_conc + n_disc) * 100
    ax.set_title(f'Direction Concordance ({pct:.0f}%)\nbinomial p={bt.pvalue:.3f} vs. chance', fontsize=10, fontweight='bold')
    ax.set_ylabel('Number of lncRNAs', fontsize=10)
    ax.legend(fontsize=8); ax.spines[['top', 'right']].set_visible(False)

    fig.suptitle('Bulk Validation of Trajectory-Associated lncRNAs (GSE221521)', fontsize=12, fontweight='bold', y=1.03)
    savefig(fig, outdir / 'Fig14_bulk_validation.png')
    print("  NOTE: Fig14 panel C (per-sample expression boxplots, originally Fig15) requires the raw")
    print("  bulk count matrix, which is not one of this script's 4 documented inputs. Run")
    print("  scripts/analysis/04_bulk_validation.py's plot_key_lncrna_boxplots() directly if needed.")


def main():
    parser = argparse.ArgumentParser(description='Regenerate final publication figures from computed results')
    parser.add_argument('--full_h5ad', required=True)
    parser.add_argument('--mono_h5ad', required=True)
    parser.add_argument('--traj_csv', required=True)
    parser.add_argument('--bulk_csv', required=True)
    parser.add_argument('--outdir', default='results/figures')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Loading monocyte h5ad...")
    adata_mono = sc.read_h5ad(args.mono_h5ad)
    print("Loading trajectory results...")
    traj_df = pd.read_csv(args.traj_csv)
    print("Loading bulk concordance...")
    bulk_df = pd.read_csv(args.bulk_csv)

    fig9_pseudotime(adata_mono, outdir)
    fig10_heatmap_and_lollipop(adata_mono, traj_df, outdir)

    print("Loading full PBMC h5ad (for Fig13 -- this is the large file)...")
    adata_full = sc.read_h5ad(args.full_h5ad)
    fig13_celltype_specificity(adata_full, outdir)
    del adata_full

    fig14_bulk_validation(bulk_df, outdir)

    print(f"\nFigures written to: {outdir}/")


if __name__ == '__main__':
    main()
