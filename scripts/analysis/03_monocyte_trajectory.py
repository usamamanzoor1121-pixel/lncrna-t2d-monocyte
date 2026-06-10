"""
Phase 3c — Final Correct Monocyte Annotation + Trajectory
Author: Usama Manzoor

Uses leiden_0.3 clusters from the MONOCYTE SUBSET h5ad (not full PBMC).
Cluster mapping verified from marker expression:

  CD14_Classical   : 5 (n=29,264), 6 (n=6,541)        — CD14 high, S100A8 high
  CD16_NonClassical: 7 (n=5,664)                        — FCGR3A high, MS4A7 high, CDKN1C high
  Intermediate     : 1 (n=1,068), 2 (n=1,520)          — moderate CD14+FCGR3A, elevated HLA-DRA
  Exclude_DC       : 3,8,9,10,12                        — HLA-DRA very high, monocyte markers ~0
  Exclude_Ambiguous: 0,4,13,14                          — all markers low / too small

Final clean counts (biologically validated):
  CD14 Classical   : 35,805 (81.3%)   literature: 75-85% ✓
  Intermediate     :  2,588  (5.9%)   literature:  2-10% ✓
  CD16 NonClassical:  5,664 (12.9%)   literature:  5-15% ✓

Run:
  conda activate scrna
  python3 phase3c_final.py \
      --mono_h5ad ./phase2_figures/GSE268210_monocytes_annotated.h5ad \
      --full_h5ad ./phase1_output/GSE268210_phase1_full.h5ad \
      --outdir    ./phase3c_results
"""

import gc, sys, warnings, argparse, logging
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import to_rgba
from pathlib import Path
from datetime import datetime
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')
sc.settings.verbosity = 0

log_file = f"phase3c_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s: %(message)s',
    datefmt='%H:%M:%S',
    handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(log_file)]
)
log = logging.getLogger(__name__)
SEED = 42
np.random.seed(SEED)

# ── VERIFIED CLUSTER MAPPING (leiden_0.3 on monocyte subset) ─────────────────
CLUSTER_MAP = {
    '5' : 'CD14_Classical',    # CD14=1.57, S100A8=3.81
    '6' : 'CD14_Classical',    # CD14=1.11, S100A8=2.73
    '7' : 'CD16_NonClassical', # FCGR3A=2.65, MS4A7=2.03, CDKN1C=2.25
    '1' : 'Intermediate',      # CD14=0.26, FCGR3A=0.32, both moderate
    '2' : 'Intermediate',      # FCGR3A=0.45, HLA-DRA=1.56 — transitional
    # Excluded
    '0' : 'Exclude_Ambiguous',
    '3' : 'Exclude_DC',
    '4' : 'Exclude_Ambiguous',
    '8' : 'Exclude_DC',
    '9' : 'Exclude_DC',
    '10': 'Exclude_DC',
    '12': 'Exclude_DC',
    '13': 'Exclude_Ambiguous',
    '14': 'Exclude_Ambiguous',
}

SUBTYPE_ORDER  = ['CD14_Classical', 'Intermediate', 'CD16_NonClassical']
SUBTYPE_LABELS = ['CD14\nClassical', 'Intermediate', 'CD16\nNon-Classical']
SUBTYPE_COLS   = {
    'CD14_Classical'    : '#E53E3E',
    'Intermediate'      : '#D69E2E',
    'CD16_NonClassical' : '#3182CE',
}


def savefig(fig, path, dpi=150):
    fig.savefig(path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    log.info(f"  Saved: {Path(path).name}")


# ════════════════════════════════════════════════════════════════════════════
# STEP 1 — ANNOTATION
# ════════════════════════════════════════════════════════════════════════════
def annotate(adata_mono: ad.AnnData) -> ad.AnnData:
    log.info("Step 1: Cluster-based annotation (leiden_0.3, monocyte subset)...")

    leiden_col = 'leiden_0.3'
    assert leiden_col in adata_mono.obs.columns, \
        f"'{leiden_col}' not found. Available: {list(adata_mono.obs.columns)}"

    subtype = adata_mono.obs[leiden_col].astype(str).map(CLUSTER_MAP).fillna('Exclude_Ambiguous')
    adata_mono.obs['mono_subtype_raw'] = subtype

    # Keep only clean monocytes
    keep = subtype.isin(SUBTYPE_ORDER)
    adata_clean = adata_mono[keep].copy()
    adata_clean.obs['mono_subtype'] = pd.Categorical(
        adata_clean.obs['mono_subtype_raw'],
        categories=SUBTYPE_ORDER
    )

    total = adata_clean.n_obs
    log.info(f"  Clean monocytes : {total:,}")
    for st in SUBTYPE_ORDER:
        n = (adata_clean.obs['mono_subtype'] == st).sum()
        log.info(f"    {st:<25}: {n:>6,}  ({n/total*100:.1f}%)")
    log.info(f"  Excluded        : {(~keep).sum():,}  (DC contamination + ambiguous)")

    # Validate markers
    log.info("  Marker validation:")
    for st in SUBTYPE_ORDER:
        mask = adata_clean.obs['mono_subtype'] == st
        row  = {}
        for gene in ['CD14', 'FCGR3A', 'S100A8', 'HLA-DRA']:
            if gene in adata_clean.var_names:
                x = adata_clean[mask, gene].X
                e = np.array(x.todense()).flatten() if sp.issparse(x) else x.flatten()
                row[gene] = round(float(e.mean()), 3)
        log.info(f"    {st:<25}: {row}")

    return adata_clean


# ════════════════════════════════════════════════════════════════════════════
# STEP 2 — PSEUDOTIME
# ════════════════════════════════════════════════════════════════════════════
def run_trajectory(adata_mono: ad.AnnData, outdir: Path) -> ad.AnnData:
    log.info("Step 2: Pseudotime trajectory...")

    use_rep = 'X_pca_harmony' if 'X_pca_harmony' in adata_mono.obsm else \
              'X_pca'         if 'X_pca'         in adata_mono.obsm else None

    if use_rep is None:
        log.info("  No PCA found — computing fresh PCA...")
        sc.pp.highly_variable_genes(adata_mono, flavor='seurat', n_top_genes=500)
        sc.tl.pca(adata_mono, n_comps=20, use_highly_variable=True)
        use_rep = 'X_pca'

    log.info(f"  Using: {use_rep}")
    sc.pp.neighbors(adata_mono, use_rep=use_rep, n_neighbors=15, n_pcs=20, random_state=SEED)
    sc.tl.umap(adata_mono, random_state=SEED)
    sc.tl.diffmap(adata_mono, n_comps=15)

    # Root = most CD14-high cell in CD14_Classical
    cd14_mask = adata_mono.obs['mono_subtype'] == 'CD14_Classical'
    if 'CD14' in adata_mono.var_names:
        x     = adata_mono[cd14_mask, 'CD14'].X
        cd14e = np.array(x.todense()).flatten() if sp.issparse(x) else x.flatten()
        root  = adata_mono.obs_names[cd14_mask][np.argmax(cd14e)]
    else:
        root  = adata_mono.obs_names[cd14_mask][0]
    adata_mono.uns['iroot'] = adata_mono.obs_names.get_loc(root)
    log.info(f"  Root cell: {root}")

    sc.tl.dpt(adata_mono, n_dcs=10)

    log.info("  Pseudotime per subtype:")
    for st in SUBTYPE_ORDER:
        mask = adata_mono.obs['mono_subtype'] == st
        med  = adata_mono.obs.loc[mask, 'dpt_pseudotime'].median()
        log.info(f"    {st:<25} median = {med:.4f}")

    _plot_trajectory(adata_mono, outdir)
    return adata_mono


def _plot_trajectory(adata_mono, outdir):
    umap = adata_mono.obsm['X_umap']
    pt   = adata_mono.obs['dpt_pseudotime'].values
    counts = {st: (adata_mono.obs['mono_subtype'] == st).sum() for st in SUBTYPE_ORDER}
    total  = sum(counts.values())

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(
        f'Monocyte Pseudotime Trajectory — T2D PBMCs\n'
        f'n={total:,} cells  |  '
        f"CD14={counts['CD14_Classical']:,} (81.3%)  "
        f"Inter={counts['Intermediate']:,} (5.9%)  "
        f"CD16={counts['CD16_NonClassical']:,} (12.9%)",
        fontsize=11, fontweight='bold'
    )

    # Panel A: subtype UMAP
    ax = axes[0]
    for st in SUBTYPE_ORDER:
        mask = adata_mono.obs['mono_subtype'] == st
        ax.scatter(umap[mask, 0], umap[mask, 1],
                   c=SUBTYPE_COLS[st], s=0.5, alpha=0.5, rasterized=True,
                   label=f"{st.replace('_',' ')} (n={mask.sum():,})")
    ax.set_title('Monocyte Subtypes', fontsize=11, fontweight='bold')
    ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
    ax.legend(markerscale=6, fontsize=8, loc='best', framealpha=0.8)
    ax.spines[['top','right']].set_visible(False)

    # Panel B: pseudotime UMAP
    ax = axes[1]
    sc_plot = ax.scatter(umap[:,0], umap[:,1], c=pt,
                         cmap='viridis', s=0.5, alpha=0.6, rasterized=True)
    plt.colorbar(sc_plot, ax=ax, shrink=0.8, label='Pseudotime')
    ax.set_title('Diffusion Pseudotime', fontsize=11, fontweight='bold')
    ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
    ax.spines[['top','right']].set_visible(False)

    # Panel C: violin per subtype
    ax = axes[2]
    for i, st in enumerate(SUBTYPE_ORDER):
        mask = adata_mono.obs['mono_subtype'] == st
        if mask.sum() > 0:
            vp = ax.violinplot(
                [pt[mask.values]], positions=[i],
                showmedians=True, showextrema=True
            )
            for pc in vp['bodies']:
                pc.set_facecolor(SUBTYPE_COLS[st])
                pc.set_alpha(0.75)
            med = np.median(pt[mask.values])
            ax.text(i, med + 0.015, f'{med:.3f}',
                    ha='center', fontsize=9, fontweight='bold', color='black')
    ax.set_xticks(range(3))
    ax.set_xticklabels(SUBTYPE_LABELS, fontsize=9)
    ax.set_ylabel('Pseudotime', fontsize=10)
    ax.set_title('Pseudotime per Subtype', fontsize=11, fontweight='bold')
    ax.spines[['top','right']].set_visible(False)

    plt.tight_layout()
    savefig(fig, outdir / 'Fig9_pseudotime_final.png', dpi=200)


# ════════════════════════════════════════════════════════════════════════════
# STEP 3 — TRAJECTORY lncRNAs
# ════════════════════════════════════════════════════════════════════════════
def find_trajectory_lncrnas(adata_mono, outdir) -> pd.DataFrame:
    log.info("Step 3: Trajectory-associated lncRNAs...")

    lnc_mask = adata_mono.var.get(
        'is_lncrna', pd.Series(False, index=adata_mono.var_names)
    ).astype(bool)
    lncrna_genes = adata_mono.var_names[lnc_mask].tolist()

    X_lnc = adata_mono[:, lncrna_genes].X
    if sp.issparse(X_lnc):
        det = np.array((X_lnc > 0).mean(axis=0)).flatten()
    else:
        det = (X_lnc > 0).mean(axis=0)
    usable = [g for g, d in zip(lncrna_genes, det) if d >= 0.02]
    log.info(f"  Usable lncRNAs (>2% detection): {len(usable):,}")

    pt  = adata_mono.obs['dpt_pseudotime'].values
    n_sub = min(adata_mono.n_obs, 30000)
    idx   = np.random.choice(adata_mono.n_obs, n_sub, replace=False)
    pt_sub = pt[idx]

    results = []
    for i in range(0, len(usable), 200):
        batch = usable[i:i+200]
        X_b   = adata_mono[idx, :][:, batch].X
        if sp.issparse(X_b):
            X_b = X_b.toarray()
        for j, gene in enumerate(batch):
            expr = X_b[:, j]
            if expr.std() < 1e-6:
                continue
            rho, pval = stats.spearmanr(pt_sub, expr)
            results.append({'gene': gene, 'rho': rho, 'pvalue': pval})

    res = pd.DataFrame(results)
    _, padj, _, _ = multipletests(res['pvalue'], method='fdr_bh')
    res['padj']    = padj
    res['abs_rho'] = res['rho'].abs()
    res = res.sort_values('abs_rho', ascending=False)

    sig = res[(res['padj'] < 0.05) & (res['abs_rho'] >= 0.15)]
    log.info(f"  Significant (|rho|>=0.15, FDR<0.05): {len(sig)}")
    for _, row in sig.head(10).iterrows():
        d = '↑ toward CD16' if row['rho'] > 0 else '↑ toward CD14'
        log.info(f"    {row['gene']:<20} rho={row['rho']:+.3f}  padj={row['padj']:.2e}  {d}")

    res.to_csv(outdir / 'trajectory_lncrna_final.csv', index=False)
    _plot_heatmap(adata_mono, sig.head(30), pt, outdir)
    _plot_lollipop(sig, outdir)
    return res


def _plot_heatmap(adata_mono, traj_df, pseudotime, outdir):
    genes = [g for g in traj_df['gene'].tolist() if g in adata_mono.var_names]
    if len(genes) < 3:
        log.warning("  Too few genes for heatmap")
        return

    pt_order = np.argsort(pseudotime)
    n_show   = min(5000, len(pt_order))
    step     = max(1, len(pt_order) // n_show)
    show_idx = pt_order[::step][:n_show]

    X = adata_mono[show_idx, :][:, genes].X
    if sp.issparse(X):
        X = X.toarray()
    X = X.T.astype(np.float32)

    win   = max(5, n_show // 50)
    X_sm  = pd.DataFrame(X).T.rolling(win, center=True, min_periods=1).mean().T.values
    rmean = X_sm.mean(axis=1, keepdims=True)
    rstd  = X_sm.std(axis=1,  keepdims=True)
    rstd[rstd < 1e-6] = 1.0
    X_z   = np.clip((X_sm - rmean) / rstd, -2.5, 2.5)

    n_genes = len(genes)
    fig  = plt.figure(figsize=(14, max(7, n_genes * 0.32 + 3)))
    gs   = gridspec.GridSpec(3, 2, height_ratios=[0.04, 0.04, 1],
                             width_ratios=[3, 1], hspace=0.05, wspace=0.08)

    # Subtype colour bar
    ax_st  = fig.add_subplot(gs[0, 0])
    st_vals = adata_mono.obs['mono_subtype'].values[show_idx]
    rgba   = np.array([[to_rgba(SUBTYPE_COLS.get(str(s), '#9E9E9E')) for s in st_vals]])
    ax_st.imshow(rgba, aspect='auto', interpolation='nearest')
    ax_st.set_xticks([]); ax_st.set_yticks([])
    ax_st.set_ylabel('Subtype', fontsize=7, rotation=0, labelpad=38, va='center')
    for st, col in SUBTYPE_COLS.items():
        ax_st.bar(0, 0, color=col, label=st.replace('_', ' '))
    ax_st.legend(loc='lower right', bbox_to_anchor=(1.0, 1.3),
                 fontsize=6, ncol=3, frameon=False)

    # Pseudotime gradient bar
    ax_pt = fig.add_subplot(gs[1, 0])
    ax_pt.imshow(pseudotime[show_idx].reshape(1, -1),
                 aspect='auto', cmap='viridis', interpolation='nearest')
    ax_pt.set_xticks([]); ax_pt.set_yticks([])
    ax_pt.set_ylabel('PT', fontsize=7, rotation=0, labelpad=20, va='center')

    # Expression heatmap
    ax_hm = fig.add_subplot(gs[2, 0])
    im    = ax_hm.imshow(X_z, aspect='auto', cmap='RdBu_r',
                         vmin=-2, vmax=2, interpolation='nearest')
    ax_hm.set_yticks(range(n_genes))
    ax_hm.set_yticklabels(genes, fontsize=8)
    ax_hm.set_xlabel('Cells ordered by pseudotime  (CD14 → Intermediate → CD16)', fontsize=9)
    ax_hm.set_xticks([])
    plt.colorbar(im, ax=ax_hm, shrink=0.4, label='Z-score', pad=0.01)

    # Spearman rho bar
    ax_rho = fig.add_subplot(gs[2, 1])
    rl     = traj_df.set_index('gene')['rho'].to_dict()
    rvals  = [rl.get(g, 0) for g in genes]
    cols   = ['#E53E3E' if r > 0 else '#3182CE' for r in rvals]
    ax_rho.barh(range(n_genes), rvals, color=cols, alpha=0.85, height=0.7)
    ax_rho.set_yticks(range(n_genes)); ax_rho.set_yticklabels([])
    ax_rho.axvline(x=0, color='black', linewidth=1.0)
    ax_rho.set_xlabel("Spearman ρ\n(+ = ↑ toward CD16)", fontsize=8)
    ax_rho.set_title('Correlation\nwith PT', fontsize=8)
    ax_rho.spines[['top', 'right']].set_visible(False)

    fig.suptitle(
        f'Trajectory-Associated lncRNAs in T2D Monocytes  (n={n_genes})\n'
        'Expression along pseudotime: CD14 Classical → Intermediate → CD16 Non-Classical',
        fontsize=11, fontweight='bold', y=1.01
    )
    savefig(fig, outdir / 'Fig10_trajectory_heatmap_final.png', dpi=180)


def _plot_lollipop(sig, outdir):
    if len(sig) == 0:
        return
    df  = sig.copy().sort_values('rho')
    fig, ax = plt.subplots(figsize=(8, max(6, len(df) * 0.35 + 2)))
    for i, (_, row) in enumerate(df.iterrows()):
        col = '#3182CE' if row['rho'] < 0 else '#E53E3E'
        ax.plot([0, row['rho']], [i, i], color=col, lw=2, alpha=0.75)
        ax.scatter(row['rho'], i, color=col, s=60, zorder=3)
        ha = 'left' if row['rho'] > 0 else 'right'
        offset = 0.005 if row['rho'] > 0 else -0.005
        ax.text(row['rho'] + offset, i, row['gene'],
                va='center', ha=ha, fontsize=8.5)
    ax.axvline(x=0, color='black', lw=1.0)
    ax.set_yticks([])
    ax.set_xlabel('Spearman ρ with pseudotime', fontsize=11)
    ax.set_title(
        f'Trajectory-Associated lncRNAs — T2D Monocytes\n'
        f'|ρ| ≥ 0.15, FDR < 0.05  |  n = {len(df)}',
        fontsize=11, fontweight='bold'
    )
    ax.text(0.97, 0.03, '↑ toward CD16', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=9, color='#E53E3E', style='italic')
    ax.text(0.03, 0.03, '↑ toward CD14', transform=ax.transAxes,
            ha='left', va='bottom', fontsize=9, color='#3182CE', style='italic')
    ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    savefig(fig, outdir / 'Fig10b_trajectory_lollipop_final.png')


# ════════════════════════════════════════════════════════════════════════════
# STEP 4 — CELL-TYPE SPECIFICITY (Fig 13 regenerated)
# ════════════════════════════════════════════════════════════════════════════
def lncrna_celltype_specificity(adata_full: ad.AnnData, outdir: Path):
    log.info("Step 4: Cell-type-specific lncRNA expression (Fig 13)...")

    lnc_mask = adata_full.var.get(
        'is_lncrna', pd.Series(False, index=adata_full.var_names)
    ).astype(bool)
    lncrna_genes = adata_full.var_names[lnc_mask].tolist()

    X_lnc = adata_full[:, lncrna_genes].X
    if sp.issparse(X_lnc):
        det = np.array((X_lnc > 0).mean(axis=0)).flatten()
    else:
        det = (X_lnc > 0).mean(axis=0)
    usable = [g for g, d in zip(lncrna_genes, det) if d >= 0.01]
    log.info(f"  Usable lncRNAs (>1% detection): {len(usable):,}")

    celltypes = ['T_cell', 'NK', 'Monocyte', 'B_cell', 'DC', 'Platelet']
    mean_expr = {}
    for ct in celltypes:
        mask = adata_full.obs['broad_celltype'] == ct
        if mask.sum() < 50:
            continue
        X_ct = adata_full[mask, :][:, usable].X
        if sp.issparse(X_ct):
            mean_expr[ct] = np.array(X_ct.mean(axis=0)).flatten()
        else:
            mean_expr[ct] = X_ct.mean(axis=0)
        log.info(f"  {ct}: {mask.sum():,} cells")

    mean_df = pd.DataFrame(mean_expr, index=usable)
    mean_df.to_csv(outdir / 'lncrna_celltype_mean_expression.csv')

    # Top 40 most variable
    var_genes = mean_df.var(axis=1).nlargest(40).index.tolist()
    plot_data = mean_df.loc[var_genes]
    row_mean  = plot_data.mean(axis=1)
    row_std   = plot_data.std(axis=1).replace(0, 1)
    plot_z    = plot_data.sub(row_mean, axis=0).div(row_std, axis=0)

    fig, ax = plt.subplots(figsize=(10, 13))
    im = ax.imshow(plot_z.values, aspect='auto', cmap='RdBu_r',
                   vmin=-2, vmax=2, interpolation='nearest')
    ax.set_xticks(range(len(mean_df.columns)))
    ax.set_xticklabels(mean_df.columns, fontsize=11, fontweight='bold')
    ax.set_yticks(range(len(var_genes)))
    ax.set_yticklabels(var_genes, fontsize=8)
    for x in np.arange(-0.5, len(mean_df.columns), 1):
        ax.axvline(x, color='white', lw=1.5)
    cbar = plt.colorbar(im, ax=ax, shrink=0.35, pad=0.02)
    cbar.set_label('Z-score (across cell types)', fontsize=9)
    ax.set_title(
        'Cell-Type-Specific lncRNA Expression in T2D PBMCs\n'
        f'(top 40 variable lncRNAs  |  {len(usable)} total with >1% detection)',
        fontsize=12, fontweight='bold', pad=12
    )
    plt.tight_layout()
    savefig(fig, outdir / 'Fig13_lncrna_celltype_specificity.png')

    log.info("  Top 3 cell-type-specific lncRNAs:")
    for ct in mean_df.columns:
        top3 = plot_z[ct].nlargest(3).index.tolist()
        log.info(f"    {ct:<12}: {top3}")


# ════════════════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(description='Phase 3c: Final annotation')
    parser.add_argument('--mono_h5ad', required=True,
                        help='./phase2_figures/GSE268210_monocytes_annotated.h5ad')
    parser.add_argument('--full_h5ad', required=True,
                        help='./phase1_output/GSE268210_phase1_full.h5ad')
    parser.add_argument('--outdir', default='./phase3c_results')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("=" * 65)
    log.info("PHASE 3c — FINAL CORRECT ANNOTATION + TRAJECTORY")
    log.info(f"Start: {datetime.now()}")
    log.info("=" * 65)

    # Load monocyte subset
    log.info(f"\nLoading monocyte subset: {args.mono_h5ad}")
    adata_mono = sc.read_h5ad(args.mono_h5ad)
    log.info(f"  {adata_mono.n_obs:,} cells × {adata_mono.n_vars:,} genes")

    # Load full PBMC
    log.info(f"Loading full PBMC: {args.full_h5ad}")
    adata_full = sc.read_h5ad(args.full_h5ad)
    log.info(f"  {adata_full.n_obs:,} cells × {adata_full.n_vars:,} genes")

    # Step 1
    adata_mono = annotate(adata_mono)

    # Step 2
    adata_mono = run_trajectory(adata_mono, outdir)

    # Step 3
    traj_df = find_trajectory_lncrnas(adata_mono, outdir)

    # Step 4 — Fig 13
    lncrna_celltype_specificity(adata_full, outdir)
    del adata_full; gc.collect()

    # Save annotated monocyte object
    out_path = outdir / 'GSE268210_monocytes_final.h5ad'
    adata_mono.write_h5ad(out_path, compression='gzip')
    log.info(f"\nSaved: {out_path}")

    # Summary
    sig = traj_df[(traj_df['padj'] < 0.05) & (traj_df['abs_rho'] >= 0.15)]
    total = adata_mono.n_obs
    lines = [
        "", "=" * 65,
        "PHASE 3c — FINAL SUMMARY",
        f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "=" * 65, "",
        "MONOCYTE SUBTYPES (final, biologically validated)",
    ]
    for st in SUBTYPE_ORDER:
        n = (adata_mono.obs['mono_subtype'] == st).sum()
        lines.append(f"  {st:<25}: {n:>6,}  ({n/total*100:.1f}%)")
    lines += [
        "", "PSEUDOTIME MEDIANS",
    ]
    for st in SUBTYPE_ORDER:
        mask = adata_mono.obs['mono_subtype'] == st
        med  = adata_mono.obs.loc[mask, 'dpt_pseudotime'].median()
        lines.append(f"  {st:<25}: {med:.4f}")
    lines += [
        "", f"TRAJECTORY lncRNAs: {len(sig)}  (|rho|>=0.15, FDR<0.05)",
    ]
    for _, r in sig.head(14).iterrows():
        d = '↑CD16' if r['rho'] > 0 else '↑CD14'
        lines.append(f"  {r['gene']:<20} rho={r['rho']:+.3f}  padj={r['padj']:.2e}  {d}")
    lines += [
        "",
        "FIGURES PRODUCED",
        "  Fig9_pseudotime_final.png",
        "  Fig10_trajectory_heatmap_final.png",
        "  Fig10b_trajectory_lollipop_final.png",
        "  Fig13_lncrna_celltype_specificity.png",
        "",
        "USE THESE FILES FOR THE PAPER — discard all previous phase3/3b outputs",
        "=" * 65,
    ]
    report = '\n'.join(lines)
    print(report)
    with open(outdir / 'phase3c_summary.txt', 'w') as f:
        f.write(report)

    log.info(f"\nPhase 3c complete. Log: {log_file}")


if __name__ == '__main__':
    main()
