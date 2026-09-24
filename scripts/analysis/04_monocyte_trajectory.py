"""
04_monocyte_trajectory.py
Phase 4 -- Monocyte subtype annotation and pseudotime trajectory.

Uses leiden_0.3 clusters from the monocyte subset (not the full PBMC object).
Cluster-to-subtype mapping is derived from marker expression and validated in
docs/methods_detail.md:

  CD14_Classical   : 5 (n=29,264), 6 (n=6,541)         -- CD14 high, S100A8 high
  CD16_NonClassical: 7 (n=5,664)                          -- FCGR3A high, MS4A7 high, CDKN1C high
  Intermediate     : 1 (n=1,068), 2 (n=1,520)           -- moderate CD14+FCGR3A, elevated HLA-DRA
  Excluded         : dendritic-cell-contaminated and ambiguous clusters

Final clean counts: CD14 Classical 35,805 (81.3%), Intermediate 2,588 (5.9%),
CD16 Non-Classical 5,664 (12.9%) -- within published flow-cytometry ranges.

Computes diffusion pseudotime and cell-type-specific lncRNA expression
(Fig9, Fig13). Trajectory-lncRNA significance testing itself is performed
separately, donor-aware, in 05_donor_aware_statistics.py.

Run:
  python3 scripts/analysis/04_monocyte_trajectory.py \
      --mono_h5ad data/processed/GSE268210_monocytes_annotated.h5ad \
      --full_h5ad data/processed/GSE268210_phase1_full.h5ad \
      --outdir results/tables
"""

import gc, sys, warnings, argparse, logging
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime

warnings.filterwarnings('ignore')
sc.settings.verbosity = 0

log_file = f"phase4_monocyte_trajectory_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s: %(message)s',
    datefmt='%H:%M:%S',
    handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(log_file)]
)
log = logging.getLogger(__name__)
SEED = 42
np.random.seed(SEED)

CLUSTER_MAP = {
    '5': 'CD14_Classical', '6': 'CD14_Classical',
    '7': 'CD16_NonClassical',
    '1': 'Intermediate', '2': 'Intermediate',
    '0': 'Exclude_Ambiguous', '3': 'Exclude_DC', '4': 'Exclude_Ambiguous',
    '8': 'Exclude_DC', '9': 'Exclude_DC', '10': 'Exclude_DC', '12': 'Exclude_DC',
    '13': 'Exclude_Ambiguous', '14': 'Exclude_Ambiguous',
}

SUBTYPE_ORDER = ['CD14_Classical', 'Intermediate', 'CD16_NonClassical']
SUBTYPE_LABELS = ['CD14\nClassical', 'Intermediate', 'CD16\nNon-Classical']
SUBTYPE_COLS = {'CD14_Classical': '#E53E3E', 'Intermediate': '#D69E2E', 'CD16_NonClassical': '#3182CE'}


def savefig(fig, path, dpi=150):
    fig.savefig(path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    log.info(f"  Saved: {Path(path).name}")


def annotate(adata_mono: ad.AnnData) -> ad.AnnData:
    log.info("Step 1: Cluster-based subtype annotation (leiden_0.3)...")
    leiden_col = 'leiden_0.3'
    assert leiden_col in adata_mono.obs.columns, f"'{leiden_col}' not found. Available: {list(adata_mono.obs.columns)}"

    subtype = adata_mono.obs[leiden_col].astype(str).map(CLUSTER_MAP).fillna('Exclude_Ambiguous')
    adata_mono.obs['mono_subtype_raw'] = subtype
    keep = subtype.isin(SUBTYPE_ORDER)
    adata_clean = adata_mono[keep].copy()
    adata_clean.obs['mono_subtype'] = pd.Categorical(adata_clean.obs['mono_subtype_raw'], categories=SUBTYPE_ORDER)

    total = adata_clean.n_obs
    log.info(f"  Clean monocytes : {total:,}")
    for st in SUBTYPE_ORDER:
        n = (adata_clean.obs['mono_subtype'] == st).sum()
        log.info(f"    {st:<25}: {n:>6,}  ({n/total*100:.1f}%)")
    return adata_clean


def run_trajectory(adata_mono: ad.AnnData, outdir: Path) -> ad.AnnData:
    log.info("Step 2: Pseudotime trajectory...")
    use_rep = 'X_pca_harmony' if 'X_pca_harmony' in adata_mono.obsm else \
              'X_pca' if 'X_pca' in adata_mono.obsm else None
    if use_rep is None:
        sc.pp.highly_variable_genes(adata_mono, flavor='seurat', n_top_genes=500)
        sc.tl.pca(adata_mono, n_comps=20, use_highly_variable=True)
        use_rep = 'X_pca'

    sc.pp.neighbors(adata_mono, use_rep=use_rep, n_neighbors=15, n_pcs=20, random_state=SEED)
    sc.tl.umap(adata_mono, random_state=SEED)
    sc.tl.diffmap(adata_mono, n_comps=15)

    cd14_mask = adata_mono.obs['mono_subtype'] == 'CD14_Classical'
    x = adata_mono[cd14_mask, 'CD14'].X
    cd14e = np.array(x.todense()).flatten() if sp.issparse(x) else x.flatten()
    root = adata_mono.obs_names[cd14_mask][np.argmax(cd14e)]
    adata_mono.uns['iroot'] = adata_mono.obs_names.get_loc(root)

    sc.tl.dpt(adata_mono, n_dcs=10)

    log.info("  Pseudotime per subtype:")
    for st in SUBTYPE_ORDER:
        mask = adata_mono.obs['mono_subtype'] == st
        med = adata_mono.obs.loc[mask, 'dpt_pseudotime'].median()
        log.info(f"    {st:<25} median = {med:.4f}")

    _plot_trajectory(adata_mono, outdir)
    return adata_mono


def _plot_trajectory(adata_mono, outdir):
    umap = adata_mono.obsm['X_umap']
    pt = adata_mono.obs['dpt_pseudotime'].values
    counts = {st: (adata_mono.obs['mono_subtype'] == st).sum() for st in SUBTYPE_ORDER}
    total = sum(counts.values())

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'Monocyte Pseudotime Trajectory -- T2D PBMCs\nn={total:,} cells', fontsize=11, fontweight='bold')

    ax = axes[0]
    for st in SUBTYPE_ORDER:
        mask = adata_mono.obs['mono_subtype'] == st
        ax.scatter(umap[mask, 0], umap[mask, 1], c=SUBTYPE_COLS[st], s=0.5, alpha=0.5, rasterized=True,
                   label=f"{st.replace('_',' ')} (n={mask.sum():,})")
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


def lncrna_celltype_specificity(adata_full: ad.AnnData, outdir: Path):
    log.info("Step 3: Cell-type-specific lncRNA expression (Fig 13)...")
    lnc_mask = adata_full.var.get('is_lncrna', pd.Series(False, index=adata_full.var_names)).astype(bool)
    lncrna_genes = adata_full.var_names[lnc_mask].tolist()

    X_lnc = adata_full[:, lncrna_genes].X
    det = np.array((X_lnc > 0).mean(axis=0)).flatten() if sp.issparse(X_lnc) else (X_lnc > 0).mean(axis=0)
    usable = [g for g, d in zip(lncrna_genes, det) if d >= 0.01]

    celltypes = ['T_cell', 'NK', 'Monocyte', 'B_cell', 'DC', 'Platelet']
    mean_expr = {}
    for ct in celltypes:
        mask = adata_full.obs['broad_celltype'] == ct
        if mask.sum() < 50:
            continue
        X_ct = adata_full[mask, :][:, usable].X
        mean_expr[ct] = np.array(X_ct.mean(axis=0)).flatten() if sp.issparse(X_ct) else X_ct.mean(axis=0)

    mean_df = pd.DataFrame(mean_expr, index=usable)
    mean_df.to_csv(outdir / 'lncrna_celltype_mean_expression.csv')

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
    savefig(fig, outdir.parent / 'figures' / 'Fig13_lncrna_celltype_specificity.png')


def main():
    parser = argparse.ArgumentParser(description='Phase 4: Monocyte subtype annotation + pseudotime')
    parser.add_argument('--mono_h5ad', required=True)
    parser.add_argument('--full_h5ad', required=True)
    parser.add_argument('--outdir', default='results/tables')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    (outdir.parent / 'figures').mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info(f"Loading monocyte subset: {args.mono_h5ad}")
    adata_mono = sc.read_h5ad(args.mono_h5ad)
    log.info(f"  {adata_mono.n_obs:,} cells x {adata_mono.n_vars:,} genes")

    adata_mono = annotate(adata_mono)
    adata_mono = run_trajectory(adata_mono, outdir.parent / 'figures')

    log.info(f"Loading full PBMC: {args.full_h5ad}")
    adata_full = sc.read_h5ad(args.full_h5ad)
    lncrna_celltype_specificity(adata_full, outdir)
    del adata_full; gc.collect()

    out_path = outdir.parent.parent / 'data' / 'processed' / 'GSE268210_monocytes_final.h5ad'
    adata_mono.write_h5ad(out_path, compression='gzip')
    log.info(f"Saved: {out_path}")
    log.info("Next step: scripts/analysis/05_donor_aware_statistics.py")


if __name__ == '__main__':
    main()
