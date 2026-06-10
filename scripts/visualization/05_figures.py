"""
Phase 2 — QC Visualisation + Cell Type Verification
Project: GSE268210 T2D PBMC scRNA-seq lncRNA Analysis
Author: Usama Manzoor

Produces 7 publication-ready figure panels:
  Fig1_QC_violins.png          — QC metrics per sample
  Fig2_UMAP_overview.png       — UMAP coloured by sample/celltype/batch
  Fig3_marker_dotplot.png      — canonical marker expression per cluster
  Fig4_monocyte_subcluster.png — CD14/CD16/intermediate monocyte resolution
  Fig5_celltype_proportions.png— cell type composition per sample
  Fig6_lncrna_overview.png     — lncRNA detection rates across cell types
  Fig7_QC_summary_table.png    — per-sample QC statistics table

Run:
  conda activate scrna
  python3 phase2_visualize.py --h5ad ./phase1_output/GSE268210_phase1_full.h5ad
"""

import os
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from pathlib import Path

warnings.filterwarnings('ignore')
sc.settings.verbosity = 0

# ── COLOUR PALETTE ────────────────────────────────────────────────────────────
CELL_COLOURS = {
    'T_cell'  : '#2196F3',
    'NK'      : '#4CAF50',
    'Monocyte': '#FF5722',
    'B_cell'  : '#9C27B0',
    'DC'      : '#FF9800',
    'Platelet': '#F44336',
    'HSPC'    : '#00BCD4',
    'Other'   : '#9E9E9E',
}

SAMPLE_COLOURS = [
    '#1F77B4','#FF7F0E','#2CA02C','#D62728','#9467BD',
    '#8C564B','#E377C2','#7F7F7F','#BCBD22',
]

MONO_COLOURS = {
    'CD14_Classical'   : '#E53E3E',
    'CD16_NonClassical': '#3182CE',
    'Intermediate'     : '#D69E2E',
    'Other_Mono'       : '#A0AEC0',
}

# ── CANONICAL MARKERS FOR DOT PLOT ───────────────────────────────────────────
MARKER_GENES = {
    'CD14 Mono'  : ['CD14', 'LYZ', 'S100A8', 'VCAN'],
    'CD16 Mono'  : ['FCGR3A', 'MS4A7', 'CDKN1C'],
    'Inter Mono' : ['CD14', 'FCGR3A', 'HLA-DRA'],
    'CD4 T'      : ['CD3D', 'CD4', 'IL7R', 'CCR7'],
    'CD8 T'      : ['CD3D', 'CD8A', 'GZMK'],
    'NK'         : ['NCAM1', 'NKG7', 'GNLY'],
    'B cell'     : ['CD79A', 'MS4A1'],
    'pDC'        : ['LILRA4', 'CLEC4C'],
    'cDC'        : ['CD1C', 'FCER1A'],
    'Platelet'   : ['PPBP', 'PF4'],
}

# ── HELPER ────────────────────────────────────────────────────────────────────
def savefig(fig, path, dpi=150):
    fig.savefig(path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: {path}")


# ── FIG 1: QC VIOLIN PLOTS ───────────────────────────────────────────────────
def plot_qc_violins(adata, outdir):
    print("Fig 1: QC violin plots...")

    samples   = adata.obs['sample_id'].cat.categories.tolist()
    n_samples = len(samples)
    metrics   = ['n_genes_by_counts', 'total_counts', 'pct_counts_mito']
    labels    = ['Genes per cell', 'UMI counts', 'Mito %']

    fig, axes = plt.subplots(3, 1, figsize=(max(10, n_samples * 1.2), 10))
    fig.suptitle('QC Metrics per Sample — Post-Filtering', fontsize=14,
                 fontweight='bold', y=1.01)

    for ax, metric, label in zip(axes, metrics, labels):
        data_per_sample = [
            adata.obs.loc[adata.obs['sample_id'] == s, metric].values
            for s in samples
        ]
        parts = ax.violinplot(data_per_sample, positions=range(n_samples),
                              showmedians=True, showextrema=False)

        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(SAMPLE_COLOURS[i % len(SAMPLE_COLOURS)])
            pc.set_alpha(0.7)
        parts['cmedians'].set_color('black')
        parts['cmedians'].set_linewidth(2)

        ax.set_xticks(range(n_samples))
        ax.set_xticklabels(samples, rotation=45, ha='right', fontsize=9)
        ax.set_ylabel(label, fontsize=10)
        ax.spines[['top', 'right']].set_visible(False)

        # Add median text
        for i, d in enumerate(data_per_sample):
            ax.text(i, np.median(d), f'{np.median(d):.0f}',
                    ha='center', va='bottom', fontsize=7, color='black')

    plt.tight_layout()
    savefig(fig, outdir / 'Fig1_QC_violins.png')


# ── FIG 2: UMAP OVERVIEW ─────────────────────────────────────────────────────
def plot_umap_overview(adata, outdir):
    print("Fig 2: UMAP overview...")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle('UMAP Overview — 209,289 T2D PBMCs (9 patients)',
                 fontsize=13, fontweight='bold')

    # Panel A: coloured by cell type
    ax = axes[0]
    broad = adata.obs['broad_celltype'].values
    umap  = adata.obsm['X_umap']
    for ct, col in CELL_COLOURS.items():
        mask = broad == ct
        if mask.sum() > 0:
            ax.scatter(umap[mask, 0], umap[mask, 1],
                       c=col, s=0.3, alpha=0.4, rasterized=True, label=ct)
    ax.set_title('Cell Type', fontsize=11, fontweight='bold')
    ax.set_xlabel('UMAP 1', fontsize=9); ax.set_ylabel('UMAP 2', fontsize=9)
    ax.tick_params(labelsize=7)
    ax.spines[['top', 'right']].set_visible(False)
    legend = ax.legend(markerscale=5, fontsize=7, loc='lower left',
                       framealpha=0.8, ncol=1)

    # Panel B: coloured by sample
    ax = axes[1]
    samples = adata.obs['sample_id'].cat.categories.tolist()
    for i, s in enumerate(samples):
        mask = adata.obs['sample_id'] == s
        ax.scatter(umap[mask, 0], umap[mask, 1],
                   c=SAMPLE_COLOURS[i % len(SAMPLE_COLOURS)],
                   s=0.3, alpha=0.3, rasterized=True, label=s)
    ax.set_title('Sample (Batch)', fontsize=11, fontweight='bold')
    ax.set_xlabel('UMAP 1', fontsize=9); ax.set_ylabel('UMAP 2', fontsize=9)
    ax.tick_params(labelsize=7)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(markerscale=5, fontsize=7, loc='lower left', framealpha=0.8)

    # Panel C: coloured by leiden cluster
    ax = axes[2]
    leiden_col = 'leiden_0.5' if 'leiden_0.5' in adata.obs.columns else \
                 [c for c in adata.obs.columns if c.startswith('leiden')][0]
    clusters = adata.obs[leiden_col].astype(str)
    n_clusters = clusters.nunique()
    cmap = plt.cm.get_cmap('tab20', n_clusters)
    for i, cl in enumerate(sorted(clusters.unique(), key=int)):
        mask = clusters == cl
        ax.scatter(umap[mask, 0], umap[mask, 1],
                   c=[cmap(i)], s=0.3, alpha=0.4, rasterized=True, label=cl)
        # Label cluster centroid
        cx, cy = umap[mask, 0].mean(), umap[mask, 1].mean()
        ax.text(cx, cy, cl, fontsize=7, ha='center', va='center',
                fontweight='bold', color='black')
    ax.set_title(f'Leiden Clusters (res=0.5, n={n_clusters})',
                 fontsize=11, fontweight='bold')
    ax.set_xlabel('UMAP 1', fontsize=9); ax.set_ylabel('UMAP 2', fontsize=9)
    ax.tick_params(labelsize=7)
    ax.spines[['top', 'right']].set_visible(False)

    plt.tight_layout()
    savefig(fig, outdir / 'Fig2_UMAP_overview.png', dpi=200)


# ── FIG 3: MARKER GENE DOT PLOT ──────────────────────────────────────────────
def plot_marker_dotplot(adata, outdir):
    print("Fig 3: Marker gene dot plot...")

    # Use leiden clusters as groupby
    leiden_col = 'leiden_0.5' if 'leiden_0.5' in adata.obs.columns else \
                 [c for c in adata.obs.columns if c.startswith('leiden')][0]

    # Flatten marker dict — keep only genes present in adata
    all_markers = []
    marker_labels = []
    for ct, genes in MARKER_GENES.items():
        present = [g for g in genes if g in adata.var_names]
        all_markers.extend(present)
        marker_labels.extend([ct] * len(present))

    # Remove duplicates while preserving order
    seen = set()
    unique_markers = []
    for g in all_markers:
        if g not in seen:
            unique_markers.append(g)
            seen.add(g)

    present_markers = [g for g in unique_markers if g in adata.var_names]

    if len(present_markers) < 3:
        print("  WARNING: Too few marker genes found — skipping dot plot")
        return

    fig, ax = plt.subplots(figsize=(len(present_markers) * 0.6 + 2, 6))

    sc.pl.dotplot(
        adata,
        var_names=present_markers,
        groupby=leiden_col,
        ax=ax,
        show=False,
        colorbar_title='Mean\nexpression',
        size_title='Fraction\nexpressing',
        standard_scale='var',
    )

    ax.set_title('Marker Gene Expression per Cluster', fontsize=12,
                 fontweight='bold', pad=15)
    ax.set_xlabel('Marker genes', fontsize=10)
    ax.set_ylabel('Leiden cluster', fontsize=10)

    plt.tight_layout()
    savefig(fig, outdir / 'Fig3_marker_dotplot.png')


# ── FIG 4: MONOCYTE SUBCLUSTERING ────────────────────────────────────────────
def plot_monocyte_subcluster(adata, outdir):
    print("Fig 4: Monocyte subclustering...")

    # Load monocyte subset
    mono_path = outdir.parent / 'phase1_output' / 'GSE268210_phase1_monocytes.h5ad'
    if not mono_path.exists():
        # Try same dir as full h5ad
        mono_path = outdir / 'GSE268210_phase1_monocytes.h5ad'
    if not mono_path.exists():
        mono_mask = adata.obs['broad_celltype'] == 'Monocyte'
        adata_mono = adata[mono_mask].copy()
        print(f"  Using monocytes from full adata: {adata_mono.n_obs:,} cells")
    else:
        adata_mono = sc.read_h5ad(mono_path)
        print(f"  Loaded monocyte subset: {adata_mono.n_obs:,} cells")

    # Re-cluster monocytes at higher resolution for subtype resolution
    print("  Subclustering monocytes...")
    if 'X_pca_harmony' in adata_mono.obsm:
        use_rep = 'X_pca_harmony'
    elif 'X_pca' in adata_mono.obsm:
        use_rep = 'X_pca'
    else:
        # Need to recompute PCA for monocyte subset
        sc.pp.highly_variable_genes(adata_mono, flavor='seurat', n_top_genes=500)
        sc.tl.pca(adata_mono, n_comps=20, use_highly_variable=True)
        use_rep = 'X_pca'

    sc.pp.neighbors(adata_mono, use_rep=use_rep, n_neighbors=15, n_pcs=20)
    sc.tl.umap(adata_mono, random_state=42)

    # Leiden at high resolution to separate CD14/CD16/intermediate
    try:
        sc.tl.leiden(adata_mono, resolution=0.8, key_added='mono_leiden',
                     flavor='igraph', n_iterations=2, directed=False)
    except TypeError:
        sc.tl.leiden(adata_mono, resolution=0.8, key_added='mono_leiden')

    # Score monocyte subtypes
    cd14_markers  = [g for g in ['CD14', 'S100A8', 'S100A9', 'VCAN', 'FCN1', 'LYZ']
                     if g in adata_mono.var_names]
    cd16_markers  = [g for g in ['FCGR3A', 'MS4A7', 'CDKN1C', 'LST1', 'VMO1']
                     if g in adata_mono.var_names]
    inter_markers = [g for g in ['CD14', 'FCGR3A', 'HLA-DRA', 'HLA-DRB1', 'CD74']
                     if g in adata_mono.var_names]

    if len(cd14_markers) >= 2:
        sc.tl.score_genes(adata_mono, cd14_markers,  score_name='score_CD14',  random_state=42)
    if len(cd16_markers) >= 2:
        sc.tl.score_genes(adata_mono, cd16_markers,  score_name='score_CD16',  random_state=42)
    if len(inter_markers) >= 2:
        sc.tl.score_genes(adata_mono, inter_markers, score_name='score_Inter', random_state=42)

    # Assign subtype per cell
    score_cols = [c for c in ['score_CD14', 'score_CD16', 'score_Inter']
                  if c in adata_mono.obs.columns]
    if len(score_cols) >= 2:
        subtype_map = {'score_CD14': 'CD14_Classical',
                       'score_CD16': 'CD16_NonClassical',
                       'score_Inter': 'Intermediate'}
        adata_mono.obs['mono_subtype'] = (
            adata_mono.obs[score_cols].idxmax(axis=1).map(subtype_map).fillna('Other_Mono')
        )
    else:
        adata_mono.obs['mono_subtype'] = 'CD14_Classical'

    # Subtype counts
    subtype_counts = adata_mono.obs['mono_subtype'].value_counts()
    print(f"  Monocyte subtypes:")
    for st, cnt in subtype_counts.items():
        print(f"    {st}: {cnt:,} ({cnt/adata_mono.n_obs*100:.1f}%)")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'Monocyte Subtype Resolution — {adata_mono.n_obs:,} cells',
                 fontsize=13, fontweight='bold')

    umap = adata_mono.obsm['X_umap']

    # Panel A: subtype
    ax = axes[0]
    for st, col in MONO_COLOURS.items():
        mask = adata_mono.obs['mono_subtype'] == st
        if mask.sum() > 0:
            ax.scatter(umap[mask, 0], umap[mask, 1], c=col, s=1, alpha=0.5,
                       rasterized=True, label=f"{st} (n={mask.sum():,})")
    ax.set_title('Monocyte Subtypes', fontsize=11, fontweight='bold')
    ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
    ax.legend(markerscale=4, fontsize=8, loc='best', framealpha=0.8)
    ax.spines[['top', 'right']].set_visible(False)

    # Panel B: CD14 expression
    ax = axes[1]
    if 'CD14' in adata_mono.var_names:
        cd14_expr = np.array(adata_mono[:, 'CD14'].X.todense()).flatten()
        sc_plot = ax.scatter(umap[:, 0], umap[:, 1], c=cd14_expr,
                             cmap='Reds', s=0.5, alpha=0.6, rasterized=True,
                             vmin=0, vmax=np.percentile(cd14_expr, 99))
        plt.colorbar(sc_plot, ax=ax, shrink=0.8, label='Log-norm expression')
    ax.set_title('CD14 Expression', fontsize=11, fontweight='bold')
    ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
    ax.spines[['top', 'right']].set_visible(False)

    # Panel C: FCGR3A (CD16) expression
    ax = axes[2]
    if 'FCGR3A' in adata_mono.var_names:
        cd16_expr = np.array(adata_mono[:, 'FCGR3A'].X.todense()).flatten()
        sc_plot = ax.scatter(umap[:, 0], umap[:, 1], c=cd16_expr,
                             cmap='Blues', s=0.5, alpha=0.6, rasterized=True,
                             vmin=0, vmax=np.percentile(cd16_expr, 99))
        plt.colorbar(sc_plot, ax=ax, shrink=0.8, label='Log-norm expression')
    ax.set_title('FCGR3A (CD16) Expression', fontsize=11, fontweight='bold')
    ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
    ax.spines[['top', 'right']].set_visible(False)

    plt.tight_layout()
    savefig(fig, outdir / 'Fig4_monocyte_subcluster.png', dpi=200)

    # Save annotated monocyte subset for Phase 3
    mono_out = outdir / 'GSE268210_monocytes_annotated.h5ad'
    adata_mono.write_h5ad(mono_out, compression='gzip')
    print(f"  Saved annotated monocytes: {mono_out}")

    return adata_mono


# ── FIG 5: CELL TYPE PROPORTIONS ─────────────────────────────────────────────
def plot_celltype_proportions(adata, outdir):
    print("Fig 5: Cell type proportions...")

    samples   = adata.obs['sample_id'].cat.categories.tolist()
    celltypes = list(CELL_COLOURS.keys())

    # Build proportion matrix
    prop_matrix = pd.DataFrame(index=samples, columns=celltypes, dtype=float)
    for s in samples:
        s_obs   = adata.obs[adata.obs['sample_id'] == s]
        n_total = len(s_obs)
        for ct in celltypes:
            prop_matrix.loc[s, ct] = (s_obs['broad_celltype'] == ct).sum() / n_total * 100
    prop_matrix = prop_matrix.fillna(0)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle('Cell Type Composition per Patient', fontsize=13, fontweight='bold')

    # Panel A: stacked bar
    ax = axes[0]
    bottom = np.zeros(len(samples))
    x      = np.arange(len(samples))
    for ct in celltypes:
        vals = prop_matrix[ct].values.astype(float)
        if vals.sum() > 0:
            ax.bar(x, vals, bottom=bottom,
                   color=CELL_COLOURS.get(ct, '#9E9E9E'), label=ct, width=0.7)
            bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Cell type proportion (%)', fontsize=10)
    ax.set_title('Stacked Proportions', fontsize=11)
    ax.legend(loc='upper right', fontsize=8, bbox_to_anchor=(1.15, 1))
    ax.spines[['top', 'right']].set_visible(False)

    # Panel B: heatmap of proportions
    ax = axes[1]
    im = ax.imshow(prop_matrix.T.values.astype(float),
                   aspect='auto', cmap='YlOrRd', vmin=0)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha='right', fontsize=9)
    non_zero_ct = [ct for ct in celltypes if prop_matrix[ct].sum() > 0]
    ax.set_yticks(range(len(celltypes)))
    ax.set_yticklabels(celltypes, fontsize=9)
    plt.colorbar(im, ax=ax, shrink=0.8, label='% of cells')
    ax.set_title('Proportion Heatmap', fontsize=11)

    # Add values to heatmap cells
    for i, ct in enumerate(celltypes):
        for j, s in enumerate(samples):
            val = prop_matrix.loc[s, ct]
            if val > 3:
                ax.text(j, i, f'{val:.0f}', ha='center', va='center',
                        fontsize=7, color='black' if val < 50 else 'white')

    plt.tight_layout()
    savefig(fig, outdir / 'Fig5_celltype_proportions.png')


# ── FIG 6: lncRNA DETECTION OVERVIEW ─────────────────────────────────────────
def plot_lncrna_overview(adata, outdir):
    print("Fig 6: lncRNA overview...")

    lncrna_mask = adata.var.get('is_lncrna', pd.Series(False, index=adata.var_names)).astype(bool)
    n_lncrna    = lncrna_mask.sum()

    if n_lncrna == 0:
        print("  No lncRNA annotations found — skipping")
        return

    adata_lnc = adata[:, lncrna_mask]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(f'lncRNA Detection Overview — {n_lncrna:,} lncRNA genes',
                 fontsize=13, fontweight='bold')

    # Panel A: lncRNA genes per cell (violin by cell type)
    ax = axes[0]
    lncrna_counts = np.array((adata_lnc.X > 0).sum(axis=1)).flatten()
    adata.obs['n_lncrna_detected'] = lncrna_counts

    celltypes  = [ct for ct in CELL_COLOURS.keys()
                  if (adata.obs['broad_celltype'] == ct).sum() > 100]
    data_parts = [adata.obs.loc[adata.obs['broad_celltype'] == ct,
                                'n_lncrna_detected'].values for ct in celltypes]
    parts = ax.violinplot(data_parts, positions=range(len(celltypes)),
                          showmedians=True, showextrema=False)
    for i, pc in enumerate(parts['bodies']):
        ct = celltypes[i]
        pc.set_facecolor(CELL_COLOURS.get(ct, '#9E9E9E'))
        pc.set_alpha(0.7)
    parts['cmedians'].set_color('black')
    ax.set_xticks(range(len(celltypes)))
    ax.set_xticklabels(celltypes, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('lncRNA genes detected per cell', fontsize=9)
    ax.set_title('lncRNA Detection per Cell Type', fontsize=10, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)

    # Panel B: Top 20 most detected lncRNAs
    ax = axes[1]
    detection_rate = np.array((adata_lnc.X > 0).mean(axis=0)).flatten() * 100
    top20_idx  = np.argsort(detection_rate)[-20:][::-1]
    top20_genes = adata_lnc.var_names[top20_idx]
    top20_rates = detection_rate[top20_idx]

    colors = ['#E53E3E' if r > 20 else '#3182CE' if r > 10 else '#718096'
              for r in top20_rates]
    bars = ax.barh(range(len(top20_genes)), top20_rates, color=colors, alpha=0.8)
    ax.set_yticks(range(len(top20_genes)))
    ax.set_yticklabels(top20_genes, fontsize=8)
    ax.set_xlabel('% cells expressing', fontsize=9)
    ax.set_title('Top 20 Detected lncRNAs\n(across all cells)', fontsize=10, fontweight='bold')
    ax.axvline(x=5, color='red', linestyle='--', alpha=0.5, linewidth=1, label='>5% threshold')
    ax.legend(fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)

    for bar, rate in zip(bars, top20_rates):
        ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height()/2,
                f'{rate:.1f}%', va='center', fontsize=7)

    # Panel C: lncRNA detection rate histogram
    ax = axes[2]
    ax.hist(detection_rate, bins=50, color='#2196F3', alpha=0.7, edgecolor='white')
    ax.axvline(x=5, color='red', linestyle='--', linewidth=2,
               label=f'>5% filter\n(n={int((detection_rate>5).sum())} lncRNAs)')
    ax.set_xlabel('Detection rate (% cells)', fontsize=9)
    ax.set_ylabel('Number of lncRNA genes', fontsize=9)
    ax.set_title('lncRNA Detection Rate Distribution', fontsize=10, fontweight='bold')
    ax.legend(fontsize=9)
    n_pass = int((detection_rate > 5).sum())
    ax.text(0.98, 0.95, f'Pass >5% filter:\n{n_pass} lncRNAs',
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    ax.spines[['top', 'right']].set_visible(False)

    plt.tight_layout()
    savefig(fig, outdir / 'Fig6_lncrna_overview.png')

    # Print key stats
    print(f"  Total lncRNA genes: {n_lncrna:,}")
    print(f"  lncRNAs >5% detection: {n_pass:,} → these are your usable lncRNAs for Phase 3")
    print(f"  Top 5 lncRNAs: {list(top20_genes[:5])}")

    return n_pass


# ── FIG 7: QC SUMMARY TABLE ──────────────────────────────────────────────────
def plot_qc_table(adata, outdir, qc_csv=None):
    print("Fig 7: QC summary table...")

    # Build table from adata or CSV
    if qc_csv and Path(qc_csv).exists():
        df = pd.read_csv(qc_csv)
    else:
        rows = []
        for s in adata.obs['sample_id'].cat.categories:
            s_obs = adata.obs[adata.obs['sample_id'] == s]
            rows.append({
                'Sample'         : s,
                'Cells (post-QC)': f"{len(s_obs):,}",
                'Med. genes/cell': f"{s_obs['n_genes_by_counts'].median():.0f}",
                'Med. UMI/cell'  : f"{s_obs['total_counts'].median():.0f}",
                'Med. mito%'     : f"{s_obs['pct_counts_mito'].median():.1f}%",
                'Monocytes'      : f"{(s_obs['broad_celltype']=='Monocyte').sum():,}",
                'T cells'        : f"{(s_obs['broad_celltype']=='T_cell').sum():,}",
                'NK cells'       : f"{(s_obs['broad_celltype']=='NK').sum():,}",
            })
        df = pd.DataFrame(rows)

    # Always build from adata directly (avoids CSV column name mismatch)
    rows = []
    for s in adata.obs['sample_id'].cat.categories:
        s_obs = adata.obs[adata.obs['sample_id'] == s]
        rows.append({
            'Sample'          : s,
            'Cells (post-QC)' : f"{len(s_obs):,}",
            'Med. genes/cell' : f"{s_obs['n_genes_by_counts'].median():.0f}",
            'Med. UMI/cell'   : f"{s_obs['total_counts'].median():.0f}",
            'Med. mito%'      : f"{s_obs['pct_counts_mito'].median():.1f}%",
            'Monocytes'       : f"{(s_obs['broad_celltype']=='Monocyte').sum():,}",
            'T cells'         : f"{(s_obs['broad_celltype']=='T_cell').sum():,}",
            'NK cells'        : f"{(s_obs['broad_celltype']=='NK').sum():,}",
        })
    df = pd.DataFrame(rows)

    display_cols = list(df.columns)
    df_display   = df[display_cols]

    if len(display_cols) == 0:
        print("  WARNING: No columns to display — skipping table")
        return

    fig, ax = plt.subplots(figsize=(16, max(4, len(df) * 0.6 + 2)))
    ax.axis('off')

    table = ax.table(
        cellText  = df_display.values,
        colLabels = df_display.columns.tolist(),
        cellLoc   = 'center',
        loc       = 'center',
        colWidths = [0.12] * len(display_cols),
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.2, 1.8)

    # Style header
    for j in range(len(display_cols)):
        table[0, j].set_facecolor('#1A3A5C')
        table[0, j].set_text_props(color='white', fontweight='bold')

    # Alternating row colours
    for i in range(1, len(df) + 1):
        for j in range(len(display_cols)):
            table[i, j].set_facecolor('#F0F4F8' if i % 2 == 0 else 'white')

    ax.set_title('Per-Sample QC Statistics — GSE268210 T2D PBMC scRNA-seq',
                 fontsize=12, fontweight='bold', pad=20)

    plt.tight_layout()
    savefig(fig, outdir / 'Fig7_QC_summary_table.png')


# ── MAIN ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description='Phase 2: QC Visualisation')
    parser.add_argument('--h5ad',    required=True,
                        help='Path to GSE268210_phase1_full.h5ad')
    parser.add_argument('--outdir',  default='./phase2_figures',
                        help='Output directory for figures')
    parser.add_argument('--qc_csv',  default=None,
                        help='Path to phase1_qc_report.csv (optional)')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("PHASE 2 — QC VISUALISATION")
    print(f"Loading: {args.h5ad}")
    print("=" * 60)

    adata = sc.read_h5ad(args.h5ad)
    print(f"Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"obs columns: {list(adata.obs.columns)}")
    print(f"obsm keys: {list(adata.obsm.keys())}")
    print(f"var columns: {list(adata.var.columns)}")

    # Verify UMAP exists
    if 'X_umap' not in adata.obsm:
        print("ERROR: X_umap not found. Re-run Phase 1.")
        sys.exit(1)

    # Verify broad_celltype exists
    if 'broad_celltype' not in adata.obs.columns:
        print("WARNING: broad_celltype not found — using leiden_0.5 as proxy")
        leiden_col = [c for c in adata.obs.columns if c.startswith('leiden')][0]
        adata.obs['broad_celltype'] = adata.obs[leiden_col]

    print()
    plot_qc_violins(adata, outdir)
    plot_umap_overview(adata, outdir)
    plot_marker_dotplot(adata, outdir)
    adata_mono = plot_monocyte_subcluster(adata, outdir)
    plot_celltype_proportions(adata, outdir)
    n_usable_lncrna = plot_lncrna_overview(adata, outdir)
    plot_qc_table(adata, outdir, args.qc_csv)

    print()
    print("=" * 60)
    print("PHASE 2 COMPLETE")
    print(f"Figures saved to: {outdir}/")
    print()
    print("REVIEW CHECKLIST before proceeding to Phase 3:")
    print("  [ ] Fig2: UMAP shows clean cell type separation")
    print("  [ ] Fig2: No batch effects visible (samples mix well)")
    print("  [ ] Fig3: Marker genes match expected cell types")
    print("  [ ] Fig4: Monocyte subtypes (CD14/CD16/intermediate) resolved")
    print("  [ ] Fig5: Cell type proportions consistent across samples")
    print("  [ ] Fig6: lncRNA detection — confirm usable count > 100")
    if n_usable_lncrna:
        print(f"  ✓  {n_usable_lncrna} lncRNAs pass >5% detection threshold")
    print()
    print("If all checks pass → proceed to Phase 3 (lncRNA analysis)")
    print("=" * 60)


if __name__ == '__main__':
    main()
