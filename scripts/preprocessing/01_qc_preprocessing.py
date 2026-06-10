"""
Phase 1 — Pre-Processing Pipeline
GSE268210 | T2D PBMC scRNA-seq | lncRNA Analysis Project
Author: Usama Manzoor

Tuned for:
  - scanpy 1.11.5 / anndata 0.11.4
  - 11GB RAM (processes samples sequentially, not all at once)
  - 10 T2D GEX samples (DM_1 through DM_10)
  - Files in flat directory (no subdirectories)

Run:
  conda activate scrna
  python3 phase1_preprocess.py --data_dir /mnt/d/Cancer\ research/SC_RNA_Seq/files/

Outputs (in ./phase1_output/):
  GSE268210_phase1_full.h5ad        — all cells, all genes
  GSE268210_phase1_lncrna_only.h5ad — all cells, lncRNA genes only
  GSE268210_phase1_monocytes.h5ad   — monocyte subset only
  phase1_qc_report.csv              — per-sample QC stats
  phase1_summary.txt                — final counts
"""

import os
import sys
import gc
import gzip
import zlib
import argparse
import warnings
import logging
import re
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
from pathlib import Path
from datetime import datetime

warnings.filterwarnings('ignore')

# ── LOGGING ───────────────────────────────────────────────────────────────────
log_file = f"phase1_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s: %(message)s',
    datefmt='%H:%M:%S',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_file)
    ]
)
log = logging.getLogger(__name__)

sc.settings.verbosity = 0   # quiet scanpy internal messages

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
CONFIG = {
    # QC thresholds — PBMC blood-appropriate
    'min_genes'      : 200,
    'max_genes'      : 5000,
    'min_counts'     : 500,
    'max_pct_mito'   : 20.0,
    'min_cells'      : 10,

    # PCA / clustering
    'n_pcs'          : 30,
    'n_neighbors'    : 15,
    'leiden_resolutions': [0.3, 0.5, 0.8],

    # Reproducibility
    'seed'           : 42,

    # Memory: process this many samples before writing intermediate checkpoint
    'checkpoint_every': 5,
}

# ── KNOWN LNCRNA NAME PATTERNS ────────────────────────────────────────────────
# Partial detection (~70% recall) — use when GTF file not available
LNCRNA_PATTERN = re.compile(
    r'^(LINC\d|AL\d{6}|AC\d{6}|SNHG|NEAT\d|MALAT|HOTAIR|XIST|'
    r'KCNQ\dOT|MEG\d|MIR\d+HG|DLEU|MIAT|NORAD|PVT1|FTX|'
    r'HAGLR|DANCR|CASC|GUARDIN|TERRA|GAS\d|PURPL|DIGIT|'
    r'THUMPD3-AS|LINC|RP\d+\.\d+)',
    re.IGNORECASE
)

# PBMC canonical marker genes
PBMC_MARKERS = {
    'CD14_Mono'    : ['CD14', 'LYZ', 'S100A8', 'S100A9', 'VCAN', 'FCN1'],
    'CD16_Mono'    : ['FCGR3A', 'MS4A7', 'CDKN1C', 'LST1'],
    'Inter_Mono'   : ['CD14', 'FCGR3A', 'HLA-DRA', 'CD74'],
    'CD4_T'        : ['CD3D', 'CD4', 'IL7R', 'CCR7'],
    'CD8_T'        : ['CD3D', 'CD8A', 'CD8B', 'GZMK'],
    'NK'           : ['NCAM1', 'NKG7', 'GNLY', 'KLRF1'],
    'B_cell'       : ['CD79A', 'MS4A1', 'CD79B'],
    'pDC'          : ['LILRA4', 'IL3RA', 'CLEC4C'],
    'cDC'          : ['CD1C', 'FCER1A', 'CLEC10A'],
    'Platelet'     : ['PPBP', 'PF4', 'GP1BA'],
    'HSPC'         : ['CD34', 'SPINK2'],
}

BROAD_MAP = {
    'CD14_Mono': 'Monocyte', 'CD16_Mono': 'Monocyte', 'Inter_Mono': 'Monocyte',
    'CD4_T': 'T_cell', 'CD8_T': 'T_cell',
    'NK': 'NK', 'B_cell': 'B_cell',
    'pDC': 'DC', 'cDC': 'DC',
    'Platelet': 'Platelet', 'HSPC': 'HSPC',
}


# ── STEP 1: BUILD SAMPLE MANIFEST ─────────────────────────────────────────────
def build_manifest(data_dir: Path) -> pd.DataFrame:
    """Find all T2D GEX samples in the flat directory."""
    barcodes = sorted(data_dir.glob("GSM*_GEX_barcodes.tsv.gz"))

    if not barcodes:
        log.error(f"No GEX barcode files found in {data_dir}")
        log.error("Expected pattern: GSM*_GEX_barcodes.tsv.gz")
        sys.exit(1)

    rows = []
    for bf in barcodes:
        stem   = bf.name.replace('_barcodes.tsv.gz', '')  # GSM8287977_PBMC_DM_1_GEX
        gsm_id = stem.split('_')[0]

        # Extract batch number
        m = re.search(r'_(DM|ND|CTRL|HC|control)_(\d+)_GEX', stem, re.I)
        batch  = int(m.group(2)) if m else 0
        group  = 'T2D' if m and m.group(1).upper() == 'DM' else 'Control'

        feat_f   = data_dir / bf.name.replace('barcodes', 'features')
        matrix_f = data_dir / bf.name.replace('barcodes.tsv.gz', 'matrix.mtx.gz')
        meta_f   = data_dir / f"{gsm_id}_PBMC_{m.group(1) if m else 'DM'}_{batch}_meta_data.csv.gz"

        # Verify matrix exists and is not corrupted (check file size)
        mat_size_mb = matrix_f.stat().st_size / 1e6 if matrix_f.exists() else 0

        rows.append({
            'sample_id'   : f"{group}_{batch:02d}",
            'gsm_id'      : gsm_id,
            'group'       : group,
            'batch'       : batch,
            'barcodes'    : str(bf),
            'features'    : str(feat_f),
            'matrix'      : str(matrix_f),
            'meta'        : str(meta_f) if meta_f.exists() else None,
            'matrix_mb'   : round(mat_size_mb, 1),
            'files_ok'    : bf.exists() and feat_f.exists() and matrix_f.exists() and mat_size_mb > 10,
        })

    df = pd.DataFrame(rows).sort_values('batch').reset_index(drop=True)

    log.info(f"Sample manifest ({len(df)} samples):")
    for _, r in df.iterrows():
        status = '✓' if r['files_ok'] else '✗ SKIP'
        log.info(f"  {status} {r['sample_id']} ({r['gsm_id']}) | matrix {r['matrix_mb']}MB")

    bad = df[~df['files_ok']]
    if len(bad):
        log.warning(f"  {len(bad)} samples will be skipped (missing/corrupted files)")

    return df[df['files_ok']].reset_index(drop=True)


# ── STEP 2: LOAD ONE SAMPLE ───────────────────────────────────────────────────
def load_sample(row: pd.Series) -> ad.AnnData:
    """Load a single sample's 10x count matrix."""
    log.info(f"  Loading matrix...")

    # Load features (genes)
    features = pd.read_csv(
        row['features'], sep='\t', header=None,
        names=['ensembl_id', 'gene_name', 'feature_type'],
        compression='gzip'
    )
    # Keep Gene Expression features only
    features = features[features['feature_type'] == 'Gene Expression'].reset_index(drop=True)

    # Load barcodes
    barcodes = pd.read_csv(row['barcodes'], header=None, compression='gzip')[0].tolist()

    # Load sparse matrix — try gzip first, fall back to direct read
    try:
        with gzip.open(row['matrix'], 'rb') as f:
            mat = scipy.io.mmread(f).tocsc().astype(np.float32)
    except (OSError, zlib.error) as e:
        log.warning(f"  gzip read failed ({e}), trying direct scipy read...")
        try:
            mat = scipy.io.mmread(row['matrix']).tocsc().astype(np.float32)
        except Exception as e2:
            raise RuntimeError(
                f"Matrix file corrupted: {row['matrix']}\n"
                f"Re-download with:\n"
                f"  rm {row['matrix']}\n"
                f"  aria2c -x 8 -s 8 'https://ftp.ncbi.nlm.nih.gov/geo/samples/"
                f"GSM8287nnn/{row['gsm_id']}/suppl/{Path(row['matrix']).name}'"
            ) from e2

    # mat is genes x cells — transpose to cells x genes
    # Subset to Gene Expression features only (by row index)
    gex_idx = features.index.tolist()
    mat = mat[gex_idx, :].T  # cells x genes

    # Build AnnData
    adata = ad.AnnData(
        X   = sp.csr_matrix(mat),
        obs = pd.DataFrame(index=barcodes),
        var = features[['ensembl_id', 'gene_name']].set_index('gene_name'),
    )
    adata.var_names_make_unique()

    # Prefix barcodes with sample ID to ensure uniqueness after merging
    adata.obs_names = [f"{row['sample_id']}_{bc}" for bc in barcodes]

    # Add sample metadata
    adata.obs['sample_id'] = row['sample_id']
    adata.obs['gsm_id']    = row['gsm_id']
    adata.obs['group']     = row['group']
    adata.obs['batch']     = str(row['batch'])

    # Load cell-level metadata if available (cluster labels from original paper)
    if row['meta'] and Path(row['meta']).exists():
        try:
            meta_df = pd.read_csv(row['meta'], index_col=0, compression='gzip')
            # Align by barcode suffix
            bc_map = {f"{row['sample_id']}_{bc}": bc for bc in meta_df.index}
            for col in ['celltype', 'cell_type', 'Celltype', 'cluster']:
                if col in meta_df.columns:
                    adata.obs['original_celltype'] = (
                        adata.obs_names.map(
                            {f"{row['sample_id']}_{k}": v
                             for k, v in meta_df[col].items()}
                        )
                    )
                    log.info(f"    Loaded original cell type labels ({col})")
                    break
        except Exception as e:
            log.warning(f"    Could not load metadata: {e}")

    log.info(f"    Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    return adata


# ── STEP 3: QC METRICS + FILTER ───────────────────────────────────────────────
def run_qc(adata: ad.AnnData, sample_id: str) -> ad.AnnData:
    """Calculate QC metrics, plot, filter."""
    n_before = adata.n_obs

    # Annotate gene types
    adata.var['mito'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mito'],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    # Apply filters
    keep = (
        (adata.obs['n_genes_by_counts'] >= CONFIG['min_genes']) &
        (adata.obs['n_genes_by_counts'] <= CONFIG['max_genes']) &
        (adata.obs['total_counts']       >= CONFIG['min_counts']) &
        (adata.obs['pct_counts_mito']    <= CONFIG['max_pct_mito'])
    )
    adata = adata[keep].copy()
    sc.pp.filter_genes(adata, min_cells=CONFIG['min_cells'])

    n_after = adata.n_obs
    pct_kept = n_after / n_before * 100
    log.info(f"    QC: {n_before:,} → {n_after:,} cells ({pct_kept:.0f}% kept)")
    log.info(f"    Genes remaining: {adata.n_vars:,}")
    log.info(f"    Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
    log.info(f"    Median UMI/cell:   {adata.obs['total_counts'].median():.0f}")
    log.info(f"    Median mito%:      {adata.obs['pct_counts_mito'].median():.1f}%")

    return adata


# ── STEP 4: DOUBLET DETECTION ─────────────────────────────────────────────────
def detect_doublets(adata: ad.AnnData, sample_id: str) -> ad.AnnData:
    """Detect doublets per sample using Scrublet via scanpy."""
    try:
        # scanpy 1.10+ has built-in scrublet wrapper
        sc.pp.scrublet(adata, random_state=CONFIG['seed'])
        n_doublets = adata.obs['predicted_doublet'].sum()
        pct = n_doublets / adata.n_obs * 100
        log.info(f"    Doublets: {n_doublets:,} ({pct:.1f}%) removed")
        adata = adata[~adata.obs['predicted_doublet']].copy()
    except Exception as e:
        log.warning(f"    Scrublet failed ({e}) — skipping doublet removal")
        adata.obs['doublet_score'] = 0.0
        adata.obs['predicted_doublet'] = False
    return adata


# ── STEP 5: NORMALISE ─────────────────────────────────────────────────────────
def normalize_sample(adata: ad.AnnData) -> ad.AnnData:
    """Store raw counts, normalise, log-transform."""
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['log_norm'] = adata.X.copy()
    return adata


# ── STEP 6: lncRNA ANNOTATION ─────────────────────────────────────────────────
def annotate_lncrna(adata: ad.AnnData, lncrna_ids_file: str = None) -> ad.AnnData:
    """
    Add is_lncrna boolean column to adata.var.
    Method A: cross-reference Gencode v32 GTF-derived Ensembl ID list (preferred)
    Method B: gene name patterns (fallback, ~70% recall)
    """
    if lncrna_ids_file and Path(lncrna_ids_file).exists():
        log.info("  Annotating lncRNAs from Gencode v32 ID list...")
        with open(lncrna_ids_file) as f:
            lncrna_set = set(l.strip().split('.')[0] for l in f if l.strip())

        if 'ensembl_id' in adata.var.columns:
            base_ids = adata.var['ensembl_id'].str.split('.').str[0]
            adata.var['is_lncrna'] = base_ids.isin(lncrna_set)
        else:
            adata.var['is_lncrna'] = adata.var_names.map(
                lambda x: bool(LNCRNA_PATTERN.match(str(x)))
            )
    else:
        log.warning("  No lncRNA IDs file — using gene name patterns (~70% recall)")
        log.warning("  To improve: generate lncrna_gencode_v32_ids.txt (see setup_environment.sh)")
        adata.var['is_lncrna'] = adata.var_names.map(
            lambda x: bool(LNCRNA_PATTERN.match(str(x)))
        )

    n_lnc = int(adata.var['is_lncrna'].sum())
    n_cod = int((~adata.var['is_lncrna'] & ~adata.var.get('mito', False)).sum())
    log.info(f"  lncRNA genes: {n_lnc:,} | Protein-coding (approx): {n_cod:,}")
    return adata


# ── STEP 7: MERGE ALL SAMPLES ─────────────────────────────────────────────────
def merge_samples(adatas: list) -> ad.AnnData:
    """Concatenate all per-sample AnnData objects."""
    log.info(f"Merging {len(adatas)} samples...")

    combined = ad.concat(
        adatas,
        join='inner',       # keep only genes present in ALL samples (safer for 11GB RAM)
        merge='first',
        label='sample',
        keys=[a.obs['sample_id'].iloc[0] for a in adatas],
    )
    combined.obs_names_make_unique()

    # Restore var metadata lost during concat (mito, ribo, is_lncrna)
    ref_var = adatas[0].var[['ensembl_id', 'mito', 'ribo', 'is_lncrna']].copy()
    common_genes = combined.var_names.intersection(ref_var.index)
    for col in ['ensembl_id', 'mito', 'ribo', 'is_lncrna']:
        if col in ref_var.columns:
            combined.var[col] = ref_var.loc[combined.var_names, col]

    log.info(f"  Combined: {combined.n_obs:,} cells × {combined.n_vars:,} genes")
    log.info(f"  Group counts:\n{combined.obs['group'].value_counts().to_string()}")
    log.info(f"  Sample counts:\n{combined.obs['sample_id'].value_counts().sort_index().to_string()}")

    return combined


# ── STEP 8: HVG SELECTION ─────────────────────────────────────────────────────
def select_hvgs(adata: ad.AnnData, n_top: int = 2000) -> ad.AnnData:
    """
    Select highly variable genes from protein-coding genes only.
    Memory-efficient: builds lean AnnData from counts layer, no full matrix copy.
    """
    log.info(f"Selecting top {n_top} highly variable genes...")

    exclude = (
        adata.var.get("mito",      pd.Series(False, index=adata.var_names)).astype(bool) |
        adata.var.get("ribo",      pd.Series(False, index=adata.var_names)).astype(bool) |
        adata.var.get("is_lncrna", pd.Series(False, index=adata.var_names)).astype(bool)
    )
    coding_genes = adata.var_names[~exclude]
    log.info(f"  Protein-coding genes for HVG selection: {len(coding_genes):,}")

    # Build lean AnnData from counts layer (avoids copying log-norm matrix)
    X_counts = adata[:, coding_genes].layers["counts"]
    tiny = ad.AnnData(
        X   = X_counts,
        obs = adata.obs[["sample_id"]].copy(),
        var = pd.DataFrame(index=coding_genes),
    )
    sc.pp.normalize_total(tiny, target_sum=1e4)
    sc.pp.log1p(tiny)
    sc.pp.highly_variable_genes(
        tiny, flavor="seurat", n_top_genes=n_top,
        batch_key="sample_id", inplace=True
    )
    hvg_names = tiny.var_names[tiny.var["highly_variable"]]
    adata.var["highly_variable"] = False
    adata.var.loc[hvg_names, "highly_variable"] = True

    log.info(f"  HVGs selected: {int(adata.var['highly_variable'].sum()):,}")
    del tiny, X_counts
    gc.collect()
    return adata


# ── STEP 9: PCA + HARMONY ─────────────────────────────────────────────────────
def run_pca_harmony(adata: ad.AnnData) -> ad.AnnData:
    """
    Memory-efficient PCA via IncrementalPCA + Harmony batch correction.
    Peak RAM per batch: 5000 cells x 2000 HVGs x 4 bytes = ~40 MB.
    """
    import harmonypy as hm
    from sklearn.decomposition import IncrementalPCA

    n_pcs    = CONFIG["n_pcs"]
    n_cells  = adata.n_obs
    hvg_mask = adata.var["highly_variable"].values.astype(bool)
    n_hvg    = int(hvg_mask.sum())
    log.info(f"IncrementalPCA + Harmony | {n_cells:,} cells x {n_hvg:,} HVGs...")

    # Extract HVG sparse matrix and clip values (never densify full matrix)
    X_hvg = adata.X[:, hvg_mask].copy()
    X_hvg.data = np.clip(X_hvg.data, 0, 10)

    batch_size = 5000
    log.info(f"  Fitting IncrementalPCA (batch={batch_size})...")
    ipca = IncrementalPCA(n_components=n_pcs, batch_size=batch_size)
    for start in range(0, n_cells, batch_size):
        end = min(start + batch_size, n_cells)
        ipca.partial_fit(X_hvg[start:end].toarray().astype(np.float32))
        if start % 40000 == 0:
            log.info(f"    fit {end:,}/{n_cells:,}")

    log.info("  Transforming...")
    X_pca = np.zeros((n_cells, n_pcs), dtype=np.float32)
    for start in range(0, n_cells, batch_size):
        end = min(start + batch_size, n_cells)
        X_pca[start:end] = ipca.transform(
            X_hvg[start:end].toarray().astype(np.float32)
        )

    adata.obsm["X_pca"] = X_pca
    del X_hvg, X_pca
    gc.collect()
    log.info(f"  PCA done: {adata.obsm['X_pca'].shape}")

    log.info("  Running Harmony...")
    ho = hm.run_harmony(
        adata.obsm["X_pca"], adata.obs,
        vars_use=["sample_id"],
        max_iter_harmony=30,
        random_state=CONFIG["seed"],
        verbose=False
    )
    # ho.Z_corr shape is (n_pcs, n_cells) — need (n_cells, n_pcs)
    Z = ho.Z_corr.T if ho.Z_corr.shape[0] != adata.n_obs else ho.Z_corr
    adata.obsm["X_pca_harmony"] = Z.astype(np.float32)
    log.info(f"  Harmony done: {adata.obsm['X_pca_harmony'].shape}")
    return adata


# ── STEP 10: NEIGHBOURS + UMAP + LEIDEN ──────────────────────────────────────
def run_clustering(adata: ad.AnnData) -> ad.AnnData:
    """Build graph, UMAP, multi-resolution Leiden clustering."""
    log.info("Building neighbour graph...")
    sc.pp.neighbors(
        adata,
        n_neighbors=CONFIG['n_neighbors'],
        n_pcs=CONFIG['n_pcs'],
        use_rep='X_pca_harmony',
        random_state=CONFIG['seed']
    )

    log.info("Computing UMAP...")
    sc.tl.umap(adata, random_state=CONFIG['seed'])

    log.info("Leiden clustering (3 resolutions)...")
    for res in CONFIG['leiden_resolutions']:
        try:
            # scanpy 1.10+ supports flavor='igraph' without leidenalg
            sc.tl.leiden(
                adata,
                resolution=res,
                key_added=f'leiden_{res}',
                flavor='igraph',
                n_iterations=2,
                directed=False,
                random_state=CONFIG['seed']
            )
        except TypeError:
            # Fallback for older scanpy
            sc.tl.leiden(
                adata,
                resolution=res,
                key_added=f'leiden_{res}',
                random_state=CONFIG['seed']
            )
        n_clusters = adata.obs[f'leiden_{res}'].nunique()
        log.info(f"  res={res}: {n_clusters} clusters")

    return adata


# ── STEP 11: CELL TYPE ANNOTATION ────────────────────────────────────────────
def annotate_celltypes(adata: ad.AnnData) -> ad.AnnData:
    """Score cell types using canonical markers."""
    log.info("Scoring cell types...")

    for ct, markers in PBMC_MARKERS.items():
        present = [m for m in markers if m in adata.var_names]
        if len(present) >= 2:
            sc.tl.score_genes(
                adata,
                gene_list=present,
                score_name=f'score_{ct}',
                random_state=CONFIG['seed']
            )

    # Assign cell type as highest-scoring label
    score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
    if score_cols:
        adata.obs['predicted_celltype'] = (
            adata.obs[score_cols].idxmax(axis=1).str.replace('score_', '', regex=False)
        )
        adata.obs['broad_celltype'] = (
            adata.obs['predicted_celltype'].map(BROAD_MAP).fillna('Other')
        )
        log.info("  Cell type distribution:")
        for ct, cnt in adata.obs['broad_celltype'].value_counts().items():
            pct = cnt / adata.n_obs * 100
            log.info(f"    {ct:<15}: {cnt:>6,} ({pct:.1f}%)")

    return adata


# ── STEP 12: SAVE OUTPUTS ────────────────────────────────────────────────────
def save_outputs(adata: ad.AnnData, output_dir: Path, qc_stats: list):
    """Write h5ad files and summary report."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Full AnnData
    path_full = output_dir / 'GSE268210_phase1_full.h5ad'
    log.info(f"Saving full AnnData ({adata.n_obs:,} cells × {adata.n_vars:,} genes)...")
    adata.write_h5ad(path_full, compression='gzip')
    size_gb = path_full.stat().st_size / 1e9
    log.info(f"  → {path_full} ({size_gb:.2f} GB)")

    # 2. lncRNA-only AnnData
    if 'is_lncrna' in adata.var.columns:
        lnc_mask = adata.var['is_lncrna'].astype(bool)
        if lnc_mask.sum() > 0:
            adata_lnc = adata[:, lnc_mask].copy()
            path_lnc = output_dir / 'GSE268210_phase1_lncrna_only.h5ad'
            adata_lnc.write_h5ad(path_lnc, compression='gzip')
            log.info(f"  → {path_lnc} ({lnc_mask.sum():,} lncRNA genes)")
            del adata_lnc
            gc.collect()

    # 3. Monocyte subset
    if 'broad_celltype' in adata.obs.columns:
        mono_mask = adata.obs['broad_celltype'] == 'Monocyte'
        if mono_mask.sum() > 100:
            adata_mono = adata[mono_mask].copy()
            path_mono = output_dir / 'GSE268210_phase1_monocytes.h5ad'
            adata_mono.write_h5ad(path_mono, compression='gzip')
            log.info(f"  → {path_mono} ({mono_mask.sum():,} monocytes)")
            del adata_mono
            gc.collect()

    # 4. QC report
    qc_df = pd.DataFrame(qc_stats)
    qc_df.to_csv(output_dir / 'phase1_qc_report.csv', index=False)

    # 5. Summary
    summary = [
        "=" * 60,
        "PHASE 1 COMPLETE — SUMMARY",
        "=" * 60,
        f"Total cells      : {adata.n_obs:,}",
        f"Total genes      : {adata.n_vars:,}",
        f"lncRNA genes     : {int(adata.var.get('is_lncrna', pd.Series(False)).sum()):,}",
        f"HVGs             : {int(adata.var.get('highly_variable', pd.Series(False)).sum()):,}",
        f"Clusters (0.5)   : {adata.obs.get('leiden_0.5', pd.Series(['?'])).nunique()}",
        "",
        "Cell type distribution:",
    ]
    if 'broad_celltype' in adata.obs.columns:
        for ct, cnt in adata.obs['broad_celltype'].value_counts().items():
            summary.append(f"  {ct:<15}: {cnt:,}")
    summary += [
        "",
        "Output files:",
        f"  {path_full}",
        f"  {output_dir}/GSE268210_phase1_lncrna_only.h5ad",
        f"  {output_dir}/GSE268210_phase1_monocytes.h5ad",
        f"  {output_dir}/phase1_qc_report.csv",
        "",
        "Next step: python3 phase2_visualize.py",
        "=" * 60,
    ]
    summary_text = '\n'.join(summary)
    print('\n' + summary_text)
    with open(output_dir / 'phase1_summary.txt', 'w') as f:
        f.write(summary_text)


# ── MAIN ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description='Phase 1: GSE268210 scRNA-seq pre-processing')
    parser.add_argument('--data_dir',    required=True,
                        help='Directory containing GSM*_GEX_*.gz files')
    parser.add_argument('--output_dir',  default='./phase1_output')
    parser.add_argument('--lncrna_ids',  default='./lncrna_gencode_v32_ids.txt',
                        help='Gencode v32 lncRNA Ensembl ID list (optional but recommended)')
    parser.add_argument('--test',        action='store_true',
                        help='Process only 2 samples (fast test run)')
    parser.add_argument('--skip_doublets', action='store_true')
    parser.add_argument('--skip_samples', default='',
                        help='Comma-separated sample_ids to skip e.g. T2D_01,T2D_08')
    args = parser.parse_args()

    data_dir   = Path(args.data_dir)
    output_dir = Path(args.output_dir)

    log.info("=" * 60)
    log.info("PHASE 1 — scRNA-seq PRE-PROCESSING")
    log.info(f"Dataset  : GSE268210 (T2D PBMC)")
    log.info(f"Data dir : {data_dir}")
    log.info(f"Output   : {output_dir}")
    log.info(f"Start    : {datetime.now()}")
    log.info("=" * 60)

    # Build manifest
    manifest = build_manifest(data_dir)
    if args.test:
        manifest = manifest.head(2)
        log.info("[TEST MODE] Processing 2 samples only")

    # Per-sample processing
    adatas    = []
    qc_stats  = []

    skip_ids = set(s.strip() for s in args.skip_samples.split(',') if s.strip())
    if skip_ids:
        log.info(f"Skipping samples: {skip_ids}")
        manifest = manifest[~manifest['sample_id'].isin(skip_ids)].reset_index(drop=True)

    for i, (_, row) in enumerate(manifest.iterrows(), 1):
        log.info(f"\n[{i}/{len(manifest)}] Processing {row['sample_id']} ({row['gsm_id']})")

        try:
            adata = load_sample(row)
        except (RuntimeError, Exception) as e:
            log.error(f"  SKIPPING {row['sample_id']} — load failed: {e}")
            log.error(f"  Fix: re-download matrix, then re-run with --skip_samples for other samples")
            continue

        # QC
        n_raw = adata.n_obs
        adata = run_qc(adata, row['sample_id'])

        # Doublets
        if not args.skip_doublets:
            adata = detect_doublets(adata, row['sample_id'])

        # Normalise
        adata = normalize_sample(adata)

        # lncRNA annotation (per sample, before merge)
        adata = annotate_lncrna(adata, args.lncrna_ids)

        # Store QC stats
        qc_stats.append({
            'sample_id'    : row['sample_id'],
            'gsm_id'       : row['gsm_id'],
            'group'        : row['group'],
            'n_raw'        : n_raw,
            'n_after_qc'   : adata.n_obs,
            'n_genes'      : adata.n_vars,
            'median_genes' : float(adata.obs['n_genes_by_counts'].median()),
            'median_umi'   : float(adata.obs['total_counts'].median()),
            'median_mito'  : float(adata.obs['pct_counts_mito'].median()),
        })

        adatas.append(adata)

        # Save per-sample checkpoint (allows resume if killed)
        ckpt_path = output_dir / f"checkpoint_{row['sample_id']}.h5ad"
        adata.write_h5ad(ckpt_path, compression='gzip')
        log.info(f"    Checkpoint: {ckpt_path.name}")
        gc.collect()

    if not adatas:
        log.error("No samples processed. Check data_dir path.")
        sys.exit(1)

    # Merge all samples
    adata = merge_samples(adatas)
    del adatas
    gc.collect()

    # Save merged checkpoint before expensive PCA (safety net)
    merged_ckpt = output_dir / 'checkpoint_merged.h5ad'
    log.info("Saving merged checkpoint (safety net before PCA)...")
    adata.write_h5ad(merged_ckpt, compression='gzip')
    log.info(f"  Saved: {merged_ckpt}")

    # lncRNA annotation on merged
    log.info("Re-annotating lncRNAs on merged object...")
    adata = annotate_lncrna(adata, args.lncrna_ids)

    # HVG selection
    adata = select_hvgs(adata, n_top=2000)

    # PCA + Harmony
    adata = run_pca_harmony(adata)

    # Clustering
    adata = run_clustering(adata)

    # Cell type annotation
    adata = annotate_celltypes(adata)

    # Save final outputs
    save_outputs(adata, output_dir, qc_stats)

    # Clean up checkpoints on success
    for ckpt in output_dir.glob("checkpoint_*.h5ad"):
        ckpt.unlink()
    log.info("Checkpoints cleaned up.")

    log.info(f"\nPhase 1 complete. Log: {log_file}")


if __name__ == '__main__':
    main()
