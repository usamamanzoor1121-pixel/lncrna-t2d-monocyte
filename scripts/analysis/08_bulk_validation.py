"""
07_bulk_validation.py
Phase 7 — Bulk Validation of Trajectory lncRNAs
Dataset: GSE221521 (whole blood RNA-seq: T2D + Pre-DM + Control)

Validates the 21 donor-aware-significant, GENCODE-confirmed trajectory
lncRNAs (from 04_donor_aware_statistics.py) against an independent bulk
RNA-seq cohort spanning three disease stages.

Steps:
  1. Download GSE221521 from GEO
  2. Normalize counts
  3. Differential expression: T2D vs Control + Pre-DM vs Control
  4. Direction concordance with scRNA-seq trajectory findings
  5. Generate validation figures

Run:
  pip install GEOparse pydeseq2
  python3 scripts/analysis/07_bulk_validation.py --outdir results/tables

Expected runtime: ~45-90 minutes (download-dependent)
"""

import os
import gc
import sys
import gzip
import warnings
import argparse
import logging
import subprocess
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
from datetime import datetime
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')

log_file = f"phase4_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s: %(message)s',
    datefmt='%H:%M:%S',
    handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(log_file)]
)
log = logging.getLogger(__name__)

# ── DONOR-AWARE TRAJECTORY lncRNAs (from 04_donor_aware_statistics.py) ────────
# 21 GENCODE-confirmed, donor-aware-significant loci. AC119396.1 and GAS7 were
# excluded (GENCODE v32: protein-coding, not lncRNA); SNHG12 did not reach the
# |rho|>=0.15 threshold under donor-aware statistics and is also excluded.
TRAJECTORY_LNCRNAS = {
    # gene: (donor_aware_rho, direction)
    'MALAT1'    : (+0.371, '↑CD16'),
    'AC020916.1': (-0.354, '↑CD14'),
    'AC104809.2': (+0.331, '↑CD16'),
    'AC020656.1': (-0.335, '↑CD14'),
    'NEAT1'     : (-0.259, '↑CD14'),
    'AC007952.4': (-0.236, '↑CD14'),
    'LINC00861' : (+0.225, '↑CD16'),
    'AL133415.1': (-0.199, '↑CD14'),
    'AC020651.2': (+0.208, '↑CD16'),
    'LINC00937' : (-0.206, '↑CD14'),
    'LINC02432' : (+0.204, '↑CD16'),
    'AC064805.1': (+0.197, '↑CD16'),
    'AL139246.5': (+0.193, '↑CD16'),
    'SNHG1'     : (+0.189, '↑CD16'),
    'LINC02345' : (+0.185, '↑CD16'),
    'LINC02384' : (+0.183, '↑CD16'),
    'AC253572.2': (-0.167, '↑CD14'),
    'SNHG8'     : (+0.165, '↑CD16'),
    'LINC02773' : (+0.164, '↑CD16'),
    'LINC01578' : (+0.159, '↑CD16'),
    'AC243960.1': (+0.156, '↑CD16'),
}

TRAJ_GENES = list(TRAJECTORY_LNCRNAS.keys())


# ════════════════════════════════════════════════════════════════════════════
# STEP 1 — DOWNLOAD GSE221521
# ════════════════════════════════════════════════════════════════════════════
def download_gse221521(outdir: Path) -> Path:
    """Download GSE221521 count matrix and metadata from GEO."""
    log.info("Step 1: Downloading GSE221521...")

    data_dir = outdir / 'GSE221521_data'
    data_dir.mkdir(exist_ok=True)

    # Check if already downloaded
    matrix_file = data_dir / 'GSE221521_counts.csv.gz'
    meta_file   = data_dir / 'GSE221521_metadata.csv'

    if matrix_file.exists() and meta_file.exists():
        log.info("  Already downloaded — skipping")
        return data_dir

    # Download count matrix
    urls = [
        # Try multiple possible file names
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221521/suppl/GSE221521_raw_counts.csv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221521/suppl/GSE221521_count_matrix.csv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221521/suppl/GSE221521_RAW.tar",
    ]

    downloaded = False
    for url in urls:
        fname = data_dir / url.split('/')[-1]
        log.info(f"  Trying: {url}")
        ret = subprocess.run(
            ['wget', '-q', '--tries=3', '-O', str(fname), url],
            capture_output=True
        )
        if ret.returncode == 0 and fname.stat().st_size > 1000:
            log.info(f"  Downloaded: {fname.name} ({fname.stat().st_size/1e6:.1f} MB)")
            downloaded = True
            break
        else:
            if fname.exists():
                fname.unlink()

    if not downloaded:
        log.error("  Direct download failed. Trying GEOparse...")
        try:
            import GEOparse
            gse = GEOparse.get_GEO("GSE221521", destdir=str(data_dir), silent=True)
            log.info("  GEOparse fetch complete")
        except Exception as e:
            log.error(f"  GEOparse also failed: {e}")
            log.error("  Manual download required:")
            log.error("  1. Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221521")
            log.error("  2. Download supplementary files")
            log.error(f"  3. Place in: {data_dir}/")
            sys.exit(1)

    # Download series matrix for metadata
    series_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221521/matrix/GSE221521_series_matrix.txt.gz"
    series_file = data_dir / "GSE221521_series_matrix.txt.gz"
    if not series_file.exists():
        log.info("  Downloading series matrix (metadata)...")
        subprocess.run(
            ['wget', '-q', '--tries=3', '-O', str(series_file), series_url],
            capture_output=True
        )

    return data_dir


# ════════════════════════════════════════════════════════════════════════════
# STEP 2 — LOAD AND PARSE DATA
# ════════════════════════════════════════════════════════════════════════════
def load_gse221521(data_dir: Path):
    """
    Load count matrix and metadata from GSE221521.
    Returns: counts_df (genes x samples), meta_df (samples x metadata)
    """
    log.info("Step 2: Loading GSE221521 data...")

    # Find count matrix file — handles Excel, CSV, TSV formats
    # GSE221521 uses: GSE221521_gene_expression.xls.gz (Excel with FPKM + counts)
    all_files = [f for f in data_dir.iterdir()
                 if 'series_matrix' not in f.name and f.name != 'GSE221521_family.soft.gz']

    if not all_files:
        log.error(f"  No data files found in {data_dir}")
        sys.exit(1)

    log.info(f"  Files found: {[f.name for f in all_files]}")

    # Pick the data file
    count_file = max(all_files, key=lambda f: f.stat().st_size)
    log.info(f"  Loading: {count_file.name} ({count_file.stat().st_size/1e6:.1f} MB)")

    counts_df = None

    # GSE221521: file extension is .xls but it's actually TSV (tab-separated)
    # Content starts with: gene_id\t... — pandas read_csv with sep='\t'
    if '.xls' in count_file.name or counts_df is None:
        log.info("  Reading as TSV (GEO files often have misleading .xls extension)...")
        try:
            # Read FULL TSV first (keep all columns including gene_name)
            full_df = pd.read_csv(count_file, index_col=0, sep='\t', compression='gzip')
            log.info(f"  Shape: {full_df.shape[0]:,} genes × {full_df.shape[1]:,} columns")
            log.info(f"  First 5 columns : {list(full_df.columns[:5])}")
            log.info(f"  Last 5 columns  : {list(full_df.columns[-5:])}")
            log.info(f"  First 3 index   : {list(full_df.index[:3])}")

            # ── Step A: remap index Ensembl → gene symbol BEFORE filtering ──
            if 'gene_name' in full_df.columns:
                log.info("  Remapping: Ensembl ID → gene symbol via gene_name column")
                gene_names = full_df['gene_name'].values
                full_df.index = gene_names
                full_df = full_df.drop(columns=['gene_name'])
                log.info(f"  Index now: {list(full_df.index[:3])}")
            else:
                log.warning("  No gene_name column found — index stays as Ensembl IDs")

            # Drop other annotation columns
            anno_cols = [c for c in full_df.columns
                         if c in ['description','gene_type','locus','chr',
                                  'start','end','strand','length','gene_id']]
            if anno_cols:
                full_df = full_df.drop(columns=anno_cols)

            # ── Step B: select count columns ──
            col_names  = [str(c).lower() for c in full_df.columns]
            count_cols = [c for c, cn in zip(full_df.columns, col_names) if 'count' in cn]
            fpkm_cols  = [c for c, cn in zip(full_df.columns, col_names) if 'fpkm' in cn]

            if count_cols:
                log.info(f"  Selecting {len(count_cols)} count columns")
                counts_df = full_df[count_cols].copy()
            elif fpkm_cols:
                non_fpkm = [c for c in full_df.columns if c not in fpkm_cols]
                numeric  = [c for c in non_fpkm
                            if pd.api.types.is_numeric_dtype(full_df[c])]
                log.info(f"  Selecting {len(numeric)} non-FPKM numeric columns")
                counts_df = full_df[numeric].copy() if numeric else full_df.copy()
            else:
                counts_df = full_df.copy()

            # ── Step C: strip _count suffix from column names → clean sample IDs ──
            counts_df.columns = [str(c).replace('_count','').replace('_FPKM','')
                                  for c in counts_df.columns]

            # Remove duplicate gene names (keep first)
            counts_df = counts_df[~counts_df.index.duplicated(keep='first')]

            log.info(f"  Final: {counts_df.shape[0]:,} genes × {counts_df.shape[1]:,} samples")
            log.info(f"  Sample IDs (first 3): {list(counts_df.columns[:3])}")
            log.info(f"  Gene names (first 3): {list(counts_df.index[:3])}")

            # Verify trajectory lncRNAs now findable
            traj_present = [g for g in TRAJ_GENES if g in counts_df.index]
            log.info(f"  Trajectory lncRNAs found in bulk: {len(traj_present)}/{len(TRAJ_GENES)}")
            if traj_present:
                log.info(f"  Examples: {traj_present[:5]}")

        except Exception as e:
            log.error(f"  TSV load failed: {e}")
            import traceback
            log.error(traceback.format_exc())
            sys.exit(1)

    # Other formats
    if counts_df is None:
        try:
            counts_df = pd.read_csv(count_file, index_col=0, compression='gzip')
        except Exception:
            try:
                counts_df = pd.read_csv(count_file, index_col=0, sep='\t', compression='gzip')
            except Exception as e:
                log.error(f"  Failed to load count matrix: {e}")
                sys.exit(1)

    log.info(f"  Count matrix: {counts_df.shape[0]:,} genes × {counts_df.shape[1]:,} samples")

    # Parse series matrix for metadata
    series_file = data_dir / "GSE221521_series_matrix.txt.gz"
    meta_df = _parse_series_matrix(series_file) if series_file.exists() else None

    if meta_df is None or len(meta_df) == 0:
        log.warning("  Could not parse metadata from series matrix")
        meta_df = pd.DataFrame(
            {'condition': ['Unknown'] * counts_df.shape[1]},
            index=counts_df.columns
        )
    else:
        log.info(f"  Raw metadata: {len(meta_df)} samples")
        log.info(f"  Conditions raw: {meta_df['condition'].value_counts().to_dict()}")

    # ── Re-align metadata to count matrix column names ──
    # Count columns after stripping _count suffix are sample IDs like RNA1, R_JS001
    # Series matrix sample titles contain condition info
    # Try to match by parsing the soft/series file more carefully
    series_file = data_dir / "GSE221521_series_matrix.txt.gz"
    if series_file.exists():
        meta_df = _parse_series_matrix_detailed(series_file, list(counts_df.columns))

    log.info(f"  Metadata: {len(meta_df)} samples")
    log.info(f"  Conditions: {meta_df['condition'].value_counts().to_dict()}")

    return counts_df, meta_df


def _parse_series_matrix(series_file: Path) -> pd.DataFrame:
    """Extract sample metadata from GEO series matrix."""
    try:
        with gzip.open(series_file, 'rt') as f:
            lines = f.readlines()
    except Exception:
        return None

    samples = []
    titles  = []
    chars   = {}

    for line in lines:
        line = line.strip()
        if line.startswith('!Sample_geo_accession'):
            samples = line.split('\t')[1:]
            samples = [s.strip('"') for s in samples]
        elif line.startswith('!Sample_title'):
            titles = line.split('\t')[1:]
            titles = [t.strip('"') for t in titles]
        elif line.startswith('!Sample_characteristics_ch1'):
            parts = line.split('\t')[1:]
            parts = [p.strip('"') for p in parts]
            for i, p in enumerate(parts):
                if ':' in p:
                    key, val = p.split(':', 1)
                    key = key.strip().lower().replace(' ', '_')
                    if i not in chars:
                        chars[i] = {}
                    chars[i][key] = val.strip()

    if not samples:
        return None

    meta = pd.DataFrame(index=samples)
    if titles:
        meta['title'] = titles[:len(samples)]

    # Determine condition from title or characteristics
    conditions = []
    for i, s in enumerate(samples):
        title = titles[i].lower() if titles and i < len(titles) else ''
        char_dict = chars.get(i, {})
        char_str  = ' '.join(char_dict.values()).lower()

        combined = title + ' ' + char_str
        if any(x in combined for x in ['type 2', 't2d', 'diabete', 'dm2', 't2dm']):
            if any(x in combined for x in ['pre', 'prediab', 'igt', 'ifg', 'impaired']):
                conditions.append('PreDM')
            else:
                conditions.append('T2D')
        elif any(x in combined for x in ['healthy', 'control', 'normal', 'ctrl', 'nd ']):
            conditions.append('Control')
        elif any(x in combined for x in ['pre', 'prediab', 'igt', 'ifg', 'impaired']):
            conditions.append('PreDM')
        else:
            conditions.append('Unknown')

    meta['condition'] = conditions

    # Align with count matrix columns if needed
    return meta


# ════════════════════════════════════════════════════════════════════════════
# STEP 3 — DIFFERENTIAL EXPRESSION (pydeseq2)
# ════════════════════════════════════════════════════════════════════════════

def _parse_series_matrix_detailed(series_file: Path, count_sample_ids: list) -> pd.DataFrame:
    """
    Parse GSE221521 series matrix.
    Title format: "leukocytes, DM group RNA1"
                  "leukocytes, DR group RNA102"
                  "leukocytes, Control group RNA109"
    Count columns: RNA1, RNA102, RNA109 (after stripping _count suffix)
    Mapping: DM = T2D | DR = PreDM | Control = Control
    """
    try:
        with gzip.open(series_file, 'rt') as f:
            lines = f.readlines()
    except Exception:
        log.warning("  Could not read series matrix")
        return pd.DataFrame({'condition': ['Unknown'] * len(count_sample_ids)},
                            index=count_sample_ids)

    titles = []
    for line in lines:
        line = line.strip()
        if line.startswith('!Sample_title'):
            titles = [t.strip('"') for t in line.split('\t')[1:]]
            break

    if not titles:
        return pd.DataFrame({'condition': ['Unknown'] * len(count_sample_ids)},
                            index=count_sample_ids)

    # Build mapping: sample_id (e.g. RNA1) → condition
    # Title: "leukocytes, DM group RNA1" → RNA1 → T2D
    sample_to_condition = {}
    for title in titles:
        # Extract sample ID (last token after last space)
        parts  = title.strip().split()
        sid    = parts[-1] if parts else ''
        t_low  = title.lower()
        if ' dm ' in t_low or t_low.startswith('dm') or ', dm ' in t_low:
            cond = 'T2D'
        elif ' dr ' in t_low or t_low.startswith('dr') or ', dr ' in t_low:
            cond = 'PreDM'
        elif 'control' in t_low:
            cond = 'Control'
        else:
            cond = 'Unknown'
        sample_to_condition[sid] = cond

    log.info(f"  Title-based condition map built: "
             f"{pd.Series(list(sample_to_condition.values())).value_counts().to_dict()}")
    log.info(f"  Example mappings: {dict(list(sample_to_condition.items())[:4])}")

    # Map count column IDs to conditions
    # count_sample_ids are already stripped of _count suffix: RNA1, RNA102...
    conditions = []
    for col in count_sample_ids:
        cond = sample_to_condition.get(col, 'Unknown')
        conditions.append(cond)

    result = pd.DataFrame({'condition': conditions}, index=count_sample_ids)
    log.info(f"  Final condition counts: "
             f"{result['condition'].value_counts().to_dict()}")
    return result


def run_deseq2_validation(counts_df: pd.DataFrame, meta_df: pd.DataFrame,
                          outdir: Path) -> dict:
    """
    Run DESeq2 for T2D vs Control and Pre-DM vs Control.
    Focus on the 24 trajectory lncRNAs.
    """
    log.info("Step 3: DESeq2 differential expression...")

    # Filter to trajectory lncRNAs present in bulk data
    present = [g for g in TRAJ_GENES if g in counts_df.index]
    missing = [g for g in TRAJ_GENES if g not in counts_df.index]
    log.info(f"  Trajectory lncRNAs in bulk data: {len(present)}/{len(TRAJ_GENES)}")
    if missing:
        log.info(f"  Missing from bulk: {missing}")

    if len(present) == 0:
        log.error("  No trajectory lncRNAs found in bulk count matrix")
        log.error("  Check gene naming — bulk may use Ensembl IDs")
        # Try to match by partial name
        bulk_genes = counts_df.index.tolist()
        log.info("  First 10 gene names in bulk:")
        for g in bulk_genes[:10]:
            log.info(f"    {g}")
        return {}

    results = {}
    contrasts = [
        ('T2D',   'Control'),
        ('PreDM', 'Control'),
    ]

    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds  import DeseqStats
        use_pydeseq2 = True
    except ImportError:
        log.warning("  pydeseq2 not installed — using edgeR-equivalent (NB GLM) fallback")
        use_pydeseq2 = False

    for test_group, ref_group in contrasts:
        log.info(f"  Contrast: {test_group} vs {ref_group}")

        # Get samples for this contrast
        mask = meta_df['condition'].isin([test_group, ref_group])
        if mask.sum() < 4:
            log.warning(f"  Insufficient samples for {test_group} vs {ref_group} — skipping")
            continue

        sub_meta   = meta_df[mask].copy()
        sub_counts = counts_df.loc[present, sub_meta.index].T.astype(int)
        sub_meta['condition_bin'] = (sub_meta['condition'] == test_group).astype(int)

        # Filter low-count genes
        keep = (sub_counts > 5).mean(axis=0) >= 0.2
        sub_counts = sub_counts.loc[:, keep]
        log.info(f"  Samples: {len(sub_meta)}  |  Genes: {len(sub_counts.columns)}")

        if use_pydeseq2:
            try:
                dds = DeseqDataSet(
                    counts   = sub_counts,
                    metadata = sub_meta[['condition_bin']],
                    design_factors = ['condition_bin'],
                    refit_cooks    = True,
                    quiet          = True,
                )
                dds.deseq2()
                stat = DeseqStats(dds, contrast=['condition_bin','1','0'], quiet=True)
                stat.summary()
                res = stat.results_df.copy()
                res['gene']     = res.index
                res['contrast'] = f"{test_group}_vs_{ref_group}"
                results[f"{test_group}_vs_{ref_group}"] = res
                sig = res[(res['padj'] < 0.05) & (res['log2FoldChange'].abs() > 0.3)]
                log.info(f"  Significant: {len(sig)}")
            except Exception as e:
                log.error(f"  DESeq2 failed: {e} — falling back to t-test")
                use_pydeseq2 = False

        if not use_pydeseq2:
            # Fallback: Welch t-test on log-normalised counts
            rows = []
            log_counts = np.log2(sub_counts + 1)
            g1_idx = sub_meta[sub_meta['condition'] == ref_group].index
            g2_idx = sub_meta[sub_meta['condition'] == test_group].index
            g1 = log_counts.loc[g1_idx]
            g2 = log_counts.loc[g2_idx]

            for gene in sub_counts.columns:
                v1 = g1[gene].values
                v2 = g2[gene].values
                if v1.std() + v2.std() < 1e-6:
                    continue
                t, pval = stats.ttest_ind(v2, v1, equal_var=False)
                lfc = v2.mean() - v1.mean()
                rows.append({'gene': gene, 'log2FoldChange': lfc, 'pvalue': pval,
                             'contrast': f"{test_group}_vs_{ref_group}"})

            if rows:
                res = pd.DataFrame(rows)
                _, padj, _, _ = multipletests(res['pvalue'], method='fdr_bh')
                res['padj'] = padj
                results[f"{test_group}_vs_{ref_group}"] = res
                sig = res[(res['padj'] < 0.05) & (res['log2FoldChange'].abs() > 0.3)]
                log.info(f"  Significant (t-test fallback): {len(sig)}")

    # Save results
    if results:
        all_res = pd.concat(results.values(), ignore_index=True)
        all_res.to_csv(outdir / 'bulk_validation_DE_results.csv', index=False)
        log.info(f"  Saved: bulk_validation_DE_results.csv")

    return results


# ════════════════════════════════════════════════════════════════════════════
# STEP 4 — CONCORDANCE ANALYSIS
# ════════════════════════════════════════════════════════════════════════════
def concordance_analysis(results: dict, outdir: Path) -> pd.DataFrame:
    """
    Check direction concordance between scRNA-seq trajectory and bulk DE.
    A lncRNA is concordant if:
      - ↑CD14 in trajectory AND upregulated in T2D vs Control in bulk
        (CD14 classical = inflammatory = T2D-associated)
      - ↑CD16 in trajectory AND downregulated in T2D vs Control in bulk
        (CD16 non-classical = patrolling = less inflammatory)
    """
    log.info("Step 4: Concordance analysis...")

    if not results:
        log.warning("  No DE results — skipping concordance")
        return pd.DataFrame()

    t2d_res = results.get('T2D_vs_Control', pd.DataFrame())
    if len(t2d_res) == 0:
        log.warning("  No T2D vs Control results")
        return pd.DataFrame()

    t2d_res = t2d_res.set_index('gene')

    concordance_rows = []
    for gene, (rho, direction) in TRAJECTORY_LNCRNAS.items():
        if gene not in t2d_res.index:
            continue

        bulk_lfc  = float(t2d_res.loc[gene, 'log2FoldChange'])
        bulk_padj = float(t2d_res.loc[gene, 'padj']) if 'padj' in t2d_res.columns else 1.0

        # Concordance logic:
        # ↑CD14 (high in inflammatory/classical) → expect higher in T2D → bulk LFC > 0
        # ↑CD16 (high in patrolling/non-classical) → expect lower in T2D → bulk LFC < 0
        if direction == '↑CD14':
            expected_bulk_direction = '+'   # should be upregulated in T2D
            concordant = bulk_lfc > 0
        else:  # ↑CD16
            expected_bulk_direction = '-'   # should be downregulated in T2D
            concordant = bulk_lfc < 0

        concordance_rows.append({
            'gene'             : gene,
            'scrna_rho'        : rho,
            'scrna_direction'  : direction,
            'bulk_lfc'         : bulk_lfc,
            'bulk_padj'        : bulk_padj,
            'expected_direction': expected_bulk_direction,
            'concordant'       : concordant,
            'bulk_sig'         : bulk_padj < 0.05,
        })

    conc_df = pd.DataFrame(concordance_rows)
    conc_df.to_csv(outdir / 'bulk_concordance.csv', index=False)

    n_conc       = conc_df['concordant'].sum()
    n_conc_sig   = (conc_df['concordant'] & conc_df['bulk_sig']).sum()
    n_total      = len(conc_df)
    pct          = n_conc / n_total * 100 if n_total > 0 else 0

    # Binomial test: is concordance > 50% (chance)?
    from scipy.stats import binomtest
    binom_p = binomtest(int(n_conc), n_total, 0.5, alternative='greater').pvalue

    log.info(f"  Concordant lncRNAs: {n_conc}/{n_total} ({pct:.0f}%)")
    log.info(f"  Concordant AND bulk-significant: {n_conc_sig}")
    log.info(f"  Binomial test (vs chance): p = {binom_p:.4f}")

    log.info("\n  Per-gene concordance:")
    log.info(f"  {'Gene':<20} {'scRNA rho':>10} {'Bulk LFC':>10} {'Bulk padj':>12} {'Concordant':>12}")
    log.info("  " + "-" * 65)
    for _, row in conc_df.sort_values('scrna_rho').iterrows():
        mark = '✓' if row['concordant'] else '✗'
        sig  = '*' if row['bulk_sig'] else ''
        log.info(f"  {row['gene']:<20} {row['scrna_rho']:>+10.3f} "
                 f"{row['bulk_lfc']:>+10.3f} {row['bulk_padj']:>12.4f} "
                 f"  {mark}{sig}")

    return conc_df


# ════════════════════════════════════════════════════════════════════════════
# STEP 5 — VALIDATION FIGURES
# ════════════════════════════════════════════════════════════════════════════
def plot_validation(conc_df: pd.DataFrame, results: dict, outdir: Path):
    """Generate Fig14: validation panel."""
    log.info("Step 5: Generating validation figures...")

    if len(conc_df) == 0:
        log.warning("  No concordance data — skipping plots")
        return

    fig = plt.figure(figsize=(16, 12))
    gs  = gridspec.GridSpec(2, 2, hspace=0.4, wspace=0.35)

    # ── Panel A: scRNA rho vs bulk LFC scatter ─────────────────────────────
    ax = fig.add_subplot(gs[0, 0])
    conc = conc_df[conc_df['concordant']]
    disc = conc_df[~conc_df['concordant']]

    ax.scatter(conc['scrna_rho'], conc['bulk_lfc'],
               c='#2ECC71', s=60, alpha=0.8, zorder=3, label='Concordant')
    ax.scatter(disc['scrna_rho'], disc['bulk_lfc'],
               c='#E74C3C', s=60, alpha=0.8, zorder=3, label='Discordant')

    # Label top genes
    for _, row in conc_df.nlargest(6, 'bulk_lfc').iterrows():
        ax.annotate(row['gene'],
                    xy=(row['scrna_rho'], row['bulk_lfc']),
                    xytext=(4, 4), textcoords='offset points', fontsize=7)
    for _, row in conc_df.nsmallest(4, 'bulk_lfc').iterrows():
        ax.annotate(row['gene'],
                    xy=(row['scrna_rho'], row['bulk_lfc']),
                    xytext=(4, -10), textcoords='offset points', fontsize=7)

    # Reference lines
    ax.axhline(0, color='gray', lw=0.8, ls='--', alpha=0.5)
    ax.axvline(0, color='gray', lw=0.8, ls='--', alpha=0.5)

    # Pearson correlation
    if len(conc_df) > 3:
        r, p = stats.pearsonr(conc_df['scrna_rho'], conc_df['bulk_lfc'])
        ax.text(0.05, 0.95, f'r = {r:.3f}\np = {p:.3f}',
                transform=ax.transAxes, va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    ax.set_xlabel('scRNA-seq Spearman ρ\n(trajectory correlation)', fontsize=10)
    ax.set_ylabel('Bulk log₂FC\n(T2D vs Control)', fontsize=10)
    ax.set_title('scRNA-seq Trajectory vs\nBulk T2D Differential Expression',
                 fontsize=10, fontweight='bold')
    ax.legend(fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)

    # ── Panel B: concordance bar ───────────────────────────────────────────
    ax = fig.add_subplot(gs[0, 1])
    n_conc = conc_df['concordant'].sum()
    n_disc = (~conc_df['concordant']).sum()
    n_conc_sig = (conc_df['concordant'] & conc_df['bulk_sig']).sum()

    bars = ax.bar(['Concordant', 'Discordant'],
                  [n_conc, n_disc],
                  color=['#2ECC71', '#E74C3C'], alpha=0.85, width=0.5)
    ax.bar(['Concordant'], [n_conc_sig],
           color='#27AE60', alpha=1.0, width=0.5,
           label=f'FDR<0.05 in bulk\n(n={n_conc_sig})')

    for bar, val in zip(bars, [n_conc, n_disc]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                str(val), ha='center', fontsize=12, fontweight='bold')

    pct = n_conc / (n_conc + n_disc) * 100 if (n_conc + n_disc) > 0 else 0
    ax.axhline(y=(n_conc + n_disc) / 2, color='gray',
               ls='--', lw=1.5, alpha=0.7, label='50% chance')
    ax.set_title(f'Direction Concordance\n({pct:.0f}% concordant)',
                 fontsize=10, fontweight='bold')
    ax.set_ylabel('Number of lncRNAs', fontsize=10)
    ax.legend(fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)

    # ── Panel C: lncRNA-level concordance heatmap ──────────────────────────
    ax = fig.add_subplot(gs[1, :])
    plot_df = conc_df.set_index('gene')[['scrna_rho', 'bulk_lfc']].copy()
    plot_df = plot_df.sort_values('scrna_rho')

    x = np.arange(len(plot_df))
    width = 0.35

    # Normalise both to [-1, 1] for visual comparison
    scrna_norm = plot_df['scrna_rho'] / plot_df['scrna_rho'].abs().max()
    bulk_norm  = plot_df['bulk_lfc']  / (plot_df['bulk_lfc'].abs().max() + 1e-6)

    bars1 = ax.bar(x - width/2, scrna_norm, width,
                   color=['#E74C3C' if v > 0 else '#3498DB' for v in scrna_norm],
                   alpha=0.8, label='scRNA-seq ρ (normalised)')
    bars2 = ax.bar(x + width/2, bulk_norm, width,
                   color=['#E74C3C' if v > 0 else '#3498DB' for v in bulk_norm],
                   alpha=0.4, label='Bulk LFC (normalised)', edgecolor='black', lw=0.5)

    # Mark concordant pairs
    for i, (_, row) in enumerate(plot_df.iterrows()):
        if conc_df[conc_df['gene'] == row.name]['concordant'].values[0]:
            ax.text(i, max(abs(scrna_norm.iloc[i]), abs(bulk_norm.iloc[i])) + 0.05,
                    '✓', ha='center', fontsize=8, color='#27AE60')

    ax.axhline(0, color='black', lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(plot_df.index, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Normalised effect size', fontsize=10)
    ax.set_title(
        'Per-Gene Concordance: scRNA-seq Trajectory vs Bulk T2D DE\n'
        '(✓ = concordant direction)',
        fontsize=10, fontweight='bold'
    )
    ax.legend(fontsize=9, loc='upper left')
    ax.spines[['top', 'right']].set_visible(False)

    fig.suptitle(
        'Bulk Validation of Trajectory-Associated lncRNAs\n'
        'GSE221521 (n=1,190 whole blood RNA-seq: T2D + Pre-DM + Control)',
        fontsize=12, fontweight='bold', y=1.01
    )

    out_path = outdir / 'Fig14_bulk_validation.png'
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    log.info(f"  Saved: Fig14_bulk_validation.png")


# ════════════════════════════════════════════════════════════════════════════
# STEP 6 — BOXPLOTS FOR KEY lncRNAs
# ════════════════════════════════════════════════════════════════════════════
def plot_key_lncrna_boxplots(counts_df, meta_df, outdir):
    """Boxplot expression of top lncRNAs across T2D / Pre-DM / Control."""
    log.info("Step 6: Key lncRNA expression boxplots...")

    key_genes = ['NEAT1', 'MALAT1', 'AC020916.1', 'AC020656.1',
                 'AC007952.4', 'LINC00861']
    present   = [g for g in key_genes if g in counts_df.index]
    if not present:
        log.warning("  Key genes not found in bulk — skipping boxplots")
        return

    order      = ['Control', 'PreDM', 'T2D']
    colours    = {'Control': '#3498DB', 'PreDM': '#F39C12', 'T2D': '#E74C3C'}
    log_counts = np.log2(counts_df.loc[present] + 1)

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle(
        'Key Trajectory lncRNA Expression in Bulk Blood\n'
        'GSE221521 (T2D, Pre-DM, Control)',
        fontsize=12, fontweight='bold'
    )

    for ax, gene in zip(axes.flatten(), present[:6]):
        data_by_group = []
        labels_shown  = []
        for cond in order:
            samples = meta_df[meta_df['condition'] == cond].index
            samples = [s for s in samples if s in log_counts.columns]
            if len(samples) < 3:
                continue
            data_by_group.append(log_counts.loc[gene, samples].values)
            labels_shown.append(f"{cond}\n(n={len(samples)})")

        if len(data_by_group) < 2:
            ax.text(0.5, 0.5, f'{gene}\nInsufficient data',
                    ha='center', va='center', transform=ax.transAxes)
            continue

        vp = ax.violinplot(data_by_group, positions=range(len(data_by_group)),
                           showmedians=True, showextrema=False)
        for i, (pc, lbl) in enumerate(zip(vp['bodies'], labels_shown)):
            cond_key = lbl.split('\n')[0]
            pc.set_facecolor(colours.get(cond_key, '#9E9E9E'))
            pc.set_alpha(0.75)

        # Add p-value annotations (Wilcoxon)
        if len(data_by_group) >= 2:
            try:
                _, p = stats.mannwhitneyu(data_by_group[0], data_by_group[-1],
                                          alternative='two-sided')
                p_str = f'p={p:.3f}' if p > 0.001 else f'p={p:.2e}'
                y_max = max([d.max() for d in data_by_group]) + 0.3
                ax.plot([0, len(data_by_group)-1], [y_max, y_max],
                        color='black', lw=1.2)
                ax.text((len(data_by_group)-1)/2, y_max + 0.05,
                        p_str, ha='center', fontsize=8)
            except Exception:
                pass

        ax.set_xticks(range(len(labels_shown)))
        ax.set_xticklabels(labels_shown, fontsize=8)
        ax.set_ylabel('log₂(counts+1)', fontsize=9)
        direction = TRAJECTORY_LNCRNAS.get(gene, (0, ''))[1]
        ax.set_title(f'{gene}\n({direction})', fontsize=10, fontweight='bold')
        ax.spines[['top', 'right']].set_visible(False)

    plt.tight_layout()
    out_path = outdir / 'Fig15_key_lncrna_boxplots.png'
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    log.info(f"  Saved: Fig15_key_lncrna_boxplots.png")


# ════════════════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(description='Phase 4: Bulk Validation')
    parser.add_argument('--outdir',    default='./phase4_results')
    parser.add_argument('--data_dir',  default=None,
                        help='Pre-downloaded GSE221521 data directory (skips download)')
    parser.add_argument('--counts',    default=None,
                        help='Path to counts CSV/TSV directly')
    parser.add_argument('--meta',      default=None,
                        help='Path to metadata CSV directly')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("=" * 65)
    log.info("PHASE 4 — BULK VALIDATION")
    log.info(f"Start: {datetime.now()}")
    log.info("=" * 65)

    # Step 1: Download or use provided data
    if args.counts and Path(args.counts).exists():
        log.info(f"Using provided counts: {args.counts}")
        counts_df = pd.read_csv(args.counts, index_col=0)
        if args.meta and Path(args.meta).exists():
            meta_df = pd.read_csv(args.meta, index_col=0)
        else:
            meta_df = pd.DataFrame(
                {'condition': ['Unknown'] * counts_df.shape[1]},
                index=counts_df.columns
            )
    elif args.data_dir and Path(args.data_dir).exists():
        counts_df, meta_df = load_gse221521(Path(args.data_dir))
    else:
        data_dir  = download_gse221521(outdir)
        counts_df, meta_df = load_gse221521(data_dir)

    log.info(f"Data shape: {counts_df.shape}")
    log.info(f"Conditions: {meta_df['condition'].value_counts().to_dict()}")

    # Step 2: Run DESeq2
    results = run_deseq2_validation(counts_df, meta_df, outdir)

    # Step 3: Concordance
    conc_df = concordance_analysis(results, outdir)

    # Step 4: Figures
    if len(conc_df) > 0:
        plot_validation(conc_df, results, outdir)
        plot_key_lncrna_boxplots(counts_df, meta_df, outdir)

    # Summary
    if len(conc_df) > 0:
        n_conc = conc_df['concordant'].sum()
        n_sig  = (conc_df['concordant'] & conc_df['bulk_sig']).sum()
        pct    = n_conc / len(conc_df) * 100
        summary = [
            "", "=" * 65,
            "PHASE 4 — BULK VALIDATION SUMMARY",
            f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
            "=" * 65,
            f"Bulk dataset: GSE221521 (n={counts_df.shape[1]})",
            f"Trajectory lncRNAs tested: {len(conc_df)}/{len(TRAJ_GENES)}",
            f"Concordant direction: {n_conc}/{len(conc_df)} ({pct:.0f}%)",
            f"Concordant AND bulk-significant: {n_sig}",
            "",
            "CONCORDANT lncRNAs:",
        ]
        for _, row in conc_df[conc_df['concordant']].sort_values('bulk_padj').iterrows():
            sig = '*' if row['bulk_sig'] else ''
            summary.append(f"  {row['gene']:<20} scRNA ρ={row['scrna_rho']:+.3f}  "
                           f"bulk LFC={row['bulk_lfc']:+.3f}  {sig}")
        summary += ["", "FILES:", "  bulk_validation_DE_results.csv",
                    "  bulk_concordance.csv", "  Fig14_bulk_validation.png",
                    "  Fig15_key_lncrna_boxplots.png", "=" * 65]
        report = '\n'.join(summary)
        print(report)
        with open(outdir / 'phase4_summary.txt', 'w') as f:
            f.write(report)

    log.info(f"\nPhase 4 complete. Log: {log_file}")


if __name__ == '__main__':
    main()
