# scRNA-seq lncRNA Analysis of T2D Monocytes

<div align="center">

![Python](https://img.shields.io/badge/Python-3.10-blue?logo=python)
![Scanpy](https://img.shields.io/badge/Scanpy-1.11.5-green)
![License](https://img.shields.io/badge/License-MIT-yellow)
![Status](https://img.shields.io/badge/Status-Preprint-orange)
![GEO](https://img.shields.io/badge/GEO-GSE268210-red)

**Single-Cell Transcriptomic Dissection of lncRNA Regulatory Programmes Along the Monocyte Activation Continuum in Type 2 Diabetes**

*Usama Manzoor — JSMU Diagnostic Laboratory & Blood Bank, Karachi, Pakistan*

</div>

---

## Overview

This repository contains the complete analytical pipeline for identifying and validating **long non-coding RNA (lncRNA) regulatory programmes** along the CD14 Classical → Intermediate → CD16 Non-Classical monocyte activation continuum in type 2 diabetes (T2D), using single-cell RNA sequencing.

### Key Findings

| Finding | Result |
|---------|--------|
| Trajectory-associated lncRNAs | **24** (Spearman \|ρ\| ≥ 0.15, FDR < 0.05) |
| Strongest CD14-enriched candidate | **AC020656.1** (ρ = −0.328) |
| Bulk validation concordance | **r = −0.630, p = 0.003** |
| AC020656.1 T2D upregulation | **log₂FC = +0.75, FDR = 0.017** |
| Progressive T2D gradient | **p = 8.2 × 10⁻⁴** (Control → Pre-DM → T2D) |
| Cell-type-specific lncRNAs | **506** detectable across 6 PBMC subsets |

---

## Repository Structure

```
.
├── README.md
├── LICENSE
├── CITATION.cff
│
├── envs/
│   ├── environment.yml          # Conda environment (Python 3.10)
│   └── requirements.txt         # pip requirements
│
├── data/
│   ├── raw/                     # Raw data (not tracked — see Data Availability)
│   │   └── .gitkeep
│   ├── processed/               # Intermediate h5ad checkpoints
│   │   └── .gitkeep
│   └── bulk_validation/         # GSE221521 bulk data
│       └── .gitkeep
│
├── scripts/
│   ├── preprocessing/
│   │   ├── 01_qc_preprocessing.py     # Phase 1: QC, normalisation, batch correction
│   │   └── 02_cell_annotation.py      # Phase 2: Cell type annotation + UMAP
│   │
│   ├── analysis/
│   │   ├── 03_monocyte_trajectory.py  # Phase 3c: Pseudotime + lncRNA trajectory
│   │   └── 04_bulk_validation.py      # Phase 4: Bulk RNA-seq validation
│   │
│   ├── visualization/
│   │   └── 05_figures.py              # Regenerate all publication figures
│   │
│   └── utils/
│       ├── lncrna_annotation.py       # lncRNA detection + GENCODE mapping
│       └── pseudobulk.py              # Pseudobulk aggregation utilities
│
├── results/
│   ├── figures/                 # Publication figures (PNG, 150+ dpi)
│   ├── tables/                  # CSV result tables
│   ├── phase3c/                 # Trajectory analysis outputs
│   └── phase4/                  # Bulk validation outputs
│
└── docs/
    ├── methods_detail.md        # Extended methods notes
    └── figure_guide.md          # Figure descriptions + interpretation
```

---

## Datasets

| Dataset | Description | n | Access |
|---------|-------------|---|--------|
| [GSE268210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268210) | T2D PBMC scRNA-seq (10x 5' v2) | 9 patients | Public |
| [GSE221521](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221521) | Whole blood RNA-seq (T2D/Pre-DM/Control) | 193 samples | Public |

---

## Quick Start

### 1. Clone and set up environment

```bash
git clone https://github.com/usamamanzoor1121-pixel/lncrna-t2d-monocyte.git
cd lncrna-t2d-monocyte

# Create conda environment
conda env create -f envs/environment.yml
conda activate scrna
```

### 2. Download data

```bash
# scRNA-seq data (GSE268210, ~1.4 GB)
cd data/raw
aria2c -x 8 -s 8 \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268210/suppl/GSE268210_RAW.tar"
tar -xf GSE268210_RAW.tar

# Bulk validation data (GSE221521, ~56 MB)
cd ../bulk_validation
aria2c -x 8 -s 8 \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221521/suppl/GSE221521_gene_expression.xls.gz"
```

### 3. Run pipeline

```bash
# Step 1: QC + pre-processing (45-90 min)
python scripts/preprocessing/01_qc_preprocessing.py \
    --data_dir data/raw \
    --output_dir data/processed \
    --skip_samples T2D_01

# Step 2: Cell type annotation + visualisation (10-15 min)
python scripts/preprocessing/02_cell_annotation.py \
    --h5ad data/processed/GSE268210_phase1_full.h5ad \
    --outdir results/figures

# Step 3: Monocyte trajectory + lncRNA analysis (25-35 min)
python scripts/analysis/03_monocyte_trajectory.py \
    --mono_h5ad results/figures/GSE268210_monocytes_annotated.h5ad \
    --full_h5ad data/processed/GSE268210_phase1_full.h5ad \
    --outdir results/phase3c

# Step 4: Bulk validation (45-90 min)
python scripts/analysis/04_bulk_validation.py \
    --data_dir data/bulk_validation \
    --outdir results/phase4
```

### 4. Regenerate figures only

If you have the processed h5ad files, regenerate all publication figures:

```bash
python scripts/visualization/05_figures.py \
    --full_h5ad data/processed/GSE268210_phase1_full.h5ad \
    --mono_h5ad results/phase3c/GSE268210_monocytes_final.h5ad \
    --traj_csv  results/phase3c/trajectory_lncrna_final.csv \
    --bulk_csv  results/phase4/bulk_concordance.csv \
    --outdir    results/figures
```

---

## Key Results

### Monocyte Subtype Proportions

| Subtype | Clusters (leiden_0.3) | n | % | Literature range |
|---------|----------------------|---|---|-----------------|
| CD14 Classical | 5, 6 | 35,805 | 81.3% | 75–85% |
| Intermediate | 1, 2 | 2,588 | 5.9% | 2–10% |
| CD16 Non-Classical | 7 | 5,664 | 12.9% | 5–15% |

### Top Trajectory lncRNAs

| Gene | ρ | FDR | Direction | Bulk LFC | Bulk FDR |
|------|---|-----|-----------|----------|----------|
| AC020916.1 | −0.343 | <2.2e-308 | ↑ CD14 | +0.21 | 0.70 |
| AC020656.1 | −0.328 | <2.2e-308 | ↑ CD14 | +0.75 | **0.017** |
| NEAT1 | −0.272 | <2.2e-308 | ↑ CD14 | +0.26 | 0.16 |
| MALAT1 | +0.379 | <2.2e-308 | ↑ CD16 | +0.10 | 0.71 |
| LINC00861 | +0.237 | <2.2e-308 | ↑ CD16 | +0.11 | 0.70 |

Full table: [`results/tables/trajectory_lncrna_final.csv`](results/tables/trajectory_lncrna_final.csv)

---

## Figures

| Figure | File | Description |
|--------|------|-------------|
| Fig 1 | `Fig1_QC_violins.png` | QC metrics per patient |
| Fig 2 | `Fig2_UMAP_overview.png` | PBMC atlas + batch correction |
| Fig 3 | `Fig3_marker_dotplot.png` | Canonical marker validation |
| Fig 4 | `Fig9_pseudotime_final.png` | Monocyte trajectory |
| Fig 5 | `Fig10_trajectory_heatmap_final.png` | Trajectory lncRNA heatmap |
| Fig 6 | `Fig10b_trajectory_lollipop_final.png` | lncRNA-pseudotime correlations |
| Fig 7 | `Fig13_lncrna_celltype_specificity.png` | PBMC lncRNA atlas |
| Fig 8 | `Fig14_bulk_validation.png` | Bulk validation concordance |
| Fig 9 | `Fig15_key_lncrna_boxplots.png` | AC020656.1 + NEAT1 expression |

---

## Computational Environment

| Tool | Version | Purpose |
|------|---------|---------|
| Python | 3.10.20 | Core language |
| Scanpy | 1.11.5 | scRNA-seq analysis |
| AnnData | 0.11.4 | Data structures |
| Harmony | 0.0.9 | Batch correction |
| Scrublet (via Scanpy) | — | Doublet detection |
| NumPy | 2.4.4 | Numerical computing |
| SciPy | 1.17.1 | Statistical tests |
| pandas | 2.3.3 | Data manipulation |
| matplotlib | 3.8 | Visualisation |
| seaborn | 0.13 | Statistical graphics |
| pydeseq2 | 0.4.9 | Differential expression |
| statsmodels | 0.14 | Multiple testing correction |

---

## Citation

If you use this code or findings, please cite:

```bibtex
@article{manzoor2026lncrna,
  title   = {Single-Cell Transcriptomic Dissection of lncRNA Regulatory Programmes
             Along the Monocyte Activation Continuum in Type 2 Diabetes},
  author  = {Manzoor, Usama},
  journal = {Briefings in Bioinformatics},
  year    = {2026},
  note    = {Under review}
}
```

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

---

## Contact

**Usama Manzoor**  
JSMU Diagnostic Laboratory & Blood Bank  
Jinnah Sindh Medical University, Karachi, Pakistan  
📧 usama.manzoor1121@gmail.com  
🐙 [@usamamanzoor1121-pixel](https://github.com/usamamanzoor1121-pixel)
