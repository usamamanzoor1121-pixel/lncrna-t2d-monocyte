# scRNA-seq lncRNA Analysis of T2D Monocytes

<div align="center">

![Python](https://img.shields.io/badge/Python-3.10-blue?logo=python)
![Scanpy](https://img.shields.io/badge/Scanpy-1.11.5-green)
![License](https://img.shields.io/badge/License-MIT-yellow)
![Status](https://img.shields.io/badge/Status-Preprint-orange)
![GEO](https://img.shields.io/badge/GEO-GSE268210-red)

**Donor-Aware Single-Cell Analysis Identifies AC020656.1, a Locus Antisense to LYZ, as a Candidate Monocyte-State and Type 2 Diabetes-Associated Transcript**

*Usama Manzoor — JSMU Diagnostic Laboratory & Blood Bank, Karachi, Pakistan*

</div>

---

## Overview

This repository contains the complete analytical pipeline for identifying and validating **long non-coding RNA (lncRNA) regulatory loci** across monocyte states in type 2 diabetes (T2D), using single-cell RNA sequencing (GSE268210) with independent bulk RNA-seq validation (GSE221521).

lncRNA identification uses a direct GENCODE v32 cross-reference by Ensembl ID. Trajectory-lncRNA significance is assessed with a donor-aware statistic — per-donor Spearman correlations combined via a DerSimonian-Laird random-effects Fisher-z meta-analysis with a Hartung-Knapp-Sidik-Jonkman small-sample correction — because cells from the same donor are not independent observations. Every major statistic is independently cross-validated in a second, separately-implemented R pipeline (`scripts/r_analysis/`).

### Key Findings

| Finding | Result |
|---------|--------|
| Trajectory-associated lncRNAs (donor-aware, GENCODE-confirmed) | **21 of 22** GENCODE-confirmed candidates (see `results/tables/trajectory_donor_aware_results.csv`); a broader unbiased screen of 481 loci finds **28** significant loci total |
| Strongest CD14-enriched candidate | **AC020656.1** (donor-aware ρ = −0.335, p = 2.7×10⁻⁸; unanimous direction in 9/9 donors; leave-one-donor-out robust; cross-validated to 4 decimal places in an independent R implementation) |
| Bulk validation — continuous metric | **r = −0.647, p = 0.0037** (18 GENCODE-confirmed loci) |
| Bulk validation — categorical concordance | 55.6% (10/18); binomial p = 0.41 — not distinguishable from chance on its own; the continuous metric above is what carries statistical weight |
| AC020656.1 T2D upregulation | log₂FC = +0.75–0.77 (two independent statistical implementations), FDR = 0.017 (targeted panel) / 0.28 (genome-wide) |
| Progressive T2D gradient | p = 8.2 × 10⁻⁴ (Control → Pre-DM → T2D) |
| Monocyte subtype topology (PAGA) | CD14↔CD16 connectivity (0.044) > 3x stronger than CD14↔Intermediate (0.012) — a shared-origin, divergent-branch structure, not a linear continuum |
| **AC020656.1–LYZ genomic relationship** | AC020656.1's entire annotated locus (both exons) lies within LYZ's terminal exon, antisense strand. No BAM/FASTQ data exist for this dataset, so strand-of-origin cannot be verified; AC020656.1 is reported as a **candidate** disease-associated locus, not a confirmed independently-regulated lncRNA. See `scripts/r_analysis/ac0206561_lyz/` and the manuscript's Discussion/Limitations for the full investigation. |

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
│   ├── processed/               # Intermediate h5ad checkpoints
│   ├── bulk_validation/         # GSE221521 bulk data
│   └── reference/                # GENCODE v32 GTF
│
├── scripts/
│   ├── preprocessing/
│   │   ├── 01_qc_preprocessing.py       # QC, normalisation, batch correction
│   │   └── 02_cell_annotation.py        # Cell type annotation + UMAP
│   │
│   ├── annotation/
│   │   └── 03_gencode_lncrna_annotation.py  # GENCODE v32 lncRNA cross-reference
│   │
│   ├── analysis/
│   │   ├── 04_monocyte_trajectory.py    # Monocyte subtype annotation + pseudotime
│   │   ├── 05_donor_aware_statistics.py # Donor-aware trajectory-lncRNA significance
│   │   ├── 06_paga_topology.py          # Branch-aware topology (PAGA)
│   │   ├── 07_ac020656_validation.py    # AC020656.1 robustness analysis
│   │   └── 08_bulk_validation.py        # Bulk RNA-seq validation (GSE221521)
│   │
│   ├── visualization/
│   │   └── 09_figures.py                # Regenerate final publication figures
│   │
│   ├── utils/
│   │   ├── lncrna_annotation.py         # Pattern-based lncRNA screen (fast first-pass heuristic)
│   │   └── pseudobulk.py                # Pseudobulk aggregation utilities
│   │
│   └── r_analysis/                      # Independent R validation layer
│       ├── 00_export_for_R.py / 01_export_bulk_for_R.py
│       ├── 02_pseudobulk_qc.R
│       ├── 03_differential_state_edgeR.R
│       ├── 04_ac020656_validation.R
│       ├── 05_neat1_validation.R
│       ├── 06_bulk_limma.R
│       ├── 07_bulk_deconvolution_nnls.R
│       ├── 08_pathway_fgsea.R
│       ├── 09_genomic_characterization.R
│       ├── 10_coexpression_screen.R
│       ├── run_all.R
│       └── ac0206561_lyz/               # AC020656.1–LYZ locus investigation
│           ├── 01_locus_overlap.R
│           ├── 02_locus_map_figure.R
│           ├── 03_correlation_ratio_models.R
│           ├── 04_disease_lyz_adjustment.R
│           └── 05_cell_level_independence.py
│
├── results/
│   ├── figures/                 # Publication figures (PNG, 150+ dpi)
│   ├── tables/                  # CSV result tables
│   └── supplementary_figures/   # R-validation and locus-investigation figures
│
├── manuscript/
│   ├── manuscript.md            # Authoritative manuscript text (source of truth)
│   └── manuscript.docx          # Submission/editing version, with figures and tables
│                                 # embedded; regenerate after editing manuscript.md via:
│                                 #   python3 scripts/utils/build_manuscript_docx.py
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

# GENCODE v32 reference (~42 MB)
cd ../reference
wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz"
```

### 3. Run the pipeline

```bash
# Step 1: QC + pre-processing (45-90 min)
python scripts/preprocessing/01_qc_preprocessing.py \
    --data_dir data/raw --output_dir data/processed --skip_samples T2D_01

# Step 2: Cell type annotation + monocyte subclustering (10-15 min)
python scripts/preprocessing/02_cell_annotation.py \
    --h5ad data/processed/GSE268210_phase1_full.h5ad --outdir results/figures

# Step 3: GENCODE v32 lncRNA annotation (~2 min)
python scripts/annotation/03_gencode_lncrna_annotation.py \
    --gtf data/reference/gencode.v32.annotation.gtf.gz \
    --h5ad data/processed/GSE268210_monocytes_final.h5ad \
    --outdir results/tables

# Step 4: Monocyte subtype annotation + pseudotime (15-20 min)
python scripts/analysis/04_monocyte_trajectory.py \
    --mono_h5ad results/figures/GSE268210_monocytes_annotated.h5ad \
    --full_h5ad data/processed/GSE268210_phase1_full.h5ad \
    --outdir results/tables

# Step 5: Donor-aware trajectory-lncRNA statistics (5-10 min)
python scripts/analysis/05_donor_aware_statistics.py \
    --mono_h5ad data/processed/GSE268210_monocytes_final.h5ad \
    --lncrna_table results/tables/lncrna_annotation_gencode_v32.csv \
    --outdir results/tables

# Step 6: PAGA topology (2-5 min)
python scripts/analysis/06_paga_topology.py \
    --mono_h5ad data/processed/GSE268210_monocytes_final.h5ad \
    --outdir results/tables --figdir results/figures

# Step 7: AC020656.1 robustness analysis (2-5 min)
python scripts/analysis/07_ac020656_validation.py \
    --mono_h5ad data/processed/GSE268210_monocytes_final.h5ad --outdir results/tables

# Step 8: Bulk validation (45-90 min)
python scripts/analysis/08_bulk_validation.py \
    --data_dir data/bulk_validation --outdir results/tables
```

### 4. Regenerate figures only

```bash
python scripts/visualization/09_figures.py \
    --full_h5ad data/processed/GSE268210_phase1_full.h5ad \
    --mono_h5ad data/processed/GSE268210_monocytes_final.h5ad \
    --traj_csv  results/tables/trajectory_donor_aware_results.csv \
    --bulk_csv  results/tables/bulk_concordance.csv \
    --outdir    results/figures
```

### 5. Independent R validation layer (optional)

```bash
cd scripts/r_analysis
Rscript run_all.R
```

See `scripts/r_analysis/` for the R environment setup (edgeR, limma, metafor, fgsea, nnls).

---

## Monocyte Subtype Proportions

| Subtype | Clusters (leiden_0.3) | n | % | Literature range |
|---------|----------------------|---|---|-----------------|
| CD14 Classical | 5, 6 | 35,805 | 81.3% | 75–85% |
| Intermediate | 1, 2 | 2,588 | 5.9% | 2–10% |
| CD16 Non-Classical | 7 | 5,664 | 12.9% | 5–15% |

## Top Trajectory lncRNAs (donor-aware statistics)

| Gene | Donor-aware ρ | FDR | Direction | Bulk LFC | Bulk FDR |
|------|---|-----|-----------|----------|----------|
| AC020916.1 | −0.354 | <10⁻⁵ | ↑ CD14 | +0.21 | 0.70 |
| AC020656.1 | −0.335 | <10⁻⁵ | ↑ CD14 | +0.75 | **0.017** |
| NEAT1 | −0.259 | <10⁻⁴ | ↑ CD14 | +0.26 | 0.16 |
| MALAT1 | +0.371 | <10⁻⁴ | ↑ CD16 | +0.10 | 0.71 |
| LINC00861 | +0.225 | <10⁻⁵ | ↑ CD16 | +0.11 | 0.70 |

Full table: [`results/tables/trajectory_donor_aware_results.csv`](results/tables/trajectory_donor_aware_results.csv)

**AC020656.1** additionally: 100% of its annotated locus lies within LYZ's terminal exon on the antisense strand; its disease association survives adjustment for LYZ expression (30.8% attenuation, p=0.00086) and for estimated blood cell composition (15.4% attenuation, p=0.0018). See `scripts/r_analysis/ac0206561_lyz/` for the complete locus investigation and `manuscript/manuscript.md` (Section 3.10, Discussion 4.5) for the full discussion of what this does and does not establish.

---

## Figures

| Figure | File | Description |
|--------|------|-------------|
| Fig 1 | `Fig1_QC_violins.png` | QC metrics per patient |
| Fig 2 | `Fig2_UMAP_overview.png` | PBMC atlas + batch correction |
| Fig 3 | `Fig3_marker_dotplot.png` | Canonical marker validation |
| Fig 4 | `Fig9_pseudotime_final.png` | Monocyte pseudotime trajectory |
| Fig 5 | `Fig_paga_topology.png` | PAGA connectivity — monocyte-state topology |
| Fig 6 | `Fig10_trajectory_heatmap_final.png` | Trajectory lncRNA heatmap |
| Fig 7 | `Fig10b_trajectory_lollipop_final.png` | lncRNA-pseudotime correlations |
| Fig 8 | `Fig13_lncrna_celltype_specificity.png` | PBMC lncRNA atlas |
| Fig 9 | `Fig14_bulk_validation.png` | Bulk validation concordance |
| Fig 10 | `Fig_AC0206561_LYZ_locus.png` | AC020656.1/LYZ genomic locus |

See `results/supplementary_figures/` for the R-validation and locus-investigation supplementary figures, and `docs/figure_guide.md` for the complete list.

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
| statsmodels | 0.14 | Multiple testing correction |
| R | 4.3.3 | Independent validation layer |
| edgeR / limma / metafor / fgsea | Bioconductor/CRAN | Donor-aware differential state, meta-analysis, pathway enrichment |

---

## Citation

If you use this code or findings, please cite:

```bibtex
@article{manzoor2026lncrna,
  title   = {Donor-Aware Single-Cell Analysis Identifies AC020656.1, a Locus Antisense
             to LYZ, as a Candidate Monocyte-State and Type 2 Diabetes-Associated Transcript},
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
