# Independent R-Based Validation Layer

An independent statistical layer, using R libraries with no code shared with the Python implementation, that cross-validates the headline findings and extends them where R provides a genuine advantage (edgeR/limma differential-state modeling, `metafor` meta-analysis, `fgsea` pathway enrichment, direct GTF-based genomic characterization).

## Installation

R can be installed system-wide (`conda`/`mamba`) or user-locally via [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) if you don't have admin rights:

```bash
micromamba create -y -n renv -c conda-forge -c bioconda \
  r-base=4.3 r-essentials bioconductor-edger bioconductor-limma \
  r-metafor r-nnls r-data.table bioconductor-fgsea r-msigdbr r-ggplot2 r-jsonlite r-ggrepel
micromamba activate renv
```

Full package list: see the packages listed above; versions are not pinned beyond conda-forge/bioconda's current release at time of installation.

**Packages deliberately not used**: `Seurat`, `SingleCellExperiment`, `scuttle`, `muscat` — all cell-level data needed for R was already exported to flat CSVs by the Python export scripts, so pulling in these heavier Bioconductor infrastructure packages adds install time for no analytical benefit. `WGCNA` itself is intentionally not used at all — see `09_coexpression_screen.R`'s decision memo: with only 9 donors (27 donor×subtype pseudobulk samples), the sample size does not support a stable co-expression network.

## Input requirements

- The GENCODE-annotated lncRNA table from the Python pipeline (`results/tables/lncrna_annotation_gencode_v32.csv`).
- The processed `.h5ad` checkpoints under `data/processed/` (read-only).
- Network access for the one-time GENCODE v32 GTF download (see the main README's Quick Start) and for `msigdbr`'s Hallmark gene-set download on first `fgsea` run.

## Execution

```bash
cd scripts/r_analysis
# 1. Export fresh data from Python (only needed once, or if the upstream h5ad changes)
python3 00_export_for_R.py
python3 01_export_bulk_for_R.py
# 2. Run the full R pipeline
Rscript run_all.R
```

Or run individual numbered scripts directly for iterative work.

## Expected runtime

~5-10 minutes total (dominated by the pseudobulk h5ad export in Python, ~1-2 min, and `fgsea`'s first-run Hallmark download, ~1 min). Individual R scripts run in seconds given the exported flat-file inputs.

## Outputs

- `results/tables/` — differential-state, bulk validation, pathway, and sensitivity-analysis CSVs, organized by theme.
- `results/supplementary_figures/` — forest plots, volcano plots, pathway enrichment bar charts, composition-fraction boxplots.
- `ac0206561_lyz/` — the AC020656.1–LYZ genomic locus investigation (locus overlap, correlation models, disease/LYZ adjustment).

## Known limitations

- No ambient-RNA (SoupX) correction — no raw/unfiltered 10x matrices exist for this GEO series.
- No independent reference-based cell-type annotation (SingleR/CellTypist) — deferred.
- No WGCNA network — sample size inadequate, declined by design (see `09_coexpression_screen.R`).
- No RNA-Chrom/HiMoRNA database query — interactive web resources, no bulk API used; manual lookup recommended.
- AC020656.1 is antisense to and nearly fully overlaps LYZ (ρ≈0.94 pseudobulk co-expression); the available data cannot distinguish genuine independent transcription from antisense read bleed-through. See `ac0206561_lyz/` and the manuscript's Discussion for the full treatment.
