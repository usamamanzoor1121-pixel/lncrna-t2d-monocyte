#!/usr/bin/env Rscript
# Single entry point for the independent R validation layer.
# Runs the full pipeline in order, failing loudly (stop(), not silent skip)
# if a required input is missing. Optional / expensive steps are explicitly
# marked and can be skipped via flags below.
#
# Usage:
#   cd scripts/r_analysis && Rscript run_all.R
# Requires the `renv` conda/micromamba environment described in README.md.
# Expected total runtime: ~5-10 minutes (excluding the one-time GENCODE GTF
# download, cached locally after the first run).

suppressMessages(library(data.table))

root <- "."
data_dir <- file.path(root, "data")

check_input <- function(path, step) {
  if (!file.exists(path)) {
    stop(sprintf("[run_all.R] Missing required input for step '%s': %s\n", step, path),
         "Run the corresponding export step (00_export_for_R.py / 01_export_bulk_for_R.py) first.")
  }
}

cat("=== R validation layer: run_all.R ===\n\n")

cat("[1/12] Checking required data exports exist...\n")
required_data <- c("pseudobulk_counts_donor_subtype.csv", "pseudobulk_meta_donor_subtype.csv",
                    "gene_annotation.csv", "cell_level_key_genes.csv", "bulk_counts.csv",
                    "bulk_meta.csv", "celltype_signature_matrix.csv")
for (f in required_data) check_input(file.path(data_dir, f), "data export")
cat("  OK -- all data exports present. (If any of these are stale, re-run the Python export\n")
cat("  scripts 00_export_for_R.py and 01_export_bulk_for_R.py before continuing.)\n\n")

steps <- list(
  list(script = "02_pseudobulk_qc.R", desc = "Pseudobulk QC", required = TRUE),
  list(script = "03_differential_state_edgeR.R", desc = "Donor-aware differential state (edgeR)", required = TRUE),
  list(script = "04_ac020656_validation.R", desc = "AC020656.1 central case study", required = TRUE),
  list(script = "05_neat1_validation.R", desc = "NEAT1 validation", required = TRUE),
  list(script = "06_bulk_limma.R", desc = "Bulk reanalysis (limma-voom)", required = TRUE),
  list(script = "07_bulk_deconvolution_nnls.R", desc = "Bulk composition sensitivity (NNLS)", required = FALSE,
       note = "Optional: marker-based, approximate; depends on 06_bulk_limma.R."),
  list(script = "08_pathway_fgsea.R", desc = "Pathway enrichment (fgsea)", required = FALSE,
       note = "Optional: downloads MSigDB Hallmark gene sets on first run (network required); depends on 03_differential_state_edgeR.R output."),
  list(script = "09_coexpression_screen.R", desc = "WGCNA decision + correlation screen", required = TRUE),
  list(script = "10_remaining_figures.R", desc = "Bulk disease-stage and volcano figures", required = FALSE),
  list(script = "ac0206561_lyz/01_locus_overlap.R", desc = "AC020656.1-LYZ locus overlap quantification", required = TRUE,
       note = "Requires the GENCODE v32 GTF at data/reference/gencode.v32.annotation.gtf.gz (see README Quick Start)."),
  list(script = "ac0206561_lyz/02_locus_map_figure.R", desc = "AC020656.1-LYZ locus diagram", required = FALSE),
  list(script = "ac0206561_lyz/03_correlation_ratio_models.R", desc = "AC020656.1-LYZ correlation/ratio models", required = TRUE),
  list(script = "ac0206561_lyz/04_disease_lyz_adjustment.R", desc = "AC020656.1 disease association, LYZ-adjusted", required = TRUE)
)

log_dir <- file.path(root, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(steps)) {
  s <- steps[[i]]
  cat(sprintf("[%d/%d] %s (%s)...\n", i + 1, length(steps) + 1, s$desc, s$script))
  if (!is.null(s$note)) cat("  NOTE:", s$note, "\n")
  script_path <- file.path(root, s$script)
  if (!file.exists(script_path)) {
    msg <- sprintf("Script not found: %s", script_path)
    if (s$required) stop(msg) else { cat("  SKIPPED (optional, script missing):", msg, "\n"); next }
  }
  result <- tryCatch({
    source(script_path, echo = FALSE)
    "OK"
  }, error = function(e) {
    if (s$required) {
      stop(sprintf("[run_all.R] REQUIRED step '%s' failed: %s", s$script, conditionMessage(e)))
    }
    cat("  OPTIONAL step failed (continuing):", conditionMessage(e), "\n")
    "FAILED (optional)"
  })
  cat(sprintf("  -> %s\n\n", result))
}

cat("=== run_all.R complete ===\n")
cat("Results in: ../../results/tables and ../../results/supplementary_figures\n")
cat("\nSteps intentionally NOT automated here (data-unavailability blocks, not just slow):\n")
cat("SoupX (no raw/unfiltered matrices exist for this dataset), SingleR/CellTypist independent\n")
cat("annotation (deferred), RNA-Chrom/HiMoRNA manual database lookup.\n")
