#!/usr/bin/env Rscript
# R_ANALYSIS Phase 5 -- NEAT1 independent R analysis (single-cell side).
# Disease-stage association for NEAT1 is a BULK-only question here (no
# disease-group variable exists in the single-cell data (all donors are T2D),
# section "disease-group contrast: No"); this script covers what the
# single-cell data CAN support: subtype distribution, donor variability,
# detection, and pseudotime association. Disease-stage evidence for NEAT1
# comes from the bulk reanalysis script (06_bulk_limma.R).

suppressMessages({
  library(data.table)
  library(metafor)
  library(ggplot2)
})

data_dir <- "scripts/r_analysis/data"
out_dir  <- "results/tables/neat1"
fig_dir  <- "results/supplementary_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cell <- fread(file.path(data_dir, "cell_level_key_genes.csv"))
GENE <- "NEAT1"

det_by_subtype <- cell[, .(n_cells = .N, detection_rate = mean(get(paste0(GENE, "_detected"))),
                            mean_lognorm = mean(get(paste0(GENE, "_lognorm")))), by = mono_subtype]
fwrite(det_by_subtype, file.path(out_dir, "NEAT1_detection_by_subtype.csv"))
cat("=== NEAT1 detection/expression by subtype ===\n"); print(det_by_subtype)

per_donor <- cell[, .(n_cells = .N, detection_rate = mean(get(paste0(GENE, "_detected"))),
                       mean_lognorm = mean(get(paste0(GENE, "_lognorm")))), by = sample_id]
fwrite(per_donor, file.path(out_dir, "NEAT1_per_donor_summary.csv"))
cat("\n=== NEAT1 per-donor variability ===\n"); print(per_donor)
cat(sprintf("Coefficient of variation across donors (mean_lognorm): %.1f%%\n",
            100 * sd(per_donor$mean_lognorm) / mean(per_donor$mean_lognorm)))

p <- ggplot(cell, aes(x = mono_subtype, y = get(paste0(GENE, "_lognorm")), fill = mono_subtype)) +
  geom_violin(scale = "width") +
  labs(title = "NEAT1 expression by monocyte subtype", y = "log-norm expression", x = NULL) +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(fig_dir, "FigD_NEAT1_subtype_profile.png"), p, width = 6, height = 5, dpi = 150)

## per-donor Spearman rho vs pseudotime -> metafor meta-analysis (same donor-aware framework as AC020656.1)
per_donor_rho <- cell[, {
  ct <- cor.test(dpt_pseudotime, get(paste0(GENE, "_lognorm")), method = "spearman", exact = FALSE)
  .(rho = unname(ct$estimate), p = ct$p.value, n = .N)
}, by = sample_id]
per_donor_rho[, yi := atanh(rho)]; per_donor_rho[, vi := 1 / (n - 3)]
res_re <- rma(yi = yi, vi = vi, data = per_donor_rho, method = "DL", test = "knha")
combined_rho <- tanh(res_re$beta[1])
cat(sprintf("\nNEAT1 donor-aware combined rho = %.4f, p = %.3e, I^2 = %.1f%%\n",
            combined_rho, res_re$pval, res_re$I2))
fwrite(per_donor_rho, file.path(out_dir, "NEAT1_per_donor_rho.csv"))
sink(file.path(out_dir, "NEAT1_metafor_summary.txt")); print(summary(res_re)); sink()

cat("\nDone. Outputs in:", out_dir, "\n")
cat("NOTE: NEAT1's Pre-DM expression peak is a BULK-cohort finding; see results/bulk_validation/ for that test.\n")
