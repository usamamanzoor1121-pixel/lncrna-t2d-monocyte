#!/usr/bin/env Rscript
# R_ANALYSIS Phase 7 -- bulk cell-composition sensitivity via marker-based
# NNLS deconvolution (Section 13). This is explicitly a lightweight,
# marker-based approach using the project's OWN scRNA-derived signature --
# not a validated external reference like CIBERSORTx, and its precision
# should not be overstated. The question it answers: does AC020656.1's
# T2D-vs-Control bulk association survive adjustment for estimated monocyte
# fraction, or could it be fully explained by a shift in cell composition?

suppressMessages({
  library(data.table)
  library(nnls)
  library(ggplot2)
})

data_dir <- "scripts/r_analysis/data"
out_dir  <- "results/tables"
fig_dir  <- "results/supplementary_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sig <- fread(file.path(data_dir, "celltype_signature_matrix.csv"))
gene_col <- names(sig)[1]
sig_mat <- as.matrix(sig[, -1, with = FALSE]); rownames(sig_mat) <- sig[[gene_col]]
# signature is in log-norm scRNA units; exponentiate to a linear-ish scale for a marker-ratio-style decomposition
sig_lin <- expm1(sig_mat)

bulk <- fread(file.path(data_dir, "bulk_counts.csv"))
bulk_gene_col <- names(bulk)[1]
bulk_mat <- as.matrix(bulk[, -1, with = FALSE]); rownames(bulk_mat) <- bulk[[bulk_gene_col]]
meta <- fread(file.path(data_dir, "bulk_meta.csv"), header = TRUE); setnames(meta, 1, "sample_id")

common_genes <- intersect(rownames(sig_lin), rownames(bulk_mat))
cat(sprintf("Marker genes available in bulk data: %d / %d\n", length(common_genes), nrow(sig_lin)))
sig_use <- sig_lin[common_genes, , drop = FALSE]

# CPM-normalize bulk for comparability with the (already normalized-scale) signature
libsize <- colSums(bulk_mat)
bulk_cpm <- sweep(bulk_mat, 2, libsize, "/") * 1e6
bulk_use <- bulk_cpm[common_genes, , drop = FALSE]

fractions <- matrix(NA, nrow = ncol(bulk_use), ncol = ncol(sig_use),
                     dimnames = list(colnames(bulk_use), colnames(sig_use)))
for (j in seq_len(ncol(bulk_use))) {
  fit <- nnls(sig_use, bulk_use[, j])
  fractions[j, ] <- fit$x / sum(fit$x)  # normalize to sum to 1
}
frac_dt <- as.data.table(fractions, keep.rownames = "sample_id")
frac_dt <- merge(frac_dt, meta, by = "sample_id")
fwrite(frac_dt, file.path(out_dir, "bulk_estimated_cellfractions_nnls.csv"))

cat("\n=== Estimated cell fractions by condition (mean) ===\n")
print(frac_dt[, lapply(.SD, mean), by = condition, .SDcols = colnames(sig_use)])

p <- ggplot(melt(frac_dt, id.vars = c("sample_id", "condition"), measure.vars = colnames(sig_use)),
            aes(x = condition, y = value, fill = condition)) +
  geom_boxplot() + facet_wrap(~variable, scales = "free_y") +
  labs(title = "NNLS-estimated cell-type fractions in bulk blood (marker-based, approximate)",
       y = "Estimated fraction") + theme_minimal()
ggsave(file.path(fig_dir, "FigH_bulk_composition_fractions.png"), p, width = 10, height = 6, dpi = 150)

## Does AC020656.1's T2D association survive adjustment for estimated Monocyte fraction?
if ("AC020656.1" %in% rownames(bulk_cpm) && "Monocyte" %in% colnames(fractions)) {
  ac_expr <- log2(bulk_cpm["AC020656.1", ] + 1)
  df <- data.table(sample_id = names(ac_expr), ac_expr = ac_expr)
  df <- merge(df, frac_dt[, .(sample_id, condition, Monocyte)], by = "sample_id")
  df <- df[condition %in% c("Control", "T2D")]
  df$condition <- factor(df$condition, levels = c("Control", "T2D"))

  m_unadj <- lm(ac_expr ~ condition, data = df)
  m_adj   <- lm(ac_expr ~ condition + Monocyte, data = df)

  cat("\n=== AC020656.1 T2D vs Control: unadjusted vs monocyte-fraction-adjusted ===\n")
  cat("Unadjusted model:\n"); print(summary(m_unadj)$coefficients)
  cat("\nAdjusted for estimated monocyte fraction:\n"); print(summary(m_adj)$coefficients)

  res_compare <- data.table(
    model = c("unadjusted", "adjusted_for_monocyte_fraction"),
    beta_T2D = c(coef(m_unadj)["conditionT2D"], coef(m_adj)["conditionT2D"]),
    se = c(summary(m_unadj)$coefficients["conditionT2D", "Std. Error"],
           summary(m_adj)$coefficients["conditionT2D", "Std. Error"]),
    p_value = c(summary(m_unadj)$coefficients["conditionT2D", "Pr(>|t|)"],
                summary(m_adj)$coefficients["conditionT2D", "Pr(>|t|)"])
  )
  fwrite(res_compare, file.path(out_dir, "AC0206561_composition_adjustment.csv"))
  print(res_compare)
  pct_change <- 100 * (res_compare$beta_T2D[2] - res_compare$beta_T2D[1]) / res_compare$beta_T2D[1]
  cat(sprintf("\nEffect size change after composition adjustment: %.1f%%\n", pct_change))
  cat(ifelse(abs(pct_change) < 25,
             "-> Association is STABLE after composition adjustment (supports a monocyte-specific transcriptional signal).\n",
             "-> Association changes MATERIALLY after composition adjustment (composition may partly explain the bulk signal).\n"))
}

cat("\nDone. Outputs in:", out_dir, "\n")
cat("CAVEAT: this is a marker-based approximate deconvolution using the project's own scRNA data as reference,\n")
cat("not a validated external method (CIBERSORTx etc. were not available in this offline session). Treat fractions\n")
cat("as directional estimates, not precise cell-count measurements.\n")
