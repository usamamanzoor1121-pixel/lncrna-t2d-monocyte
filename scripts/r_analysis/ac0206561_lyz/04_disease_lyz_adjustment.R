#!/usr/bin/env Rscript
# AC0206561_LYZ Phase 4 -- Sections 15, 19, 20: LYZ-adjusted disease
# association (bulk, the only place a disease-status variable exists),
# side-by-side LYZ vs AC020656.1 disease comparison, and detection-threshold
# sensitivity for the CD14-specificity claim.
#
# IMPORTANT SCOPE NOTE: all 9
# single-cell donors are T2D -- there is NO disease-status variable in the
# single-cell data. "Within-CD14 disease association" (mandate Section 17)
# is therefore NOT TESTABLE at single-cell resolution with this dataset.
# The closest available substitute is the bulk composition-adjustment
# analysis already performed (adjusting for estimated monocyte fraction),
# which is reported again here for completeness, NOT re-labeled as a
# single-cell within-CD14 test.

suppressMessages({ library(data.table); library(ggplot2) })

data_dir <- "scripts/r_analysis/data"
out_dir  <- "results/tables/ac0206561_lyz"
fig_dir  <- "results/supplementary_figures"

## ---- Section 15/20: AC020656.1 ~ disease_status  vs  ~ disease_status + LYZ (bulk) ----
counts <- fread(file.path(data_dir, "bulk_counts.csv"))
gene_names <- counts[[1]]
counts_mat <- as.matrix(counts[, -1, with = FALSE]); rownames(counts_mat) <- gene_names
counts_mat <- counts_mat[!duplicated(rownames(counts_mat)), ]
libsize <- colSums(counts_mat)
logcpm <- log2(sweep(counts_mat, 2, libsize, "/") * 1e6 + 1)

meta <- fread(file.path(data_dir, "bulk_meta.csv"), header = TRUE); setnames(meta, 1, "sample_id")
meta$condition <- factor(meta$condition, levels = c("Control", "PreDM", "T2D"))

df <- data.table(sample_id = colnames(logcpm), condition = meta$condition,
                  ac = logcpm["AC020656.1", ], lyz = logcpm["LYZ", ])
df_tc <- df[condition %in% c("Control", "T2D")]
df_tc$condition <- factor(df_tc$condition, levels = c("Control", "T2D"))

m_unadj <- lm(ac ~ condition, data = df_tc)
m_adj   <- lm(ac ~ condition + lyz, data = df_tc)

cat("=== AC020656.1 ~ disease_status (T2D vs Control), bulk, unadjusted ===\n")
print(summary(m_unadj)$coefficients); print(confint(m_unadj))
cat("\n=== AC020656.1 ~ disease_status + LYZ ===\n")
print(summary(m_adj)$coefficients); print(confint(m_adj))

beta_unadj <- coef(m_unadj)["conditionT2D"]; beta_adj <- coef(m_adj)["conditionT2D"]
pct_atten <- 100 * (beta_unadj - beta_adj) / beta_unadj
cat(sprintf("\nUnadjusted beta = %.4f | LYZ-adjusted beta = %.4f | attenuation = %.1f%%\n",
            beta_unadj, beta_adj, pct_atten))
cat("INTERPRETATION (per mandate Section 16): if the association survives, this means\n")
cat("'the disease association is not fully explained by measured LYZ expression' --\n")
cat("it does NOT mean 'AC020656.1 is proven to be independently transcribed.'\n")

comparison <- data.table(
  model = c("unadjusted", "LYZ-adjusted"),
  beta_T2D = c(beta_unadj, beta_adj),
  ci_low = c(confint(m_unadj)["conditionT2D", 1], confint(m_adj)["conditionT2D", 1]),
  ci_high = c(confint(m_unadj)["conditionT2D", 2], confint(m_adj)["conditionT2D", 2]),
  p_value = c(summary(m_unadj)$coefficients["conditionT2D", 4], summary(m_adj)$coefficients["conditionT2D", 4])
)
fwrite(comparison, file.path(out_dir, "04_disease_association_LYZ_adjusted.csv"))

## ---- Section 20: side-by-side AC020656.1 vs LYZ disease-stage comparison ----
plot_df <- rbindlist(list(
  cbind(df, gene = "AC020656.1", expr = df$ac),
  cbind(df, gene = "LYZ", expr = df$lyz)
), use.names = TRUE, fill = TRUE)
p <- ggplot(plot_df, aes(x = condition, y = expr, fill = condition)) +
  geom_violin(scale = "width", alpha = 0.7) + geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
  facet_wrap(~gene, scales = "free_y") +
  labs(title = "AC020656.1 vs LYZ across disease stage (bulk, independent limma model already in bulk_validation/)",
       y = "log2(CPM+1)", x = NULL) +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(fig_dir, "Figure_AC0206561_vs_LYZ_disease_stage.png"), p, width = 9, height = 5, dpi = 150)

## ---- Section 19: detection-threshold sensitivity for CD14-specificity ----
cell <- fread(file.path(data_dir, "cell_level_key_genes.csv"))
# "detection" is currently defined as raw count > 0; test whether requiring
# >=2 or >=3 counts changes the qualitative CD14-specificity conclusion
sens <- rbindlist(lapply(c(1, 2, 3, 5), function(thresh) {
  cell[, .(threshold = thresh, n_cells = .N,
           detection_rate = mean(AC020656.1_counts >= thresh)), by = mono_subtype]
}))
sens_wide <- dcast(sens, mono_subtype + n_cells ~ threshold, value.var = "detection_rate")
setnames(sens_wide, as.character(c(1, 2, 3, 5)), paste0("detection_at_ge", c(1, 2, 3, 5)))
cat("\n=== Detection-threshold sensitivity for AC020656.1 CD14-specificity ===\n")
print(sens_wide)
fwrite(sens_wide, file.path(out_dir, "04_detection_threshold_sensitivity.csv"))

cd14_ratio <- sens_wide[mono_subtype == "CD14_Classical", detection_at_ge1] /
              sens_wide[mono_subtype == "Intermediate", detection_at_ge1]
cat(sprintf("\nCD14:Intermediate detection ratio at threshold>=1: %.1fx\n", cd14_ratio))
for (t in c(2, 3, 5)) {
  col <- paste0("detection_at_ge", t)
  ratio <- sens_wide[mono_subtype == "CD14_Classical", get(col)] / max(sens_wide[mono_subtype == "Intermediate", get(col)], 1e-6)
  cat(sprintf("CD14:Intermediate detection ratio at threshold>=%d: %.1fx\n", t, ratio))
}
cat("-> If this ratio stays large across thresholds, the CD14-specificity finding is not an artifact\n")
cat("   of the specific 'count>0' detection definition.\n")

cat("\nDone. Outputs in:", out_dir, "\n")
