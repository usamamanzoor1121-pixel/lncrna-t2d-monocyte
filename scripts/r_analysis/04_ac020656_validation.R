#!/usr/bin/env Rscript
# R_ANALYSIS Phase 4 -- AC020656.1 central case study, independently in R.
# Uses `metafor` (a different implementation than the hand-written
# DerSimonian-Laird+HKSJ implementation in Python) as genuine
# cross-validation of the meta-analytic combination, plus donor-level and
# pseudobulk-based subtype comparisons that never pool cells across donors.

suppressMessages({
  library(data.table)
  library(metafor)
  library(ggplot2)
})

data_dir <- "scripts/r_analysis/data"
out_dir  <- "results/tables/ac020656"
fig_dir  <- "results/supplementary_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cell <- fread(file.path(data_dir, "cell_level_key_genes.csv"))
GENE <- "AC020656.1"

## ---- 8.1 / 8.3: donor-level expression + detection rate ----
per_donor <- cell[, .(
  n_cells = .N,
  detection_rate = mean(get(paste0(GENE, "_detected"))),
  mean_lognorm = mean(get(paste0(GENE, "_lognorm"))),
  mean_counts = mean(get(paste0(GENE, "_counts")))
), by = sample_id]
fwrite(per_donor, file.path(out_dir, "AC0206561_per_donor_summary.csv"))
cat("=== Per-donor AC020656.1 summary ===\n"); print(per_donor)

p_donor <- ggplot(cell, aes(x = sample_id, y = get(paste0(GENE, "_lognorm")), fill = mono_subtype)) +
  geom_violin(scale = "width", alpha = 0.7) +
  labs(title = "AC020656.1 expression per donor and subtype", y = "log-norm expression", x = "Donor") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "FigA_AC0206561_donor_subtype.png"), p_donor, width = 10, height = 5, dpi = 150)

## ---- 8.3: detection rate by subtype, exact denominator/definition stated ----
det_by_subtype <- cell[, .(
  n_cells = .N,
  n_detected = sum(get(paste0(GENE, "_detected"))),
  detection_rate = mean(get(paste0(GENE, "_detected")))
), by = mono_subtype]
fwrite(det_by_subtype, file.path(out_dir, "AC0206561_detection_by_subtype.csv"))
cat("\n=== Detection rate by subtype (denominator = all cells of that subtype, definition = raw count > 0) ===\n")
print(det_by_subtype)

## ---- 8.2 / 8.4 / 8.5: per-donor Spearman rho vs pseudotime, metafor meta-analysis ----
per_donor_rho <- cell[, {
  ct <- cor.test(dpt_pseudotime, get(paste0(GENE, "_lognorm")), method = "spearman", exact = FALSE)
  .(rho = unname(ct$estimate), p = ct$p.value, n = .N)
}, by = sample_id]
fwrite(per_donor_rho, file.path(out_dir, "AC0206561_per_donor_rho.csv"))
cat("\n=== Per-donor Spearman rho vs pseudotime ===\n"); print(per_donor_rho)

# Fisher z transform + metafor random-effects meta-analysis (independent library)
per_donor_rho[, yi := atanh(rho)]
per_donor_rho[, vi := 1 / (n - 3)]
res_re <- rma(yi = yi, vi = vi, data = per_donor_rho, method = "DL", test = "knha")  # knha = Hartung-Knapp, independent implementation
cat("\n=== metafor random-effects meta-analysis (DerSimonian-Laird tau^2, Hartung-Knapp test) ===\n")
print(summary(res_re))
combined_rho <- tanh(res_re$beta[1])
cat(sprintf("\nCombined rho (back-transformed) = %.4f,  p = %.3e,  I^2 = %.1f%%\n",
            combined_rho, res_re$pval, res_re$I2))

sink(file.path(out_dir, "AC0206561_metafor_summary.txt"))
print(summary(res_re))
cat(sprintf("\nCombined rho (back-transformed) = %.4f\n", combined_rho))
sink()

png(file.path(fig_dir, "FigB_AC0206561_forest_plot.png"), width = 900, height = 700, res = 130)
forest(res_re, slab = per_donor_rho$sample_id, xlab = "Fisher z (AC020656.1 vs pseudotime)",
       main = "AC020656.1 -- per-donor effect (metafor random-effects meta-analysis)")
dev.off()

## ---- leave-one-donor-out, independently in R via metafor ----
loo_rows <- list()
for (i in seq_len(nrow(per_donor_rho))) {
  sub <- per_donor_rho[-i]
  r <- rma(yi = yi, vi = vi, data = sub, method = "DL", test = "knha")
  loo_rows[[i]] <- data.table(
    excluded_donor = per_donor_rho$sample_id[i],
    combined_rho = tanh(r$beta[1]), p_value = r$pval, I2 = r$I2
  )
}
loo_dt <- rbindlist(loo_rows)
fwrite(loo_dt, file.path(out_dir, "AC0206561_leave_one_donor_out_R.csv"))
cat("\n=== Leave-one-donor-out (R / metafor) ===\n"); print(loo_dt)
cat(sprintf("\nrho range: %.4f  |  worst-case p: %.3e\n",
            max(loo_dt$combined_rho) - min(loo_dt$combined_rho), max(loo_dt$p_value)))

## ---- subtype comparison using donor-level pseudobulk means (not pooled cells) ----
donor_subtype_means <- cell[, .(mean_expr = mean(get(paste0(GENE, "_lognorm")))), by = .(sample_id, mono_subtype)]
donor_subtype_wide <- dcast(donor_subtype_means, sample_id ~ mono_subtype, value.var = "mean_expr")
fwrite(donor_subtype_wide, file.path(out_dir, "AC0206561_donor_subtype_means.csv"))

wt_cd14_inter <- wilcox.test(donor_subtype_wide$CD14_Classical, donor_subtype_wide$Intermediate, paired = TRUE)
wt_cd14_cd16  <- wilcox.test(donor_subtype_wide$CD14_Classical, donor_subtype_wide$CD16_NonClassical, paired = TRUE)
wt_inter_cd16 <- wilcox.test(donor_subtype_wide$Intermediate, donor_subtype_wide$CD16_NonClassical, paired = TRUE)

cat("\n=== Paired (by donor) Wilcoxon signed-rank tests on donor-level mean expression (n=9 donors each) ===\n")
cat(sprintf("CD14 vs Intermediate: V=%s p=%.4f\n", wt_cd14_inter$statistic, wt_cd14_inter$p.value))
cat(sprintf("CD14 vs CD16:         V=%s p=%.4f\n", wt_cd14_cd16$statistic, wt_cd14_cd16$p.value))
cat(sprintf("Intermediate vs CD16: V=%s p=%.4f\n", wt_inter_cd16$statistic, wt_inter_cd16$p.value))

subtype_test_results <- data.table(
  comparison = c("CD14_vs_Intermediate", "CD14_vs_CD16", "Intermediate_vs_CD16"),
  test = "paired Wilcoxon signed-rank (n=9 donors)",
  statistic = c(wt_cd14_inter$statistic, wt_cd14_cd16$statistic, wt_inter_cd16$statistic),
  p_value = c(wt_cd14_inter$p.value, wt_cd14_cd16$p.value, wt_inter_cd16$p.value)
)
subtype_test_results[, fdr := p.adjust(p_value, method = "BH")]
fwrite(subtype_test_results, file.path(out_dir, "AC0206561_subtype_comparison_tests.csv"))
print(subtype_test_results)

cat("\nDone. Outputs in:", out_dir, "\n")
