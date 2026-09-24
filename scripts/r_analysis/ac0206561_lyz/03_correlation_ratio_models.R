#!/usr/bin/env Rscript
# AC0206561_LYZ Phase 3 -- Sections 12-15: correlation, ratio, and
# LYZ-adjusted disease-association models, done properly (donor/subtype-level,
# not cell-pooled; explicit models; CIs where the model supports them).

suppressMessages({ library(data.table); library(ggplot2) })

data_dir <- "scripts/r_analysis/data"
out_dir  <- "results/tables/ac0206561_lyz"
fig_dir  <- "results/supplementary_figures"

## ---- pseudobulk: primary evidence (n=27 donor x subtype samples) ----
pb <- fread(file.path(data_dir, "pseudobulk_counts_donor_subtype.csv"))
pb_genes <- pb[[1]]; pb_mat <- as.matrix(pb[, -1, with = FALSE]); rownames(pb_mat) <- pb_genes
cpm <- sweep(pb_mat, 2, colSums(pb_mat), "/") * 1e6
logcpm <- log2(cpm + 1)
meta <- fread(file.path(data_dir, "pseudobulk_meta_donor_subtype.csv"))

df <- data.table(pseudobulk_id = colnames(logcpm), donor = meta$donor, subtype = meta$subtype,
                  ac = logcpm["AC020656.1", ], lyz = logcpm["LYZ", ])

## Section 12: Figure 1 -- AC020656.1 vs LYZ, donor/subtype-level (primary evidence)
p1 <- ggplot(df, aes(x = lyz, y = ac, color = subtype, shape = donor)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_shape_manual(values = 1:9) +
  geom_smooth(aes(group = 1), method = "lm", color = "black", linetype = "dashed", se = TRUE) +
  labs(title = "AC020656.1 vs LYZ -- donor x subtype pseudobulk (n=27, the primary evidence unit)",
       x = "LYZ log2(CPM+1)", y = "AC020656.1 log2(CPM+1)") +
  theme_minimal()
ggsave(file.path(fig_dir, "Figure1_AC0206561_vs_LYZ_donor_subtype.png"), p1, width = 8, height = 6, dpi = 150)

overall_pearson <- cor.test(df$ac, df$lyz, method = "pearson")
overall_spearman <- cor.test(df$ac, df$lyz, method = "spearman", exact = FALSE)
cat(sprintf("Overall (n=27): Pearson r=%.3f (p=%.3e) | Spearman rho=%.3f (p=%.3e)\n",
            overall_pearson$estimate, overall_pearson$p.value, overall_spearman$estimate, overall_spearman$p.value))

within_subtype <- df[, {
  pear <- cor.test(ac, lyz, method = "pearson")
  spear <- cor.test(ac, lyz, method = "spearman", exact = FALSE)
  .(n = .N, pearson_r = pear$estimate, pearson_p = pear$p.value,
    spearman_rho = spear$estimate, spearman_p = spear$p.value)
}, by = subtype]
cat("\n=== Within-subtype AC020656.1-LYZ correlation (n=9 donors per subtype) ===\n")
print(within_subtype)
fwrite(within_subtype, file.path(out_dir, "03_within_subtype_correlation.csv"))

## Section 13: does AC020656.1 contain information beyond LYZ abundance?
m0 <- lm(ac ~ lyz, data = df)
m1 <- lm(ac ~ lyz + subtype, data = df)
cat("\n=== AC020656.1 ~ LYZ ===\n"); print(summary(m0)$coefficients)
cat(sprintf("R^2 = %.3f\n", summary(m0)$r.squared))
cat("\n=== AC020656.1 ~ LYZ + subtype ===\n"); print(summary(m1)$coefficients)
cat(sprintf("R^2 = %.3f (vs %.3f without subtype)\n", summary(m1)$r.squared, summary(m0)$r.squared))
av <- anova(m0, m1)
cat("\nANOVA (does subtype add explanatory power beyond LYZ alone?):\n"); print(av)
cat(ifelse(av$`Pr(>F)`[2] < 0.05,
    "-> Subtype adds SIGNIFICANT explanatory power beyond LYZ abundance alone: AC020656.1 is not simply\n   a linear function of LYZ; some subtype-specific structure exists in the AC020656.1 signal that LYZ\n   alone does not capture. This does not prove independent transcription (subtype-specific antisense\n   noise from LYZ's own subtype-specific promoter activity is equally consistent with this result).\n",
    "-> Subtype does NOT add significant explanatory power beyond LYZ abundance: consistent with (but does\n   not prove) AC020656.1 being essentially a function of LYZ abundance alone.\n"))

sink(file.path(out_dir, "03_lm_ac_vs_lyz.txt"))
cat("=== AC020656.1 ~ LYZ ===\n"); print(summary(m0))
cat("\n=== AC020656.1 ~ LYZ + subtype ===\n"); print(summary(m1))
cat("\n=== ANOVA comparison ===\n"); print(av)
sink()

## Section 14: ratio analysis (verify existing p=0.29 claim independently, with model + CI)
df[, log_ratio := ac - lyz]  # log2(AC/LYZ) since both already log2
kt <- kruskal.test(log_ratio ~ subtype, data = df)
lm_ratio <- lm(log_ratio ~ subtype, data = df)
cat(sprintf("\n=== Ratio analysis: log2(AC020656.1/LYZ) by subtype ===\n"))
cat(sprintf("Kruskal-Wallis: p = %.4f (previously reported: p=0.29)\n", kt$p.value))
print(summary(lm_ratio)$coefficients)
ci <- confint(lm_ratio)
cat("95% CI on subtype coefficients:\n"); print(ci)
fwrite(as.data.table(df)[, .(pseudobulk_id, donor, subtype, ac, lyz, log_ratio)],
       file.path(out_dir, "03_ac_lyz_ratio_by_sample.csv"))

## also test log_ratio by donor (does the ratio vary more across donors than across subtypes?)
kt_donor <- kruskal.test(log_ratio ~ donor, data = df)
cat(sprintf("\nKruskal-Wallis (log-ratio by DONOR instead of subtype): p = %.4f\n", kt_donor$p.value))

cat("\nDone (Sections 12-14). Outputs in:", out_dir, "\n")
