#!/usr/bin/env Rscript
# R_ANALYSIS Phase 12 -- Figure C (AC020656.1/NEAT1 bulk disease-stage) and
# Figure E (CD14 vs CD16 pseudobulk differential-state volcano).

suppressMessages({ library(data.table); library(ggplot2); library(ggrepel) })

data_dir <- "scripts/r_analysis/data"
de_dir   <- "results/tables/differential_state"
fig_dir  <- "results/supplementary_figures"

## Figure C: bulk disease-stage expression for AC020656.1 and NEAT1
counts <- fread(file.path(data_dir, "bulk_counts.csv"))
gene_names <- counts[[1]]
counts_mat <- as.matrix(counts[, -1, with = FALSE]); rownames(counts_mat) <- gene_names
counts_mat <- counts_mat[!duplicated(rownames(counts_mat)), ]
libsize <- colSums(counts_mat)
cpm <- sweep(counts_mat, 2, libsize, "/") * 1e6
logcpm <- log2(cpm + 1)

meta <- fread(file.path(data_dir, "bulk_meta.csv"), header = TRUE); setnames(meta, 1, "sample_id")
meta$condition <- factor(meta$condition, levels = c("Control", "PreDM", "T2D"))

plot_df <- rbindlist(lapply(c("AC020656.1", "NEAT1"), function(g) {
  data.table(sample_id = colnames(logcpm), gene = g, expr = logcpm[g, ], condition = meta$condition)
}))

p <- ggplot(plot_df, aes(x = condition, y = expr, fill = condition)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  facet_wrap(~gene, scales = "free_y") +
  labs(title = "Bulk blood expression across disease stage (independent limma-voom reanalysis)",
       y = "log2(CPM+1)", x = NULL) +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(fig_dir, "FigC_bulk_disease_stage.png"), p, width = 9, height = 5, dpi = 150)
cat("Saved FigC_bulk_disease_stage.png\n")

## Figure E: CD14 vs CD16 pseudobulk differential-state volcano
de <- fread(file.path(de_dir, "edgeR_CD14_vs_CD16.csv"))
de[, sig := FDR < 0.05 & abs(logFC) > 1]
de[, label := ifelse(gene %in% c("AC020656.1", "NEAT1", "MALAT1", "LYZ", "CD14", "FCGR3A", "S100A8"), gene, NA)]

p2 <- ggplot(de, aes(x = logFC, y = -log10(PValue), color = sig)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c(`TRUE` = "#E53E3E", `FALSE` = "grey70")) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20, na.rm = TRUE) +
  labs(title = "Pseudobulk differential state: CD14 Classical vs CD16 Non-Classical (edgeR, donor-paired)",
       subtitle = "Positive logFC = higher in CD16; negative = higher in CD14",
       x = "log2 fold change", y = "-log10(P value)") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(fig_dir, "FigE_CD14_vs_CD16_volcano.png"), p2, width = 8, height = 6, dpi = 150)
cat("Saved FigE_CD14_vs_CD16_volcano.png\n")
