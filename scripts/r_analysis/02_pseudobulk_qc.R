#!/usr/bin/env Rscript
# R_ANALYSIS Phase 2 — donor x subtype pseudobulk QC
# Independent R-side sanity check on the pseudobulk counts exported from
# Python (00_export_for_R.py performed only data reshaping, no modeling).

suppressMessages({
  library(data.table)
  library(ggplot2)
})

data_dir <- "scripts/r_analysis/data"
out_dir  <- "results/tables/pseudobulk"
fig_dir  <- "results/supplementary_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

counts <- fread(file.path(data_dir, "pseudobulk_counts_donor_subtype.csv"))
gene_names <- counts[[1]]
counts_mat <- as.matrix(counts[, -1, with = FALSE])
rownames(counts_mat) <- gene_names

meta <- fread(file.path(data_dir, "pseudobulk_meta_donor_subtype.csv"))
stopifnot(all(colnames(counts_mat) == meta$pseudobulk_id))

cat("=== Pseudobulk matrix ===\n")
cat(sprintf("Genes: %d   Samples: %d\n", nrow(counts_mat), ncol(counts_mat)))
cat("Samples per subtype:\n"); print(table(meta$subtype))
cat("Samples per donor:\n"); print(table(meta$donor))

meta$total_counts <- colSums(counts_mat)
meta$n_genes_detected <- colSums(counts_mat > 0)

fwrite(meta, file.path(out_dir, "pseudobulk_qc_summary.csv"))

p1 <- ggplot(meta, aes(x = pseudobulk_id, y = library_size, fill = subtype)) +
  geom_col() + coord_flip() +
  labs(title = "Pseudobulk library size (donor x subtype)", x = NULL, y = "Total raw counts") +
  theme_minimal()
ggsave(file.path(fig_dir, "pseudobulk_library_sizes.png"), p1, width = 8, height = 7, dpi = 150)

p2 <- ggplot(meta, aes(x = subtype, y = n_cells, color = subtype)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, size = 2) +
  labs(title = "Cells contributing per donor x subtype pseudobulk sample", y = "n_cells", x = NULL) +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(fig_dir, "pseudobulk_ncells_per_sample.png"), p2, width = 6, height = 5, dpi = 150)

# Outlier flag: library size or gene count > 3 MAD from median (within subtype)
meta[, lib_z := abs(library_size - median(library_size)) / (1.4826 * mad(library_size)), by = subtype]
meta[, outlier_flag := lib_z > 3]
cat("\nPotential library-size outliers (>3 MAD within subtype):\n")
print(meta[outlier_flag == TRUE, .(pseudobulk_id, donor, subtype, library_size, n_cells)])

fwrite(meta, file.path(out_dir, "pseudobulk_qc_summary.csv"))
cat("\nWrote:", file.path(out_dir, "pseudobulk_qc_summary.csv"), "\n")
cat("Figures in:", fig_dir, "\n")
