#!/usr/bin/env Rscript
# R_ANALYSIS Phase 10 -- WGCNA feasibility decision (Section 10).
#
# DECISION: full WGCNA network construction is NOT performed on this dataset.
# Rationale (written before looking at any co-expression numbers, to avoid
# post-hoc rationalization):
#   - The only biologically defensible pseudobulk unit is donor x subtype,
#     giving at most 27 samples total, or 9 samples if a network is
#     restricted to a single subtype (the natural choice for a
#     monocyte-subtype-specific module).
#   - Standard WGCNA guidance (Langfelder & Horvath 2008 and subsequent
#     practice) recommends a minimum of ~15-20 samples for stable
#     topological overlap estimation, with 50+ preferred; correlations
#     estimated from n=9-27 are known to be unstable, and soft-threshold
#     power selection is unreliable at this n (scale-free topology fit
#     typically cannot be achieved).
#   - Running WGCNA anyway would produce a network whose module assignments
#     could not be distinguished from sampling noise -- exactly the kind of
#     "manufactured" result from an inadequate design.
#
# INSTEAD: a much more modest, explicitly-labeled correlation screen is run
# to answer the specific, answerable question -- "is AC020656.1 isolated, or
# does it correlate with a coherent set of other genes?" -- without claiming
# network/module structure that the sample size cannot support.

suppressMessages({
  library(data.table)
  library(ggplot2)
})

data_dir <- "scripts/r_analysis/data"
out_dir  <- "results/tables/wgcna"
fig_dir  <- "results/supplementary_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

decision_memo <- "
WGCNA DECISION MEMO
====================
Decision: WGCNA network construction was NOT performed.
Reason: only 9 donors (27 donor x subtype pseudobulk samples, or 9 within a
single subtype) are available -- below the sample size WGCNA methodology
requires for stable module detection (commonly cited guidance: >=15-20
samples minimum, 50+ preferred). Forcing a network at this n would produce
unstable, non-reproducible modules that could not be distinguished from
noise, which the donor-aware statistical framework used throughout this project
exists to avoid -- the same 'do not manufacture
significance from an inadequate design' principle applies here).
Alternative performed: a simple pseudobulk correlation screen for
AC020656.1 and NEAT1 (Spearman rho against all protein-coding genes, ranked
and reported with sample size n=27 clearly stated), explicitly NOT
interpreted as network/module membership.
"
writeLines(decision_memo, file.path(out_dir, "WGCNA_decision_memo.txt"))
cat(decision_memo)

gene_anno <- fread(file.path(data_dir, "gene_annotation.csv"))
coding_genes <- gene_anno[is_lncrna_gencode == FALSE, gene_symbol]

pb <- fread(file.path(data_dir, "pseudobulk_counts_donor_subtype.csv"))
pb_genes <- pb[[1]]; pb_mat <- as.matrix(pb[, -1, with = FALSE]); rownames(pb_mat) <- pb_genes
cpm <- sweep(pb_mat, 2, colSums(pb_mat), "/") * 1e6
logcpm <- log2(cpm + 1)
# keep reasonably expressed, variable genes only (avoids spurious correlations from near-zero genes)
gene_var <- apply(logcpm, 1, var); gene_mean <- rowMeans(logcpm)
keep_genes <- names(gene_var)[gene_var > quantile(gene_var, 0.5) & gene_mean > 1]
coding_keep <- intersect(coding_genes, keep_genes)
cat(sprintf("\nGenes entering the correlation screen: %d (variance- and expression-filtered, protein-coding only)\n", length(coding_keep)))
cat(sprintf("Pseudobulk samples (n): %d\n", ncol(logcpm)))

for (target_gene in c("AC020656.1", "NEAT1")) {
  if (!(target_gene %in% rownames(logcpm))) next
  target_expr <- logcpm[target_gene, ]
  cors <- sapply(coding_keep, function(g) {
    cor.test(logcpm[g, ], target_expr, method = "spearman", exact = FALSE)$estimate
  })
  ps <- sapply(coding_keep, function(g) {
    cor.test(logcpm[g, ], target_expr, method = "spearman", exact = FALSE)$p.value
  })
  res <- data.table(gene = coding_keep, rho = cors, p_value = ps)
  res[, fdr := p.adjust(p_value, method = "BH")]
  res <- res[order(-abs(rho))]
  fname <- file.path(out_dir, sprintf("%s_correlation_screen_n%d.csv", gsub("\\.", "", target_gene), ncol(logcpm)))
  fwrite(res, fname)

  n_strong <- sum(abs(res$rho) > 0.7 & res$fdr < 0.05)
  cat(sprintf("\n%s: %d genes with |rho|>0.7 AND FDR<0.05 (of %d tested, n=%d samples)\n",
              target_gene, n_strong, nrow(res), ncol(logcpm)))
  print(head(res, 10))

  p <- ggplot(res, aes(x = rho)) + geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    labs(title = sprintf("%s pseudobulk correlation screen (n=%d samples, NOT a validated network)",
                          target_gene, ncol(logcpm)),
         x = "Spearman rho", y = "Number of genes") + theme_minimal()
  ggsave(file.path(fig_dir, sprintf("Fig_%s_correlation_distribution.png", gsub("\\.", "", target_gene))),
         p, width = 7, height = 5, dpi = 150)
}

cat("\nDone. Outputs in:", out_dir, "\n")
cat("REMINDER: results above are a correlation screen at n=27 pseudobulk samples, not a WGCNA\n")
cat("module or validated co-expression network. Do not describe genes here as 'module members'.\n")
