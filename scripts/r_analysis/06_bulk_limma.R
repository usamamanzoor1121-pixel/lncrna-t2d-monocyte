#!/usr/bin/env Rscript
# R_ANALYSIS Phase 6 -- independent bulk RNA-seq reanalysis (GSE221521) in R.
# The Python pipeline used Welch's t-test (and pydeseq2 where available); this
# uses limma-voom, a genuinely different statistical engine, on the SAME raw
# counts+metadata export, as an independent cross-check -- not a copy of the
# Python numbers.

suppressMessages({
  library(edgeR)
  library(limma)
  library(data.table)
})

data_dir <- "scripts/r_analysis/data"
out_dir  <- "results/tables/bulk_validation"
fig_dir  <- "results/supplementary_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- fread(file.path(data_dir, "bulk_counts.csv"))
gene_names <- counts[[1]]
counts_mat <- as.matrix(counts[, -1, with = FALSE]); rownames(counts_mat) <- gene_names
counts_mat <- counts_mat[!duplicated(rownames(counts_mat)), ]
storage.mode(counts_mat) <- "integer"

meta <- fread(file.path(data_dir, "bulk_meta.csv"), header = TRUE)
setnames(meta, 1, "sample_id")
stopifnot(all(colnames(counts_mat) == meta$sample_id))
meta$condition <- factor(meta$condition, levels = c("Control", "PreDM", "T2D"))
cat("Sample sizes:\n"); print(table(meta$condition))

dge <- DGEList(counts = counts_mat, group = meta$condition)
keep <- filterByExpr(dge, group = meta$condition)
cat(sprintf("Genes retained after filterByExpr: %d / %d\n", sum(keep), length(keep)))
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

design <- model.matrix(~ 0 + condition, data = meta)
colnames(design) <- levels(meta$condition)
v <- voom(dge, design)
fit <- lmFit(v, design)

contrasts <- makeContrasts(
  T2D_vs_Control = T2D - Control,
  PreDM_vs_Control = PreDM - Control,
  T2D_vs_PreDM = T2D - PreDM,
  levels = design
)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

run_contrast <- function(coef_name) {
  tt <- topTable(fit2, coef = coef_name, n = Inf, sort.by = "P", confint = TRUE)
  tt$gene <- rownames(tt)
  fwrite(tt, file.path(out_dir, sprintf("limma_%s.csv", coef_name)))
  n_sig <- sum(tt$adj.P.Val < 0.05)
  cat(sprintf("%s: %d genes FDR<0.05 (of %d tested)\n", coef_name, n_sig, nrow(tt)))
  invisible(tt)
}
res_t2d <- run_contrast("T2D_vs_Control")
res_predm <- run_contrast("PreDM_vs_Control")
res_t2d_predm <- run_contrast("T2D_vs_PreDM")

cat("\n=== Key genes: AC020656.1, NEAT1, MALAT1 ===\n")
key_results <- list()
for (g in c("AC020656.1", "NEAT1", "MALAT1")) {
  cat(sprintf("\n--- %s ---\n", g))
  for (pair in list(list(res_t2d, "T2D_vs_Control"), list(res_predm, "PreDM_vs_Control"), list(res_t2d_predm, "T2D_vs_PreDM"))) {
    tt <- pair[[1]]; label <- pair[[2]]
    row <- tt[tt$gene == g, ]
    if (nrow(row) > 0) {
      cat(sprintf("  %s: logFC=%.3f  CI=[%.3f, %.3f]  P=%.4f  FDR=%.4f  AveExpr=%.2f\n",
                  label, row$logFC, row$CI.L, row$CI.R, row$P.Value, row$adj.P.Val, row$AveExpr))
      key_results[[paste(g, label)]] <- data.table(gene = g, contrast = label, logFC = row$logFC,
                                                     CI.L = row$CI.L, CI.R = row$CI.R,
                                                     P.Value = row$P.Value, FDR = row$adj.P.Val, n_samples = nrow(meta))
    }
  }
}
fwrite(rbindlist(key_results), file.path(out_dir, "limma_key_genes_summary.csv"))

## progressive gradient test for AC020656.1 (ordered: Control < PreDM < T2D), independent of Python's Mann-Whitney
if ("AC020656.1" %in% rownames(dge)) {
  expr <- v$E["AC020656.1", ]
  ord_meta <- data.table(condition = meta$condition, expr = expr)
  ord_meta[, ord_code := as.numeric(condition)]
  kt <- cor.test(ord_meta$ord_code, ord_meta$expr, method = "spearman", exact = FALSE)
  cat(sprintf("\nAC020656.1 progressive-gradient test (Spearman, ordered Control<PreDM<T2D): rho=%.3f, p=%.3e\n",
              kt$estimate, kt$p.value))
  jt <- kruskal.test(expr ~ condition, data = ord_meta)
  cat(sprintf("Kruskal-Wallis across 3 groups: p=%.3e\n", jt$p.value))
}

cat("\nDone. Outputs in:", out_dir, "\n")
