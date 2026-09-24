#!/usr/bin/env Rscript
# R_ANALYSIS Phase 3 — donor-aware differential-state analysis between
# monocyte subtypes using edgeR on donor x subtype pseudobulk (raw counts).
# Design: paired by donor (9 donors, each contributing all 3 subtypes) ->
# a proper repeated-measures / paired design, NOT 44,057 pseudo-independent
# cells. This is the R-independent counterpart to the Python donor-aware
# trajectory statistic -- a different question (categorical subtype contrast
# vs. continuous pseudotime correlation) using a different, purpose-built
# donor-aware tool (edgeR's paired GLM) rather than a re-implementation of
# the Python method.

suppressMessages({
  library(edgeR)
  library(data.table)
})

data_dir <- "scripts/r_analysis/data"
out_dir  <- "results/tables/differential_state"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- fread(file.path(data_dir, "pseudobulk_counts_donor_subtype.csv"))
gene_names <- counts[[1]]
counts_mat <- as.matrix(counts[, -1, with = FALSE]); rownames(counts_mat) <- gene_names
meta <- fread(file.path(data_dir, "pseudobulk_meta_donor_subtype.csv"))
gene_anno <- fread(file.path(data_dir, "gene_annotation.csv"))
stopifnot(all(colnames(counts_mat) == meta$pseudobulk_id))

meta$subtype <- factor(meta$subtype, levels = c("CD14_Classical", "Intermediate", "CD16_NonClassical"))
meta$donor <- factor(meta$donor)

# filter low-count genes (standard edgeR practice; group-aware)
keep <- filterByExpr(counts_mat, group = meta$subtype)
cat(sprintf("Genes retained after filterByExpr: %d / %d\n", sum(keep), length(keep)))
counts_f <- counts_mat[keep, ]

dge <- DGEList(counts = counts_f, group = meta$subtype)
dge <- calcNormFactors(dge)

# paired design: donor + subtype (donor as blocking factor -> donor-aware by construction)
design <- model.matrix(~ donor + subtype, data = meta)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

run_contrast <- function(coef_name, label) {
  qlf <- glmQLFTest(fit, coef = coef_name)
  tt <- topTags(qlf, n = Inf)$table
  tt$gene <- rownames(tt)
  tt <- merge(tt, gene_anno[, .(gene_symbol, is_lncrna_gencode)], by.x = "gene", by.y = "gene_symbol", all.x = TRUE)
  tt <- tt[order(tt$PValue), ]
  fwrite(tt, file.path(out_dir, sprintf("edgeR_%s.csv", label)))
  n_sig <- sum(tt$FDR < 0.05, na.rm = TRUE)
  cat(sprintf("%s: %d genes FDR<0.05 (of %d tested)\n", label, n_sig, nrow(tt)))
  invisible(tt)
}

cat("\n=== Design coefficients ===\n"); print(colnames(design))

res_inter <- run_contrast("subtypeIntermediate", "CD14_vs_Intermediate")
res_cd16  <- run_contrast("subtypeCD16_NonClassical", "CD14_vs_CD16")

# Intermediate vs CD16 needs a releveled design
meta2 <- copy(meta); meta2$subtype <- relevel(meta2$subtype, ref = "Intermediate")
design2 <- model.matrix(~ donor + subtype, data = meta2)
dge2 <- DGEList(counts = counts_f, group = meta2$subtype); dge2 <- calcNormFactors(dge2)
dge2 <- estimateDisp(dge2, design2); fit2 <- glmQLFit(dge2, design2)
qlf2 <- glmQLFTest(fit2, coef = "subtypeCD16_NonClassical")
tt2 <- topTags(qlf2, n = Inf)$table; tt2$gene <- rownames(tt2)
tt2 <- merge(tt2, gene_anno[, .(gene_symbol, is_lncrna_gencode)], by.x = "gene", by.y = "gene_symbol", all.x = TRUE)
tt2 <- tt2[order(tt2$PValue), ]
fwrite(tt2, file.path(out_dir, "edgeR_Intermediate_vs_CD16.csv"))
cat(sprintf("Intermediate_vs_CD16: %d genes FDR<0.05 (of %d tested)\n", sum(tt2$FDR < 0.05, na.rm = TRUE), nrow(tt2)))

# highlight AC020656.1 and NEAT1 specifically in each contrast
for (g in c("AC020656.1", "NEAT1", "MALAT1")) {
  cat(sprintf("\n--- %s ---\n", g))
  for (df_name in list(list(res_inter, "CD14 vs Intermediate"), list(res_cd16, "CD14 vs CD16"))) {
    row <- df_name[[1]][df_name[[1]]$gene == g, ]
    if (nrow(row) > 0) {
      cat(sprintf("  %s: logFC=%.3f  PValue=%.3e  FDR=%.3e\n", df_name[[2]], row$logFC[1], row$PValue[1], row$FDR[1]))
    }
  }
  row3 <- tt2[tt2$gene == g, ]
  if (nrow(row3) > 0) cat(sprintf("  Intermediate vs CD16: logFC=%.3f  PValue=%.3e  FDR=%.3e\n", row3$logFC[1], row3$PValue[1], row3$FDR[1]))
}

cat("\nDone. Outputs in:", out_dir, "\n")
