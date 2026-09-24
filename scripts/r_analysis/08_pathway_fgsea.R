#!/usr/bin/env Rscript
# R_ANALYSIS Phase 8 -- pathway enrichment (Section 11).
# lncRNAs (AC020656.1, NEAT1) have no canonical pathway annotation of their
# own. We therefore run fgsea on PROTEIN-CODING gene rankings that are
# ASSOCIATED with each lncRNA (by pseudobulk correlation) or that define the
# CD14-vs-CD16 differential-state signature already computed in R
# (03_differential_state_edgeR.R). Results are explicitly labeled
# "associated program", never "pathway regulated by the lncRNA".

suppressMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

data_dir <- "scripts/r_analysis/data"
de_dir   <- "results/tables/differential_state"
out_dir  <- "results/tables/pathways"
fig_dir  <- "results/supplementary_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gene_anno <- fread(file.path(data_dir, "gene_annotation.csv"))
coding_genes <- gene_anno[is_lncrna_gencode == FALSE, gene_symbol]

hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)
cat(sprintf("Loaded %d Hallmark gene sets\n", length(pathways)))

run_fgsea_on_ranks <- function(ranks, label) {
  ranks <- ranks[!is.na(ranks) & is.finite(ranks)]
  ranks <- ranks[!duplicated(names(ranks))]
  set.seed(42)
  res <- fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500, eps = 0)
  res <- res[order(pval)]
  fwrite(res[, .(pathway, pval, padj, ES, NES, size)], file.path(out_dir, sprintf("fgsea_%s.csv", label)))
  n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
  cat(sprintf("%s: %d Hallmark sets FDR<0.05 (of %d tested)\n", label, n_sig, nrow(res)))
  if (n_sig > 0) {
    top <- head(res[padj < 0.05][order(-abs(NES))], 10)
    p <- ggplot(top, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
      geom_col() + coord_flip() +
      labs(title = sprintf("Top enriched Hallmark sets: %s", label), x = NULL, y = "NES") +
      theme_minimal() + theme(legend.position = "none")
    ggsave(file.path(fig_dir, sprintf("FigG_fgsea_%s.png", label)), p, width = 8, height = 5, dpi = 150)
  }
  invisible(res)
}

## A/B: AC020656.1- and NEAT1-associated programs (pseudobulk correlation ranking, coding genes only)
pb <- fread(file.path(data_dir, "pseudobulk_counts_donor_subtype.csv"))
pb_genes <- pb[[1]]; pb_mat <- as.matrix(pb[, -1, with = FALSE]); rownames(pb_mat) <- pb_genes
cpm <- sweep(pb_mat, 2, colSums(pb_mat), "/") * 1e6
logcpm <- log2(cpm + 1)

for (target_gene in c("AC020656.1", "NEAT1")) {
  if (!(target_gene %in% rownames(logcpm))) next
  target_expr <- logcpm[target_gene, ]
  coding_in_mat <- intersect(coding_genes, rownames(logcpm))
  cors <- apply(logcpm[coding_in_mat, , drop = FALSE], 1, function(x) {
    if (sd(x) < 1e-6) return(NA)
    cor(x, target_expr, method = "spearman")
  })
  cors <- cors[!is.na(cors)]
  fwrite(data.table(gene = names(cors), rho = cors)[order(-rho)],
         file.path(out_dir, sprintf("%s_associated_gene_correlations.csv", gsub("\\.", "", target_gene))))
  run_fgsea_on_ranks(cors, sprintf("%s_associated_program", gsub("\\.", "", target_gene)))
}

## C: CD14 vs CD16 differential-state signature (from edgeR results already computed)
cd14_cd16_file <- file.path(de_dir, "edgeR_CD14_vs_CD16.csv")
if (file.exists(cd14_cd16_file)) {
  de <- fread(cd14_cd16_file)
  de <- de[is_lncrna_gencode == FALSE]  # protein-coding only, per Section 11 instruction
  ranks <- setNames(de$logFC, de$gene)
  run_fgsea_on_ranks(ranks, "CD14_vs_CD16_state")
}

cat("\nDone. Outputs in:", out_dir, "\n")
cat("REMINDER: these are gene programs statistically ASSOCIATED with lncRNA expression or subtype state,\n")
cat("not pathways demonstrated to be regulated BY the lncRNA. No mechanism is claimed.\n")
