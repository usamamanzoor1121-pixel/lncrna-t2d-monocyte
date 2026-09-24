#!/usr/bin/env Rscript
# AC0206561_LYZ Phase 1 -- precise locus/overlap quantification (Sections 4-6).
# Uses the same local GENCODE v32 GTF as the prior genomic-characterization
# script, but computes exon-level and gene-body-level overlap explicitly and
# separately, rather than repeating the single "87% overlap" headline number
# without defining what it means.

suppressMessages(library(data.table))

GTF <- "data/reference/gencode.v32.annotation.gtf.gz"
out_dir <- "results/tables/ac0206561_lyz"
rep_dir <- "scripts/r_analysis/AC0206561_LYZ/reports"

extract_attr <- function(attr_str, key) {
  pat <- sprintf('%s "([^"]+)"', key)
  pos <- regexpr(pat, attr_str)
  raw_matches <- regmatches(attr_str, pos)
  out <- rep(NA_character_, length(attr_str))
  out[pos > 0] <- gsub(sprintf('%s "|"', key), "", raw_matches)
  out
}

cat("Reading GENCODE v32 GTF...\n")
gtf <- fread(cmd = paste("zcat", shQuote(GTF)), sep = "\t", header = FALSE, quote = "",
             col.names = c("chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"))
gtf[, gene_name := extract_attr(attributes, "gene_name")]
gtf[, gene_id := sub("\\..*", "", extract_attr(attributes, "gene_id"))]
gtf[, transcript_id := extract_attr(attributes, "transcript_id")]
gtf[, gene_type := extract_attr(attributes, "gene_type")]
gtf[, transcript_type := extract_attr(attributes, "transcript_type")]

get_gene_block <- function(symbol) {
  g <- gtf[feature == "gene" & gene_name == symbol]
  tx <- gtf[feature == "transcript" & gene_name == symbol]
  ex <- gtf[feature == "exon" & gene_name == symbol]
  list(gene = g, transcripts = tx, exons = ex)
}

ac <- get_gene_block("AC020656.1")
lyz <- get_gene_block("LYZ")

cat("\n=== AC020656.1 ===\n"); print(ac$gene[, .(chrom, start, end, strand, gene_id, gene_type)])
cat(sprintf("Transcripts: %d | Exons (all transcripts): %d\n", nrow(ac$transcripts), nrow(ac$exons)))
cat("\n=== LYZ ===\n"); print(lyz$gene[, .(chrom, start, end, strand, gene_id, gene_type)])
cat(sprintf("Transcripts: %d | Exons (all transcripts): %d\n", nrow(lyz$transcripts), nrow(lyz$exons)))

## ---- gene-body overlap ----
ac_start <- ac$gene$start[1]; ac_end <- ac$gene$end[1]
lyz_start <- lyz$gene$start[1]; lyz_end <- lyz$gene$end[1]
ac_len <- ac_end - ac_start + 1
lyz_len <- lyz_end - lyz_start + 1

overlap_start <- max(ac_start, lyz_start)
overlap_end <- min(ac_end, lyz_end)
overlap_len <- max(0, overlap_end - overlap_start + 1)

pct_ac_overlapping <- 100 * overlap_len / ac_len
pct_lyz_overlapping <- 100 * overlap_len / lyz_len

cat("\n=== GENE-BODY OVERLAP (whole annotated gene span, not exon-specific) ===\n")
cat(sprintf("AC020656.1 total length         : %d bp (chr12:%d-%d)\n", ac_len, ac_start, ac_end))
cat(sprintf("LYZ total length                : %d bp (chr12:%d-%d)\n", lyz_len, lyz_start, lyz_end))
cat(sprintf("Overlapping interval             : %d bp (chr12:%d-%d)\n", overlap_len, overlap_start, overlap_end))
cat(sprintf("%% of AC020656.1 gene body overlapping LYZ gene body : %.1f%%\n", pct_ac_overlapping))
cat(sprintf("%% of LYZ gene body overlapping AC020656.1 gene body : %.1f%%\n", pct_lyz_overlapping))
cat(sprintf("AC020656.1 unique (non-overlapping) length : %d bp\n", ac_len - overlap_len))
cat(sprintf("LYZ unique (non-overlapping) length         : %d bp\n", lyz_len - overlap_len))

## ---- exon-level overlap (the biologically relevant number for a single-exon
##      transcript being sequenced, since only exonic sequence is captured
##      by polyA-primed/TSO-based 10x chemistry) ----
ac_exons <- ac$exons[, .(start, end)]
lyz_exons <- unique(lyz$exons[, .(start, end)])  # LYZ may have multiple transcripts; dedupe identical exon coords

exon_overlap_bp <- function(a_start, a_end, exon_table) {
  ov <- pmax(0, pmin(a_end, exon_table$end) - pmax(a_start, exon_table$start) + 1)
  sum(ov[ov > 0])
}

ac_total_exonic_len <- sum(ac_exons$end - ac_exons$start + 1)
ac_exon_overlap_with_lyz_exons <- sum(sapply(seq_len(nrow(ac_exons)), function(i)
  exon_overlap_bp(ac_exons$start[i], ac_exons$end[i], lyz_exons)))

cat("\n=== EXON-LEVEL OVERLAP (AC020656.1 exonic bp vs. any LYZ exon) ===\n")
cat(sprintf("AC020656.1 total exonic length            : %d bp\n", ac_total_exonic_len))
cat(sprintf("AC020656.1 exonic bp overlapping a LYZ exon : %d bp\n", ac_exon_overlap_with_lyz_exons))
cat(sprintf("%% of AC020656.1 EXONIC sequence overlapping a LYZ EXON : %.1f%%\n",
            100 * ac_exon_overlap_with_lyz_exons / ac_total_exonic_len))
cat(sprintf("AC020656.1 exonic sequence with NO overlapping LYZ exon : %d bp\n",
            ac_total_exonic_len - ac_exon_overlap_with_lyz_exons))

## is the AC020656.1 exon inside a LYZ EXON specifically, or only inside a LYZ INTRON?
lyz_introns <- data.table()
if (nrow(lyz$transcripts) > 0) {
  # build introns per transcript from its own exons, then check
  for (txid in unique(lyz$exons$transcript_id)) {
    tx_ex <- lyz$exons[transcript_id == txid][order(start)]
    if (nrow(tx_ex) > 1) {
      for (i in seq_len(nrow(tx_ex) - 1)) {
        lyz_introns <- rbind(lyz_introns, data.table(start = tx_ex$end[i] + 1, end = tx_ex$start[i + 1] - 1))
      }
    }
  }
}
ac_overlap_with_lyz_introns <- if (nrow(lyz_introns) > 0) {
  sum(sapply(seq_len(nrow(ac_exons)), function(i) exon_overlap_bp(ac_exons$start[i], ac_exons$end[i], lyz_introns)))
} else 0
cat(sprintf("\nAC020656.1 exonic bp overlapping a LYZ INTRON (not exon) : %d bp\n", ac_overlap_with_lyz_introns))
cat("-> If AC020656.1's exon overlaps a LYZ EXON (not just an intron), reads mapping\n")
cat("   there are exonic for BOTH genes -- this is the scenario where 10x's exon-based\n")
cat("   counting is most prone to strand-dependent (not sequence-dependent) disambiguation,\n")
cat("   i.e. correct gene assignment depends entirely on the read's strand of origin being\n")
cat("   preserved and correctly interpreted, not on sequence uniqueness.\n")

summary_dt <- data.table(
  metric = c("AC020656.1_total_length_bp", "AC020656.1_gene_body_overlap_bp", "AC020656.1_gene_body_overlap_pct",
             "LYZ_total_length_bp", "LYZ_gene_body_overlap_bp", "LYZ_gene_body_overlap_pct",
             "AC020656.1_total_exonic_bp", "AC020656.1_exonic_overlap_with_LYZ_exon_bp",
             "AC020656.1_exonic_overlap_with_LYZ_exon_pct", "AC020656.1_exonic_overlap_with_LYZ_intron_bp"),
  value = c(ac_len, overlap_len, round(pct_ac_overlapping, 1),
            lyz_len, overlap_len, round(pct_lyz_overlapping, 1),
            ac_total_exonic_len, ac_exon_overlap_with_lyz_exons,
            round(100 * ac_exon_overlap_with_lyz_exons / ac_total_exonic_len, 1), ac_overlap_with_lyz_introns)
)
fwrite(summary_dt, file.path(out_dir, "01_locus_overlap_quantification.csv"))

fwrite(ac$gene[, .(chrom, start, end, strand, gene_id, gene_name, gene_type)],
       file.path(out_dir, "01_locus_annotation.tsv"), sep = "\t")
fwrite(rbind(
  cbind(gene = "AC020656.1", ac$transcripts[, .(transcript_id, transcript_type, start, end, strand)]),
  cbind(gene = "LYZ", lyz$transcripts[, .(transcript_id, transcript_type, start, end, strand)])
), file.path(out_dir, "01_locus_annotation.tsv"), sep = "\t", append = TRUE)

cat("\nDone. Outputs in:", out_dir, "\n")
