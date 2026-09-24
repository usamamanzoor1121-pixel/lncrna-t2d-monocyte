#!/usr/bin/env Rscript
# AC0206561_LYZ Phase 2 -- publication-style locus diagram (Section 5), using
# EXACT exon coordinates read fresh from the GENCODE v32 GTF (not schematic
# placeholders) since accurate representation of the overlap is the entire
# point of this figure.

suppressMessages({ library(data.table); library(ggplot2) })

GTF <- "data/reference/gencode.v32.annotation.gtf.gz"
fig_dir <- "results/supplementary_figures"

extract_attr <- function(attr_str, key) {
  pat <- sprintf('%s "([^"]+)"', key)
  pos <- regexpr(pat, attr_str)
  raw_matches <- regmatches(attr_str, pos)
  out <- rep(NA_character_, length(attr_str))
  out[pos > 0] <- gsub(sprintf('%s "|"', key), "", raw_matches)
  out
}

gtf <- fread(cmd = paste("zcat", shQuote(GTF)), sep = "\t", header = FALSE, quote = "",
             col.names = c("chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"))
gtf[, gene_name := extract_attr(attributes, "gene_name")]
gtf[, transcript_id := extract_attr(attributes, "transcript_id")]

lyz_exons <- unique(gtf[feature == "exon" & gene_name == "LYZ", .(start, end, transcript_id)])
ac_exons <- unique(gtf[feature == "exon" & gene_name == "AC020656.1", .(start, end, transcript_id)])
lyz_gene <- gtf[feature == "gene" & gene_name == "LYZ"]
ac_gene <- gtf[feature == "gene" & gene_name == "AC020656.1"]

cat("LYZ exons (all transcripts):\n"); print(lyz_exons)
cat("AC020656.1 exons:\n"); print(ac_exons)

# stack LYZ's 3 transcripts on separate rows so overlapping exons across
# transcripts remain visually distinguishable, rather than merging them
lyz_tx_ids <- unique(lyz_exons$transcript_id)
lyz_exons[, row := match(transcript_id, lyz_tx_ids) + 1]  # rows 2,3,4 for LYZ's transcripts
ac_row <- 1

p <- ggplot() +
  geom_segment(data = data.frame(tx = lyz_tx_ids, row = seq_along(lyz_tx_ids) + 1),
               aes(x = lyz_gene$start, xend = lyz_gene$end, y = row, yend = row), color = "grey50", linewidth = 0.6) +
  geom_rect(data = lyz_exons, aes(xmin = start, xmax = end, ymin = row - 0.15, ymax = row + 0.15), fill = "#3182CE") +
  geom_segment(aes(x = ac_gene$start, xend = ac_gene$end, y = ac_row, yend = ac_row), color = "grey50", linewidth = 0.6) +
  geom_rect(data = ac_exons, aes(xmin = start, xmax = end, ymin = ac_row - 0.15, ymax = ac_row + 0.15), fill = "#E53E3E") +
  annotate("rect", xmin = ac_gene$start, xmax = ac_gene$end, ymin = 0.5, ymax = max(lyz_exons$row) + 0.5,
           fill = "grey50", alpha = 0.12) +
  annotate("text", x = (ac_gene$start + ac_gene$end) / 2, y = max(lyz_exons$row) + 0.7,
           label = sprintf("AC020656.1 (%d bp) fully nested within LYZ's\nterminal exon region -- 100%% gene-body overlap",
                            ac_gene$end - ac_gene$start + 1),
           size = 3, vjust = 0) +
  scale_y_continuous(breaks = c(ac_row, seq_along(lyz_tx_ids) + 1),
                      labels = c("AC020656.1\n(- strand)", paste0("LYZ ", lyz_tx_ids, "\n(+ strand)")),
                      limits = c(0.3, max(lyz_exons$row) + 1.2)) +
  labs(title = sprintf("AC020656.1 / LYZ locus (chr12:%d-%d, GRCh38/GENCODE v32)", lyz_gene$start, lyz_gene$end),
       subtitle = "Exact exon coordinates from GENCODE v32; AC020656.1 exon sits within LYZ's 3' exonic region, opposite strand",
       x = "Genomic position (chr12)", y = NULL) +
  theme_minimal() + theme(panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "Figure_AC0206561_LYZ_locus.png"), p, width = 9, height = 5, dpi = 200)
ggsave(file.path(fig_dir, "Figure_AC0206561_LYZ_locus.pdf"), p, width = 9, height = 5)
cat("\nSaved exact-coordinate locus figure (PNG+PDF) to", fig_dir, "\n")
