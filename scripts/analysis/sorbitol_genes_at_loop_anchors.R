# ##############################################################################
# filename:    sorbitol_genes_at_loop_anchors.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Verifies that most sorbitol-responsive DEGs are NOT at gained/lost
#              loop anchors; barplot of % DEGs at vs not at differential anchors
# ##############################################################################

# Libraries ----
library(GenomicRanges)
library(InteractionSet)
library(DESeq2)
library(tidyverse)
library(mariner)

# Parameters ----
dds_rds        <- "data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds"
diff_loops_rds <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
output_pdf     <- "plots/sorbitol_genes_at_loop_anchors.pdf"
padj_cutoff    <- 0.05
lfc_cutoff     <- 2

# Data import ----
dds <- readRDS(dds_rds)

diff_loopCounts <- readRDS(diff_loops_rds) |>
  interactions() |>
  keepStandardChromosomes(pruning.mode = "coarse")

# Analysis ----
gainedLoops <- diff_loopCounts[mcols(diff_loopCounts)$padj < 0.1 &
                                 mcols(diff_loopCounts)$log2FoldChange > 0]
lostLoops   <- diff_loopCounts[mcols(diff_loopCounts)$padj < 0.1 &
                                 mcols(diff_loopCounts)$log2FoldChange < 0]

extract_unique_anchors <- function(gi) {
  c(anchors(gi, "first"), anchors(gi, "second")) |> unique()
}

gained_anchors <- extract_unique_anchors(gainedLoops)
lost_anchors   <- extract_unique_anchors(lostLoops)
diff_anchors   <- c(gained_anchors, lost_anchors) |> unique()

dge_gr <- results(dds, format = "GRanges")
mcols(dge_gr) <- rowData(dds)

dge_gr_prom <- promoters(dge_gr, upstream = 2000, downstream = 2000) |>
  keepStandardChromosomes(pruning.mode = "coarse")
seqlevelsStyle(dge_gr_prom) <- "UCSC"

time_coefficients <- resultsNames(dds) |> str_subset("Time")

extract_genes <- function(direction = c("both", "up", "down")) {
  direction <- match.arg(direction)
  lapply(time_coefficients, function(coef) {
    res <- results(dds, name = coef) |> as.data.frame()
    if (direction == "up")   res <- filter(res, padj < padj_cutoff & log2FoldChange >  lfc_cutoff)
    if (direction == "down") res <- filter(res, padj < padj_cutoff & log2FoldChange < -lfc_cutoff)
    if (direction == "both") res <- filter(res, padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff)
    rownames(res)
  }) |> unlist() |> unique()
}

all_degs      <- extract_genes("both")
upreg_genes   <- extract_genes("up")
downreg_genes <- extract_genes("down")

all_deg_prom    <- dge_gr_prom[names(dge_gr_prom) %in% all_degs]
upreg_prom      <- dge_gr_prom[names(dge_gr_prom) %in% upreg_genes]
downreg_prom    <- dge_gr_prom[names(dge_gr_prom) %in% downreg_genes]

deg_at_gained    <- subsetByOverlaps(all_deg_prom, gained_anchors)
deg_at_lost      <- subsetByOverlaps(all_deg_prom, lost_anchors)
deg_at_any_diff  <- subsetByOverlaps(all_deg_prom, diff_anchors)
deg_NOT_at_diff  <- all_deg_prom[countOverlaps(all_deg_prom, diff_anchors) == 0]

upreg_at_gained   <- subsetByOverlaps(upreg_prom,   gained_anchors)
upreg_at_lost     <- subsetByOverlaps(upreg_prom,   lost_anchors)
downreg_at_gained <- subsetByOverlaps(downreg_prom, gained_anchors)
downreg_at_lost   <- subsetByOverlaps(downreg_prom, lost_anchors)

pct <- function(n, d) round(100 * n / d, 1)

results_summary <- tibble(
  category   = c("Total DEGs",
                 "DEGs at gained anchors", "DEGs at lost anchors",
                 "DEGs at ANY differential anchor", "DEGs NOT at differential anchors",
                 NA,
                 "Upregulated genes total",
                 "Upregulated at gained anchors", "Upregulated at lost anchors",
                 NA,
                 "Downregulated genes total",
                 "Downregulated at gained anchors", "Downregulated at lost anchors"),
  count      = c(length(all_degs),
                 length(deg_at_gained), length(deg_at_lost),
                 length(deg_at_any_diff), length(deg_NOT_at_diff),
                 NA,
                 length(upreg_genes), length(upreg_at_gained), length(upreg_at_lost),
                 NA,
                 length(downreg_genes), length(downreg_at_gained), length(downreg_at_lost)),
  percentage = c(100,
                 pct(length(deg_at_gained), length(all_degs)),
                 pct(length(deg_at_lost),   length(all_degs)),
                 pct(length(deg_at_any_diff), length(all_degs)),
                 pct(length(deg_NOT_at_diff), length(all_degs)),
                 NA,
                 100,
                 pct(length(upreg_at_gained), length(upreg_genes)),
                 pct(length(upreg_at_lost),   length(upreg_genes)),
                 NA,
                 100,
                 pct(length(downreg_at_gained), length(downreg_genes)),
                 pct(length(downreg_at_lost),   length(downreg_genes)))
)

print(results_summary)

# Visualization ----
plot_data <- tibble(
  location   = factor(c("At anchors", "Not at anchors"),
                      levels = c("At anchors", "Not at anchors")),
  count      = c(length(deg_at_any_diff), length(deg_NOT_at_diff)),
  percentage = c(pct(length(deg_at_any_diff), length(all_degs)),
                 pct(length(deg_NOT_at_diff), length(all_degs)))
)

p <- ggplot(plot_data, aes(x = location, y = percentage, fill = location)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("n = %d\n%.1f%%", count, percentage)),
            vjust = -0.3, size = 5) +
  scale_fill_manual(values = c("At anchors" = "#7570B3", "Not at anchors" = "#E7298A")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Sorbitol-responsive genes at differential loop anchors",
       x = NULL, y = "Percentage of all DEGs (%)") +
  theme_minimal() +
  theme(
    legend.position    = "none",
    axis.text.x        = element_text(size = 12, color = "black"),
    axis.text.y        = element_text(size = 11, color = "black"),
    axis.title.y       = element_text(size = 12),
    plot.title         = element_text(size = 14, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.background   = element_rect(fill = "white", color = NA)
  )

# Save outputs ----
ggsave(output_pdf, p, width = 8, height = 6)

sessionInfo()
