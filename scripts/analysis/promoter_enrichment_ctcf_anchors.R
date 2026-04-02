# ##############################################################################
# filename:    promoter_enrichment_ctcf_anchors.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Stacked barplot of promoter overlap at CTCF-bound gained vs lost
#              loop anchors; unique anchors, seqlevel harmonization, Fisher test
# ##############################################################################

# Libraries ----
library(InteractionSet)
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(plyranges)
library(GenomeInfoDb)

# Parameters ----
diff_loops_rds  <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
ctcf_peaks_file <- "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
output_pdf      <- "plots/promoter_enrichment_ctcf_anchors.pdf"

# Data import ----
loops <- readRDS(diff_loops_rds) |> interactions()

gainedLoops <- loops[loops$padj < 0.1 & loops$log2FoldChange > 0]
lostLoops   <- loops[loops$padj < 0.1 & loops$log2FoldChange < 0]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
ctcf <- plyranges::read_narrowpeaks(ctcf_peaks_file)

# Analysis ----
get_loop_anchors_labeled <- function(loops, loop_type) {
  all_anchors <- c(anchors(loops, "first"), anchors(loops, "second"))
  mcols(all_anchors)$loop_type <- loop_type
  all_anchors
}

analyze_promoter_enrichment <- function(anchors) {
  promoter_regions  <- promoters(genes(txdb))
  promoter_overlaps <- countOverlaps(anchors, promoter_regions) > 0
  total             <- length(anchors)
  n_promoter        <- sum(promoter_overlaps)
  list(total                = total,
       promoter_count        = n_promoter,
       promoter_fraction     = n_promoter / total,
       non_promoter_fraction = 1 - n_promoter / total)
}

gained_anchors <- unique(get_loop_anchors_labeled(gainedLoops, "Gained"))
lost_anchors   <- unique(get_loop_anchors_labeled(lostLoops,   "Lost"))

seqlevelsStyle(gained_anchors) <- "UCSC"
seqlevelsStyle(lost_anchors)   <- "UCSC"
seqlevelsStyle(ctcf)           <- "UCSC"

gained_ctcf <- gained_anchors |> join_overlap_inner(ctcf) |> unique()
lost_ctcf   <- lost_anchors   |> join_overlap_inner(ctcf) |> unique()

cat("Unique anchors:\n")
cat("  Gained (all):      ", length(gained_anchors), "\n")
cat("  Lost (all):        ", length(lost_anchors),   "\n")
cat("  Gained CTCF-bound: ", length(gained_ctcf),    "\n")
cat("  Lost CTCF-bound:   ", length(lost_ctcf),      "\n\n")

promoter_results <- list(
  Gained = analyze_promoter_enrichment(gained_ctcf),
  Lost   = analyze_promoter_enrichment(lost_ctcf)
)

contingency_table <- matrix(
  c(promoter_results$Gained$promoter_count,
    promoter_results$Gained$total - promoter_results$Gained$promoter_count,
    promoter_results$Lost$promoter_count,
    promoter_results$Lost$total   - promoter_results$Lost$promoter_count),
  nrow = 2,
  dimnames = list(c("Gained", "Lost"), c("Promoter", "Non-promoter"))
)

print(fisher.test(contingency_table))

# Visualization ----
promoter_viz_data <- data.frame(
  loop_type             = c("Gained", "Lost"),
  promoter_fraction     = c(promoter_results$Gained$promoter_fraction,
                            promoter_results$Lost$promoter_fraction),
  non_promoter_fraction = c(promoter_results$Gained$non_promoter_fraction,
                            promoter_results$Lost$non_promoter_fraction)
) |>
  pivot_longer(cols = c(promoter_fraction, non_promoter_fraction),
               names_to = "category", values_to = "fraction") |>
  mutate(category = factor(category,
                           levels = c("non_promoter_fraction", "promoter_fraction"))) |>
  group_by(loop_type) |>
  mutate(pos        = cumsum(fraction) - 0.5 * fraction,
         percentage = scales::percent(fraction, accuracy = 0.1),
         fill_group = interaction(category, loop_type))

theme_loops <- theme_minimal() +
  theme(
    axis.text.x        = element_text(angle = 0, hjust = 0.5, size = 12),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "none"
  )

p <- ggplot(promoter_viz_data, aes(x = loop_type, y = fraction, fill = fill_group)) +
  geom_bar(position = "stack", stat = "identity", width = 0.7) +
  geom_text(aes(y = pos, label = percentage), color = "black", size = 3.5) +
  scale_fill_manual(values = c(
    non_promoter_fraction.Gained = "lightgrey",
    non_promoter_fraction.Lost   = "lightgrey",
    promoter_fraction.Gained     = "#F8766D",
    promoter_fraction.Lost       = "#619CFF"
  )) +
  scale_y_continuous(labels = scales::percent_format(suffix = ""), limits = c(0, 1)) +
  theme_loops +
  labs(x = "", y = "% overlapping a promoter",
       title = "Promoter Overlap at CTCF-Bound Loop Anchors", fill = "")

# Save outputs ----
pdf(output_pdf, width = 6, height = 5)
print(p)
dev.off()

sessionInfo()
