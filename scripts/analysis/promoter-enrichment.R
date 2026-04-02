# ##############################################################################
# filename:    promoter-enrichment.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Stacked barplot of the fraction of gained, lost, and static loop
#              anchors overlapping promoters; per-loop-type color coding
# ##############################################################################

# Libraries ----
library(InteractionSet)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(plyranges)
library(ggplot2)
library(scales)
library(BiocGenerics)

# Parameters ----
diff_loops_rds <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
output_pdf     <- "plots/promoter-enrichment.pdf"

# Data import ----
noDroso_loops <- readRDS(diff_loops_rds) |>
  interactions() |>
  as.data.frame() |>
  as_ginteractions()

mcols(noDroso_loops)$loop_type <- case_when(
  mcols(noDroso_loops)$padj < 0.05 & mcols(noDroso_loops)$log2FoldChange >  1 ~ "gained",
  mcols(noDroso_loops)$padj < 0.05 & mcols(noDroso_loops)$log2FoldChange < -1 ~ "lost",
  mcols(noDroso_loops)$padj > 0.05 ~ "static",
  TRUE ~ "other"
)

mcols(noDroso_loops)$loop_size <- pairdist(noDroso_loops)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene |> keepStandardChromosomes()
seqlevels(noDroso_loops) <- seqlevels(txdb)
seqinfo(noDroso_loops)   <- seqinfo(txdb)

gainedLoops <- noDroso_loops |> subset(loop_type == "gained")
lostLoops   <- noDroso_loops |> subset(loop_type == "lost")
staticLoops <- noDroso_loops |> subset(loop_type == "static")

genes            <- genes(txdb)
promoter_regions <- promoters(genes)

# Analysis ----
analyze_promoter_enrichment <- function(loops) {
  anchors_first <- GRanges(
    seqnames = Rle(seqnames1(loops)),
    ranges   = IRanges(start = start(anchors(loops, type = "first")),
                       end   = end(anchors(loops,   type = "first")))
  )
  anchors_second <- GRanges(
    seqnames = Rle(seqnames2(loops)),
    ranges   = IRanges(start = start(anchors(loops, type = "second")),
                       end   = end(anchors(loops,   type = "second")))
  )
  all_anchors       <- c(anchors_first, anchors_second)
  promoter_overlaps <- countOverlaps(all_anchors, promoter_regions) > 0
  total             <- length(all_anchors)
  n_promoter        <- sum(promoter_overlaps)
  list(total             = total,
       promoter_count    = n_promoter,
       promoter_fraction = n_promoter / total)
}

promoter_results <- list(
  gained = analyze_promoter_enrichment(gainedLoops),
  lost   = analyze_promoter_enrichment(lostLoops),
  static = analyze_promoter_enrichment(staticLoops)
)

promoter_viz_data <- data.frame(
  loop_type            = c("Gained", "Lost", "Static"),
  promoter_fraction    = c(promoter_results$gained$promoter_fraction,
                           promoter_results$lost$promoter_fraction,
                           promoter_results$static$promoter_fraction),
  non_promoter_fraction = c(1 - promoter_results$gained$promoter_fraction,
                            1 - promoter_results$lost$promoter_fraction,
                            1 - promoter_results$static$promoter_fraction)
) |>
  pivot_longer(cols = c(promoter_fraction, non_promoter_fraction),
               names_to = "category", values_to = "fraction")

promoter_viz_data$category <- factor(promoter_viz_data$category,
                                     levels = c("non_promoter_fraction", "promoter_fraction"))

promoter_viz_data <- promoter_viz_data |>
  group_by(loop_type) |>
  mutate(pos        = cumsum(fraction) - 0.5 * fraction,
         percentage = scales::percent(fraction, accuracy = 0.1),
         fill_group = interaction(category, loop_type))

# Visualization ----
theme_loops <- theme_minimal() +
  theme(
    axis.text.x        = element_text(angle = 0, hjust = 0.55, size = 12),
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
    non_promoter_fraction.Static = "lightgrey",
    promoter_fraction.Gained     = "#F8766D",
    promoter_fraction.Lost       = "#619CFF",
    promoter_fraction.Static     = "#999999"
  )) +
  scale_y_continuous(labels = scales::percent_format(suffix = ""), limits = c(0, 1)) +
  theme_loops +
  labs(x = "", y = "% of loop anchors overlapping a promoter", fill = "")

# Save outputs ----
pdf(output_pdf)
print(p)
dev.off()

sessionInfo()
