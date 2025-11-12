## Promoter enrichment at loop anchors

library(InteractionSet)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(plyranges)
library(ggplot2)
library(scales)
library(BiocGenerics)

## Load Hi-C differential loops
noDroso_loops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions() |> 
  as.data.frame() |>
  as_ginteractions()

## Add loop type and size metadata
mcols(noDroso_loops)$loop_type <- case_when(
  mcols(noDroso_loops)$padj < 0.05 & mcols(noDroso_loops)$log2FoldChange > 1 ~ "gained",
  mcols(noDroso_loops)$padj < 0.05 & mcols(noDroso_loops)$log2FoldChange < -1 ~ "lost",
  mcols(noDroso_loops)$padj > 0.05 ~ "static",
  is.character("NA") ~ "other")

mcols(noDroso_loops)$loop_size <- pairdist(noDroso_loops)

## Add seqinfo
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene |> 
  keepStandardChromosomes()
seqlevels(noDroso_loops) <- seqlevels(txdb)
seqinfo(noDroso_loops) <- seqinfo(txdb)

## Subset loops by type
gainedLoops <- noDroso_loops |> 
  subset(loop_type == "gained")
lostLoops <- noDroso_loops |> 
  subset(loop_type == "lost")
staticLoops <- noDroso_loops |> 
  subset(loop_type == "static")

## Get genes and promoters
genes <- genes(txdb)
promoter_regions <- promoters(genes)

## Helper function for promoter enrichment for loops
analyze_promoter_enrichment <- function(loops) {
  # Get all anchors
  anchors_first <- GRanges(
    seqnames = Rle(seqnames1(loops)),
    ranges = IRanges(
      start = start(anchors(loops, type = "first")),
      end = end(anchors(loops, type = "first"))
    )
  )
  
  anchors_second <- GRanges(
    seqnames = Rle(seqnames2(loops)),
    ranges = IRanges(
      start = start(anchors(loops, type = "second")),
      end = end(anchors(loops, type = "second"))
    )
  )
  
  # Combine all anchors
  all_anchors <- c(anchors_first, anchors_second)
  
  # Find overlaps with promoters
  promoter_overlaps <- countOverlaps(all_anchors, promoter_regions) > 0
  
  # Calculate statistics
  total_anchors <- length(all_anchors)
  promoter_anchors <- sum(promoter_overlaps)
  promoter_fraction <- promoter_anchors / total_anchors
  
  return(list(
    total = total_anchors,
    promoter_count = promoter_anchors,
    promoter_fraction = promoter_fraction
  ))
}

## Analyze promoter enrichment for each category
promoter_results <- list(
  gained = analyze_promoter_enrichment(gainedLoops),
  lost = analyze_promoter_enrichment(lostLoops),
  static = analyze_promoter_enrichment(staticLoops)
)

## Create data frame for viz
promoter_viz_data <- data.frame(
  loop_type = c("Gained", "Lost", "Static"),
  promoter_fraction = c(
    promoter_results$gained$promoter_fraction,
    promoter_results$lost$promoter_fraction,
    promoter_results$static$promoter_fraction
  ),
  non_promoter_fraction = c(
    1 - promoter_results$gained$promoter_fraction,
    1 - promoter_results$lost$promoter_fraction,
    1 - promoter_results$static$promoter_fraction
  )
) |>
  pivot_longer(
    cols = c(promoter_fraction, non_promoter_fraction),
    names_to = "category",
    values_to = "fraction"
  )

## Modify factor levels to control stacking order
promoter_viz_data$category <- factor(promoter_viz_data$category, 
                                     levels = c("non_promoter_fraction", "promoter_fraction"))

## Create theme for consistent plotting
theme_loops <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.55,
                               size = 12),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 12, face = "bold")
  )

## Add percentage labels for plotting
promoter_viz_data <- promoter_viz_data |>
  group_by(loop_type) |>
  mutate(
    pos = cumsum(fraction) - 0.5 * fraction,
    percentage = scales::percent(fraction, accuracy = 0.1)
  )

## Create promoter enrichment plot
pdf("plots/promoter-enrichment.pdf")

# Create interaction between category and loop_type for custom coloring
promoter_viz_data$fill_group <- interaction(promoter_viz_data$category, promoter_viz_data$loop_type)

ggplot(promoter_viz_data, aes(x = loop_type, y = fraction, fill = fill_group)) +
  geom_bar(position = "stack", stat = "identity", width = 0.7) +
  geom_text(aes(y = pos, label = percentage), 
            color = "black", size = 3.5) +
  scale_fill_manual(
    values = c(
      non_promoter_fraction.Gained = "lightgrey",
      non_promoter_fraction.Lost = "lightgrey",
      non_promoter_fraction.Static = "lightgrey",
      promoter_fraction.Gained = "#F8766D",    # Red for Gained
      promoter_fraction.Lost = "#619CFF",      # Blue for Lost
      promoter_fraction.Static = "#999999"     # Neutral gray for Static
    )
  ) +
  theme_loops +
  labs(
    x = "",
    y = "% of loop anchors overlapping a promoter",
    title = "",
    fill = ""
  ) +
  # Use scales::percent_format() with the suffix parameter set to empty string
  scale_y_continuous(labels = scales::percent_format(suffix = ""), limits = c(0, 1))

dev.off()
