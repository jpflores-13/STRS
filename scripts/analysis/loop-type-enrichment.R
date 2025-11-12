# enhancer-promoter loop analysis and visualization
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
  TRUE ~ "other")

mcols(noDroso_loops)$loop_size <- pairdist(noDroso_loops)

## Add seqinfo
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene |> keepStandardChromosomes()
seqlevels(noDroso_loops) <- seqlevels(txdb)
seqinfo(noDroso_loops) <- seqinfo(txdb)

## Subset loops by type
gainedLoops <- noDroso_loops |> subset(loop_type == "gained")
lostLoops <- noDroso_loops |> subset(loop_type == "lost")
staticLoops <- noDroso_loops |> subset(loop_type == "static")

## Get genes and promoters
genes <- genes(txdb)
promoter_regions <- promoters(genes, upstream=2000, downstream=500)

## Load and process H3K27ac peaks
k27ac_peaks <- list.files("data/processed/cutntag/output/subsamples/output/peaks/",
                          full.names = TRUE,
                          pattern = ".narrowPeak") |> 
  str_subset("H3K27ac") |> 
  lapply(read_narrowpeaks) |> 
  lapply(keepStandardChromosomes, pruning.mode = "coarse")

## Combine peaks 
k27ac_consensus <- GenomicRanges::reduce(unlist(GRangesList(k27ac_peaks)))

## Identify enhancers (H3K27ac peaks not overlapping promoters)
overlapping_peaks <- subsetByOverlaps(k27ac_consensus, promoter_regions)
enhancers <- GenomicRanges::setdiff(k27ac_consensus, overlapping_peaks)

## Modified function to analyze all loop types
analyze_loop_types <- function(loops, promoters, enhancers) {
  if (length(loops) == 0) {
    return(list(
      total = 0,
      ep_count = 0,
      pp_count = 0,
      ee_count = 0,
      other_count = 0,
      ep_fraction = 0,
      pp_fraction = 0,
      ee_fraction = 0,
      other_fraction = 0
    ))
  }
  
  # Convert anchors to GRanges objects
  anchors1 <- anchors(loops, "first")
  anchors2 <- anchors(loops, "second")
  
  # Create matrices for overlap results
  p1 <- overlapsAny(anchors1, promoters)
  p2 <- overlapsAny(anchors2, promoters)
  e1 <- overlapsAny(anchors1, enhancers)
  e2 <- overlapsAny(anchors2, enhancers)
  
  # Each index represents a loop
  all_indices <- seq_len(length(loops))
  
  # Create mutually exclusive categories
  # First identify EP loops (promoter at one end, enhancer at other, but not both promoters or both enhancers)
  ep_loops <- which(
    ((p1 & e2) | (e1 & p2)) &  # Must have one promoter and one enhancer
      !(p1 & p2) &               # Must not be PP
      !(e1 & e2)                 # Must not be EE
  )
  
  # Then PP loops (must not be already categorized as EP)
  pp_loops <- which(
    (p1 & p2) &               # Both ends must be promoters
      !(e1 | e2)                # Neither end should be an enhancer
  )
  
  # Then EE loops (must not be already categorized as EP or PP)
  ee_loops <- which(
    (e1 & e2) &               # Both ends must be enhancers
      !(p1 | p2)                # Neither end should be a promoter
  )
  
  # Verify no overlap between categories
  stopifnot(
    length(intersect(ep_loops, pp_loops)) == 0,
    length(intersect(ep_loops, ee_loops)) == 0,
    length(intersect(pp_loops, ee_loops)) == 0
  )
  
  # Other loops are those that don't fit any category
  categorized <- unique(c(ep_loops, pp_loops, ee_loops))
  other_loops <- setdiff(all_indices, categorized)
  
  # Calculate counts
  total_loops <- length(loops)
  ep_count <- length(ep_loops)
  pp_count <- length(pp_loops)
  ee_count <- length(ee_loops)
  other_count <- length(other_loops)
  
  # Calculate fractions
  ep_fraction <- ep_count / total_loops
  pp_fraction <- pp_count / total_loops
  ee_fraction <- ee_count / total_loops
  other_fraction <- other_count / total_loops
  
  return(list(
    total = total_loops,
    ep_count = ep_count,
    pp_count = pp_count,
    ee_count = ee_count,
    other_count = other_count,
    ep_fraction = ep_fraction,
    pp_fraction = pp_fraction,
    ee_fraction = ee_fraction,
    other_fraction = other_fraction
  ))
}

## Analyze loops for each category
results <- list(
  gained = analyze_loop_types(gainedLoops, promoter_regions, enhancers),
  lost = analyze_loop_types(lostLoops, promoter_regions, enhancers),
  static = analyze_loop_types(staticLoops, promoter_regions, enhancers)
)

## Create data frame for visualization
viz_data <- data.frame(
  loop_type = rep(c("Gained", "Lost", "Static"), each = 4),
  category = rep(c("E-P Loops", "P-P Loops", "E-E Loops", "Other Loops"), times = 3),
  fraction = c(
    results$gained$ep_fraction, results$gained$pp_fraction, 
    results$gained$ee_fraction, results$gained$other_fraction,
    results$lost$ep_fraction, results$lost$pp_fraction, 
    results$lost$ee_fraction, results$lost$other_fraction,
    results$static$ep_fraction, results$static$pp_fraction, 
    results$static$ee_fraction, results$static$other_fraction
  )
)

## Modify factor levels to control stacking order
viz_data$category <- factor(viz_data$category, 
                            levels = c("Other Loops", "E-E Loops", "P-P Loops", "E-P Loops"))
viz_data$loop_type <- factor(viz_data$loop_type,
                             levels = c("Gained", "Lost", "Static"))

## Add percentage labels for plotting
viz_data <- viz_data |>
  group_by(loop_type) |>
  mutate(
    pos = cumsum(fraction) - 0.5 * fraction,
    percentage = scales::percent(fraction, accuracy = 0.1)
  )

## Create theme for consistent plotting
theme_loops <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

## Create stacked bar plot
pdf("plots/loop-type-enrichment.pdf", width = 8, height = 6)

ggplot(viz_data, aes(x = loop_type, y = fraction, fill = category)) +
  geom_bar(position = "stack", stat = "identity", width = 0.7) +
  geom_text(aes(y = pos, label = percentage), 
            color = "black", size = 3.5) +
  scale_fill_manual(
    values = c("lightgrey", "#FDB462", "#80B1D3", "#33A02C"),
    labels = c("Other Loops", "E-E Loops", "P-P Loops", "E-P Loops")
  ) +
  theme_loops +
  labs(
    x = "",
    y = "",
    title = "",
    fill = "Loop Type"
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0))

dev.off()