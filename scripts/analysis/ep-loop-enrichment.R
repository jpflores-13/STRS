## Enhancer-promoter loop analysis and visualization

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
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene |> keepStandardChromosomes()
seqlevels(noDroso_loops) <- seqlevels(txdb)
seqinfo(noDroso_loops) <- seqinfo(txdb)

## Subset loops
gainedLoops <- noDroso_loops |>
  subset(loop_type == "gained")
lostLoops <- noDroso_loops |>
  subset(loop_type == "lost")
staticLoops <- noDroso_loops |>
  subset(loop_type == "static")

# Get genes and promoters
genes <- genes(txdb)
promoter_regions <- promoters(genes)

# Load and process H3K27ac peaks
k27ac_peaks <- list.files("data/processed/cutntag/output/subsamples/output/peaks/",
                          full.names = TRUE,
                          pattern = ".narrowPeak") |> 
  str_subset("H3K27ac") |> 
  lapply(read_narrowpeaks) |> 
  lapply(keepStandardChromosomes, pruning.mode = "coarse")

# Combine peaks from replicates
k27ac_consensus <- GenomicRanges::reduce(unlist(GRangesList(k27ac_peaks)))

# Identify enhancers (H3K27ac peaks not overlapping promoters)
overlapping_peaks <- subsetByOverlaps(k27ac_consensus, promoter_regions)
enhancers <- GenomicRanges::setdiff(k27ac_consensus, overlapping_peaks)

# Function to analyze E-P loops
analyze_ep_loops <- function(loops, promoters, enhancers) {
  # Find E-P connections and get unique loops
  ep_loops <- loops |> 
    linkOverlaps(promoters, enhancers) |>
    as.data.frame() |>  
    dplyr::select(query) |>
    dplyr::distinct()
  
  # Calculate statistics
  total_loops <- length(loops)
  ep_count <- nrow(ep_loops)
  ep_fraction <- ep_count / total_loops
  
  return(list(
    total = total_loops,
    ep_count = ep_count,
    ep_fraction = ep_fraction
  ))
}

## Analyze E-P loops for each category
ep_results <- list(
  gained = analyze_ep_loops(gainedLoops, promoter_regions, enhancers),
  lost = analyze_ep_loops(lostLoops, promoter_regions, enhancers),
  static = analyze_ep_loops(staticLoops, promoter_regions, enhancers)
)

## Create data frame for viz
ep_viz_data <- data.frame(
  loop_type = c("Gained", "Lost", "Static"),
  ep_fraction = c(
    ep_results$gained$ep_fraction,
    ep_results$lost$ep_fraction,
    ep_results$static$ep_fraction
  ),
  non_ep_fraction = c(
    1 - ep_results$gained$ep_fraction,
    1 - ep_results$lost$ep_fraction,
    1 - ep_results$static$ep_fraction
  )
) |>
  pivot_longer(
    cols = c(ep_fraction, non_ep_fraction),
    names_to = "category",
    values_to = "fraction"
  )

## Modify factor levels to control stacking order
ep_viz_data$category <- factor(ep_viz_data$category, 
                               levels = c("non_ep_fraction", "ep_fraction"))

## Create theme for consistent plotting
theme_loops <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.55,
                               size = 12),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    plot.title = element_text(size = 12, face = "bold")
  )

## Add percentage labels for plotting
ep_viz_data <- ep_viz_data |>
  group_by(loop_type) |>
  mutate(
    pos = cumsum(fraction) - 0.5 * fraction,
    percentage = scales::percent(fraction, accuracy = 0.1)
  )

## Create E-P loop enrichment plot
pdf("plots/ep-loop-enrichment.pdf")

ggplot(ep_viz_data, aes(x = loop_type, y = fraction, fill = category)) +
  geom_bar(position = "stack", stat = "identity", width = 0.7) +
  geom_text(aes(y = pos, label = percentage), 
            color = "black", size = 3.5) +
  scale_fill_manual(
    values = c("lightgrey", "#33A02C"),
    labels = c("Other Loops", "E-P Loops")
  ) +
  theme_loops +
  labs(
    x = "",
    y = "",
    title = "",
    subtitle = "",
    fill = "Loop Type"
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1))

dev.off()