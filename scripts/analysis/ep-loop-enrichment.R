# ##############################################################################
# filename:    ep-loop-enrichment.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Stacked barplot showing fraction of gained, lost, and static
#              loops that are enhancer-promoter (E-P) loops; enhancers defined
#              as H3K27ac peaks not overlapping promoters
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
diff_loops_rds    <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
cutntag_peaks_dir <- "data/processed/cutntag/output/peaks/"
output_pdf        <- "plots/ep-loop-enrichment.pdf"

# Data import ----
noDroso_loops <- readRDS(diff_loops_rds) |>
  interactions() |>
  as.data.frame() |>
  as_ginteractions()

mcols(noDroso_loops)$loop_type <- case_when(
  mcols(noDroso_loops)$padj <= 0.1 & mcols(noDroso_loops)$log2FoldChange > 0 ~ "gained",
  mcols(noDroso_loops)$padj <= 0.1 & mcols(noDroso_loops)$log2FoldChange < 0 ~ "lost",
  mcols(noDroso_loops)$padj >  0.1 ~ "static",
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

k27ac_peaks <- list.files(cutntag_peaks_dir, full.names = TRUE, pattern = ".narrowPeak") |>
  stringr::str_subset("H3K27ac") |>
  lapply(plyranges::read_narrowpeaks) |>
  lapply(keepStandardChromosomes, pruning.mode = "coarse")

k27ac_consensus  <- GenomicRanges::reduce(unlist(GRangesList(k27ac_peaks)))
overlapping_peaks <- subsetByOverlaps(k27ac_consensus, promoter_regions)
enhancers         <- GenomicRanges::setdiff(k27ac_consensus, overlapping_peaks)

# Analysis ----

analyze_ep_loops <- function(loops, promoters, enhancers) {
  ep_loops <- loops |>
    linkOverlaps(promoters, enhancers) |>
    as.data.frame() |>
    dplyr::select(query) |>
    dplyr::distinct()

  total_loops <- length(loops)
  ep_count    <- nrow(ep_loops)

  list(total = total_loops, ep_count = ep_count,
       ep_fraction = ep_count / total_loops)
}

ep_results <- list(
  gained = analyze_ep_loops(gainedLoops, promoter_regions, enhancers),
  lost   = analyze_ep_loops(lostLoops,   promoter_regions, enhancers),
  static = analyze_ep_loops(staticLoops, promoter_regions, enhancers)
)

ep_viz_data <- data.frame(
  loop_type    = c("Gained", "Lost", "Static"),
  ep_fraction  = c(ep_results$gained$ep_fraction,
                   ep_results$lost$ep_fraction,
                   ep_results$static$ep_fraction),
  non_ep_fraction = c(1 - ep_results$gained$ep_fraction,
                      1 - ep_results$lost$ep_fraction,
                      1 - ep_results$static$ep_fraction)
) |>
  pivot_longer(cols = c(ep_fraction, non_ep_fraction),
               names_to = "category", values_to = "fraction")

ep_viz_data$category <- factor(ep_viz_data$category,
                               levels = c("non_ep_fraction", "ep_fraction"))

ep_viz_data <- ep_viz_data |>
  group_by(loop_type) |>
  mutate(pos        = cumsum(fraction) - 0.5 * fraction,
         percentage = scales::percent(fraction, accuracy = 0.1))

# Visualization ----
theme_loops <- theme_minimal() +
  theme(axis.text.x   = element_text(angle = 0, hjust = 0.55, size = 12),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position    = "top",
        plot.title         = element_text(size = 12, face = "bold"))

p <- ggplot(ep_viz_data, aes(x = loop_type, y = fraction, fill = category)) +
  geom_bar(position = "stack", stat = "identity", width = 0.7) +
  geom_text(aes(y = pos, label = percentage), color = "black", size = 3.5) +
  scale_fill_manual(values = c("lightgrey", "#807DBA"),
                    labels = c("Other Loops", "E-P Loops")) +
  theme_loops +
  labs(x = "", y = "", fill = "Loop Type") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1))

# Save outputs ----
pdf(output_pdf)
print(p)
dev.off()

sessionInfo()
