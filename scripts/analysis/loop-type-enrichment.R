# ##############################################################################
# filename:    loop-type-enrichment.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Stacked barplot of loop-type fractions (E-P, P-P, E-E, Other)
#              for gained, lost, and static loops; enhancers defined as H3K27ac
#              peaks not overlapping promoters (±2 kb/500 bp window)
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
cutntag_peaks_dir <- "data/processed/cutntag/output/subsamples/output/peaks/"
output_pdf        <- "plots/loop-type-enrichment.pdf"

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
promoter_regions <- promoters(genes, upstream = 2000, downstream = 500)

k27ac_peaks <- list.files(cutntag_peaks_dir, full.names = TRUE,
                           pattern = ".narrowPeak") |>
  stringr::str_subset("H3K27ac") |>
  lapply(plyranges::read_narrowpeaks) |>
  lapply(keepStandardChromosomes, pruning.mode = "coarse")

k27ac_consensus   <- GenomicRanges::reduce(unlist(GRangesList(k27ac_peaks)))
overlapping_peaks <- subsetByOverlaps(k27ac_consensus, promoter_regions)
enhancers         <- GenomicRanges::setdiff(k27ac_consensus, overlapping_peaks)

# Analysis ----
analyze_loop_types <- function(loops, promoters, enhancers) {
  if (length(loops) == 0) {
    return(list(total = 0, ep_count = 0, pp_count = 0, ee_count = 0,
                other_count = 0, ep_fraction = 0, pp_fraction = 0,
                ee_fraction = 0, other_fraction = 0))
  }

  anchors1 <- anchors(loops, "first")
  anchors2 <- anchors(loops, "second")

  p1 <- overlapsAny(anchors1, promoters)
  p2 <- overlapsAny(anchors2, promoters)
  e1 <- overlapsAny(anchors1, enhancers)
  e2 <- overlapsAny(anchors2, enhancers)

  ep_loops <- which(((p1 & e2) | (e1 & p2)) & !(p1 & p2) & !(e1 & e2))
  pp_loops <- which((p1 & p2) & !(e1 | e2))
  ee_loops <- which((e1 & e2) & !(p1 | p2))

  stopifnot(
    length(intersect(ep_loops, pp_loops)) == 0,
    length(intersect(ep_loops, ee_loops)) == 0,
    length(intersect(pp_loops, ee_loops)) == 0
  )

  other_loops <- setdiff(seq_len(length(loops)),
                         unique(c(ep_loops, pp_loops, ee_loops)))

  total <- length(loops)
  list(total       = total,
       ep_count    = length(ep_loops),    pp_count    = length(pp_loops),
       ee_count    = length(ee_loops),    other_count = length(other_loops),
       ep_fraction = length(ep_loops)    / total,
       pp_fraction = length(pp_loops)    / total,
       ee_fraction = length(ee_loops)    / total,
       other_fraction = length(other_loops) / total)
}

loop_results <- list(
  gained = analyze_loop_types(gainedLoops, promoter_regions, enhancers),
  lost   = analyze_loop_types(lostLoops,   promoter_regions, enhancers),
  static = analyze_loop_types(staticLoops, promoter_regions, enhancers)
)

viz_data <- data.frame(
  loop_type = rep(c("Gained", "Lost", "Static"), each = 4),
  category  = rep(c("E-P Loops", "P-P Loops", "E-E Loops", "Other Loops"), times = 3),
  fraction  = c(
    loop_results$gained$ep_fraction, loop_results$gained$pp_fraction,
    loop_results$gained$ee_fraction, loop_results$gained$other_fraction,
    loop_results$lost$ep_fraction,   loop_results$lost$pp_fraction,
    loop_results$lost$ee_fraction,   loop_results$lost$other_fraction,
    loop_results$static$ep_fraction, loop_results$static$pp_fraction,
    loop_results$static$ee_fraction, loop_results$static$other_fraction
  )
)

viz_data$category  <- factor(viz_data$category,
                              levels = c("Other Loops", "E-E Loops", "P-P Loops", "E-P Loops"))
viz_data$loop_type <- factor(viz_data$loop_type, levels = c("Gained", "Lost", "Static"))

viz_data <- viz_data |>
  group_by(loop_type) |>
  mutate(pos        = cumsum(fraction) - 0.5 * fraction,
         percentage = scales::percent(fraction, accuracy = 0.1))

# Visualization ----
theme_loops <- theme_minimal() +
  theme(
    axis.text.x        = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y        = element_text(size = 10),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "top",
    legend.text        = element_text(size = 10),
    legend.title       = element_text(size = 12)
  )

p <- ggplot(viz_data, aes(x = loop_type, y = fraction, fill = category)) +
  geom_bar(position = "stack", stat = "identity", width = 0.7) +
  geom_text(aes(y = pos, label = percentage), color = "black", size = 3.5) +
  scale_fill_manual(values = c("lightgrey", "#FDB462", "#80B1D3", "#33A02C"),
                    labels = c("Other Loops", "E-E Loops", "P-P Loops", "E-P Loops")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0)) +
  theme_loops +
  labs(x = "", y = "", fill = "Loop Type")

# Save outputs ----
pdf(output_pdf, width = 8, height = 6)
print(p)
dev.off()

sessionInfo()
