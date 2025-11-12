## Promoter enrichment at CTCF-bound loop anchors
## Standalone script - loads data, analyzes, and plots
## Minimal edits: unique anchors, restrict to CTCF-bound, harmonize seqlevels

library(InteractionSet)
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(plyranges)
library(GenomeInfoDb)   # <- added for seqlevelsStyle()

# Load data ---------------------------------------------------------------

## Load loops
loops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions()

## Gained and lost loops
gainedLoops <- loops[loops$padj < 0.1 & loops$log2FoldChange > 0]
lostLoops   <- loops[loops$padj < 0.1 & loops$log2FoldChange < 0]

## Load TxDb
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## NEW: Load CTCF peaks (for “CTCF-bound” restriction)
ctcf <- plyranges::read_narrowpeaks(
  "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
)

# Helper functions --------------------------------------------------------

get_loop_anchors_labeled <- function(loops, loop_type) {
  first_anchors  <- anchors(loops, "first")
  second_anchors <- anchors(loops, "second")
  all_anchors <- c(first_anchors, second_anchors)
  mcols(all_anchors)$loop_type <- loop_type
  return(all_anchors)
}

## unchanged — but we’ll pass CTCF-bound anchors into this
analyze_promoter_enrichment <- function(anchors) {
  genes <- genes(txdb)
  promoter_regions <- promoters(genes)  # default: TSS -2000/+200
  promoter_overlaps <- countOverlaps(anchors, promoter_regions) > 0
  total_anchors <- length(anchors)
  promoter_anchors <- sum(promoter_overlaps)
  promoter_fraction <- promoter_anchors / total_anchors
  
  return(list(
    total = total_anchors,
    promoter_count = promoter_anchors,
    promoter_fraction = promoter_fraction,
    non_promoter_fraction = 1 - promoter_fraction
  ))
}

# Analysis ----------------------------------------------------------------

## Original anchor extraction
gained_anchors <- get_loop_anchors_labeled(gainedLoops, "Gained")
lost_anchors   <- get_loop_anchors_labeled(lostLoops,  "Lost")

## NEW: make anchors unique (avoid counting the same site multiple times)
gained_anchors <- unique(gained_anchors)
lost_anchors   <- unique(lost_anchors)

## NEW: harmonize chromosome naming to be safe (chr1 vs 1 issues)
seqlevelsStyle(gained_anchors) <- "UCSC"
seqlevelsStyle(lost_anchors)   <- "UCSC"
seqlevelsStyle(ctcf)           <- "UCSC"

## NEW: restrict to CTCF-bound anchors (intersection with CTCF peaks)
gained_ctcf <- gained_anchors %>% join_overlap_inner(ctcf) %>% unique()
lost_ctcf   <- lost_anchors   %>% join_overlap_inner(ctcf)   %>% unique()

cat("Counts (unique anchors):\n")
cat("  Gained (all anchors): ", length(gained_anchors), "\n")
cat("  Lost (all anchors):   ", length(lost_anchors),   "\n")
cat("  Gained CTCF-bound:    ", length(gained_ctcf),    "\n")
cat("  Lost CTCF-bound:      ", length(lost_ctcf),      "\n\n")

## Run enrichment on the CTCF-bound sets
promoter_results <- list(
  Gained = analyze_promoter_enrichment(gained_ctcf),
  Lost   = analyze_promoter_enrichment(lost_ctcf)
)

## Fisher's exact test
contingency_table <- matrix(
  c(promoter_results$Gained$promoter_count, 
    promoter_results$Gained$total - promoter_results$Gained$promoter_count,
    promoter_results$Lost$promoter_count,
    promoter_results$Lost$total - promoter_results$Lost$promoter_count),
  nrow = 2,
  dimnames = list(c("Gained", "Lost"), c("Promoter", "Non-promoter"))
)

fisher_test <- fisher.test(contingency_table)
print(fisher_test)  # <- NEW: show p-value, OR, CI in console

# Prepare visualization data ----------------------------------------------

promoter_viz_data <- data.frame(
  loop_type = c("Gained", "Lost"),
  promoter_fraction = c(
    promoter_results$Gained$promoter_fraction,
    promoter_results$Lost$promoter_fraction
  ),
  non_promoter_fraction = c(
    promoter_results$Gained$non_promoter_fraction,
    promoter_results$Lost$non_promoter_fraction
  )
) |>
  pivot_longer(
    cols = c(promoter_fraction, non_promoter_fraction),
    names_to = "category",
    values_to = "fraction"
  )

promoter_viz_data$category <- factor(
  promoter_viz_data$category, 
  levels = c("non_promoter_fraction", "promoter_fraction")
)

promoter_viz_data <- promoter_viz_data |>
  group_by(loop_type) |>
  mutate(
    pos = cumsum(fraction) - 0.5 * fraction,
    percentage = scales::percent(fraction, accuracy = 0.1)
  )

promoter_viz_data$fill_group <- interaction(
  promoter_viz_data$category, 
  promoter_viz_data$loop_type
)

# Theme -------------------------------------------------------------------

theme_loops <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )

# Create plot -------------------------------------------------------------

pdf("plots/promoter_enrichment_ctcf_anchors.pdf", width = 6, height = 5)

ggplot(promoter_viz_data, aes(x = loop_type, y = fraction, fill = fill_group)) +
  geom_bar(position = "stack", stat = "identity", width = 0.7) +
  geom_text(
    aes(y = pos, label = percentage), 
    color = "black", 
    size = 3.5
  ) +
  scale_fill_manual(
    values = c(
      non_promoter_fraction.Gained = "lightgrey",
      non_promoter_fraction.Lost = "lightgrey",
      promoter_fraction.Gained = "#F8766D",
      promoter_fraction.Lost = "#619CFF"
    )
  ) +
  theme_loops +
  labs(
    x = "",
    y = "% overlapping a promoter",  # <- updated label
    title = "Promoter Overlap at CTCF-Bound Loop Anchors",       # <- updated title
    subtitle = "",                                               # keep empty per original
    fill = ""
  ) +
  scale_y_continuous(
    labels = scales::percent_format(suffix = ""), 
    limits = c(0, 1)
  )

dev.off()

cat("Plot saved to: plots/promoter_enrichment_ctcf_anchors.pdf\n")
