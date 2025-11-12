## CTCF log2FC (Sorbitol/Control)
## Standalone script - loads data, analyzes, and plots

library(InteractionSet)
library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(plyranges)
library(bamsignals)

# Load data ---------------------------------------------------------------

## Load loops
loops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |>
  interactions()

gainedLoops <- loops[loops$padj < 0.1 & loops$log2FoldChange > 0]
lostLoops   <- loops[loops$padj < 0.1 & loops$log2FoldChange < 0]

## Load CTCF peaks
ctcf_control <- plyranges::read_narrowpeaks(
  "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
)

## Load CTCF BAM files
ctcf_control_bam  <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_nodups_sorted.bam"
ctcf_sorbitol_bam <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h_nodups_sorted.bam"

# Helper function ---------------------------------------------------------

get_loop_anchors_labeled <- function(loops, loop_type) {
  first_anchors  <- anchors(loops, "first")
  second_anchors <- anchors(loops, "second")
  all_anchors <- c(first_anchors, second_anchors)
  mcols(all_anchors)$loop_type <- loop_type
  mcols(all_anchors)$anchor_id <- c(
    paste0("loop", seq_along(first_anchors), "_anchor1"),
    paste0("loop", seq_along(second_anchors), "_anchor2")
  )
  all_anchors
}

# Analysis ----------------------------------------------------------------

gained_anchors <- get_loop_anchors_labeled(gainedLoops, "Gained")
lost_anchors   <- get_loop_anchors_labeled(lostLoops,  "Lost")
all_anchors    <- c(gained_anchors, lost_anchors)

## Find CTCF peaks at anchors
ctcf_at_anchors <- subsetByOverlaps(ctcf_control, all_anchors)

overlaps <- findOverlaps(ctcf_at_anchors, all_anchors)
ctcf_at_anchors <- ctcf_at_anchors[queryHits(overlaps)]
mcols(ctcf_at_anchors)$loop_type <- mcols(all_anchors)$loop_type[subjectHits(overlaps)]
mcols(ctcf_at_anchors)$anchor_id <- mcols(all_anchors)$anchor_id[subjectHits(overlaps)]

## Calculate CTCF signal (counts over peak intervals)
control_signal  <- bamCount(ctcf_control_bam,  ctcf_at_anchors, paired.end = "midpoint")
sorbitol_signal <- bamCount(ctcf_sorbitol_bam, ctcf_at_anchors, paired.end = "midpoint")

mcols(ctcf_at_anchors)$control_signal  <- control_signal
mcols(ctcf_at_anchors)$sorbitol_signal <- sorbitol_signal
mcols(ctcf_at_anchors)$log2FC <- log2((sorbitol_signal + 1) / (control_signal + 1))

## Statistical test
gained_lfc <- mcols(ctcf_at_anchors)$log2FC[mcols(ctcf_at_anchors)$loop_type == "Gained"]
lost_lfc   <- mcols(ctcf_at_anchors)$log2FC[mcols(ctcf_at_anchors)$loop_type == "Lost"]
wilcox_lfc <- wilcox.test(gained_lfc, lost_lfc)

## Prepare data for plotting
results_df <- as.data.frame(mcols(ctcf_at_anchors))

## Compute dynamic y-limits to avoid clipping
y_rng <- range(results_df$log2FC, na.rm = TRUE)
y_pad <- diff(y_rng)
if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.5
ylim_use <- c(y_rng[1] - 0.05 * y_pad, y_rng[2] + 0.05 * y_pad)

# Create plot -------------------------------------------------------------

pdf("plots/ctcf_log2fc_comparison.pdf", width = 6, height = 5)

ggplot(results_df, aes(x = loop_type, y = log2FC, fill = loop_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.5) +
  scale_fill_manual(values = c("Gained" = "#F8766D", "Lost" = "#619CFF")) +
  labs(
    title = "CTCF occupancy change at loop anchors",
    x = "Loop Type",
    y = "log\u2082FC (sorbitol / control)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
  ) +
  coord_cartesian(ylim = ylim_use)

dev.off()

cat("Plot saved to: plots/ctcf_log2fc_comparison.pdf\n")
