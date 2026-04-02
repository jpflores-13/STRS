# ##############################################################################
# filename:    hic_timecourse_loop_enrichment.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Calculate loop enrichment scores for gained and lost loops
#              across the Hi-C timecourse (8 timepoints); uses center pixel
#              foreground and top-left/bottom-right corner background
# ##############################################################################

# Libraries ----
library(mariner)
library(tidyverse)
library(SummarizedExperiment)

# Parameters ----
hic_dir              <- "data/raw/hic/250212_hicTimecourseMerge/output/subsampled_hic/"
diff_loops_rds       <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
enrichment_mat_dir   <- "data/processed/hic/enrichment_matrices"
output_all_rds       <- "data/processed/hic/loop_enrichment_all_original_bg.rds"
output_summary_rds   <- "data/processed/hic/loop_enrichment_summary_original_bg.rds"
buffer               <- 10
bin_size             <- 10e3

# Data import ----
timepoints      <- c("0h", "10m", "30m", "1h", "3h", "6h", "12h", "24h")
timepoint_hours <- c(0, 10/60, 30/60, 1, 3, 6, 12, 24)
names(timepoint_hours) <- timepoints

## Load HiC files
tc_hic <- list.files(hic_dir,
                     recursive = TRUE,
                     full.names = TRUE,
                     pattern = "STRS_HEK293_WT") |>
  str_subset("_inter_30.hic")

## Sort by timepoint order
hic_files <- tc_hic[order(factor(
  str_extract(tc_hic, "\\d+[hm]"),
  levels = timepoints
))]

## Load differential loop calls
diffLoops <- readRDS(diff_loops_rds)

gainedLoops <- diffLoops[which(rowData(diffLoops)$padj < 0.1 &
                                 rowData(diffLoops)$log2FoldChange > 0)]
gainedLoops <- interactions(gainedLoops)

lostLoops <- diffLoops[which(rowData(diffLoops)$padj < 0.1 &
                               rowData(diffLoops)$log2FoldChange < 0)]
lostLoops <- interactions(lostLoops)

# Analysis ----

## Calculate loop enrichment for a set of loops with center pixel foreground
## and top-left + bottom-right corner background
calculate_enrichment <- function(loops, hic_files, buffer = 10,
                                 bin_size = 10e3, prefix = "") {
  regions <- pixelsToMatrices(loops, buffer = buffer)
  regions <- removeShortPairs(regions)

  matrices <- pullHicMatrices(
    x = regions,
    files = hic_files,
    binSize = bin_size,
    norm = "NONE",
    matrix = "observed",
    blockSize = 200e6
  )

  dir.create(enrichment_mat_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(matrices, file = file.path(enrichment_mat_dir,
                                     paste0(prefix, "_all_matrices.rds")))

  fg <- selectCenterPixel(mhDist = 0, buffer = buffer)
  bg <- selectTopLeft(n = 5, buffer = buffer) +
    selectBottomRight(n = 5, buffer = buffer)

  scores <- calcLoopEnrichment(x = matrices, fg = fg, bg = bg)

  df <- as.data.frame(scores) |>
    tidyr::pivot_longer(cols = tidyr::everything(),
                        names_to = "file",
                        values_to = "score")

  df$timepoint <- str_extract(df$file, "\\d+[hm]")
  df$timepoint <- factor(df$timepoint, levels = timepoints)
  df$hours     <- timepoint_hours[as.character(df$timepoint)]

  return(df)
}

gained_enrichment <- calculate_enrichment(
  loops = gainedLoops,
  hic_files = hic_files,
  buffer = buffer,
  bin_size = bin_size,
  prefix = "gained_original_bg"
)
gained_enrichment$loop_type <- "Gained"

lost_enrichment <- calculate_enrichment(
  loops = lostLoops,
  hic_files = hic_files,
  buffer = buffer,
  bin_size = bin_size,
  prefix = "lost_original_bg"
)
lost_enrichment$loop_type <- "Lost"

all_enrichment <- bind_rows(gained_enrichment, lost_enrichment)

## Summary statistics
enrichment_summary <- all_enrichment |>
  group_by(loop_type, timepoint, hours) |>
  summarize(
    mean_score   = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score     = sd(score, na.rm = TRUE),
    se_score     = sd_score / sqrt(length(score[!is.na(score)])),
    lower_ci     = mean_score - 1.96 * se_score,
    upper_ci     = mean_score + 1.96 * se_score,
    n            = length(score[!is.na(score)]),
    .groups = "drop"
  )

# Save outputs ----
saveRDS(all_enrichment,     output_all_rds)
saveRDS(enrichment_summary, output_summary_rds)

sessionInfo()
