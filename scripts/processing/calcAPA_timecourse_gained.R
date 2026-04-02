# ##############################################################################
# filename:    calcAPA_timecourse_gained.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Calculate normalized APA matrices for Hi-C timecourse gained
#              loops across 8 timepoints (0h, 10m, 30m, 1h, 3h, 6h, 12h, 24h)
# ##############################################################################

# Libraries ----
library(mariner)
library(plotgardener)
library(glue)
library(RColorBrewer)
library(tidyverse)
library(InteractionSet)

# Parameters ----
hic_dir        <- "data/raw/hic/250212_hicTimecourseMerge/output/subsampled_hic/"
diff_loops_rds <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
output_dir     <- "data/processed/hic/normalizedAPA"
buffer         <- 10
resolution     <- 10e3

# Data import ----

## Load HiC files from timecourse
tc_hic <- list.files(hic_dir,
                     recursive = TRUE,
                     full.names = TRUE,
                     pattern = "STRS_HEK293_WT") |>
  str_subset("_inter_30.hic")

## Reorder files by timepoint
tc_hic <- tc_hic[order(factor(
  str_extract(tc_hic, "\\d+[hm]"),
  levels = c("0h", "10m", "30m", "1h", "3h", "6h", "12h", "24h")
))]

## Load differential loop calls
diffLoops <- readRDS(diff_loops_rds)
gainedLoops <- diffLoops[which(rowData(diffLoops)$padj < 0.1 &
                                 rowData(diffLoops)$log2FoldChange > 0)]
gainedLoops <- interactions(gainedLoops)

# Analysis ----
nLoops <- length(gainedLoops)

## Process a single timepoint Hi-C file
process_timepoint <- function(hic_file) {
  timepoint <- str_extract(hic_file, "(control|sorbitol)_[\\d]+[hm]")

  apa_matrix <- gainedLoops |>
    pixelsToMatrices(buffer = buffer) |>
    removeShortPairs() |>
    pullHicMatrices(binSize = resolution,
                    files = hic_file,
                    half = "upper",
                    norm = "NONE",
                    matrix = "observed") |>
    aggHicMatrices(FUN = sum)

  apa_matrix_norm <- apa_matrix / nLoops

  output_file <- file.path(output_dir,
                           glue("tc_gainedLoops_{timepoint}_normalized.rds"))
  saveRDS(apa_matrix_norm, file = output_file)

  message(glue("Processed and saved APA matrix for {timepoint}"))

  return(apa_matrix_norm)
}

# Save outputs ----
apa_matrices <- lapply(tc_hic, process_timepoint)
names(apa_matrices) <- str_extract(tc_hic, "(control|sorbitol)_[\\d]+[hm]")

sessionInfo()
