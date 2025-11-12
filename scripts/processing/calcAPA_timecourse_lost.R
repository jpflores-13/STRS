## Create APA plots for Hi-C Timecourse

## load packages
library(mariner)
library(plotgardener)
library(glue)
library(RColorBrewer)
library(tidyverse)
library(InteractionSet)

# Load data ---------------------------------------------------------------

## Load HiC files
tc_hic <- list.files("data/raw/hic/250212_hicTimecourseMerge/output/subsampled_hic/", 
                     recursive = T,
                     full.names = T,
                     pattern = "STRS_HEK293_WT") |> 
  str_subset("_inter_30.hic")

## Reorder files
tc_hic <- tc_hic[order(factor(
  str_extract(tc_hic, "\\d+[hm]"),
  levels = c("0h", "10m", "30m", "1h", "3h", "6h", "12h", "24h")
))]

## Load differential loop calls
diffLoops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds")
lostLoops <- diffLoops[which(rowData(diffLoops)$padj < 0.1 &
                                 rowData(diffLoops)$log2FoldChange < 0)]
lostLoops <- interactions(lostLoops)

# calculate APA matrices --------------------------------------------------
## Parameters
buffer <- 10
resolution <- 10e3
nLoops <- length(lostLoops)

## Function to process single timepoint
process_timepoint <- function(hic_file) {
  # Extract timepoint and condition from filename
  timepoint <- str_extract(hic_file, "(control|sorbitol)_[\\d]+[hm]")
  
  # Calculate APA matrix for this timepoint
  apa_matrix <- lostLoops |> 
    pixelsToMatrices(buffer = buffer) |> 
    removeShortPairs() |> 
    pullHicMatrices(binSize = resolution,
                    files = hic_file,
                    half = "upper",
                    norm = "NONE",
                    matrix = "observed") |>
    aggHicMatrices(FUN = sum)
  
  # Normalize by number of loops
  apa_matrix_norm <- apa_matrix/nLoops
  
  # Save normalized APA matrix
  output_file <- glue("data/processed/hic/normalizedAPA/tc_lostLoops_{timepoint}_normalized.rds")
  saveRDS(apa_matrix_norm, file = output_file)
  
  message(glue("Processed and saved APA matrix for {timepoint}"))
  
  return(apa_matrix_norm)
}

# Process each timepoint
apa_matrices <- lapply(tc_hic, process_timepoint)

# Name the matrices list with timepoints for easier reference
names(apa_matrices) <- str_extract(tc_hic, "(control|sorbitol)_[\\d]+[hm]")