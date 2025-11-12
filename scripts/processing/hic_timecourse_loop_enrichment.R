## HiC Timecourse Loop Enrichment Analysis

## Load required libraries
library(mariner)
library(tidyverse)
library(SummarizedExperiment)

# Define timepoints
timepoints <- c("0h", "10m", "30m", "1h", "3h", "6h", "12h", "24h")
time_order <- setNames(seq_along(timepoints), timepoints)

# Convert timepoints to numeric hours for data organization
timepoint_hours <- c(0, 10/60, 30/60, 1, 3, 6, 12, 24)
names(timepoint_hours) <- timepoints

## Load HiC files
tc_hic <- list.files("data/raw/hic/250212_hicTimecourseMerge/output/subsampled_hic/", 
                     recursive = TRUE,
                     full.names = TRUE,
                     pattern = "STRS_HEK293_WT") |> 
  str_subset("_inter_30.hic")

# Sort the files based on timepoint order
hic_files <- tc_hic[order(factor(
  str_extract(tc_hic, "\\d+[hm]"),
  levels = timepoints
))]

# Load differential loop calls
diffLoops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds")

# Extract gained and lost loops
gainedLoops <- diffLoops[which(rowData(diffLoops)$padj < 0.1 & 
                                 rowData(diffLoops)$log2FoldChange > 0)]
gainedLoops <- interactions(gainedLoops)

lostLoops <- diffLoops[which(rowData(diffLoops)$padj < 0.1 & 
                               rowData(diffLoops)$log2FoldChange < 0)]
lostLoops <- interactions(lostLoops)

# Function to calculate loop enrichment for a set of loops with specified background
calculate_enrichment <- function(loops, hic_files, buffer = 10, bin_size = 10e3, prefix = "") {
  # Convert loops to matrices with buffer
  regions <- pixelsToMatrices(loops, buffer = buffer)
  
  # Remove short pairs that might bias the analysis
  regions <- removeShortPairs(regions)
  
  # Pull Hi-C matrices for each timepoint
  matrices <- pullHicMatrices(
    x = regions,
    files = hic_files,
    binSize = bin_size,
    norm = "NONE",
    matrix = "observed",
    blockSize = 200e6
  )
  
  # Create directory for saving matrices
  dir.create("data/processed/hic/enrichment_matrices", showWarnings = FALSE, recursive = TRUE)
  
  # Save the complete matrices object
  saveRDS(matrices, file = paste0("data/processed/hic/enrichment_matrices/", prefix, "_all_matrices.rds"))
  
  # Define the original background - top-left and bottom-right corners
  # This background controls well for distance-dependent decay of interaction frequency
  fg <- 
    selectCenterPixel(mhDist = 0, buffer = buffer)
  bg <- 
    selectTopLeft(n=5, buffer=buffer) +
    selectBottomRight(n=5, buffer=buffer)
  
  # Calculate loop enrichment with specified background
  # Using default: median(foreground + 1)/median(background + 1)
  scores <- calcLoopEnrichment(
    x = matrices,
    fg = fg,
    bg = bg
  )
  
  # Create dataframe with timepoints
  df <- as.data.frame(scores) |>
    tidyr::pivot_longer(cols = tidyr::everything(), 
                        names_to = "file", 
                        values_to = "score")
  
  # Extract timepoint from filename
  df$timepoint <- str_extract(df$file, "\\d+[hm]")
  
  # Ensure timepoints are in correct order
  df$timepoint <- factor(df$timepoint, levels = timepoints)
  
  # Add numeric hours for data organization
  df$hours <- timepoint_hours[as.character(df$timepoint)]
  
  return(df)
}

# Calculate enrichment for gained and lost loops
gained_enrichment <- calculate_enrichment(
  loops = gainedLoops, 
  hic_files = hic_files, 
  buffer = 10, 
  prefix = "gained_original_bg"
)
gained_enrichment$loop_type <- "Gained"

lost_enrichment <- calculate_enrichment(
  loops = lostLoops, 
  hic_files = hic_files, 
  buffer = 10, 
  prefix = "lost_original_bg"
)
lost_enrichment$loop_type <- "Lost"

# Combine data
all_enrichment <- bind_rows(gained_enrichment, lost_enrichment)

# Calculate summary statistics
enrichment_summary <- all_enrichment |>
  group_by(loop_type, timepoint, hours) |>
  summarize(
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    se_score = sd_score / sqrt(length(score[!is.na(score)])),
    lower_ci = mean_score - 1.96 * se_score,
    upper_ci = mean_score + 1.96 * se_score,
    n = length(score[!is.na(score)]),
    .groups = "drop"
  )

# Save the enrichment data
saveRDS(all_enrichment, "data/processed/hic/loop_enrichment_all_original_bg.rds")
saveRDS(enrichment_summary, "data/processed/hic/loop_enrichment_summary_original_bg.rds")