## HiC Timecourse APA Matrix Enrichment Analysis
## This script calculates loop enrichment scores from pre-computed APA matrices
## using the Original Background (Top-Left + Bottom-Right) selection strategy

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

# Function to calculate loop enrichment for pre-computed APA matrices
calculate_enrichment <- function(loop_type = c("gained", "lost"), buffer = 10, prefix = "") {
  
  loop_type <- match.arg(loop_type)
  
  # Load pre-computed APA matrices instead of Hi-C files
  message("Loading pre-computed APA matrices...")
  pattern <- paste0("tc_", loop_type, "Loops.*normalized\\.rds")
  apa_files <- list.files("data/processed/hic/normalizedAPA",
                          full.names = TRUE,
                          pattern = pattern)
  
  if (length(apa_files) == 0) {
    stop("No APA matrix files found for ", loop_type, " loops")
  }
  
  # Sort the files based on timepoint order
  apa_files <- apa_files[order(sapply(apa_files, function(x) {
    time <- str_extract(x, "0h|10m|30m|1h|3h|6h|12h|24h")
    time_order[time]
  }))]
  
  message("Found ", length(apa_files), " APA matrix files for ", loop_type, " loops")
  
  # Load APA matrices - this replaces the pixelsToMatrices -> pullHicMatrices pipeline
  message("Loading APA matrices (replacing individual loop matrix extraction)...")
  apa_matrices <- vector("list", length(apa_files))
  timepoint_names <- character(length(apa_files))
  
  for (i in seq_along(apa_files)) {
    timepoint <- str_extract(apa_files[i], "0h|10m|30m|1h|3h|6h|12h|24h")
    timepoint_names[i] <- timepoint
    message("  Loading timepoint: ", timepoint)
    
    apa_matrices[[i]] <- readRDS(apa_files[i])
  }
  
  # Create directory for saving matrices
  dir.create("data/processed/hic/apa_enrichment_matrices", showWarnings = FALSE, recursive = TRUE)
  
  # Save the loaded APA matrices
  message("Saving loaded APA matrices...")
  saveRDS(apa_matrices, file = paste0("data/processed/hic/apa_enrichment_matrices/", prefix, "_apa_matrices.rds"))
  
  # Define the original background - top-left and bottom-right corners
  # This background controls well for distance-dependent decay of interaction frequency
  message("Setting up foreground and background selections...")
  fg <- 
    selectCenterPixel(mhDist = 0, buffer = buffer)
  bg <- 
    selectTopLeft(n=5, buffer=buffer) +
    selectBottomRight(n=5, buffer=buffer)
  
  # Calculate loop enrichment with specified background
  # Using default: median(foreground + 1)/median(background + 1)
  # Note: Since we're working with pre-computed APA matrices, we implement the same logic manually
  message("Calculating enrichment scores (using median calculation)...")
  
  # Manual implementation of calcLoopEnrichment logic for APA matrices
  enrichment_scores <- numeric(length(apa_matrices))
  
  for (i in seq_along(apa_matrices)) {
    matrix <- apa_matrices[[i]]
    
    # Apply foreground selection (selectCenterPixel equivalent)
    center_idx <- buffer + 1
    foreground_value <- matrix[center_idx, center_idx]
    
    # Apply background selection (selectTopLeft + selectBottomRight equivalent)
    # Top-left corner (n=5)
    top_left_values <- as.vector(matrix[1:5, 1:5])
    
    # Bottom-right corner (n=5)
    matrix_size <- nrow(matrix)
    bottom_right_start <- matrix_size - 4
    bottom_right_values <- as.vector(matrix[bottom_right_start:matrix_size, 
                                            bottom_right_start:matrix_size])
    
    # Combine background values and remove any NA values
    background_values <- c(top_left_values, bottom_right_values)
    background_values <- background_values[!is.na(background_values)]
    
    # Calculate enrichment using the same formula as calcLoopEnrichment default
    if (length(background_values) > 0) {
      enrichment_scores[i] <- (foreground_value) / (median(background_values))
    } else {
      enrichment_scores[i] <- NA
    }
  }
  
  # Create dataframe with timepoints
  message("Creating dataframe with enrichment scores...")
  df <- tibble(
    timepoint = timepoint_names,
    score = enrichment_scores
  ) |>
    # Create file column to match original script format
    mutate(file = paste0(loop_type, "_", timepoint))
  
  # Convert to long format to match original script structure
  df <- df |>
    pivot_longer(cols = score, 
                 names_to = "metric", 
                 values_to = "score") |>
    select(-metric) |>
    mutate(file = paste0(file, "_apa_matrix"))
  
  # Extract timepoint and ensure proper ordering
  df$timepoint <- str_extract(df$file, "\\d+[hm]")
  df$timepoint <- factor(df$timepoint, levels = timepoints)
  
  # Add numeric hours for data organization
  df$hours <- timepoint_hours[as.character(df$timepoint)]
  
  return(df)
}

# Calculate enrichment for gained loops using APA matrices
message("========================================")
message("Calculating enrichment for gained loops using pre-computed APA matrices...")
message("Buffer size: 10 pixels (from APA matrix dimensions)")
message("========================================")
gained_enrichment <- calculate_enrichment(
  loop_type = "gained", 
  buffer = 10, 
  prefix = "gained_apa_bg"
)
gained_enrichment$loop_type <- "Gained"

message("========================================")
message("Calculating enrichment for lost loops using pre-computed APA matrices...")
message("========================================")
lost_enrichment <- calculate_enrichment(
  loop_type = "lost", 
  buffer = 10, 
  prefix = "lost_apa_bg"
)
lost_enrichment$loop_type <- "Lost"

# Combine data
message("Combining gained and lost loop enrichment data...")
all_enrichment <- bind_rows(gained_enrichment, lost_enrichment)

# Calculate summary statistics
message("Calculating summary statistics...")
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
message("Saving enrichment data...")
saveRDS(all_enrichment, "data/processed/hic/apa_loop_enrichment_all_original_bg.rds")
saveRDS(enrichment_summary, "data/processed/hic/apa_loop_enrichment_summary_original_bg.rds")

# Print summary of what was done
message("========================================")
message("Analysis complete using pre-computed APA matrices!")
message("Method: APA matrix enrichment with Original Background selection")
message("Buffer size used: 10 pixels")
message("Background selection: Top-Left + Bottom-Right corners")
message("========================================")
message("- Raw enrichment data saved to: apa_loop_enrichment_all_original_bg.rds")
message("- Summary statistics saved to: apa_loop_enrichment_summary_original_bg.rds")
message("- APA matrices saved to: apa_enrichment_matrices/")
message("========================================")