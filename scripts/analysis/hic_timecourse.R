## HiC Timecourse - Hyperosmotic stress rewires chromatin interactions

## Load required libraries
library(plotgardener)
library(tidyverse)
library(RColorBrewer)

# Set up visualization parameters -----------------------------------------

## Create PDF and plotGardener page
pdf("plots/hic_timecourse.pdf",
    width = 8,    # Width to accommodate 8 matrices
    height = 2)   # Height for 2 rows with tight spacing

pageCreate(width = 8, height = 2,
           default.units = "inches", showGuides = F)

## Set up APA parameters (matching reference script)
apaParams <- pgParams(assembly = "hg38",
                      height = 0.75,
                      width = 0.75,
                      fontsize = 5,
                      palette = colorRampPalette(brewer.pal(n = 9, "YlGnBu")))

# Define time points and get files ---------------------------------------

timepoints <- c("0h", "10m", "30m", "1h", "3h", "6h", "12h", "24h")
time_order <- setNames(seq_along(timepoints), timepoints)

# Get and sort files
gained_files <- list.files("data/processed/hic/normalizedAPA",
                           full.names = TRUE,
                           pattern = "tc_gainedLoops.*normalized.rds")

lost_files <- list.files("data/processed/hic/normalizedAPA",
                         full.names = TRUE,
                         pattern = "tc_lostLoops.*normalized.rds")

# Sort files by time point order
sort_files <- function(files) {
  files[order(sapply(files, function(x) {
    time <- str_extract(x, "0h|10m|30m|1h|3h|6h|12h|24h")
    time_order[time]
  }))]
}

gained_files <- sort_files(gained_files)
lost_files <- sort_files(lost_files)

# Load matrices
gained_matrices <- lapply(gained_files, readRDS)
lost_matrices <- lapply(lost_files, readRDS)

# Plot matrices ----------------------------------------------------------

calculate_background_median <- function(matrix) {
  # Extract background regions (top-left and bottom-right 5x5 corners)
  top_left <- matrix[1:5, 1:5]
  bottom_right <- matrix[17:21, 17:21]
  
  # Combine values
  background_values <- c(as.vector(top_left), as.vector(bottom_right))
  
  # Remove NA values
  background_values <- background_values[!is.na(background_values)]
  
  # Calculate median or return a small positive value if all values are 0 or no values remain
  if(length(background_values) == 0 || all(background_values == 0)) {
    return(0.01)  # Small default value to avoid zrange issues
  } else {
    return(median(background_values))
  }
}



# Plot each timepoint with its own range
for(i in seq_along(timepoints)) {
  x_pos <- 0.2 + (i-1) * 0.9
  
  # Calculate background for each matrix
  gainedBg <- calculate_background_median(gained_matrices[[i]])
  lostBg <- calculate_background_median(lost_matrices[[i]])
  
  # Set arbitrary value
  arbVal <- 2.75
  
  # Plot gained loops matrix with legend
  plotMatrix(gained_matrices[[i]],
             x = x_pos,
             y = 0.2,
             zrange = (c(0, arbVal) * gainedBg),
             params = apaParams) |>
    annoHeatmapLegend(x = x_pos + 0.8,
                      y = 0.2,
                      width = 0.05,
                      height = 0.75,
                      fontcolor = "black",
                      fontsize = 5,
                      zrange = (c(0, arbVal) * gainedBg))
  
  # Plot lost loops matrix with legend
  plotMatrix(lost_matrices[[i]],
             x = x_pos,
             y = 1.0,
             zrange = (c(0, arbVal) * lostBg),
             params = apaParams) |>
    annoHeatmapLegend(x = x_pos + 0.8,
                      y = 1.0,
                      width = 0.05,
                      height = 0.75,
                      fontcolor = "black",
                      fontsize = 5,
                      zrange = (c(0, arbVal) * lostBg))
  
  # Add time point label
  plotText(label = timepoints[i],
           x = x_pos + 0.375,
           y = 0.1,
           params = apaParams,
           fontsize = 8,
           just = "center")
}

# Add row labels
plotText(label = "gained",  
         x = 0.1,         
         y = 0.6,
         fontcolor = "#F8766D",
         fontsize = 10,
         rot = 90)

plotText(label = "lost",    
         x = 0.1,         
         y = 1.4,
         fontcolor = "#619CFF",
         fontsize = 10,
         rot = 90)

dev.off()