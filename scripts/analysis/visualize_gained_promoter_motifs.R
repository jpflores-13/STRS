# visualize_gained_promoter_motifs_known_with_logos.R
# Visualize KNOWN motifs with logos

library(ggplot2)
library(cowplot)
library(dplyr)
library(rsvg)
library(png)
library(grid)

# Function to parse HOMER known results ----
parse_homer_known_results <- function(homer_dir) {
  results_file <- file.path(homer_dir, "knownResults.txt")
  
  if (!file.exists(results_file)) {
    stop(paste("Cannot find", results_file))
  }
  
  # Read HOMER known results (tab-delimited)
  results <- read.delim(results_file, 
                        header = TRUE, 
                        stringsAsFactors = FALSE)
  
  # Calculate -log10(p-value) for plotting
  results$neg_log_pval <- -log10(results$P.value)
  
  # Extract motif name (remove extra info in parentheses)
  results$motif_short <- gsub("/.*", "", results$Motif.Name)
  results$motif_short <- gsub("\\(.*", "", results$motif_short)
  
  return(results)
}

# Function to get known motif logo file paths ----
get_known_logo_files <- function(homer_dir, motif_names) {
  knownResults_dir <- file.path(homer_dir, "knownResults")
  
  logo_files <- character(length(motif_names))
  
  for (i in seq_along(motif_names)) {
    # HOMER saves known motifs as known#.logo.svg
    pattern <- paste0("^known", i, "\\.logo\\.svg$")
    logo_file <- list.files(knownResults_dir, pattern = pattern, full.names = TRUE)
    
    if (length(logo_file) > 0) {
      logo_files[i] <- logo_file[1]
    } else {
      logo_files[i] <- NA
    }
  }
  
  return(logo_files)
}

# Function to convert SVG to raster for plotting ----
svg_to_raster <- function(svg_file) {
  if (is.na(svg_file) || !file.exists(svg_file)) {
    return(NULL)
  }
  
  # Create temporary PNG file
  temp_png <- tempfile(fileext = ".png")
  
  # Convert SVG to PNG using rsvg
  rsvg::rsvg_png(svg_file, temp_png, width = 1200, height = 300)
  
  # Read the PNG
  img <- png::readPNG(temp_png)
  
  # Clean up temp file
  unlink(temp_png)
  
  return(img)
}

# Function to create a combined logo column plot ----
create_logo_column <- function(logo_files) {
  n_logos <- length(logo_files)
  
  # Create a blank plot with the right dimensions
  p <- ggplot() + 
    xlim(0, 1) + 
    ylim(0.5, n_logos + 0.5) +  # Add 0.5 padding for better centering
    theme_void()
  
  # Add each logo as an annotation_raster - centered on each integer position
  for (i in seq_along(logo_files)) {
    if (!is.na(logo_files[i]) && file.exists(logo_files[i])) {
      tryCatch({
        logo_raster <- svg_to_raster(logo_files[i])
        
        if (!is.null(logo_raster)) {
          # Center each logo at position i (with 0.4 height for better fit)
          p <- p + annotation_raster(logo_raster,
                                     xmin = 0, xmax = 1,
                                     ymin = i - 0.4, ymax = i + 0.4,
                                     interpolate = TRUE)
        }
      }, error = function(e) {
        message(paste("Could not read logo:", logo_files[i]))
      })
    }
  }
  
  return(p)
}

# Function to create motif visualization ----
create_known_motif_plot <- function(homer_dir, top_n = 15, condition_name = "") {
  
  # Parse results
  results <- parse_homer_known_results(homer_dir)
  
  # Check if any results exist
  if (nrow(results) == 0) {
    stop("No motifs found in knownResults.txt")
  }
  
  # Summary stats
  cat("\n========================================\n")
  cat("HOMER Known Motif Results Summary\n")
  cat("========================================\n")
  cat("Total motifs tested:", nrow(results), "\n")
  cat("Motifs with P < 1e-10:", sum(results$P.value < 1e-10, na.rm = TRUE), "\n")
  cat("Motifs with P < 1e-5:", sum(results$P.value < 1e-5, na.rm = TRUE), "\n")
  cat("Motifs with P < 0.001:", sum(results$P.value < 0.001, na.rm = TRUE), "\n")
  cat("Motifs with P < 0.05:", sum(results$P.value < 0.05, na.rm = TRUE), "\n\n")
  
  # Take top N motifs
  n_available <- min(top_n, nrow(results))
  top_motifs <- results[1:n_available, ]
  
  # Get logo files
  logo_files <- get_known_logo_files(homer_dir, 1:n_available)
  
  # Reverse order for plotting (bottom to top)
  logo_files_reversed <- rev(logo_files)
  
  # Create factor with reverse order for plotting
  top_motifs$motif_factor <- factor(top_motifs$motif_short, 
                                    levels = rev(top_motifs$motif_short))
  
  # Add numeric position for alignment
  top_motifs$y_pos <- seq(n_available, 1, by = -1)
  
  # Create barplot with NAVY BLUE and WHITE text
  bar_plot <- ggplot(top_motifs, aes(x = y_pos, y = neg_log_pval)) +
    geom_col(fill = "#3182bd", color = "#08519c", linewidth = 0.3) +
    geom_text(aes(y = 0, label = motif_short), 
              hjust = 0, 
              nudge_y = 0.3,
              size = 3.5,
              color = "white",
              fontface = "bold") +
    geom_text(aes(label = sprintf("%.1f", neg_log_pval)), 
              hjust = -0.2, 
              size = 2.5,
              color = "white") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(expand = expansion(add = c(0.5, 0.5))) +  # Match logo padding
    labs(x = NULL, 
         y = "-log10(P-value)",
         title = condition_name) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.margin = margin(5, 10, 5, 5)
    )
  
  # Create combined logo column
  logo_column <- create_logo_column(logo_files_reversed)
  
  # Combine logo column with barplot using cowplot
  combined <- plot_grid(logo_column, bar_plot, 
                        ncol = 2, 
                        rel_widths = c(1.2, 2.5),
                        align = "h",
                        axis = "tb")
  
  return(combined)
}

# Main execution ----------------------------------------------------------
homer_dir <- "data/processed/cutntag/homer_motifs/gained_anchor_promoters"
output_dir <- "plots"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Check what files exist
cat("Checking HOMER output files...\n")
cat("knownResults.txt exists:", 
    file.exists(file.path(homer_dir, "knownResults.txt")), "\n")
cat("homerResults.txt exists:", 
    file.exists(file.path(homer_dir, "homerResults.txt")), "\n\n")

if (!file.exists(file.path(homer_dir, "knownResults.txt"))) {
  stop("knownResults.txt not found! Has HOMER finished running?")
}

# Generate plot with known motifs
message("Generating gained anchor promoter KNOWN motif plot with logos...")
gained_promoter_plot <- create_known_motif_plot(
  homer_dir = homer_dir,
  top_n = 15,
  condition_name = "Promoters at Gained Loop Anchors"
)

# Save
ggsave(file.path(output_dir, "gained_anchor_promoters_known_motifs_top15_with_logos.pdf"),
       gained_promoter_plot,
       width = 12,
       height = 8)

message("✓ Visualization complete!")
message("Plot saved to: plots/gained_anchor_promoters_known_motifs_top15_with_logos.pdf")

# Print top 15 motifs with details
cat("\n========================================\n")
cat("Top 15 Enriched Known Motifs\n")
cat("========================================\n")
results <- parse_homer_known_results(homer_dir)
print(results[1:min(15, nrow(results)), 
              c("Motif.Name", "P.value", "q.value..Benjamini.", 
                "X..of.Target.Sequences.with.Motif", 
                "X..of.Background.Sequences.with.Motif")])