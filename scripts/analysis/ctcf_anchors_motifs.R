# HOMER De Novo Motif Visualization: Retained vs Lost CTCF Peaks
# Compares de novo motif enrichment between retained and lost CTCF peaks
# Shows top 5 de novo motifs for each condition

# Load required libraries
library(ggplot2)
library(cowplot)
library(rsvg)
library(png)
library(grid)
library(dplyr)

# Function to parse HOMER de novo results ----
parse_homer_results <- function(homer_dir) {
  results_file <- file.path(homer_dir, "homerResults.txt")
  
  if (!file.exists(results_file)) {
    stop(paste("Cannot find", results_file))
  }
  
  # Read HOMER de novo results (tab-delimited)
  results <- read.delim(results_file, 
                        header = TRUE, 
                        stringsAsFactors = FALSE)
  
  # Calculate -log10(p-value) for plotting
  results$neg_log_pval <- -log10(results$P.value)
  
  return(results)
}

# Function to get de novo motif logo file paths ----
get_logo_files <- function(homer_dir, n_motifs) {
  homerResults_dir <- file.path(homer_dir, "homerResults")
  
  # HOMER saves de novo logos as SVG files named motif#.logo.svg
  logo_files <- list.files(homerResults_dir, 
                           pattern = "^motif\\d+\\.logo\\.svg$", 
                           full.names = TRUE)
  
  if (length(logo_files) == 0) {
    warning(paste("No logo files found in", homerResults_dir))
    return(rep(NA, n_motifs))
  }
  
  # Extract numeric indices and sort
  logo_indices <- as.numeric(gsub(".*motif(\\d+)\\.logo\\.svg", "\\1", basename(logo_files)))
  logo_files <- logo_files[order(logo_indices)]
  
  # Return logos matching the number requested
  return(logo_files[1:min(n_motifs, length(logo_files))])
}

# Function to convert SVG to raster for plotting ----
svg_to_raster <- function(svg_file) {
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
    ylim(0, n_logos) +
    theme_void()
  
  # Add each logo as an annotation_raster
  for (i in seq_along(logo_files)) {
    if (!is.na(logo_files[i]) && file.exists(logo_files[i])) {
      tryCatch({
        logo_raster <- svg_to_raster(logo_files[i])
        
        # Position from bottom to top (i-1 to i)
        p <- p + annotation_raster(logo_raster,
                                   xmin = 0, xmax = 1,
                                   ymin = i - 1, ymax = i,
                                   interpolate = TRUE)
      }, error = function(e) {
        message(paste("Could not read logo:", logo_files[i]))
      })
    }
  }
  
  return(p)
}

# Function to create motif visualization for one condition ----
create_motif_panel <- function(homer_dir, top_n = 5, condition_name = "") {
  
  # Parse results
  results <- parse_homer_results(homer_dir)
  
  # Take only as many motifs as are available (may be less than top_n)
  n_available <- min(top_n, nrow(results))
  top_motifs <- results[1:n_available, ]
  
  # Get logo files - only request as many as we have motifs
  logo_files <- get_logo_files(homer_dir, n_available)
  
  # Make sure logo_files matches the number of motifs
  if (length(logo_files) > nrow(top_motifs)) {
    logo_files <- logo_files[1:nrow(top_motifs)]
  } else if (length(logo_files) < nrow(top_motifs)) {
    # Pad with NAs if we have fewer logos than motifs
    logo_files <- c(logo_files, rep(NA, nrow(top_motifs) - length(logo_files)))
  }
  
  # Add logo paths to data
  top_motifs$logo_file <- logo_files
  
  # Extract motif name (first part before parenthesis)
  top_motifs$motif_short <- gsub("\\(.*", "", top_motifs$Motif.Name)
  
  # Create factor with reverse order for plotting
  top_motifs$motif_factor <- factor(top_motifs$motif_short, 
                                    levels = rev(top_motifs$motif_short))
  
  # Reverse logo files to match plotting order (bottom to top)
  logo_files_reversed <- rev(logo_files)
  
  # Create barplot
  bar_plot <- ggplot(top_motifs, aes(x = motif_factor, y = neg_log_pval)) +
    geom_col(fill = "#3182bd", color = "#08519c", linewidth = 0.3) +
    geom_text(aes(label = motif_short), 
              hjust = 1.1, 
              size = 3,
              color = "white",
              fontface = "bold") +
    geom_text(aes(label = sprintf("%.1f", neg_log_pval)), 
              hjust = -0.2, 
              size = 2.5) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = NULL, 
         y = "-log10(P-value)",
         title = condition_name) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.margin = margin(5, 10, 5, 5)
    )
  
  # Create combined logo column
  logo_column <- create_logo_column(logo_files_reversed)
  
  # Combine logo column with barplot using cowplot
  combined <- plot_grid(logo_column, bar_plot, 
                        ncol = 2, 
                        rel_widths = c(1.5, 2),
                        align = "h")
  
  return(combined)
}

# Main plotting function ----
plot_ctcf_retained_lost <- function(base_dir, output_file, top_n = 5) {
  
  # Set up paths
  retained_dir <- file.path(base_dir, "ctcf_retained_peaks")
  lost_dir <- file.path(base_dir, "ctcf_lost_peaks")
  
  # Create panels
  retained_panel <- create_motif_panel(retained_dir, top_n, "Retained CTCF Peaks")
  lost_panel <- create_motif_panel(lost_dir, top_n, "Lost CTCF Peaks")
  
  # Combine side by side
  combined_panels <- plot_grid(retained_panel, lost_panel, ncol = 2)
  
  # Add title
  title <- ggdraw() + 
    draw_label("De Novo Motif Enrichment: Retained vs Lost CTCF Peaks", 
               fontface = "bold", size = 16)
  
  final_plot <- plot_grid(title, combined_panels, 
                          ncol = 1, 
                          rel_heights = c(0.05, 1))
  
  # Save
  ggsave(output_file, final_plot, width = 14, height = 8)
  message(paste("CTCF retained vs lost plot saved to:", output_file))
}

# Main execution ----
base_dir <- "/work/users/j/p/jpflores/projects/STRS/data/processed/cutntag/homer_motifs"
output_dir <- "plots"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Generate plot with top 5 de novo motifs
message("Generating CTCF retained vs lost de novo motif comparison plot...")
plot_ctcf_retained_lost(
  base_dir = base_dir,
  output_file = file.path(output_dir, "ctcf_retained_lost_denovo_motifs_top5.pdf"),
  top_n = 5
)

message("✓ Visualization complete!")
