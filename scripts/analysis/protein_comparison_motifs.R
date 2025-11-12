# HOMER Motif Visualization: Protein Comparison
# Compares motif enrichment for CTCF, RAD21, H3K27ac, and YAP1
# across control vs sorbitol treatment

# Load required libraries
library(ggplot2)
library(cowplot)
library(rsvg)
library(png)
library(grid)
library(dplyr)

# Function to parse HOMER results ----
parse_homer_results <- function(homer_dir) {
  results_file <- file.path(homer_dir, "knownResults.txt")
  
  if (!file.exists(results_file)) {
    stop(paste("Cannot find", results_file))
  }
  
  # Read HOMER results (tab-delimited)
  results <- read.delim(results_file, 
                        header = TRUE, 
                        stringsAsFactors = FALSE)
  
  # Calculate -log10(p-value) for plotting
  results$neg_log_pval <- -log10(results$P.value)
  
  return(results)
}

# Function to get logo file paths ----
get_logo_files <- function(homer_dir, n_motifs) {
  knownResults_dir <- file.path(homer_dir, "knownResults")
  
  # HOMER saves logos as SVG files named known#.logo.svg
  logo_files <- list.files(knownResults_dir, 
                           pattern = "^known\\d+\\.logo\\.svg$", 
                           full.names = TRUE)
  
  if (length(logo_files) == 0) {
    warning(paste("No logo files found in", knownResults_dir))
    return(rep(NA, n_motifs))
  }
  
  # Extract numeric indices and sort
  logo_indices <- as.numeric(gsub(".*known(\\d+)\\.logo\\.svg", "\\1", basename(logo_files)))
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
create_motif_panel <- function(homer_dir, top_n = 15, condition_name = "") {
  
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
  
  # Create barplot with motif names on y-axis
  bar_plot <- ggplot(top_motifs, aes(x = motif_factor, y = neg_log_pval)) +
    geom_col(fill = "#3182bd", color = "#08519c", linewidth = 0.3) +
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
      axis.text.y = element_text(size = 9, hjust = 1, face = "bold"),
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
                        rel_widths = c(1.5, 2.5),
                        align = "h")
  
  return(combined)
}

# Main plotting function ----
plot_protein_comparison <- function(base_dir, output_file, top_n = 15) {
  
  proteins <- c("CTCF", "RAD21", "H3K27ac", "YAP1")
  
  # Create all panels
  all_rows <- list()
  
  for (protein in proteins) {
    protein_lower <- tolower(protein)
    
    # Control panel
    control_dir <- file.path(base_dir, "top500", paste0(protein_lower, "_control"))
    if (dir.exists(control_dir)) {
      control_panel <- create_motif_panel(control_dir, top_n, "Control")
    } else {
      control_panel <- ggplot() + theme_void()
    }
    
    # Sorbitol panel
    sorbitol_dir <- file.path(base_dir, "top500", paste0(protein_lower, "_sorbitol"))
    if (dir.exists(sorbitol_dir)) {
      sorbitol_panel <- create_motif_panel(sorbitol_dir, top_n, "Sorbitol")
    } else {
      sorbitol_panel <- ggplot() + theme_void()
    }
    
    # Combine control and sorbitol for this protein
    protein_row <- plot_grid(control_panel, sorbitol_panel, ncol = 2)
    
    # Add protein label
    protein_label <- ggdraw() + 
      draw_label(protein, fontface = "bold", size = 14, x = 0.05, hjust = 0)
    
    protein_row_with_label <- plot_grid(protein_label, protein_row, 
                                        ncol = 1, 
                                        rel_heights = c(0.05, 1))
    
    all_rows[[protein]] <- protein_row_with_label
  }
  
  # Stack all proteins vertically
  combined_rows <- plot_grid(plotlist = all_rows, ncol = 1)
  
  # Add main title
  title <- ggdraw() + 
    draw_label("Top 15 Enriched Motifs: Control vs Sorbitol Treatment", 
               fontface = "bold", size = 16)
  
  final_plot <- plot_grid(title, combined_rows, 
                          ncol = 1, 
                          rel_heights = c(0.03, 1))
  
  # Save
  ggsave(output_file, final_plot, width = 16, height = 20)
  message(paste("Protein comparison plot saved to:", output_file))
}

# Main execution ----
base_dir <- "/work/users/j/p/jpflores/projects/STRS/data/processed/cutntag/homer_motifs"
output_dir <- "plots"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Generate plot
message("Generating protein comparison plot...")
plot_protein_comparison(
  base_dir = base_dir,
  output_file = file.path(output_dir, "protein_comparison_motifs.pdf"),
  top_n = 15
)

message("✓ Visualization complete!")
