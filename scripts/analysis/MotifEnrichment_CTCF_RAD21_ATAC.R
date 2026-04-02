## Multi-Panel Motif Analysis Figure

pdf("plots/MotifEnrichment_CTCF_RAD21_ATAC.pdf",
    width = 11,
    height = 14)

library(plotgardener)
library(ggplot2)
library(cowplot)
library(rsvg)
library(png)
library(grid)
library(dplyr)

# Helper functions --------------------------------------------------------

parse_homer_known_results <- function(homer_dir) {
  results_file <- file.path(homer_dir, "knownResults.txt")
  
  if (!file.exists(results_file)) {
    stop(paste("Cannot find", results_file))
  }
  
  results <- read.delim(results_file, 
                        header = TRUE, 
                        stringsAsFactors = FALSE)
  
  results$neg_log_pval <- -log10(results$P.value)
  results$motif_short <- gsub("/.*", "", results$Motif.Name)
  results$motif_short <- gsub("\\(.*", "", results$motif_short)
  
  return(results)
}

get_known_logo_files <- function(homer_dir, motif_names) {
  knownResults_dir <- file.path(homer_dir, "knownResults")
  
  logo_files <- character(length(motif_names))
  
  for (i in seq_along(motif_names)) {
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

svg_to_raster <- function(svg_file) {
  if (is.na(svg_file) || !file.exists(svg_file)) {
    return(NULL)
  }
  
  temp_png <- tempfile(fileext = ".png")
  rsvg::rsvg_png(svg_file, temp_png, width = 1200, height = 300)
  img <- png::readPNG(temp_png)
  unlink(temp_png)
  
  return(img)
}

create_logo_column <- function(logo_files) {
  n_logos <- length(logo_files)
  
  p <- ggplot() + 
    xlim(0, 1) + 
    ylim(0.5, n_logos + 0.5) +
    theme_void()
  
  for (i in seq_along(logo_files)) {
    if (!is.na(logo_files[i]) && file.exists(logo_files[i])) {
      tryCatch({
        logo_raster <- svg_to_raster(logo_files[i])
        
        if (!is.null(logo_raster)) {
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

create_known_motif_plot <- function(homer_dir, top_n = 5) {
  results <- parse_homer_known_results(homer_dir)
  
  if (nrow(results) == 0) {
    stop("No motifs found in knownResults.txt")
  }
  
  n_available <- min(top_n, nrow(results))
  top_motifs <- results[1:n_available, ]
  
  logo_files <- get_known_logo_files(homer_dir, 1:n_available)
  logo_files_reversed <- rev(logo_files)
  
  top_motifs$motif_factor <- factor(top_motifs$motif_short, 
                                    levels = rev(top_motifs$motif_short))
  top_motifs$y_pos <- seq(n_available, 1, by = -1)
  
  # Create barplot matching Figure 4 aesthetics
  bar_plot <- ggplot(top_motifs, aes(x = y_pos, y = neg_log_pval)) +
    geom_col(fill = "#3182bd") +
    geom_text(aes(y = 0, label = motif_short), 
              hjust = 0, 
              nudge_y = 0.3,
              size = 2.8,
              color = "white",
              fontface = "bold") +
    geom_text(aes(label = sprintf("%.1f", neg_log_pval)), 
              hjust = -0.2, 
              size = 2.5,
              color = "white") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(expand = expansion(add = c(0.5, 0.5))) +
    labs(x = NULL, 
         y = "-log10(P-value)") +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 10, 5, 5)
    )
  
  logo_column <- create_logo_column(logo_files_reversed)
  
  combined <- plot_grid(logo_column, bar_plot, 
                        ncol = 2, 
                        rel_widths = c(1.2, 2.5),
                        align = "h",
                        axis = "tb")
  
  return(combined)
}

# Create page -------------------------------------------------------------

pageCreate(width = 11,
           height = 14,
           showGuides = FALSE)

# Define layout parameters ------------------------------------------------

## Panel dimensions
panel_width <- 4.8
panel_height <- 3.8

## Column positions
col1_x <- 0.5
col2_x <- col1_x + panel_width + 0.4

## Row positions
row1_y <- 0.8
row2_y <- row1_y + panel_height + 0.6
row3_y <- row2_y + panel_height + 0.6

## Label offset
label_offset_x <- -0.2
label_offset_y <- -0.3

# Define HOMER directories ------------------------------------------------

homer_base <- "data/processed/cutntag/homer_motifs"
homer_atac <- "data/processed/atac/homer_motifs"

homer_dirs <- list(
  ctcf_retained = file.path(homer_base, "ctcf_retained_peaks"),
  ctcf_lost = file.path(homer_base, "ctcf_lost_peaks"),
  rad21_retained = file.path(homer_base, "rad21_retained_peaks"),
  rad21_lost = file.path(homer_base, "rad21_lost_peaks"),
  atac_gained = file.path(homer_atac, "atac_gained_anchors"),
  atac_lost = file.path(homer_atac, "atac_lost_anchors")
)

# Row 1: CTCF Retained and Lost -------------------------------------------

## Panel A: CTCF Retained
plotText(label = "A", 
         x = col1_x + label_offset_x, 
         y = row1_y + label_offset_y, 
         fontsize = 14, 
         fontface = "bold")

plotText(label = "CTCF Retained Peaks",
         x = col1_x + panel_width/2,
         y = row1_y - 0.1,
         fontsize = 11,
         fontface = "bold",
         just = c("center", "bottom"))

if (file.exists(file.path(homer_dirs$ctcf_retained, "knownResults.txt"))) {
  ctcf_retained_plot <- create_known_motif_plot(
    homer_dir = homer_dirs$ctcf_retained,
    top_n = 5
  )
  
  plotGG(ctcf_retained_plot,
         x = col1_x,
         y = row1_y,
         width = panel_width,
         height = panel_height,
         just = c("left", "top"))
} else {
  plotText(label = "HOMER results not found",
           x = col1_x + panel_width/2,
           y = row1_y + panel_height/2,
           fontsize = 10,
           fontcolor = "red",
           just = c("center", "center"))
}

## Panel B: CTCF Lost
plotText(label = "B", 
         x = col2_x + label_offset_x, 
         y = row1_y + label_offset_y, 
         fontsize = 14, 
         fontface = "bold")

plotText(label = "CTCF Lost Peaks",
         x = col2_x + panel_width/2,
         y = row1_y - 0.1,
         fontsize = 11,
         fontface = "bold",
         just = c("center", "bottom"))

if (file.exists(file.path(homer_dirs$ctcf_lost, "knownResults.txt"))) {
  ctcf_lost_plot <- create_known_motif_plot(
    homer_dir = homer_dirs$ctcf_lost,
    top_n = 5
  )
  
  plotGG(ctcf_lost_plot,
         x = col2_x,
         y = row1_y,
         width = panel_width,
         height = panel_height,
         just = c("left", "top"))
} else {
  plotText(label = "HOMER results not found",
           x = col2_x + panel_width/2,
           y = row1_y + panel_height/2,
           fontsize = 10,
           fontcolor = "red",
           just = c("center", "center"))
}

# Row 2: RAD21 Retained and Lost ------------------------------------------

## Panel C: RAD21 Retained
plotText(label = "C", 
         x = col1_x + label_offset_x, 
         y = row2_y + label_offset_y, 
         fontsize = 14, 
         fontface = "bold")

plotText(label = "RAD21 Retained Peaks",
         x = col1_x + panel_width/2,
         y = row2_y - 0.1,
         fontsize = 11,
         fontface = "bold",
         just = c("center", "bottom"))

if (file.exists(file.path(homer_dirs$rad21_retained, "knownResults.txt"))) {
  rad21_retained_plot <- create_known_motif_plot(
    homer_dir = homer_dirs$rad21_retained,
    top_n = 5
  )
  
  plotGG(rad21_retained_plot,
         x = col1_x,
         y = row2_y,
         width = panel_width,
         height = panel_height,
         just = c("left", "top"))
} else {
  plotText(label = "HOMER results not found",
           x = col1_x + panel_width/2,
           y = row2_y + panel_height/2,
           fontsize = 10,
           fontcolor = "red",
           just = c("center", "center"))
}

## Panel D: RAD21 Lost
plotText(label = "D", 
         x = col2_x + label_offset_x, 
         y = row2_y + label_offset_y, 
         fontsize = 14, 
         fontface = "bold")

plotText(label = "RAD21 Lost Peaks",
         x = col2_x + panel_width/2,
         y = row2_y - 0.1,
         fontsize = 11,
         fontface = "bold",
         just = c("center", "bottom"))

if (file.exists(file.path(homer_dirs$rad21_lost, "knownResults.txt"))) {
  rad21_lost_plot <- create_known_motif_plot(
    homer_dir = homer_dirs$rad21_lost,
    top_n = 5
  )
  
  plotGG(rad21_lost_plot,
         x = col2_x,
         y = row2_y,
         width = panel_width,
         height = panel_height,
         just = c("left", "top"))
} else {
  plotText(label = "HOMER results not found",
           x = col2_x + panel_width/2,
           y = row2_y + panel_height/2,
           fontsize = 10,
           fontcolor = "red",
           just = c("center", "center"))
}

# Row 3: ATAC Gained and Lost Anchors -------------------------------------

## Panel E: ATAC Gained Anchors
plotText(label = "E", 
         x = col1_x + label_offset_x, 
         y = row3_y + label_offset_y, 
         fontsize = 14, 
         fontface = "bold")

plotText(label = "ATAC Peaks at Gained Loop Anchors",
         x = col1_x + panel_width/2,
         y = row3_y - 0.1,
         fontsize = 11,
         fontface = "bold",
         just = c("center", "bottom"))

if (file.exists(file.path(homer_dirs$atac_gained, "knownResults.txt"))) {
  atac_gained_plot <- create_known_motif_plot(
    homer_dir = homer_dirs$atac_gained,
    top_n = 5
  )
  
  plotGG(atac_gained_plot,
         x = col1_x,
         y = row3_y,
         width = panel_width,
         height = panel_height,
         just = c("left", "top"))
} else {
  plotText(label = "HOMER results not found",
           x = col1_x + panel_width/2,
           y = row3_y + panel_height/2,
           fontsize = 10,
           fontcolor = "red",
           just = c("center", "center"))
}

## Panel F: ATAC Lost Anchors
plotText(label = "F", 
         x = col2_x + label_offset_x, 
         y = row3_y + label_offset_y, 
         fontsize = 14, 
         fontface = "bold")

plotText(label = "ATAC Peaks at Lost Loop Anchors",
         x = col2_x + panel_width/2,
         y = row3_y - 0.1,
         fontsize = 11,
         fontface = "bold",
         just = c("center", "bottom"))

if (file.exists(file.path(homer_dirs$atac_lost, "knownResults.txt"))) {
  atac_lost_plot <- create_known_motif_plot(
    homer_dir = homer_dirs$atac_lost,
    top_n = 5
  )
  
  plotGG(atac_lost_plot,
         x = col2_x,
         y = row3_y,
         width = panel_width,
         height = panel_height,
         just = c("left", "top"))
} else {
  plotText(label = "HOMER results not found",
           x = col2_x + panel_width/2,
           y = row3_y + panel_height/2,
           fontsize = 10,
           fontcolor = "red",
           just = c("center", "center"))
}

dev.off()