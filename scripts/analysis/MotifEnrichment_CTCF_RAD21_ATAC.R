# ##############################################################################
# filename:    MotifEnrichment_CTCF_RAD21_ATAC.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Six-panel plotgardener figure of HOMER known-motif enrichment
#              for CTCF retained/lost, RAD21 retained/lost, and ATAC peaks at
#              gained/lost loop anchors; top 5 motifs per panel
# ##############################################################################

# Libraries ----
library(plotgardener)
library(ggplot2)
library(cowplot)
library(rsvg)
library(png)
library(grid)
library(dplyr)

# Parameters ----
homer_base   <- "data/processed/cutntag/homer_motifs"
homer_atac   <- "data/processed/atac/homer_motifs"
output_pdf   <- "plots/MotifEnrichment_CTCF_RAD21_ATAC.pdf"
top_n        <- 5
page_width   <- 11
page_height  <- 14
panel_width  <- 4.8
panel_height <- 3.8

# Helper functions ----
parse_homer_known_results <- function(homer_dir) {
  results_file <- file.path(homer_dir, "knownResults.txt")
  if (!file.exists(results_file)) stop("Cannot find ", results_file)

  results <- read.delim(results_file, header = TRUE, stringsAsFactors = FALSE)
  results$neg_log_pval <- -log10(results$P.value)
  results$motif_short  <- gsub("/.*",   "", results$Motif.Name)
  results$motif_short  <- gsub("\\(.*", "", results$motif_short)
  results
}

get_known_logo_files <- function(homer_dir, n) {
  knownResults_dir <- file.path(homer_dir, "knownResults")
  vapply(seq_len(n), function(i) {
    f <- list.files(knownResults_dir,
                    pattern    = paste0("^known", i, "\\.logo\\.svg$"),
                    full.names = TRUE)
    if (length(f) > 0) f[1] else NA_character_
  }, character(1))
}

svg_to_raster <- function(svg_file) {
  if (is.na(svg_file) || !file.exists(svg_file)) return(NULL)
  tmp <- tempfile(fileext = ".png")
  rsvg::rsvg_png(svg_file, tmp, width = 1200, height = 300)
  img <- png::readPNG(tmp)
  unlink(tmp)
  img
}

create_logo_column <- function(logo_files) {
  n <- length(logo_files)
  p <- ggplot() + xlim(0, 1) + ylim(0.5, n + 0.5) + theme_void()
  for (i in seq_along(logo_files)) {
    if (!is.na(logo_files[i]) && file.exists(logo_files[i])) {
      tryCatch({
        rast <- svg_to_raster(logo_files[i])
        if (!is.null(rast)) {
          p <- p + annotation_raster(rast,
                                     xmin = 0, xmax = 1,
                                     ymin = i - 0.4, ymax = i + 0.4,
                                     interpolate = TRUE)
        }
      }, error = function(e) message("Could not read logo: ", logo_files[i]))
    }
  }
  p
}

create_known_motif_plot <- function(homer_dir, top_n = 5) {
  results     <- parse_homer_known_results(homer_dir)
  n_available <- min(top_n, nrow(results))
  top_motifs  <- results[seq_len(n_available), ]

  logo_files          <- get_known_logo_files(homer_dir, n_available)
  top_motifs$motif_factor <- factor(top_motifs$motif_short,
                                    levels = rev(top_motifs$motif_short))
  top_motifs$y_pos    <- seq(n_available, 1, by = -1)

  bar_plot <- ggplot(top_motifs, aes(x = y_pos, y = neg_log_pval)) +
    geom_col(fill = "#3182bd") +
    geom_text(aes(y = 0, label = motif_short),
              hjust = 0, nudge_y = 0.3, size = 2.8,
              color = "white", fontface = "bold") +
    geom_text(aes(label = sprintf("%.1f", neg_log_pval)),
              hjust = -0.2, size = 2.5, color = "white") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(expand = expansion(add = c(0.5, 0.5))) +
    labs(x = NULL, y = "-log10(P-value)") +
    theme_minimal() +
    theme(
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      axis.text.x        = element_text(size = 8),
      axis.title.x       = element_text(size = 9),
      axis.line.x        = element_line(color = "black", linewidth = 0.5),
      axis.line.y        = element_line(color = "black", linewidth = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin        = margin(5, 10, 5, 5)
    )

  plot_grid(create_logo_column(rev(logo_files)), bar_plot,
            ncol = 2, rel_widths = c(1.2, 2.5), align = "h", axis = "tb")
}

add_panel <- function(homer_dir, label, title, x, y, pw, ph) {
  label_offset_x <- -0.2
  label_offset_y <- -0.3

  plotText(label = label,
           x = x + label_offset_x, y = y + label_offset_y,
           fontsize = 14, fontface = "bold")
  plotText(label = title,
           x = x + pw / 2, y = y - 0.1,
           fontsize = 11, fontface = "bold",
           just = c("center", "bottom"))

  if (file.exists(file.path(homer_dir, "knownResults.txt"))) {
    plotGG(create_known_motif_plot(homer_dir, top_n = top_n),
           x = x, y = y, width = pw, height = ph,
           just = c("left", "top"))
  } else {
    plotText(label = "HOMER results not found",
             x = x + pw / 2, y = y + ph / 2,
             fontsize = 10, fontcolor = "red",
             just = c("center", "center"))
  }
}

# Layout constants ----
col1_x  <- 0.5
col2_x  <- col1_x + panel_width + 0.4
row1_y  <- 0.8
row2_y  <- row1_y + panel_height + 0.6
row3_y  <- row2_y + panel_height + 0.6

homer_dirs <- list(
  ctcf_retained  = file.path(homer_base, "ctcf_retained_peaks"),
  ctcf_lost      = file.path(homer_base, "ctcf_lost_peaks"),
  rad21_retained = file.path(homer_base, "rad21_retained_peaks"),
  rad21_lost     = file.path(homer_base, "rad21_lost_peaks"),
  atac_gained    = file.path(homer_atac,  "atac_gained_anchors"),
  atac_lost      = file.path(homer_atac,  "atac_lost_anchors")
)

# Save outputs ----
pdf(output_pdf, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = FALSE)

add_panel(homer_dirs$ctcf_retained,  "A", "CTCF Retained Peaks",              col1_x, row1_y, panel_width, panel_height)
add_panel(homer_dirs$ctcf_lost,      "B", "CTCF Lost Peaks",                  col2_x, row1_y, panel_width, panel_height)
add_panel(homer_dirs$rad21_retained, "C", "RAD21 Retained Peaks",             col1_x, row2_y, panel_width, panel_height)
add_panel(homer_dirs$rad21_lost,     "D", "RAD21 Lost Peaks",                 col2_x, row2_y, panel_width, panel_height)
add_panel(homer_dirs$atac_gained,    "E", "ATAC Peaks at Gained Loop Anchors", col1_x, row3_y, panel_width, panel_height)
add_panel(homer_dirs$atac_lost,      "F", "ATAC Peaks at Lost Loop Anchors",   col2_x, row3_y, panel_width, panel_height)

dev.off()

sessionInfo()
