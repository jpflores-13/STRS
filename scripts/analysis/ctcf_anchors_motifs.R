# ##############################################################################
# filename:    ctcf_anchors_motifs.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Visualize HOMER de novo motif enrichment for retained vs lost
#              CTCF peaks; shows top 5 de novo motifs for each condition
# ##############################################################################

# Libraries ----
library(ggplot2)
library(cowplot)
library(rsvg)
library(png)
library(grid)
library(dplyr)

# Parameters ----
homer_motifs_dir <- "/work/users/j/p/jpflores/projects/STRS/data/processed/cutntag/homer_motifs"
output_dir       <- "plots"
top_n_motifs     <- 5

# Analysis ----

## Parse HOMER de novo results ----
parse_homer_results <- function(homer_dir) {
  results_file <- file.path(homer_dir, "homerResults.txt")

  if (!file.exists(results_file)) {
    stop(paste("Cannot find", results_file))
  }

  results <- read.delim(results_file, header = TRUE, stringsAsFactors = FALSE)
  results$neg_log_pval <- -log10(results$P.value)

  return(results)
}

## Get de novo motif logo file paths ----
get_logo_files <- function(homer_dir, n_motifs) {
  homerResults_dir <- file.path(homer_dir, "homerResults")

  logo_files <- list.files(homerResults_dir,
                           pattern = "^motif\\d+\\.logo\\.svg$",
                           full.names = TRUE)

  if (length(logo_files) == 0) {
    warning(paste("No logo files found in", homerResults_dir))
    return(rep(NA, n_motifs))
  }

  logo_indices <- as.numeric(gsub(".*motif(\\d+)\\.logo\\.svg", "\\1", basename(logo_files)))
  logo_files   <- logo_files[order(logo_indices)]

  return(logo_files[1:min(n_motifs, length(logo_files))])
}

## Convert SVG to raster for plotting ----
svg_to_raster <- function(svg_file) {
  temp_png <- tempfile(fileext = ".png")
  rsvg::rsvg_png(svg_file, temp_png, width = 1200, height = 300)
  img <- png::readPNG(temp_png)
  unlink(temp_png)
  return(img)
}

## Create combined logo column ----
create_logo_column <- function(logo_files) {
  n_logos <- length(logo_files)

  p <- ggplot() +
    xlim(0, 1) +
    ylim(0, n_logos) +
    theme_void()

  for (i in seq_along(logo_files)) {
    if (!is.na(logo_files[i]) && file.exists(logo_files[i])) {
      tryCatch({
        logo_raster <- svg_to_raster(logo_files[i])
        p <- p + annotation_raster(logo_raster, xmin = 0, xmax = 1,
                                   ymin = i - 1, ymax = i, interpolate = TRUE)
      }, error = function(e) {
        message(paste("Could not read logo:", logo_files[i]))
      })
    }
  }

  return(p)
}

## Create motif panel for one condition ----
create_motif_panel <- function(homer_dir, top_n = 5, condition_name = "") {

  results    <- parse_homer_results(homer_dir)
  n_available <- min(top_n, nrow(results))
  top_motifs  <- results[1:n_available, ]

  logo_files <- get_logo_files(homer_dir, n_available)

  if (length(logo_files) > nrow(top_motifs)) {
    logo_files <- logo_files[1:nrow(top_motifs)]
  } else if (length(logo_files) < nrow(top_motifs)) {
    logo_files <- c(logo_files, rep(NA, nrow(top_motifs) - length(logo_files)))
  }

  top_motifs$logo_file   <- logo_files
  top_motifs$motif_short <- gsub("\\(.*", "", top_motifs$Motif.Name)
  top_motifs$motif_factor <- factor(top_motifs$motif_short,
                                    levels = rev(top_motifs$motif_short))

  logo_files_reversed <- rev(logo_files)

  bar_plot <- ggplot(top_motifs, aes(x = motif_factor, y = neg_log_pval)) +
    geom_col(fill = "#3182bd", color = "#08519c", linewidth = 0.3) +
    geom_text(aes(label = motif_short), hjust = 1.1, size = 3,
              color = "white", fontface = "bold") +
    geom_text(aes(label = sprintf("%.1f", neg_log_pval)), hjust = -0.2, size = 2.5) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = NULL, y = "-log10(P-value)", title = condition_name) +
    theme_minimal() +
    theme(
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.margin        = margin(5, 10, 5, 5)
    )

  logo_column <- create_logo_column(logo_files_reversed)

  combined <- plot_grid(logo_column, bar_plot, ncol = 2,
                        rel_widths = c(1.5, 2), align = "h")

  return(combined)
}

## Combined retained vs lost plot ----
plot_ctcf_retained_lost <- function(base_dir, output_file, top_n = 5) {

  retained_dir <- file.path(base_dir, "ctcf_retained_peaks")
  lost_dir     <- file.path(base_dir, "ctcf_lost_peaks")

  retained_panel <- create_motif_panel(retained_dir, top_n, "Retained CTCF Peaks")
  lost_panel     <- create_motif_panel(lost_dir,     top_n, "Lost CTCF Peaks")

  combined_panels <- plot_grid(retained_panel, lost_panel, ncol = 2)

  title <- ggdraw() +
    draw_label("De Novo Motif Enrichment: Retained vs Lost CTCF Peaks",
               fontface = "bold", size = 16)

  final_plot <- plot_grid(title, combined_panels, ncol = 1,
                          rel_heights = c(0.05, 1))

  ggsave(output_file, final_plot, width = 14, height = 8)
  message(paste("CTCF retained vs lost plot saved to:", output_file))
}

# Save outputs ----
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("Generating CTCF retained vs lost de novo motif comparison plot...")
plot_ctcf_retained_lost(
  base_dir    = homer_motifs_dir,
  output_file = file.path(output_dir, "ctcf_retained_lost_denovo_motifs_top5.pdf"),
  top_n       = top_n_motifs
)

sessionInfo()
