# ##############################################################################
# filename:    visualize_gained_promoter_motifs.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Visualizes HOMER known motif enrichment results for promoters
#              at gained loop anchors; barplot of top N motifs with logo column
# ##############################################################################

# Libraries ----
library(ggplot2)
library(cowplot)
library(dplyr)
library(rsvg)
library(png)
library(grid)

# Parameters ----
homer_dir  <- "data/processed/cutntag/homer_motifs/gained_anchor_promoters"
output_pdf <- "plots/gained_anchor_promoters_known_motifs_top15_with_logos.pdf"
top_n      <- 15L

# Analysis ----
parse_homer_known_results <- function(homer_dir) {
  results_file <- file.path(homer_dir, "knownResults.txt")
  if (!file.exists(results_file))
    stop(paste("Cannot find", results_file))
  results <- read.delim(results_file, header = TRUE, stringsAsFactors = FALSE)
  results$neg_log_pval <- -log10(results$P.value)
  results$motif_short  <- gsub("/.*|\\(.*", "", results$Motif.Name)
  results
}

get_known_logo_files <- function(homer_dir, motif_indices) {
  knownResults_dir <- file.path(homer_dir, "knownResults")
  vapply(motif_indices, function(i) {
    f <- list.files(knownResults_dir,
                    pattern = paste0("^known", i, "\\.logo\\.svg$"),
                    full.names = TRUE)
    if (length(f) > 0) f[1] else NA_character_
  }, character(1))
}

svg_to_raster <- function(svg_file) {
  if (is.na(svg_file) || !file.exists(svg_file)) return(NULL)
  temp_png <- tempfile(fileext = ".png")
  rsvg::rsvg_png(svg_file, temp_png, width = 1200, height = 300)
  img <- png::readPNG(temp_png)
  unlink(temp_png)
  img
}

create_logo_column <- function(logo_files) {
  n <- length(logo_files)
  p <- ggplot() + xlim(0, 1) + ylim(0.5, n + 0.5) + theme_void()
  for (i in seq_along(logo_files)) {
    if (!is.na(logo_files[i]) && file.exists(logo_files[i])) {
      tryCatch({
        rast <- svg_to_raster(logo_files[i])
        if (!is.null(rast))
          p <- p + annotation_raster(rast,
                                     xmin = 0, xmax = 1,
                                     ymin = i - 0.4, ymax = i + 0.4,
                                     interpolate = TRUE)
      }, error = function(e) message("Could not read logo: ", logo_files[i]))
    }
  }
  p
}

create_known_motif_plot <- function(homer_dir, top_n = 15L, condition_name = "") {
  results    <- parse_homer_known_results(homer_dir)
  n_avail    <- min(top_n, nrow(results))
  top_motifs <- results[seq_len(n_avail), ]

  cat("\nHOMER Known Motif Results\n")
  cat("Total motifs tested:", nrow(results), "\n")
  cat("Motifs with P < 1e-10:", sum(results$P.value < 1e-10, na.rm = TRUE), "\n")
  cat("Motifs with P < 1e-5:",  sum(results$P.value < 1e-5,  na.rm = TRUE), "\n")
  cat("Motifs with P < 0.05:",  sum(results$P.value < 0.05,  na.rm = TRUE), "\n\n")

  logo_files          <- get_known_logo_files(homer_dir, seq_len(n_avail))
  logo_files_reversed <- rev(logo_files)

  top_motifs$motif_factor <- factor(top_motifs$motif_short,
                                    levels = rev(top_motifs$motif_short))
  top_motifs$y_pos <- seq(n_avail, 1, by = -1)

  bar_plot <- ggplot(top_motifs, aes(x = y_pos, y = neg_log_pval)) +
    geom_col(fill = "#3182bd", color = "#08519c", linewidth = 0.3) +
    geom_text(aes(y = 0, label = motif_short),
              hjust = 0, nudge_y = 0.3, size = 3.5,
              color = "white", fontface = "bold") +
    geom_text(aes(label = sprintf("%.1f", neg_log_pval)),
              hjust = -0.2, size = 2.5, color = "white") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(expand = expansion(add = c(0.5, 0.5))) +
    labs(x = NULL, y = "-log10(P-value)", title = condition_name) +
    theme_minimal() +
    theme(axis.text.y      = element_blank(),
          axis.ticks.y     = element_blank(),
          axis.line.x      = element_line(color = "black", linewidth = 0.5),
          axis.line.y      = element_line(color = "black", linewidth = 0.5),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.title         = element_text(face = "bold", hjust = 0.5, size = 14),
          plot.margin        = margin(5, 10, 5, 5))

  logo_col <- create_logo_column(logo_files_reversed)
  plot_grid(logo_col, bar_plot, ncol = 2, rel_widths = c(1.2, 2.5),
            align = "h", axis = "tb")
}

# Visualization ----
stopifnot(file.exists(file.path(homer_dir, "knownResults.txt")))

message("Generating gained anchor promoter known motif plot with logos...")
p <- create_known_motif_plot(homer_dir      = homer_dir,
                              top_n         = top_n,
                              condition_name = "Promoters at Gained Loop Anchors")

# Save outputs ----
ggsave(output_pdf, p, width = 12, height = 8)

results <- parse_homer_known_results(homer_dir)
cat("\nTop", top_n, "Enriched Known Motifs\n")
print(results[seq_len(min(top_n, nrow(results))),
              c("Motif.Name", "P.value", "q.value..Benjamini.",
                "X..of.Target.Sequences.with.Motif",
                "X..of.Background.Sequences.with.Motif")])

sessionInfo()
