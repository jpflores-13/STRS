# Supplementary Figure 6 — CUT&Tag binding scatter plots -----------------
# Author:      JP Flores
# Date:        2026-06-01
# Project:     STRS
# Description: 3-panel supplementary figure showing per-peak log-transformed
#              CUT&Tag signal for CTCF, RAD21, and YAP1 in control vs.
#              sorbitol conditions. Points above the diagonal indicate higher
#              signal in control (lost binding); below indicate gained binding.
# Input:       data/processed/cutntag/deseq2/diff_CTCF_counts.rds
#              data/processed/cutntag/deseq2/diff_RAD21_counts.rds
#              data/processed/cutntag/deseq2/diff_YAP1_counts.rds
# Output:      figures/FigureS6.pdf
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

diff_CTCF_rds  <- "data/processed/cutntag/deseq2/diff_CTCF_counts.rds"
diff_RAD21_rds <- "data/processed/cutntag/deseq2/diff_RAD21_counts.rds"
diff_YAP1_rds  <- "data/processed/cutntag/deseq2/diff_YAP1_counts.rds"
output_pdf     <- "figures/FigureS6.pdf"

page_width  <- 10   # inches — 3 panels side by side
page_height <- 3.5

ctcf_color  <- "#253494"
rad21_color <- "#41B6C4"
yap1_color  <- "#238B45"

highlight_color <- "#E91E8C"   # magenta, matching reference figure


# Libraries ---------------------------------------------------------------

library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(patchwork)


# Load data ---------------------------------------------------------------

## GRanges objects with per-replicate raw counts in mcols() and DESeq2 stats
diff_CTCF  <- readRDS(diff_CTCF_rds)
diff_RAD21 <- readRDS(diff_RAD21_rds)
diff_YAP1  <- readRDS(diff_YAP1_rds)


# Wrangle data ------------------------------------------------------------

#' Build a tidy data frame of log2-transformed mean counts per condition
#'
#' @param gr      GRanges object with per-replicate counts in mcols()
#' @param protein Character string matching the protein name in column names
#'                (e.g. "CTCF", "RAD21", "YAP1")
#'
#' @return data.frame with columns: log_control, log_sorbitol, direction
make_scatter_data <- function(gr, protein) {
  col_names <- names(mcols(gr))
  
  cont_cols <- col_names[grepl(protein, col_names) & grepl("cont", col_names)]
  sorb_cols <- col_names[grepl(protein, col_names) & grepl("sorbitol", col_names)]
  
  cont_mat <- as.matrix(as.data.frame(mcols(gr)[, cont_cols]))
  sorb_mat <- as.matrix(as.data.frame(mcols(gr)[, sorb_cols]))
  
  log_control  <- rowMeans(log2(cont_mat + 1))
  log_sorbitol <- rowMeans(log2(sorb_mat + 1))
  
  # Pull padj directly from mcols — already there from DESeq2
  padj <- mcols(gr)$padj
  
  data.frame(log_control, log_sorbitol, padj)
}


# Visualization -----------------------------------------------------------

#' Make a single binding scatter panel
#'
#' @param gr            GRanges object (output of readRDS on a diff counts RDS)
#' @param protein       Character; protein name matching mcols() column names
#' @param point_color   Hex color for highlighted off-diagonal points
#'
#' @return A ggplot object
make_binding_scatter <- function(gr, protein, point_color) {
  
  df <- make_scatter_data(gr, protein)
  axis_max <- ceiling(quantile(c(df$log_control, df$log_sorbitol), 0.999, na.rm = TRUE))
  
  df$significant <- !is.na(df$padj) & df$padj < 0.05
  
  ggplot(df, aes(x = log_sorbitol, y = log_control)) +
    # Smooth grey density raster underneath
    stat_density_2d(
      aes(fill = after_stat(density)),
      geom    = "raster",
      contour = FALSE,
      n       = 300,
      h       = c(1.5, 1.5)
    ) +
    scale_fill_gradientn(
      colors = c("white", "#d9d9d9", "#969696", "#525252"),
      guide  = "none"
    ) +
    # Non-significant points — same grey family, semi-transparent
    geom_point(
      data  = df |> dplyr::filter(!significant),
      color = "grey75",
      alpha = 0.3,
      size  = 0.3
    ) +
    # Significant points on top in protein color
    geom_point(
      data  = df |> dplyr::filter(significant),
      color = point_color,
      alpha = 0.7,
      size  = 0.4
    ) +
    geom_abline(slope = 1, intercept = 0, color = "grey60", linewidth = 0.5) +
    coord_fixed(xlim = c(0, axis_max), ylim = c(0, axis_max)) +
    scale_x_continuous(breaks = seq(0, axis_max, by = 2)) +
    scale_y_continuous(breaks = seq(0, axis_max, by = 2)) +
    labs(
      title = paste(protein, "binding"),
      x     = "Log concentration: sorbitol 1h",
      y     = "Log concentration: No stress"
    ) +
    theme_classic() +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 9, face = "bold"),
      axis.text    = element_text(size = 7),
      axis.title   = element_text(size = 8),
      axis.line    = element_line(linewidth = 0.3),
      axis.ticks   = element_line(linewidth = 0.3),
      aspect.ratio = 1,
      plot.margin  = margin(5.5, 5.5, 5.5, 5.5, "pt")
    )
}

p_ctcf  <- make_binding_scatter(diff_CTCF,  "CTCF",  point_color = ctcf_color)
p_rad21 <- make_binding_scatter(diff_RAD21, "RAD21", point_color = rad21_color)
p_yap1  <- make_binding_scatter(diff_YAP1,  "YAP1",  point_color = yap1_color)

figure_s <- p_ctcf + p_rad21 + p_yap1 +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold"))


# Save outputs ------------------------------------------------------------

pdf(output_pdf, width = page_width, height = page_height)
print(figure_s)
dev.off()


# Session info ------------------------------------------------------------

sessionInfo()