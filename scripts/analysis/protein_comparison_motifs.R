# ##############################################################################
# filename:    protein_comparison_motifs.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Multi-panel HOMER known-motif comparison for CTCF, RAD21,
#              H3K27ac, and YAP1 in control vs sorbitol conditions;
#              top 15 motifs per panel with SVG logos
# ##############################################################################

# Libraries ----
library(ggplot2)
library(cowplot)
library(rsvg)
library(png)
library(grid)
library(dplyr)

# Parameters ----
homer_base <- "data/processed/cutntag/homer_motifs"
output_pdf <- "plots/protein_comparison_motifs.pdf"
proteins   <- c("CTCF", "RAD21", "H3K27ac", "YAP1")
top_n      <- 15

# Helper functions ----
parse_homer_results <- function(homer_dir) {
  results_file <- file.path(homer_dir, "knownResults.txt")
  if (!file.exists(results_file)) stop("Cannot find ", results_file)
  results <- read.delim(results_file, header = TRUE, stringsAsFactors = FALSE)
  results$neg_log_pval <- -log10(results$P.value)
  results
}

get_logo_files <- function(homer_dir, n_motifs) {
  knownResults_dir <- file.path(homer_dir, "knownResults")
  logo_files <- list.files(knownResults_dir,
                           pattern    = "^known\\d+\\.logo\\.svg$",
                           full.names = TRUE)
  if (length(logo_files) == 0) {
    warning("No logo files found in ", knownResults_dir)
    return(rep(NA_character_, n_motifs))
  }
  logo_idx   <- as.numeric(gsub(".*known(\\d+)\\.logo\\.svg", "\\1", basename(logo_files)))
  logo_files <- logo_files[order(logo_idx)]
  logo_files[seq_len(min(n_motifs, length(logo_files)))]
}

svg_to_raster <- function(svg_file) {
  tmp <- tempfile(fileext = ".png")
  rsvg::rsvg_png(svg_file, tmp, width = 1200, height = 300)
  img <- png::readPNG(tmp)
  unlink(tmp)
  img
}

create_logo_column <- function(logo_files) {
  n <- length(logo_files)
  p <- ggplot() + xlim(0, 1) + ylim(0, n) + theme_void()
  for (i in seq_along(logo_files)) {
    if (!is.na(logo_files[i]) && file.exists(logo_files[i])) {
      tryCatch({
        rast <- svg_to_raster(logo_files[i])
        p <- p + annotation_raster(rast, xmin = 0, xmax = 1,
                                   ymin = i - 1, ymax = i, interpolate = TRUE)
      }, error = function(e) message("Could not read logo: ", logo_files[i]))
    }
  }
  p
}

create_motif_panel <- function(homer_dir, top_n = 15, condition_name = "") {
  results     <- parse_homer_results(homer_dir)
  n_available <- min(top_n, nrow(results))
  top_motifs  <- results[seq_len(n_available), ]

  logo_files <- get_logo_files(homer_dir, n_available)
  if (length(logo_files) > n_available) {
    logo_files <- logo_files[seq_len(n_available)]
  } else if (length(logo_files) < n_available) {
    logo_files <- c(logo_files, rep(NA_character_, n_available - length(logo_files)))
  }

  top_motifs$logo_file    <- logo_files
  top_motifs$motif_short  <- gsub("\\(.*", "", top_motifs$Motif.Name)
  top_motifs$motif_factor <- factor(top_motifs$motif_short,
                                    levels = rev(top_motifs$motif_short))

  bar_plot <- ggplot(top_motifs, aes(x = motif_factor, y = neg_log_pval)) +
    geom_col(fill = "#3182bd", color = "#08519c", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f", neg_log_pval)),
              hjust = -0.2, size = 2.5) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = NULL, y = "-log10(P-value)", title = condition_name) +
    theme_minimal() +
    theme(
      axis.text.y        = element_text(size = 9, hjust = 1, face = "bold"),
      axis.ticks.y       = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.margin        = margin(5, 10, 5, 5)
    )

  plot_grid(create_logo_column(rev(logo_files)), bar_plot,
            ncol = 2, rel_widths = c(1.5, 2.5), align = "h")
}

# Analysis ----
plot_protein_comparison <- function(base_dir, output_file, top_n = 15) {
  all_rows <- lapply(proteins, function(protein) {
    protein_lower <- tolower(protein)
    control_dir   <- file.path(base_dir, "top500", paste0(protein_lower, "_control"))
    sorbitol_dir  <- file.path(base_dir, "top500", paste0(protein_lower, "_sorbitol"))

    control_panel  <- if (dir.exists(control_dir))
      create_motif_panel(control_dir,  top_n, "Control")  else ggplot() + theme_void()
    sorbitol_panel <- if (dir.exists(sorbitol_dir))
      create_motif_panel(sorbitol_dir, top_n, "Sorbitol") else ggplot() + theme_void()

    protein_row <- plot_grid(control_panel, sorbitol_panel, ncol = 2)
    protein_label <- ggdraw() +
      draw_label(protein, fontface = "bold", size = 14, x = 0.05, hjust = 0)

    plot_grid(protein_label, protein_row, ncol = 1, rel_heights = c(0.05, 1))
  })
  names(all_rows) <- proteins

  title <- ggdraw() +
    draw_label(
      sprintf("Top %d Enriched Motifs: Control vs Sorbitol Treatment", top_n),
      fontface = "bold", size = 16
    )

  final_plot <- plot_grid(title, plot_grid(plotlist = all_rows, ncol = 1),
                          ncol = 1, rel_heights = c(0.03, 1))

  ggsave(output_file, final_plot, width = 16, height = 20)
  message("Saved: ", output_file)
}

# Save outputs ----
plot_protein_comparison(homer_base, output_pdf, top_n = top_n)

sessionInfo()
