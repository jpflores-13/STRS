# ##############################################################################
# filename:    eisa_scatter.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: EISA Part 4 — eisaR-style scatter plots of delta-exon vs
#              delta-intron across all timecourse timepoints; horizontal strip
#              and faceted versions; genes at gained loop anchors highlighted
# ##############################################################################

# Libraries ----
library(tidyverse)
library(ggplot2)
library(patchwork)

# Parameters ----
eisa_dir                <- "data/processed/rna/timecourse/output/EISA"
plots_dir               <- "plots"
eisa_results_rds        <- file.path(eisa_dir, "eisa_results_with_loops.rds")
output_horizontal_pdf   <- file.path(plots_dir, "eisa_scatter_horizontal_strip.pdf")
output_faceted_pdf      <- file.path(plots_dir, "eisa_scatter_faceted.pdf")
output_stats_csv        <- file.path(eisa_dir, "eisa_timecourse_summary_stats.csv")

# Data import ----
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
eisa_results <- readRDS(eisa_results_rds)
message(sprintf("Loaded %d genes", nrow(eisa_results)))

# Helper functions ----

reg_colors <- c(
  "Transcriptional"      = "#E41A1C",
  "Post-transcriptional" = "#377EB8",
  "Mixed"                = "#984EA3",
  "Not significant"      = "gray80"
)

plot_eisa_scatter <- function(eisa_df, timepoint, show_legend = TRUE) {
  plot_data <- data.frame(
    gene_id      = eisa_df$gene_id,
    symbol       = eisa_df$symbol,
    delta_intron = eisa_df[[paste0("delta_intron_", timepoint)]],
    delta_exon   = eisa_df[[paste0("delta_exon_",   timepoint)]],
    delta_posttx = eisa_df[[paste0("delta_posttx_", timepoint)]],
    loop_category = eisa_df$loop_category
  ) |>
    filter(!is.na(delta_intron) & !is.na(delta_exon)) |>
    filter(loop_category != "none") |>
    mutate(
      regulation_type = case_when(
        abs(delta_intron) > 1 & abs(delta_posttx) < 0.5 ~ "Transcriptional",
        abs(delta_intron) < 0.5 & abs(delta_posttx) > 1 ~ "Post-transcriptional",
        abs(delta_intron) > 1 & abs(delta_posttx) > 1   ~ "Mixed",
        TRUE ~ "Not significant"
      ),
      regulation_type = factor(regulation_type,
                               levels = c("Transcriptional", "Post-transcriptional",
                                          "Mixed", "Not significant")),
      at_gained = loop_category == "gained"
    )

  p <- ggplot(plot_data, aes(x = delta_intron, y = delta_exon)) +
    geom_abline(slope = 1, intercept = 0, color = "#E41A1C",
                linetype = "dashed", linewidth = 0.6) +
    geom_point(data = filter(plot_data, !at_gained),
               aes(color = regulation_type), alpha = 0.25, size = 1.2) +
    geom_point(data = filter(plot_data, at_gained),
               aes(fill = regulation_type), color = "black",
               shape = 21, size = 2, alpha = 0.7, stroke = 0.4) +
    scale_color_manual(values = reg_colors, name = "Regulation") +
    scale_fill_manual(values  = reg_colors, guide = "none") +
    coord_fixed(ratio = 1, xlim = c(-4, 4), ylim = c(-4, 4)) +
    labs(title    = timepoint,
         x = expression(paste(Delta, "intron (nascent RNA, log"[2], "FC)")),
         y = expression(paste(Delta, "exon (mature mRNA, log"[2],   "FC)"))) +
    theme_minimal(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.text  = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 9),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position  = if (show_legend) "right" else "none",
      legend.title     = element_text(size = 8, face = "bold"),
      legend.text      = element_text(size = 7),
      legend.key.size  = unit(0.4, "cm")
    )

  list(plot = p, data = plot_data)
}

# Visualization ----
timepoints <- c("1h", "3h", "6h", "9h", "12h", "24h")

## Horizontal strip ----
plots_horizontal <- list()

for (i in seq_along(timepoints)) {
  tp     <- timepoints[i]
  result <- plot_eisa_scatter(eisa_results, tp, show_legend = (i == length(timepoints)))

  plots_horizontal[[tp]] <- result$plot +
    theme(
      axis.title.y = if (i == 1) element_text(size = 9) else element_blank(),
      axis.text.y  = if (i == 1) element_text(size = 8) else element_blank()
    )
}

grid_horizontal_final <- wrap_plots(plots_horizontal, ncol = 6, nrow = 1) +
  plot_annotation(
    title    = expression(paste("EISA Timecourse: ", Delta, "exon vs ", Delta, "intron")),
    subtitle = "Red diagonal = purely transcriptional | Genes at gained loop anchors (black outline)",
    theme    = theme(plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
                     plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40",
                                                  lineheight = 1.2))
  )

ggsave(output_horizontal_pdf, grid_horizontal_final, width = 18, height = 4)
message(sprintf("Saved: %s", output_horizontal_pdf))

## Faceted version ----
prepare_long_data <- function(eisa_df) {
  tps <- c("1h", "3h", "6h", "9h", "12h", "24h")

  lapply(tps, function(tp) {
    data.frame(
      gene_id      = eisa_df$gene_id,
      symbol       = eisa_df$symbol,
      timepoint    = tp,
      delta_intron = eisa_df[[paste0("delta_intron_", tp)]],
      delta_exon   = eisa_df[[paste0("delta_exon_",   tp)]],
      delta_posttx = eisa_df[[paste0("delta_posttx_", tp)]],
      loop_category = eisa_df$loop_category
    )
  }) |>
    bind_rows() |>
    filter(!is.na(delta_intron) & !is.na(delta_exon), loop_category != "none") |>
    mutate(
      timepoint = factor(timepoint, levels = tps),
      regulation_type = case_when(
        abs(delta_intron) > 1 & abs(delta_posttx) < 0.5 ~ "Transcriptional",
        abs(delta_intron) < 0.5 & abs(delta_posttx) > 1 ~ "Post-transcriptional",
        abs(delta_intron) > 1 & abs(delta_posttx) > 1   ~ "Mixed",
        TRUE ~ "Not significant"
      ),
      at_gained = loop_category == "gained"
    )
}

facet_data <- prepare_long_data(eisa_results)

p_facet <- ggplot(facet_data, aes(x = delta_intron, y = delta_exon)) +
  geom_abline(slope = 1, intercept = 0, color = "#E41A1C",
              linetype = "dashed", linewidth = 0.5) +
  geom_point(data = filter(facet_data, !at_gained),
             aes(color = regulation_type), alpha = 0.3, size = 0.8) +
  geom_point(data = filter(facet_data,  at_gained),
             aes(fill = regulation_type), color = "black",
             shape = 21, size = 1.5, alpha = 0.7, stroke = 0.3) +
  facet_wrap(~timepoint, ncol = 3, nrow = 2) +
  scale_color_manual(values = reg_colors, name = "Regulation") +
  scale_fill_manual(values  = reg_colors, guide = "none") +
  coord_fixed(ratio = 1, xlim = c(-4, 4), ylim = c(-4, 4)) +
  labs(
    title    = "EISA Timecourse: Transcriptional vs Post-transcriptional Regulation",
    subtitle = "Genes at chromatin loop anchors | Red diagonal = purely transcriptional",
    x = expression(paste(Delta, "intron (nascent RNA, log"[2], "FC)")),
    y = expression(paste(Delta, "exon (mature mRNA, log"[2],   "FC)"))
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.background  = element_rect(fill = "white", color = NA),
    panel.grid.major  = element_line(color = "gray90", linewidth = 0.25),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.text        = element_text(size = 10, face = "bold"),
    strip.background  = element_rect(fill = "gray95", color = "black", linewidth = 0.5),
    axis.text  = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40", lineheight = 1.1),
    legend.position = "bottom"
  )

ggsave(output_faceted_pdf, p_facet, width = 11, height = 7.5)
message(sprintf("Saved: %s", output_faceted_pdf))

# Save outputs ----
stats_summary <- facet_data |>
  filter(at_gained) |>
  group_by(timepoint) |>
  summarise(
    n_genes              = n(),
    cor_pearson          = cor(delta_intron, delta_exon, method = "pearson", use = "complete.obs"),
    pct_transcriptional  = sum(regulation_type == "Transcriptional",      na.rm = TRUE) / n() * 100,
    pct_posttx           = sum(regulation_type == "Post-transcriptional", na.rm = TRUE) / n() * 100,
    pct_mixed            = sum(regulation_type == "Mixed",                na.rm = TRUE) / n() * 100,
    .groups = "drop"
  )

print(stats_summary)
write.csv(stats_summary, output_stats_csv, row.names = FALSE)

sessionInfo()
