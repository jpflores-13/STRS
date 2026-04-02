# ##############################################################################
# filename:    eisa_lineplots.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: EISA Part 3 — three-panel line plot of median delta-intron,
#              delta-exon, and delta-posttx over time at gained, lost, and
#              static loop anchor genes
# ##############################################################################

# Libraries ----
library(tidyverse)
library(ggplot2)
library(cowplot)

# Parameters ----
eisa_dir        <- "data/processed/rna/timecourse/output/EISA"
plots_dir       <- "plots"
eisa_results_rds <- file.path(eisa_dir, "eisa_results_with_loops.rds")
stats_rds        <- file.path(eisa_dir, "statistics_summary.rds")
output_pdf       <- file.path(plots_dir, "eisa_lineplots.pdf")

# Data import ----
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(eisa_results_rds)) {
  stop("EISA results not found! Run eisa_analysis.R first.")
}

eisa_results <- readRDS(eisa_results_rds)

message(sprintf("Loaded %d genes", nrow(eisa_results)))
message(sprintf("  Gained: %d", sum(eisa_results$loop_category == "gained")))
message(sprintf("  Lost: %d",   sum(eisa_results$loop_category == "lost")))
message(sprintf("  Static: %d", sum(eisa_results$loop_category == "static")))

# Analysis ----

## Prepare plot data ----
loop_colors <- c("static" = "#999999", "gained" = "#F8766D", "lost" = "#619CFF")

prepare_plot_data <- function(eisa_df, metric_prefix) {
  eisa_df |>
    filter(loop_category != "none") |>
    select(gene_id, symbol, loop_category, starts_with(metric_prefix)) |>
    pivot_longer(
      cols      = starts_with(metric_prefix),
      names_to  = "timepoint",
      values_to = "value"
    ) |>
    mutate(
      timepoint     = str_extract(timepoint, "[0-9]+h") |> str_remove("h"),
      timepoint     = factor(timepoint, levels = c("1", "3", "6", "9", "12", "24")),
      loop_category = factor(loop_category, levels = c("static", "gained", "lost"))
    )
}

data_intron  <- prepare_plot_data(eisa_results, "delta_intron_")
data_exon    <- prepare_plot_data(eisa_results, "delta_exon_")
data_posttx  <- prepare_plot_data(eisa_results, "delta_posttx_")

# Visualization ----

make_line_plot <- function(data, title, subtitle, show_ylabel = TRUE) {
  data_summary <- data |>
    group_by(timepoint, loop_category) |>
    summarise(
      mean   = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      se     = sd(value, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    ) |>
    mutate(timepoint_num = as.numeric(as.character(timepoint)))

  ggplot(data_summary, aes(x = timepoint_num, y = median, color = loop_category)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.4) +
    geom_ribbon(aes(ymin = median - se, ymax = median + se, fill = loop_category),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2, alpha = 0.9) +
    geom_point(size = 3.5, alpha = 0.9) +
    scale_color_manual(values = loop_colors) +
    scale_fill_manual(values  = loop_colors, guide = "none") +
    scale_x_continuous(breaks = c(1, 3, 6, 9, 12, 24),
                       labels = c("1", "3", "6", "9", "12", "24")) +
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    labs(title    = title,
         subtitle = subtitle,
         x        = "Hours after hyperosmotic stress",
         y        = if (show_ylabel) "log2FoldChange (treated/untreated)" else "",
         color    = "Loop category") +
    theme_minimal() +
    theme(
      panel.background   = element_rect(fill = "white", color = NA),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.major.x = element_line(color = "gray95", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      axis.line  = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.text  = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 10, face = "bold"),
      plot.title    = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = if (show_ylabel) "right" else "none",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text  = element_text(size = 8)
    )
}

p1 <- make_line_plot(data_intron,  "Transcriptional Regulation",
                     "Δintron: nascent RNA changes",          show_ylabel = TRUE)
p2 <- make_line_plot(data_exon,    "Total Regulation",
                     "Δexon: mature mRNA changes",            show_ylabel = FALSE)
p3 <- make_line_plot(data_posttx,  "Post-transcriptional Regulation",
                     "Δexon - Δintron: mRNA stability changes", show_ylabel = FALSE)

# Save outputs ----
fig_combined <- plot_grid(p1, p2, p3, ncol = 3, align = "h",
                          rel_widths = c(1.2, 1, 1))

ggsave(output_pdf, fig_combined, width = 15, height = 5, units = "in")
message(sprintf("Figure saved: %s", output_pdf))

## Print figure statistics ----
stats <- readRDS(stats_rds)

cat(sprintf("\nGenes analyzed: %d total, %d at loop anchors\n",
            nrow(eisa_results),
            sum(eisa_results$loop_category != "none")))
cat(sprintf("Correlation (12h): R = %.3f, p = %.2e\n",
            stats$correlation$estimate, stats$correlation$p.value))
cat(sprintf("Gained loop genes median Δintron (12h): %.2f\n",
            median(eisa_results$delta_intron_12h[eisa_results$loop_category == "gained"],
                   na.rm = TRUE)))

sessionInfo()
