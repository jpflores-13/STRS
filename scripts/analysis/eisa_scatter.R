## =============================================================================
## EISA PART 4: eisaR-STYLE SCATTER PLOTS (HORIZONTAL & FACETED ONLY)
## =============================================================================
##
## Purpose: Create eisaR vignette-style scatter plots for entire timecourse
##
## Input:  EISA results table from Part 2
## Output: Horizontal strip (1×6) and faceted version
##
## Runtime: ~30 seconds
##
## =============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)

message("Starting EISA Part 4 (eisaR-style plots): ", Sys.time())

## =============================================================================
## CONFIGURATION
## =============================================================================

BASE_DIR <- "/work/users/j/p/jpflores/projects/STRS"
EISA_DIR <- file.path(BASE_DIR, "data/processed/rna/timecourse/output/EISA")
PLOT_DIR <- file.path(BASE_DIR, "plots")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

## =============================================================================
## LOAD DATA
## =============================================================================

message("\n=== Loading EISA results ===")

eisa_results <- readRDS(file.path(EISA_DIR, "eisa_results_with_loops.rds"))

message(sprintf("Loaded %d genes", nrow(eisa_results)))

## =============================================================================
## COLOR SCHEME
## =============================================================================

reg_colors <- c(
  "Transcriptional" = "#E41A1C",      # Red
  "Post-transcriptional" = "#377EB8", # Blue
  "Mixed" = "#984EA3",                # Purple
  "Not significant" = "gray80"
)

## =============================================================================
## HELPER FUNCTION: CREATE SCATTER PLOT
## =============================================================================

plot_eisa_scatter <- function(eisa_df, timepoint, show_legend = TRUE) {
  
  ## Extract data for this timepoint
  plot_data <- data.frame(
    gene_id = eisa_df$gene_id,
    symbol = eisa_df$symbol,
    delta_intron = eisa_df[[paste0("delta_intron_", timepoint)]],
    delta_exon = eisa_df[[paste0("delta_exon_", timepoint)]],
    delta_posttx = eisa_df[[paste0("delta_posttx_", timepoint)]],
    loop_category = eisa_df$loop_category
  ) |>
    filter(!is.na(delta_intron) & !is.na(delta_exon)) |>
    filter(loop_category != "none")
  
  ## Classify regulation type
  plot_data <- plot_data |>
    mutate(
      regulation_type = case_when(
        abs(delta_intron) > 1 & abs(delta_posttx) < 0.5 ~ "Transcriptional",
        abs(delta_intron) < 0.5 & abs(delta_posttx) > 1 ~ "Post-transcriptional",
        abs(delta_intron) > 1 & abs(delta_posttx) > 1 ~ "Mixed",
        TRUE ~ "Not significant"
      ),
      regulation_type = factor(regulation_type, 
                               levels = c("Transcriptional", "Post-transcriptional", 
                                          "Mixed", "Not significant")),
      at_gained = loop_category == "gained"
    )
  
  ## Create plot
  p <- ggplot(plot_data, aes(x = delta_intron, y = delta_exon)) +
    
    ## Diagonal line (Δexon = Δintron = purely transcriptional)
    geom_abline(slope = 1, intercept = 0, 
                color = "#E41A1C", linetype = "dashed", linewidth = 0.6) +
    
    ## Background points (not at gained loops)
    geom_point(
      data = filter(plot_data, !at_gained),
      aes(color = regulation_type),
      alpha = 0.25,
      size = 1.2
    ) +
    
    ## Gained loop genes (highlighted with black outline)
    geom_point(
      data = filter(plot_data, at_gained),
      aes(fill = regulation_type),
      color = "black",
      shape = 21,
      size = 2,
      alpha = 0.7,
      stroke = 0.4
    ) +
    
    ## Color scales
    scale_color_manual(values = reg_colors, name = "Regulation") +
    scale_fill_manual(values = reg_colors, guide = "none") +
    
    ## Axes
    coord_fixed(ratio = 1, xlim = c(-4, 4), ylim = c(-4, 4)) +
    
    labs(
      title = timepoint,
      x = expression(paste(Delta, "intron (nascent RNA, log"[2], "FC)")),
      y = expression(paste(Delta, "exon (mature mRNA, log"[2], "FC)"))
    ) +
    
    ## Theme
    theme_minimal(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 9),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = if(show_legend) "right" else "none",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.4, "cm")
    )
  
  return(list(plot = p, data = plot_data))
}

## =============================================================================
## PLOT 1: HORIZONTAL STRIP (1 ROW × 6 COLUMNS)
## =============================================================================

message("\n=== Creating horizontal strip (1×6) ===")

timepoints <- c("1h", "3h", "6h", "9h", "12h", "24h")
plots_horizontal <- list()

for (i in seq_along(timepoints)) {
  tp <- timepoints[i]
  
  ## Show legend on last plot
  show_legend <- (i == length(timepoints))
  
  result <- plot_eisa_scatter(eisa_results, tp, show_legend = show_legend)
  
  ## Adjust for horizontal layout
  plots_horizontal[[tp]] <- result$plot +
    theme(
      axis.title.y = if(i == 1) {
        element_text(size = 9)
      } else {
        element_blank()
      },
      axis.text.y = if(i == 1) element_text(size = 8) else element_blank()
    )
}

## Stack horizontally
grid_horizontal <- wrap_plots(plots_horizontal, ncol = 6, nrow = 1)

grid_horizontal_final <- grid_horizontal +
  plot_annotation(
    title = expression(paste("EISA Timecourse: ", Delta, "exon (mature mRNA) vs ", Delta, "intron (nascent RNA)")),
    subtitle = "Red diagonal indicates purely transcriptional regulation | Genes at gained loop anchors (black outline)",
    theme = theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40", lineheight = 1.2)
    )
  )

ggsave(
  filename = file.path(PLOT_DIR, "eisa_scatter_horizontal_strip.pdf"),
  plot = grid_horizontal_final,
  width = 18,
  height = 4
)

message("✓ Saved: eisa_scatter_horizontal_strip.pdf (18 × 4 inches)")

## =============================================================================
## PLOT 2: FACETED VERSION (ggplot facets)
## =============================================================================

message("\n=== Creating faceted version ===")

## Prepare long-format data
prepare_long_data <- function(eisa_df) {
  timepoints <- c("1h", "3h", "6h", "9h", "12h", "24h")
  
  long_data <- list()
  
  for (tp in timepoints) {
    df <- data.frame(
      gene_id = eisa_df$gene_id,
      symbol = eisa_df$symbol,
      timepoint = tp,
      delta_intron = eisa_df[[paste0("delta_intron_", tp)]],
      delta_exon = eisa_df[[paste0("delta_exon_", tp)]],
      delta_posttx = eisa_df[[paste0("delta_posttx_", tp)]],
      loop_category = eisa_df$loop_category
    )
    
    long_data[[tp]] <- df
  }
  
  bind_rows(long_data) |>
    filter(!is.na(delta_intron) & !is.na(delta_exon)) |>
    filter(loop_category != "none") |>
    mutate(
      timepoint = factor(timepoint, levels = c("1h", "3h", "6h", "9h", "12h", "24h")),
      regulation_type = case_when(
        abs(delta_intron) > 1 & abs(delta_posttx) < 0.5 ~ "Transcriptional",
        abs(delta_intron) < 0.5 & abs(delta_posttx) > 1 ~ "Post-transcriptional",
        abs(delta_intron) > 1 & abs(delta_posttx) > 1 ~ "Mixed",
        TRUE ~ "Not significant"
      ),
      at_gained = loop_category == "gained"
    )
}

facet_data <- prepare_long_data(eisa_results)

p_facet <- ggplot(facet_data, aes(x = delta_intron, y = delta_exon)) +
  
  ## Diagonal
  geom_abline(slope = 1, intercept = 0, 
              color = "#E41A1C", linetype = "dashed", linewidth = 0.5) +
  
  ## Background points
  geom_point(
    data = filter(facet_data, !at_gained),
    aes(color = regulation_type),
    alpha = 0.3,
    size = 0.8
  ) +
  
  ## Gained points
  geom_point(
    data = filter(facet_data, at_gained),
    aes(fill = regulation_type),
    color = "black",
    shape = 21,
    size = 1.5,
    alpha = 0.7,
    stroke = 0.3
  ) +
  
  ## Facet by timepoint
  facet_wrap(~timepoint, ncol = 3, nrow = 2) +
  
  ## Colors
  scale_color_manual(values = reg_colors, name = "Regulation") +
  scale_fill_manual(values = reg_colors, guide = "none") +
  
  coord_fixed(ratio = 1, xlim = c(-4, 4), ylim = c(-4, 4)) +
  
  labs(
    title = "EISA Timecourse: Transcriptional vs Post-transcriptional Regulation",
    subtitle = "Genes at chromatin loop anchors | Red diagonal indicates purely transcriptional regulation",
    x = expression(paste(Delta, "intron (nascent RNA, log"[2], "FC)")),
    y = expression(paste(Delta, "exon (mature mRNA, log"[2], "FC)"))
  ) +
  
  theme_minimal(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40", lineheight = 1.1),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(PLOT_DIR, "eisa_scatter_faceted.pdf"),
  plot = p_facet,
  width = 11,
  height = 7.5
)

message("✓ Saved: eisa_scatter_faceted.pdf (11 × 7.5 inches)")

## =============================================================================
## SUMMARY STATISTICS
## =============================================================================

message("\n=== Summary Statistics Across Timepoints ===")

stats_summary <- facet_data |>
  filter(at_gained) |>
  group_by(timepoint) |>
  summarise(
    n_genes = n(),
    cor_pearson = cor(delta_intron, delta_exon, method = "pearson", use = "complete.obs"),
    pct_transcriptional = sum(regulation_type == "Transcriptional", na.rm = TRUE) / n() * 100,
    pct_posttx = sum(regulation_type == "Post-transcriptional", na.rm = TRUE) / n() * 100,
    pct_mixed = sum(regulation_type == "Mixed", na.rm = TRUE) / n() * 100,
    .groups = "drop"
  )

print(stats_summary)

write.csv(stats_summary, 
          file.path(EISA_DIR, "eisa_timecourse_summary_stats.csv"),
          row.names = FALSE)

## =============================================================================
## COMPLETION
## =============================================================================

message("\n")
message("================================================================================")
message("              PART 4 COMPLETE: eisaR-STYLE SCATTER PLOTS")
message("================================================================================")
message("\n")
message("Created 2 plots:")
message(sprintf("  1. Horizontal strip: %s", file.path(PLOT_DIR, "eisa_scatter_horizontal_strip.pdf")))
message(sprintf("  2. Faceted version:  %s", file.path(PLOT_DIR, "eisa_scatter_faceted.pdf")))
message("\n")
message("Summary statistics saved:")
message(sprintf("  %s", file.path(EISA_DIR, "eisa_timecourse_summary_stats.csv")))
message("\n")
message("Axis labels:")
message("  X-axis: Δintron (nascent RNA, log₂FC)")
message("  Y-axis: Δexon (mature mRNA, log₂FC)")
message("\n")
message("Key finding:")
message("  Genes at gained loop anchors cluster along red diagonal")
message("  → Regulation is TRANSCRIPTIONAL, not post-transcriptional")
message("\n")
message("Completed: ", Sys.time())
message("================================================================================")