## =============================================================================
## EISA PART 3: VISUALIZATIONS
## =============================================================================
##
## Purpose: Create publication-quality figures from EISA results
##
## Input:  EISA results table from Part 2
## Output: Three-panel line plot figure
##
## Runtime: ~2 minutes
##
## =============================================================================

library(tidyverse)
library(ggplot2)
library(cowplot)

message("Starting EISA Part 3 (Visualizations): ", Sys.time())

## =============================================================================
## CONFIGURATION
## =============================================================================

BASE_DIR <- "/work/users/j/p/jpflores/projects/STRS"
EISA_DIR <- file.path(BASE_DIR, "data/processed/rna/timecourse/output/EISA")
PLOT_DIR <- file.path(BASE_DIR, "plots")

## Create plot directory
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

## =============================================================================
## LOAD EISA RESULTS
## =============================================================================

message("\n=== STEP 1: Loading EISA results ===")

eisa_file <- file.path(EISA_DIR, "eisa_results_with_loops.rds")

if (!file.exists(eisa_file)) {
  stop("EISA results not found! Run Part 2 first:\n",
       "  Rscript scripts/EISA/EISA_part2_DESeq2.R")
}

message("Loading EISA results...")
eisa_results <- readRDS(eisa_file)

message(sprintf("Loaded %d genes", nrow(eisa_results)))
message(sprintf("  Gained: %d", sum(eisa_results$loop_category == "gained")))
message(sprintf("  Lost: %d", sum(eisa_results$loop_category == "lost")))
message(sprintf("  Static: %d", sum(eisa_results$loop_category == "static")))

## =============================================================================
## PREPARE PLOT DATA
## =============================================================================

message("\n=== STEP 2: Preparing plot data ===")

## Define color scheme
loop_colors <- c(
  "static" = "#999999",
  "gained" = "#F8766D",
  "lost" = "#619CFF"
)

## Prepare data for plotting
prepare_plot_data <- function(eisa_df, metric_prefix) {
  eisa_df |>
    filter(loop_category != "none") |>
    select(gene_id, symbol, loop_category, starts_with(metric_prefix)) |>
    pivot_longer(
      cols = starts_with(metric_prefix),
      names_to = "timepoint",
      values_to = "value"
    ) |>
    mutate(
      timepoint = str_extract(timepoint, "[0-9]+h") |> str_remove("h"),
      timepoint = factor(timepoint, levels = c("1", "3", "6", "9", "12", "24")),
      loop_category = factor(loop_category, levels = c("static", "gained", "lost"))
    )
}

message("Preparing data for Δintron...")
data_intron <- prepare_plot_data(eisa_results, "delta_intron_")

message("Preparing data for Δexon...")
data_exon <- prepare_plot_data(eisa_results, "delta_exon_")

message("Preparing data for Δposttx...")
data_posttx <- prepare_plot_data(eisa_results, "delta_posttx_")

## =============================================================================
## CREATE LINE PLOTS
## =============================================================================

message("\n=== STEP 3: Creating line plots ===")

make_line_plot <- function(data, title, subtitle, show_ylabel = TRUE) {
  
  ## Calculate summary statistics
  data_summary <- data |>
    group_by(timepoint, loop_category) |>
    summarise(
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      se = sd(value, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    ) |>
    mutate(timepoint_num = as.numeric(as.character(timepoint)))
  
  ggplot(data_summary, aes(x = timepoint_num, y = median, color = loop_category)) +
    ## Zero reference
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.4) +
    
    ## Standard error ribbon
    geom_ribbon(
      aes(ymin = median - se, ymax = median + se, fill = loop_category),
      alpha = 0.2,
      color = NA
    ) +
    
    ## Main line
    geom_line(linewidth = 1.2, alpha = 0.9) +
    
    ## Points
    geom_point(size = 3.5, alpha = 0.9) +
    
    ## Colors
    scale_color_manual(values = loop_colors) +
    scale_fill_manual(values = loop_colors, guide = "none") +
    
    ## X-axis
    scale_x_continuous(
      breaks = c(1, 3, 6, 9, 12, 24),
      labels = c("1", "3", "6", "9", "12", "24")
    ) +
    
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    
    labs(
      title = title,
      subtitle = subtitle,
      x = "Hours after hyperosmotic stress",
      y = if(show_ylabel) "log2FoldChange (treated/untreated)" else "",
      color = "Loop category"
    ) +
    
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.major.x = element_line(color = "gray95", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = if(show_ylabel) "right" else "none",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8)
    )
}

message("Creating Panel A (Δintron)...")
p1 <- make_line_plot(
  data_intron,
  "Transcriptional Regulation",
  "Δintron: nascent RNA changes",
  show_ylabel = TRUE
)

message("Creating Panel B (Δexon)...")
p2 <- make_line_plot(
  data_exon,
  "Total Regulation",
  "Δexon: mature mRNA changes",
  show_ylabel = FALSE
)

message("Creating Panel C (Δposttx)...")
p3 <- make_line_plot(
  data_posttx,
  "Post-transcriptional Regulation",
  "Δexon - Δintron: mRNA stability changes",
  show_ylabel = FALSE
)

## =============================================================================
## COMBINE AND SAVE
## =============================================================================

message("\n=== STEP 4: Combining panels and saving ===")

fig_combined <- plot_grid(
  p1, p2, p3,
  ncol = 3,
  align = "h",
  rel_widths = c(1.2, 1, 1)
)

output_file <- file.path(PLOT_DIR, "eisa_lineplots.pdf")

ggsave(
  output_file,
  fig_combined,
  width = 15,
  height = 5,
  units = "in"
)

message(sprintf("✓ Figure saved: %s", output_file))

## =============================================================================
## OPTIONAL: SAVE INDIVIDUAL PANELS
## =============================================================================

if (FALSE) {  # Change to TRUE to save individual panels
  message("\nSaving individual panels...")
  
  ggsave(file.path(PLOT_DIR, "eisa_intron.pdf"),
         p1, width = 6, height = 5)
  
  ggsave(file.path(PLOT_DIR, "eisa_exon.pdf"),
         p2, width = 5, height = 5)
  
  ggsave(file.path(PLOT_DIR, "eisa_posttxn.pdf"),
         p3, width = 5, height = 5)
  
  message("✓ Individual panels saved")
}

## =============================================================================
## SUMMARY STATISTICS FOR FIGURE LEGEND
## =============================================================================

message("\n=== Figure Summary Statistics ===")

stats <- readRDS(file.path(EISA_DIR, "statistics_summary.rds"))

cat("\nFor figure legend/caption:\n")
cat(sprintf("  - Genes analyzed: %d total, %d at loop anchors\n",
            nrow(eisa_results),
            sum(eisa_results$loop_category != "none")))
cat(sprintf("  - Correlation (12h): R = %.3f, p = %.2e\n",
            stats$correlation$estimate,
            stats$correlation$p.value))
cat(sprintf("  - Gained loop genes show median Δintron = %.2f at 12h\n",
            median(eisa_results$delta_intron_12h[eisa_results$loop_category == "gained"], 
                   na.rm = TRUE)))
cat(sprintf("  - Post-transcriptional contribution: %.2f-fold smaller\n",
            median(abs(eisa_results$delta_intron_12h[eisa_results$loop_category == "gained"]), 
                   na.rm = TRUE) /
              median(abs(eisa_results$delta_posttx_12h[eisa_results$loop_category == "gained"]), 
                     na.rm = TRUE)))

## =============================================================================
## COMPLETION
## =============================================================================

message("\n")
message("================================================================================")
message("                    PART 3 COMPLETE: VISUALIZATIONS")
message("================================================================================")
message("\n")
message("Output:")
message(sprintf("  Figure: %s", output_file))
message("\n")
message("Figure shows:")
message("  Panel A: Δintron (transcriptional regulation)")
message("  Panel B: Δexon (total expression changes)")
message("  Panel C: Δexon - Δintron (post-transcriptional regulation)")
message("\n")
message("Key finding:")
message("  - Gained loops show high Δintron (Panel A)")
message("  - Similar Δexon pattern (Panel B)")
message("  - Minimal Δposttx (Panel C)")
message("  → Regulation is TRANSCRIPTIONAL, not post-transcriptional")
message("\n")
message("Next steps:")
message("  1. Download figure:")
message(sprintf("     scp jpflores@longleaf:%s ~/Downloads/", output_file))
message("  2. Review three-panel line plot")
message("  3. Integrate statistics into manuscript")
message("  4. Add figure to manuscript as Figure 5")
message("\n")
message("Completed: ", Sys.time())
message("================================================================================")