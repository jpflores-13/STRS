# ##############################################################################
# filename:    rnaseqTimecourse_loopAnchors_boxplots.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Three-panel boxplots of log2FC over the RNA-seq timecourse for
#              genes at gained, lost, and static loop anchors; includes
#              consecutive-timepoint Wilcoxon tests and baseline expression
#              comparison; also outputs absolute expression timecourse panels
# ##############################################################################

# Libraries ----
library(InteractionSet)
library(mariner)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(RColorBrewer)
library(cowplot)
library(rstatix)
library(stringr)

# Parameters ----
dds_rds        <- "data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds"
diff_loops_rds <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
output_main    <- "plots/rnaseqTimecourse_loopAnchors.pdf"
output_square  <- "plots/rnaseqTimecourse_loopAnchors_square.pdf"
output_baseline <- "plots/rnaseqTimecourse_baseline_expression.pdf"
output_abexpr  <- "plots/rnaseqTimecourse_absoluteExpression_by_loopType.pdf"

loop_colors <- c("static" = "#999999", "gained" = "#F8766D", "lost" = "#619CFF")
all_timepoints <- c("0h", "1h", "3h", "6h", "9h", "12h", "24h")

# Data import ----
dge <- readRDS(dds_rds)

dge_gr      <- results(dge, format = "GRanges")
mcols(dge_gr) <- rowData(dge)
dge_gr_prom <- promoters(dge_gr) |>
  keepStandardChromosomes(pruning.mode = "coarse")
seqlevelsStyle(dge_gr_prom) <- "UCSC"

loops <- readRDS(diff_loops_rds) |>
  interactions() |>
  as.data.frame() |>
  as_ginteractions()

# Analysis ----
gainedLoops <- loops |> subset(padj < 0.1 & log2FoldChange > 0)
lostLoops   <- loops |> subset(padj < 0.1 & log2FoldChange < 0)
staticLoops <- loops |> subset(padj > 0.1)

make_hits_df <- function(loops_gi, label) {
  df       <- as.data.frame(subsetByOverlaps(dge_gr_prom, loops_gi))
  df$type  <- label
  df
}

gained_df <- make_hits_df(gainedLoops, "gained")
lost_df   <- make_hits_df(lostLoops,   "lost")
static_df <- make_hits_df(staticLoops, "static")

combined <- bind_rows(gained_df, lost_df, static_df) |>
  dplyr::select(Time_1h_vs_0h, Time_3h_vs_0h, Time_6h_vs_0h,
                Time_9h_vs_0h, Time_12h_vs_0h, Time_24h_vs_0h,
                type, symbol, gene_id) |>
  pivot_longer(cols = starts_with("Time"),
               values_to = "log2FoldChange", names_to = "timepoint") |>
  mutate(timepoint = case_when(
    timepoint == "Time_1h_vs_0h"  ~ "1",
    timepoint == "Time_3h_vs_0h"  ~ "3",
    timepoint == "Time_6h_vs_0h"  ~ "6",
    timepoint == "Time_9h_vs_0h"  ~ "9",
    timepoint == "Time_12h_vs_0h" ~ "12",
    timepoint == "Time_24h_vs_0h" ~ "24"
  ),
  timepoint = factor(timepoint, levels = c("1", "3", "6", "9", "12", "24")),
  type      = factor(type, levels = c("static", "gained", "lost")))

loop_anchor_genes <- bind_rows(gained_df, lost_df, static_df) |>
  dplyr::select(gene_id, type) |>
  dplyr::distinct()

## Baseline expression analysis ----
normalized_counts     <- counts(dge, normalized = TRUE)
baseline_idx          <- which(dge$Time == "0h")
mean_baseline         <- rowMeans(normalized_counts[, baseline_idx, drop = FALSE])

baseline_df <- data.frame(gene_id = names(mean_baseline),
                           baseline_expression = mean_baseline,
                           log2_baseline       = log2(mean_baseline + 1),
                           stringsAsFactors    = FALSE)

baseline_with_loops <- baseline_df |>
  inner_join(loop_anchor_genes, by = "gene_id") |>
  mutate(type = factor(type, levels = c("static", "gained", "lost")))

baseline_summary <- baseline_with_loops |>
  group_by(type) |>
  summarise(n_genes              = dplyr::n(),
            mean_baseline        = mean(baseline_expression),
            median_baseline      = median(baseline_expression),
            mean_log2_baseline   = mean(log2_baseline),
            median_log2_baseline = median(log2_baseline),
            q25_log2             = quantile(log2_baseline, 0.25),
            q75_log2             = quantile(log2_baseline, 0.75),
            .groups              = "drop")
print(baseline_summary)

kw_test        <- kruskal_test(baseline_with_loops, log2_baseline ~ type)
pairwise_tests <- wilcox_test(baseline_with_loops, log2_baseline ~ type,
                               p.adjust.method = "BH") |>
  add_significance("p.adj")
print(kw_test)
print(pairwise_tests)

## Absolute expression timecourse ----
timecourse_expression_df <- lapply(all_timepoints, function(tp) {
  idx         <- which(dge$Time == tp)
  mean_counts <- rowMeans(normalized_counts[, idx, drop = FALSE])
  data.frame(gene_id         = names(mean_counts),
             expression      = mean_counts,
             log2_expression = log2(mean_counts + 1),
             timepoint       = tp, stringsAsFactors = FALSE)
}) |> bind_rows()

timecourse_with_loops <- timecourse_expression_df |>
  inner_join(loop_anchor_genes, by = "gene_id") |>
  mutate(type      = factor(type, levels = c("static", "gained", "lost")),
         timepoint = factor(timepoint, levels = all_timepoints))

## Consecutive-timepoint Wilcoxon tests ----
prepare_consecutive_data <- function(data) {
  data |>
    arrange(gene_id, type, timepoint) |>
    group_by(gene_id, type) |>
    mutate(previous_expression = lag(log2FoldChange),
           previous_timepoint  = lag(timepoint),
           expression_change   = log2FoldChange - previous_expression,
           transition          = paste0(previous_timepoint, "\u2192", timepoint)) |>
    filter(!is.na(expression_change)) |>
    ungroup()
}

conduct_consecutive_tests <- function(consecutive_data) {
  lapply(c("static", "gained", "lost"), function(lt) {
    type_data <- filter(consecutive_data, type == lt)
    res <- type_data |>
      group_by(transition) |>
      rstatix::wilcox_test(expression_change ~ 1, mu = 0) |>
      mutate(median_change = type_data |> group_by(transition) |>
               summarise(med = median(expression_change, na.rm = TRUE)) |> pull(med),
             n_genes = type_data |> group_by(transition) |>
               summarise(n = dplyr::n()) |> pull(n)) |>
      ungroup()
    res$p_adjusted <- p.adjust(res$p, method = "BH")
    res$significance_stars <- case_when(
      res$p_adjusted < 0.001 ~ "***",
      res$p_adjusted < 0.01  ~ "**",
      res$p_adjusted < 0.05  ~ "*",
      TRUE ~ "ns"
    )
    res
  }) |> setNames(c("static", "gained", "lost"))
}

consecutive_data <- prepare_consecutive_data(combined)
test_results     <- conduct_consecutive_tests(consecutive_data)

# Visualization ----
theme_boxplots <- theme_minimal() +
  theme(
    panel.background   = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(color = "black", linewidth = 0.5),
    axis.ticks         = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length  = unit(0.15, "cm"),
    axis.text          = element_text(size = 9, color = "black"),
    plot.margin        = margin(20, 25, 20, 15)
  )

tp_positions <- c("1" = 1, "3" = 2, "6" = 3, "9" = 4, "12" = 5, "24" = 6)

create_panel_with_stats <- function(loop_type, show_y_axis = TRUE,
                                    show_y_title = TRUE, show_x_title = TRUE) {
  plot_data         <- filter(combined, type == loop_type)
  loop_test_results <- test_results[[loop_type]]

  p <- ggplot(plot_data, aes(x = timepoint, y = log2FoldChange)) +
    geom_boxplot(fill = "white", color = loop_colors[loop_type],
                 outlier.shape = NA, width = 0.7, linewidth = 0.5) +
    geom_hline(yintercept = 0, color = "gray30", linetype = "dashed", linewidth = 0.4) +
    stat_summary(fun = median, geom = "line", aes(group = 1),
                 color = loop_colors[loop_type], linewidth = 0.8) +
    stat_summary(fun = median, geom = "point",
                 color = loop_colors[loop_type], size = 2, shape = 18) +
    annotate("text", x = -Inf, y = Inf,
             label = str_to_title(loop_type),
             hjust = -0.1, vjust = 1.5, size = 4, fontface = "bold",
             color = loop_colors[loop_type])

  y_pos    <- 0.85
  bar_step <- 0.08

  for (i in seq_len(nrow(loop_test_results))) {
    tps      <- str_split(loop_test_results$transition[i], "\u2192")[[1]]
    x_start  <- tp_positions[tps[1]] + 0.15
    x_end    <- tp_positions[tps[2]] - 0.15
    curr_y   <- y_pos - (i - 1) * bar_step
    stars    <- loop_test_results$significance_stars[i]
    a_bar    <- if (stars != "ns") 1.0 else 0.4
    a_txt    <- if (stars != "ns") 1.0 else 0.6

    p <- p +
      annotate("segment",
               x = x_start, xend = x_end, y = curr_y, yend = curr_y,
               color = scales::alpha(loop_colors[loop_type], a_bar), linewidth = 0.6) +
      annotate("text",
               x = (tp_positions[tps[1]] + tp_positions[tps[2]]) / 2,
               y = curr_y + 0.04, label = stars,
               size      = if (stars != "ns") 4 else 3.5,
               fontface  = if (stars != "ns") "bold" else "plain",
               color     = scales::alpha(loop_colors[loop_type], a_txt))
  }

  p +
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    theme_boxplots +
    theme(
      axis.text.y  = if (show_y_axis) element_text(size = 9) else element_blank(),
      axis.title.x = if (show_x_title) element_text(size = 10, face = "bold") else element_blank(),
      axis.title.y = if (show_y_title) element_text(size = 10, face = "bold") else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line(color = "black", linewidth = 0.5) else element_blank()
    ) +
    labs(x = if (show_x_title) "Hours after hyperosmotic stress" else "",
         y = if (show_y_title) "log2FoldChange(treated/untreated)" else "")
}

create_timecourse_panel <- function(loop_type, show_y_axis = TRUE, show_y_title = TRUE) {
  plot_data <- filter(timecourse_with_loops, type == loop_type)
  ggplot(plot_data, aes(x = timepoint, y = log2_expression)) +
    geom_boxplot(fill = "white", color = loop_colors[loop_type],
                 outlier.shape = NA, width = 0.7, linewidth = 0.5) +
    stat_summary(fun = median, geom = "line", aes(group = 1),
                 color = loop_colors[loop_type], linewidth = 0.8) +
    stat_summary(fun = median, geom = "point",
                 color = loop_colors[loop_type], size = 2, shape = 18) +
    annotate("text", x = -Inf, y = Inf, label = str_to_title(loop_type),
             hjust = -0.1, vjust = 1.5, size = 4, fontface = "bold",
             color = loop_colors[loop_type]) +
    coord_cartesian(ylim = c(0, 18)) +
    theme_boxplots +
    theme(
      axis.text.y  = if (show_y_axis) element_text(size = 9) else element_blank(),
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = if (show_y_title) element_text(size = 10, face = "bold") else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line(color = "black", linewidth = 0.5) else element_blank()
    ) +
    labs(x = "Hours after hyperosmotic stress",
         y = if (show_y_title) "Expression (log2 normalized counts + 1)" else "")
}

p_static <- create_panel_with_stats("static", TRUE,  TRUE,  FALSE)
p_gained <- create_panel_with_stats("gained", FALSE, FALSE, TRUE)
p_lost   <- create_panel_with_stats("lost",   FALSE, FALSE, FALSE)

integrated_plot <- plot_grid(p_static, p_gained, p_lost, ncol = 3, align = "h") |>
  (\(cp) ggdraw(cp) + draw_text("** p < 0.01, *** p < 0.001",
                                x = 0.98, y = 0.02, hjust = 1, vjust = 0,
                                size = 9, color = "gray60"))()

baseline_plot <- ggplot(baseline_with_loops,
                        aes(x = type, y = log2_baseline, fill = type)) +
  geom_boxplot(fill = "white", aes(color = type),
               outlier.shape = NA, width = 0.7, linewidth = 0.5) +
  stat_summary(fun = median, geom = "point", aes(color = type), size = 2, shape = 18) +
  scale_fill_manual(values = loop_colors) +
  scale_color_manual(values = loop_colors) +
  theme_boxplots +
  theme(legend.position = "none") +
  labs(x = "Loop Anchor Type",
       y = "Baseline Expression (log2 normalized counts + 1)")

tc_static <- create_timecourse_panel("static", TRUE,  TRUE)
tc_gained <- create_timecourse_panel("gained", FALSE, FALSE)
tc_lost   <- create_timecourse_panel("lost",   FALSE, FALSE)
timecourse_combined_plot <- plot_grid(tc_static, tc_gained, tc_lost, ncol = 3, align = "h")

# Save outputs ----
ggsave(output_main,    integrated_plot,       width = 12, height = 4.5, units = "in")
ggsave(output_square,  integrated_plot,       width =  8, height =  8,  units = "in")
ggsave(output_baseline, baseline_plot,        width =  6, height =  5,  units = "in")
ggsave(output_abexpr,  timecourse_combined_plot, width = 12, height = 4.5, units = "in")

sessionInfo()
