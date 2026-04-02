# ##############################################################################
# filename:    loop_size_comparison.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Compares sizes of gained, lost, and static loops; Wilcoxon
#              pairwise tests with BH correction and effect sizes; produces
#              boxplot, violin, scatter, histogram, density, and cumulative
#              distribution plots
# ##############################################################################

# Libraries ----
library(InteractionSet)
library(mariner)
library(ggplot2)
library(dplyr)
library(plyranges)
library(patchwork)
library(rstatix)
library(ggpubr)

# Parameters ----
diff_loops_rds      <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
output_pdf          <- "plots/loop_size_comparison_all.pdf"
output_combined_pdf <- "plots/loop_size_comparison_combined.pdf"

# Data import ----
loops <- readRDS(diff_loops_rds) |>
  interactions() |>
  as.data.frame() |>
  as_ginteractions()

gainedLoops <- loops |> subset(padj <= 0.1 & log2FoldChange > 0)
lostLoops   <- loops |> subset(padj <= 0.1 & log2FoldChange < 0)
staticLoops <- loops |> subset(padj > 0.1)

cat("Loop counts:\n")
cat("Gained loops:", length(gainedLoops), "\n")
cat("Lost loops:",   length(lostLoops),   "\n")
cat("Static loops:", length(staticLoops), "\n")

# Analysis ----

## Calculate loop distances ----
gained_distances <- pairdist(gainedLoops, type = "mid")
lost_distances   <- pairdist(lostLoops,   type = "mid")
static_distances <- pairdist(staticLoops, type = "mid")

distance_data <- data.frame(
  distance = c(gained_distances, lost_distances, static_distances),
  category = factor(
    c(rep("Gained", length(gained_distances)),
      rep("Lost",   length(lost_distances)),
      rep("Static", length(static_distances))),
    levels = c("Static", "Lost", "Gained")
  )
)

## Summary statistics ----
summary_stats <- distance_data |>
  group_by(category) |>
  summarise(
    n      = dplyr::n(),
    mean   = mean(distance),
    median = median(distance),
    sd     = sd(distance),
    q25    = quantile(distance, 0.25),
    q75    = quantile(distance, 0.75),
    min    = min(distance),
    max    = max(distance),
    .groups = "drop"
  )

cat("\n\nSummary Statistics (distances in bp):\n")
print(summary_stats)

## Statistical testing ----
cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("STATISTICAL ANALYSIS\n")
cat(rep("=", 70), "\n", sep = "")
cat("\nApproach: Wilcoxon rank-sum test (Mann-Whitney U)\n")
cat("Multiple testing correction: Benjamini-Hochberg FDR (3 comparisons)\n\n")

stat_results <- distance_data |>
  wilcox_test(distance ~ category, p.adjust.method = "BH", detailed = TRUE) |>
  add_significance("p.adj")

stat_results <- stat_results |>
  mutate(
    effect_size_r  = abs(statistic) / sqrt(n1 + n2),
    effect_magnitude = case_when(
      effect_size_r < 0.1 ~ "negligible",
      effect_size_r < 0.3 ~ "small",
      effect_size_r < 0.5 ~ "medium",
      TRUE ~ "large"
    )
  )

cat("Pairwise Comparisons:\n")
cat(rep("-", 70), "\n", sep = "")
print(stat_results |>
        dplyr::select(group1, group2, n1, n2, statistic, p, p.adj,
                      p.adj.signif, effect_size_r, effect_magnitude))

for (i in seq_len(nrow(stat_results))) {
  cat("\n", stat_results$group1[i], " vs ", stat_results$group2[i], ":\n", sep = "")
  cat("  n1 =", stat_results$n1[i], ", n2 =", stat_results$n2[i], "\n")
  cat("  Median difference: ",
      round(median(distance_data$distance[distance_data$category == stat_results$group1[i]]) -
              median(distance_data$distance[distance_data$category == stat_results$group2[i]])),
      " bp\n", sep = "")
  cat("  p-value:", format.pval(stat_results$p[i], digits = 3), "\n")
  cat("  Adjusted p-value (BH):", format.pval(stat_results$p.adj[i], digits = 3),
      stat_results$p.adj.signif[i], "\n")
  cat("  Effect size (r):", round(stat_results$effect_size_r[i], 3),
      "(", stat_results$effect_magnitude[i], ")\n")
}

## Significance bracket positions ----
stat_for_plot <- stat_results |>
  mutate(
    y.position = max(log10(distance_data$distance / 1000)) + 0.2 + (row_number() - 1) * 0.15,
    xmin = as.numeric(factor(group1, levels = c("Static", "Lost", "Gained"))),
    xmax = as.numeric(factor(group2, levels = c("Static", "Lost", "Gained")))
  )

# Visualization ----
color_palette <- c("Gained" = "#E64B35", "Lost" = "#4DBBD5", "Static" = "grey70")

category_labels <- summary_stats |>
  mutate(label = paste0(category, "\n(n = ", format(n, big.mark = ","), ")"))

## Boxplot ----
p_boxplot <- ggplot(distance_data, aes(x = category, y = distance / 1000, fill = category)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3, outlier.size = 0.5) +
  scale_fill_manual(values = color_palette) +
  scale_x_discrete(labels = setNames(category_labels$label, category_labels$category)) +
  scale_y_log10(breaks = c(10, 50, 100, 500, 1000, 5000),
                labels = c("10", "50", "100", "500", "1000", "5000")) +
  stat_pvalue_manual(stat_for_plot, label = "p.adj.signif",
                     tip.length = 0.01, y.position = "y.position", bracket.size = 0.5) +
  labs(title = "Loop Size Distribution", x = NULL, y = "Loop Size (kb, log scale)") +
  theme_classic() +
  theme(plot.title  = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "none", aspect.ratio = 1, axis.text.x = element_text(size = 9))

## Violin plot ----
p_violin <- ggplot(distance_data, aes(x = category, y = distance / 1000, fill = category)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = color_palette) +
  scale_x_discrete(labels = setNames(category_labels$label, category_labels$category)) +
  scale_y_log10(breaks = c(10, 50, 100, 500, 1000, 5000),
                labels = c("10", "50", "100", "500", "1000", "5000")) +
  stat_pvalue_manual(stat_for_plot, label = "p.adj.signif",
                     tip.length = 0.01, y.position = "y.position", bracket.size = 0.5) +
  labs(title = "Loop Size Distribution", x = NULL, y = "Loop Size (kb, log scale)") +
  theme_classic() +
  theme(plot.title  = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "none", aspect.ratio = 1, axis.text.x = element_text(size = 9))

## Scatter with jitter ----
p_scatter <- ggplot(distance_data, aes(x = category, y = distance / 1000, color = category)) +
  geom_jitter(alpha = 0.3, size = 0.5, width = 0.2) +
  geom_boxplot(aes(fill = category), alpha = 0.3, outlier.shape = NA, color = "black") +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values  = color_palette) +
  scale_x_discrete(labels = setNames(category_labels$label, category_labels$category)) +
  scale_y_log10(breaks = c(10, 50, 100, 500, 1000, 5000),
                labels = c("10", "50", "100", "500", "1000", "5000")) +
  stat_pvalue_manual(stat_for_plot, label = "p.adj.signif",
                     tip.length = 0.01, y.position = "y.position", bracket.size = 0.5) +
  labs(title = "Loop Size Distribution", x = NULL, y = "Loop Size (kb, log scale)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "none", aspect.ratio = 1, axis.text.x = element_text(size = 9))

## Histogram ----
p_histogram <- ggplot(distance_data, aes(x = distance / 1000, fill = category)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
  scale_fill_manual(values = color_palette,
                    labels = setNames(category_labels$label, category_labels$category)) +
  scale_x_log10(breaks = c(10, 50, 100, 500, 1000, 5000),
                labels = c("10", "50", "100", "500", "1000", "5000")) +
  labs(title = "Loop Size Distribution",
       x = "Loop Size (kb, log scale)", y = "Count", fill = NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        aspect.ratio = 0.75, legend.position = "bottom", legend.text = element_text(size = 9))

## Density curves ----
p_density <- ggplot(distance_data,
                    aes(x = distance / 1000, fill = category, color = category)) +
  geom_density(alpha = 0.3, linewidth = 1) +
  scale_fill_manual(values  = color_palette,
                    labels = setNames(category_labels$label, category_labels$category)) +
  scale_color_manual(values = color_palette,
                     labels = setNames(category_labels$label, category_labels$category)) +
  scale_x_log10(breaks = c(10, 50, 100, 500, 1000, 5000),
                labels = c("10", "50", "100", "500", "1000", "5000")) +
  labs(title = "Loop Size Distribution",
       x = "Loop Size (kb, log scale)", y = "Density", fill = NULL, color = NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        aspect.ratio = 0.75, legend.position = "bottom", legend.text = element_text(size = 9))

## Cumulative distribution ----
distance_data_sorted <- distance_data |>
  arrange(category, distance) |>
  group_by(category) |>
  mutate(cumulative_prop = seq_along(distance) / dplyr::n()) |>
  ungroup()

p_cumulative <- ggplot(distance_data_sorted,
                       aes(x = distance / 1000, y = cumulative_prop, color = category)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = color_palette,
                     labels = setNames(category_labels$label, category_labels$category)) +
  scale_x_log10(breaks = c(10, 50, 100, 500, 1000, 5000),
                labels = c("10", "50", "100", "500", "1000", "5000")) +
  labs(title = "Cumulative Distribution",
       x = "Loop Size (kb, log scale)", y = "Cumulative Proportion", color = NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        aspect.ratio = 0.75, legend.position = "bottom", legend.text = element_text(size = 9))

# Save outputs ----
pdf(output_pdf, width = 10, height = 12)
print(p_boxplot)
print(p_violin)
print(p_scatter)
print(p_histogram)
print(p_density)
print(p_cumulative)
dev.off()

combined_plot <- (p_boxplot | p_violin) / (p_scatter | p_density) / (p_histogram | p_cumulative)
ggsave(output_combined_pdf, plot = combined_plot, width = 12, height = 16)

cat("\n\nAnalysis complete!\n")
cat("Individual plots:", output_pdf, "\n")
cat("Combined panel:", output_combined_pdf, "\n")

sessionInfo()
