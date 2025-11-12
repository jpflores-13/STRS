## Compare sizes of preexisting vs stress-induced loops

library(InteractionSet)
library(mariner)
library(ggplot2)
library(dplyr)
library(plyranges)
library(patchwork)
library(rstatix)
library(ggpubr)

# Data --------------------------------------------------------------------

## Load differential loops
loops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions() |> 
  as.data.frame() |> 
  as_ginteractions()

## Separate loops into gained, lost, static
gainedLoops <- loops |> 
  subset(padj <= 0.1 & log2FoldChange > 0)

lostLoops <- loops |> 
  subset(padj <= 0.1 & log2FoldChange < 0)

staticLoops <- loops |> 
  subset(padj > 0.1)

cat("Loop counts:\n")
cat("Gained loops:", length(gainedLoops), "\n")
cat("Lost loops:", length(lostLoops), "\n")
cat("Static loops:", length(staticLoops), "\n")

# Calculate Loop Distances ------------------------------------------------

## pairdist calculates the distance between anchors of each interaction
gained_distances <- pairdist(gainedLoops, type = "mid")
lost_distances <- pairdist(lostLoops, type = "mid")
static_distances <- pairdist(staticLoops, type = "mid")

## Create data frame for plotting
distance_data <- data.frame(
  distance = c(gained_distances, lost_distances, static_distances),
  category = factor(
    c(
      rep("Gained", length(gained_distances)),
      rep("Lost", length(lost_distances)),
      rep("Static", length(static_distances))
    ),
    levels = c("Static", "Lost", "Gained")
  )
)

# Summary Statistics ------------------------------------------------------

summary_stats <- distance_data |>
  group_by(category) |>
  summarise(
    n = dplyr::n(),
    mean = mean(distance),
    median = median(distance),
    sd = sd(distance),
    q25 = quantile(distance, 0.25),
    q75 = quantile(distance, 0.75),
    min = min(distance),
    max = max(distance),
    .groups = "drop"
  )

cat("\n\nSummary Statistics (distances in bp):\n")
print(summary_stats)

# Statistical Testing -----------------------------------------------------

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("STATISTICAL ANALYSIS\n")
cat(rep("=", 70), "\n", sep = "")

cat("\nApproach: Wilcoxon rank-sum test (Mann-Whitney U)\n")
cat("Rationale: Non-parametric test appropriate for non-normal distributions\n")
cat("           and unequal sample sizes. Tests whether distributions differ.\n")
cat("Multiple testing correction: Benjamini-Hochberg FDR (3 comparisons)\n\n")

## Perform pairwise Wilcoxon tests with effect sizes
stat_results <- distance_data |>
  wilcox_test(distance ~ category, 
              p.adjust.method = "BH",
              detailed = TRUE) |>
  add_significance("p.adj")

## Calculate effect sizes (r = Z / sqrt(N))
stat_results <- stat_results |>
  mutate(
    effect_size_r = abs(statistic) / sqrt(n1 + n2),
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

## Print detailed results
cat("\n\nDetailed Results:\n")
cat(rep("-", 70), "\n", sep = "")

for (i in 1:nrow(stat_results)) {
  cat("\n", stat_results$group1[i], " vs ", stat_results$group2[i], ":\n", sep = "")
  cat("  n1 = ", stat_results$n1[i], ", n2 = ", stat_results$n2[i], "\n", sep = "")
  cat("  Median difference: ", 
      round(median(distance_data$distance[distance_data$category == stat_results$group1[i]]) - 
              median(distance_data$distance[distance_data$category == stat_results$group2[i]])), 
      " bp\n", sep = "")
  cat("  p-value: ", format.pval(stat_results$p[i], digits = 3), "\n", sep = "")
  cat("  Adjusted p-value (BH): ", format.pval(stat_results$p.adj[i], digits = 3), 
      " ", stat_results$p.adj.signif[i], "\n", sep = "")
  cat("  Effect size (r): ", round(stat_results$effect_size_r[i], 3), 
      " (", stat_results$effect_magnitude[i], ")\n", sep = "")
}

# Prepare significance annotations ----------------------------------------

## Create significance bracket data for plots
stat_for_plot <- stat_results |>
  mutate(
    y.position = max(log10(distance_data$distance / 1000)) + 0.2 + (row_number() - 1) * 0.15,
    xmin = as.numeric(factor(group1, levels = c("Static", "Lost", "Gained"))),
    xmax = as.numeric(factor(group2, levels = c("Static", "Lost", "Gained")))
  )

# Visualization -----------------------------------------------------------

## Define colors
color_palette <- c(
  "Gained" = "#E64B35",
  "Lost" = "#4DBBD5",
  "Static" = "grey70"
)

## Create labels with n counts
category_labels <- summary_stats |>
  mutate(label = paste0(category, "\n(n = ", format(n, big.mark = ","), ")"))

## 1. Boxplot with significance
p_boxplot <- ggplot(distance_data, aes(x = category, y = distance / 1000, fill = category)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3, outlier.size = 0.5) +
  scale_fill_manual(values = color_palette) +
  scale_x_discrete(labels = setNames(category_labels$label, category_labels$category)) +
  scale_y_log10(
    breaks = c(10, 50, 100, 500, 1000, 5000),
    labels = c("10", "50", "100", "500", "1000", "5000")
  ) +
  stat_pvalue_manual(
    stat_for_plot,
    label = "p.adj.signif",
    tip.length = 0.01,
    y.position = "y.position",
    bracket.size = 0.5
  ) +
  labs(
    title = "Loop Size Distribution",
    x = NULL,
    y = "Loop Size (kb, log scale)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "none",
    aspect.ratio = 1,
    axis.text.x = element_text(size = 9)
  )

## 2. Violin plot with boxplot overlay and significance
p_violin <- ggplot(distance_data, aes(x = category, y = distance / 1000, fill = category)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = color_palette) +
  scale_x_discrete(labels = setNames(category_labels$label, category_labels$category)) +
  scale_y_log10(
    breaks = c(10, 50, 100, 500, 1000, 5000),
    labels = c("10", "50", "100", "500", "1000", "5000")
  ) +
  stat_pvalue_manual(
    stat_for_plot,
    label = "p.adj.signif",
    tip.length = 0.01,
    y.position = "y.position",
    bracket.size = 0.5
  ) +
  labs(
    title = "Loop Size Distribution",
    x = NULL,
    y = "Loop Size (kb, log scale)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "none",
    aspect.ratio = 1,
    axis.text.x = element_text(size = 9)
  )

## 3. Histogram
p_histogram <- ggplot(distance_data, aes(x = distance / 1000, fill = category)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
  scale_fill_manual(
    values = color_palette,
    labels = setNames(category_labels$label, category_labels$category)
  ) +
  scale_x_log10(
    breaks = c(10, 50, 100, 500, 1000, 5000),
    labels = c("10", "50", "100", "500", "1000", "5000")
  ) +
  labs(
    title = "Loop Size Distribution",
    x = "Loop Size (kb, log scale)",
    y = "Count",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    aspect.ratio = 0.75,
    legend.position = "bottom",
    legend.text = element_text(size = 9)
  )

## 4. Density curves
p_density <- ggplot(distance_data, aes(x = distance / 1000, fill = category, color = category)) +
  geom_density(alpha = 0.3, linewidth = 1) +
  scale_fill_manual(
    values = color_palette,
    labels = setNames(category_labels$label, category_labels$category)
  ) +
  scale_color_manual(
    values = color_palette,
    labels = setNames(category_labels$label, category_labels$category)
  ) +
  scale_x_log10(
    breaks = c(10, 50, 100, 500, 1000, 5000),
    labels = c("10", "50", "100", "500", "1000", "5000")
  ) +
  labs(
    title = "Loop Size Distribution",
    x = "Loop Size (kb, log scale)",
    y = "Density",
    fill = NULL,
    color = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    aspect.ratio = 0.75,
    legend.position = "bottom",
    legend.text = element_text(size = 9)
  )

## 5. Scatter plot with jitter and significance
p_scatter <- ggplot(distance_data, aes(x = category, y = distance / 1000, color = category)) +
  geom_jitter(alpha = 0.3, size = 0.5, width = 0.2) +
  geom_boxplot(aes(fill = category), alpha = 0.3, outlier.shape = NA, color = "black") +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  scale_x_discrete(labels = setNames(category_labels$label, category_labels$category)) +
  scale_y_log10(
    breaks = c(10, 50, 100, 500, 1000, 5000),
    labels = c("10", "50", "100", "500", "1000", "5000")
  ) +
  stat_pvalue_manual(
    stat_for_plot,
    label = "p.adj.signif",
    tip.length = 0.01,
    y.position = "y.position",
    bracket.size = 0.5
  ) +
  labs(
    title = "Loop Size Distribution",
    x = NULL,
    y = "Loop Size (kb, log scale)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "none",
    aspect.ratio = 1,
    axis.text.x = element_text(size = 9)
  )

## 6. Cumulative distribution plot
distance_data_sorted <- distance_data |>
  arrange(category, distance) |>
  group_by(category) |>
  mutate(cumulative_prop = seq_along(distance) / dplyr::n()) |>
  ungroup()

p_cumulative <- ggplot(distance_data_sorted, aes(x = distance / 1000, 
                                                 y = cumulative_prop, 
                                                 color = category)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = color_palette,
    labels = setNames(category_labels$label, category_labels$category)
  ) +
  scale_x_log10(
    breaks = c(10, 50, 100, 500, 1000, 5000),
    labels = c("10", "50", "100", "500", "1000", "5000")
  ) +
  labs(
    title = "Cumulative Distribution",
    x = "Loop Size (kb, log scale)",
    y = "Cumulative Proportion",
    color = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    aspect.ratio = 0.75,
    legend.position = "bottom",
    legend.text = element_text(size = 9)
  )

## Save individual plots
pdf("plots/loop_size_comparison_all.pdf", width = 10, height = 12)
print(p_boxplot)
print(p_violin)
print(p_scatter)
print(p_histogram)
print(p_density)
print(p_cumulative)
dev.off()

## Create combined plot panel
combined_plot <- (p_boxplot | p_violin) / (p_scatter | p_density) / (p_histogram | p_cumulative)

ggsave("plots/loop_size_comparison_combined.pdf",
       plot = combined_plot,
       width = 12,
       height = 16)

cat("\n\nAnalysis complete!\n")
cat("Individual plots saved to: plots/loop_size_comparison_all.pdf\n")
cat("Combined panel saved to: plots/loop_size_comparison_combined.pdf\n")
cat("Data saved to: data/processed/hic/loop_size_comparison.rds\n")
cat("Statistical results saved to: data/processed/hic/loop_size_statistical_tests.rds\n")

## Print session info
sessionInfo()