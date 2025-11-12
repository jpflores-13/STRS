library(data.table)
library(InteractionSet)
library(mariner)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(plyranges)
library(tidyr)
library(forcats)
library(ggpubr)

# Data --------------------------------------------------------------------

## Load cutntag data
cutntag <- list.files("data/processed/cutntag/output/peaks/",
                      full.names = TRUE,
                      pattern = ".narrowPeak") |> 
  lapply(read_narrowpeaks)

target <- c("CTCF", "H3K27ac", "RAD21", "YAP1")
condition <- c("control", "sorbitol")
names(cutntag) <- paste0(rep(target, each = 2), "_", condition)

## Check peak counts
lapply(cutntag, length) |> print()

## Load loops
loops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions()

## Separate loops into gained, lost, static
gainedLoops <- loops[loops$padj <= 0.1 & loops$log2FoldChange > 1]
lostLoops <- loops[loops$padj <= 0.1 & loops$log2FoldChange < -1]

cat("\nLoop counts:\n")
cat("Gained loops:", length(gainedLoops), "\n")
cat("Lost loops:", length(lostLoops), "\n")

# Helper Functions --------------------------------------------------------

## Function to convert GInteractions to GRanges
extractAnchors <- function(gi) {
  anchor1 <- anchors(gi, type = "first")
  anchor2 <- anchors(gi, type = "second")
  
  c(anchor1, anchor2) |> 
    unique()
}

## Function to calculate bootstrap confidence intervals
calculateCI <- function(anchors, peaks, n_bootstrap = 1000, conf_level = 0.95) {
  n_anchors <- length(anchors)
  
  boot_props <- replicate(n_bootstrap, {
    boot_idx <- sample(n_anchors, replace = TRUE)
    boot_anchors <- anchors[boot_idx]
    
    hasOverlap <- countOverlaps(boot_anchors, peaks) > 0
    sum(hasOverlap) / n_anchors
  })
  
  quantiles <- quantile(boot_props, probs = c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2))
  return(data.frame(
    lower = quantiles[1] * 100,
    upper = quantiles[2] * 100
  ))
}

## Function to include confidence intervals during overlap calculation
calculateOverlaps <- function(anchors, peaks_list, category) {
  overlaps <- lapply(peaks_list, function(peaks) {
    hasOverlap <- countOverlaps(anchors, peaks) > 0
    sum(hasOverlap) / length(anchors)
  })
  
  cis <- lapply(peaks_list, function(peaks) {
    calculateCI(anchors, peaks)
  })
  
  data.frame(
    Category = category,
    Target = sub("_.*", "", names(overlaps)),
    Condition = sub(".*_", "", names(overlaps)),
    Percentage = unlist(overlaps) * 100,
    Lower = sapply(cis, function(x) x$lower),
    Upper = sapply(cis, function(x) x$upper)
  )
}

# Analysis ----------------------------------------------------------------

## Extract anchors for each category
gained_anchors <- extractAnchors(gainedLoops)
lost_anchors <- extractAnchors(lostLoops)

cat("\nAnchor counts:\n")
cat("Gained loop anchors:", length(gained_anchors), "\n")
cat("Lost loop anchors:", length(lost_anchors), "\n")

## Calculate proportions for each category
results <- rbind(
  calculateOverlaps(gained_anchors, cutntag, "Gained"),
  calculateOverlaps(lost_anchors, cutntag, "Lost")
)

## Format data for plotting
plot_data <- results |>
  mutate(
    Category = factor(Category, levels = c("Lost", "Gained")),
    Target = factor(Target, levels = c("CTCF", "RAD21", "H3K27ac", "YAP1"))
  )

# CTCF Plots --------------------------------------------------------------

cat("\n", rep("=", 70), "\n", sep = "")
cat("CTCF Visualizations\n")
cat(rep("=", 70), "\n\n")

ctcf_data <- plot_data |> filter(Target == "CTCF")

## CTCF: Grouped Bar Plot
cat("=== CTCF: Grouped Bar Plot ===\n")

p_ctcf_bar <- ggplot(ctcf_data, aes(x = Category, y = Percentage, fill = Condition)) +
  geom_hline(yintercept = seq(0, 100, 25), 
             color = "gray90", linetype = "dashed") +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.9),
                width = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("control" = "#619CFF", "sorbitol" = "#F8766D")) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    aspect.ratio = 1
  ) +
  labs(y = "Percentage of anchors bound (%)",
       title = "CTCF Binding at Loop Anchors") +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.05)),
    breaks = seq(0, 100, 25)
  )

print(p_ctcf_bar)

ggsave("plots/CTCF_grouped_bar.pdf", p_ctcf_bar, width = 5, height = 5)

# RAD21 Plots -------------------------------------------------------------

cat("\n", rep("=", 70), "\n", sep = "")
cat("RAD21 Visualizations\n")
cat(rep("=", 70), "\n\n")

rad21_data <- plot_data |> filter(Target == "RAD21")

## RAD21: Grouped Bar Plot
cat("=== RAD21: Grouped Bar Plot ===\n")

p_rad21_bar <- ggplot(rad21_data, aes(x = Category, y = Percentage, fill = Condition)) +
  geom_hline(yintercept = seq(0, 100, 25), 
             color = "gray90", linetype = "dashed") +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.9),
                width = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("control" = "#619CFF", "sorbitol" = "#F8766D")) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    aspect.ratio = 1
  ) +
  labs(y = "Percentage of anchors bound (%)",
       title = "RAD21 Binding at Loop Anchors") +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.05)),
    breaks = seq(0, 100, 25)
  )

print(p_rad21_bar)

ggsave("plots/RAD21_grouped_bar.pdf", p_rad21_bar, width = 5, height = 5)
