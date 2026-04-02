# ##############################################################################
# filename:    CTCF_RAD21_grouped_barplot.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Grouped barplots of CTCF and RAD21 peak overlap at gained vs
#              lost loop anchors under control and sorbitol conditions;
#              bootstrap confidence intervals included
# ##############################################################################

# Libraries ----
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

# Parameters ----
cutntag_peaks_dir <- "data/processed/cutntag/output/peaks/"
diff_loops_rds    <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
plots_dir         <- "plots"
target            <- c("CTCF", "H3K27ac", "RAD21", "YAP1")
condition         <- c("control", "sorbitol")

# Data import ----
cutntag <- list.files(cutntag_peaks_dir, full.names = TRUE, pattern = ".narrowPeak") |>
  lapply(read_narrowpeaks)

names(cutntag) <- paste0(rep(target, each = 2), "_", condition)

lapply(cutntag, length) |> print()

loops <- readRDS(diff_loops_rds) |> interactions()

gainedLoops <- loops[loops$padj <= 0.1 & loops$log2FoldChange > 1]
lostLoops   <- loops[loops$padj <= 0.1 & loops$log2FoldChange < -1]

cat("\nLoop counts:\n")
cat("Gained loops:", length(gainedLoops), "\n")
cat("Lost loops:", length(lostLoops), "\n")

# Analysis ----

## Helper functions ----
extractAnchors <- function(gi) {
  anchor1 <- anchors(gi, type = "first")
  anchor2 <- anchors(gi, type = "second")
  c(anchor1, anchor2) |> unique()
}

calculateCI <- function(anchors, peaks, n_bootstrap = 1000, conf_level = 0.95) {
  n_anchors <- length(anchors)

  boot_props <- replicate(n_bootstrap, {
    boot_idx     <- sample(n_anchors, replace = TRUE)
    boot_anchors <- anchors[boot_idx]
    hasOverlap   <- countOverlaps(boot_anchors, peaks) > 0
    sum(hasOverlap) / n_anchors
  })

  quantiles <- quantile(boot_props, probs = c((1 - conf_level) / 2,
                                               1 - (1 - conf_level) / 2))
  return(data.frame(lower = quantiles[1] * 100, upper = quantiles[2] * 100))
}

calculateOverlaps <- function(anchors, peaks_list, category) {
  overlaps <- lapply(peaks_list, function(peaks) {
    hasOverlap <- countOverlaps(anchors, peaks) > 0
    sum(hasOverlap) / length(anchors)
  })

  cis <- lapply(peaks_list, function(peaks) calculateCI(anchors, peaks))

  data.frame(
    Category   = category,
    Target     = sub("_.*", "", names(overlaps)),
    Condition  = sub(".*_", "", names(overlaps)),
    Percentage = unlist(overlaps) * 100,
    Lower      = sapply(cis, function(x) x$lower),
    Upper      = sapply(cis, function(x) x$upper)
  )
}

gained_anchors <- extractAnchors(gainedLoops)
lost_anchors   <- extractAnchors(lostLoops)

cat("\nAnchor counts:\n")
cat("Gained loop anchors:", length(gained_anchors), "\n")
cat("Lost loop anchors:", length(lost_anchors), "\n")

results <- rbind(
  calculateOverlaps(gained_anchors, cutntag, "Gained"),
  calculateOverlaps(lost_anchors,   cutntag, "Lost")
)

plot_data <- results |>
  mutate(
    Category = factor(Category, levels = c("Lost", "Gained")),
    Target   = factor(Target, levels = c("CTCF", "RAD21", "H3K27ac", "YAP1"))
  )

# Visualization ----

## CTCF grouped barplot ----
ctcf_data <- plot_data |> filter(Target == "CTCF")

p_ctcf_bar <- ggplot(ctcf_data, aes(x = Category, y = Percentage, fill = Condition)) +
  geom_hline(yintercept = seq(0, 100, 25), color = "gray90", linetype = "dashed") +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.9), width = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("control" = "#619CFF", "sorbitol" = "#F8766D")) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), axis.title.x = element_blank(),
        legend.position = "top", aspect.ratio = 1) +
  labs(y = "Percentage of anchors bound (%)", title = "CTCF Binding at Loop Anchors") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 100, 25))

## RAD21 grouped barplot ----
rad21_data <- plot_data |> filter(Target == "RAD21")

p_rad21_bar <- ggplot(rad21_data, aes(x = Category, y = Percentage, fill = Condition)) +
  geom_hline(yintercept = seq(0, 100, 25), color = "gray90", linetype = "dashed") +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.9), width = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("control" = "#619CFF", "sorbitol" = "#F8766D")) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), axis.title.x = element_blank(),
        legend.position = "top", aspect.ratio = 1) +
  labs(y = "Percentage of anchors bound (%)", title = "RAD21 Binding at Loop Anchors") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 100, 25))

# Save outputs ----
ggsave(file.path(plots_dir, "CTCF_grouped_bar.pdf"),  p_ctcf_bar,  width = 5, height = 5)
ggsave(file.path(plots_dir, "RAD21_grouped_bar.pdf"),  p_rad21_bar, width = 5, height = 5)

sessionInfo()
