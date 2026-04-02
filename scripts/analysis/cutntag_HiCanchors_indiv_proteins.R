# ##############################################################################
# filename:    cutntag_HiCanchors_indiv_proteins.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Per-protein barplots showing percent of loop anchors bound by
#              CTCF, H3K27ac, RAD21, and YAP1 under control and sorbitol
#              conditions; bootstrap 95% CIs included
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
cutntag <- list.files(cutntag_peaks_dir, full.names = TRUE,
                      pattern = ".narrowPeak") |>
  lapply(plyranges::read_narrowpeaks)

names(cutntag) <- paste0(rep(target, each = 2), "_", condition)
lapply(cutntag, length)

loops <- readRDS(diff_loops_rds) |>
  interactions() |>
  as.data.frame() |>
  mariner::as_ginteractions()

gainedLoops <- loops[loops$padj <= 0.05 & loops$log2FoldChange >  1]
lostLoops   <- loops[loops$padj <= 0.05 & loops$log2FoldChange < -1]

# Helper functions ----

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

  quantiles <- quantile(boot_props,
                        probs = c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2))
  data.frame(lower = quantiles[1] * 100, upper = quantiles[2] * 100)
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

create_protein_plot <- function(data, protein_name) {
  protein_data <- data |> filter(Target == protein_name)

  ggplot(protein_data, aes(x = Category, y = Percentage,
                           fill = Condition, group = Condition)) +
    geom_hline(yintercept = seq(0, 100, 25), color = "gray90", linetype = "dashed") +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
                  position = position_dodge(width = 0.7), width = 0.25) +
    geom_text(aes(label = sprintf("%.1f%%", Percentage), y = Percentage + 5),
              position = position_dodge(width = 0.7), vjust = 0, size = 3) +
    scale_fill_manual(values = c("control" = "#619CFF", "sorbitol" = "#F8766D")) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.background   = element_rect(fill = "white"),
          strip.text         = element_text(face = "bold"),
          axis.title.x       = element_blank(),
          plot.title         = element_text(face = "bold", size = 14, hjust = 0.5),
          legend.position    = "bottom") +
    labs(y = "Percentage of anchors bound (%)", title = protein_name) +
    scale_y_continuous(limits = c(0, 100),
                       expand = expansion(mult = c(0, 0.15)),
                       breaks = seq(0, 100, 25))
}

# Analysis ----
gained_anchors <- extractAnchors(gainedLoops)
lost_anchors   <- extractAnchors(lostLoops)

results <- rbind(
  calculateOverlaps(gained_anchors, cutntag, "Gained"),
  calculateOverlaps(lost_anchors,   cutntag, "Lost")
)

plot_data <- results |>
  mutate(
    Category = factor(Category, levels = c("Lost", "Gained")),
    Target   = factor(Target,   levels = c("CTCF", "RAD21", "H3K27ac", "YAP1"))
  )

# Save outputs ----
for (protein in levels(plot_data$Target)) {
  protein_plot <- create_protein_plot(plot_data, protein)
  ggsave(
    filename = file.path(plots_dir, paste0("cutntag_HiCanchors_", protein, ".pdf")),
    plot     = protein_plot,
    width    = 5,
    height   = 6
  )
}

sessionInfo()
