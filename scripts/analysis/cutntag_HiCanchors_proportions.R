## Q: What % of loop anchors are bound by our proteins of interest?

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

# Print summary results
cat("\nSummary Results (% anchors bound):\n")
print(plot_data |> 
        dplyr::select(Category, Target, Condition, Percentage) |>
        pivot_wider(names_from = Condition, values_from = Percentage) |>
        arrange(Category, Target))

# Plot --------------------------------------------------------------------

pdf("plots/cutntag_HiCanchors_proportions.pdf", width = 8, height = 10)

ggplot(plot_data, aes(x = Percentage, y = Category, 
                      fill = Condition, group = Condition)) +
  geom_vline(xintercept = seq(0, 100, 25), 
             color = "gray90", linetype = "dashed") +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper),
                position = position_dodge(width = 0.9),
                width = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage),
                x = Percentage + 2),  
            position = position_dodge(width = 0.9),
            hjust = -0.75, size = 3) +
  facet_wrap(~Target, scales = "fixed", ncol = 1) +
  scale_fill_manual(values = c("control" = "#619CFF", "sorbitol" = "#F8766D")) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),  
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    axis.title.y = element_blank()
  ) +
  labs(x = "Percentage of anchors bound (%)",
       title = "Loop Anchor Binding by Chromatin Proteins") +
  scale_x_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.15)),
    breaks = seq(0, 100, 25)
  )

dev.off()

cat("\nAnalysis complete! Plot saved to: plots/cutntag_HiCanchors_proportions.pdf\n")
