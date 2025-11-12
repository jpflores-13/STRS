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

## Load cutntag data
cutntag <- list.files("data/processed/cutntag/output/peaks/",
                      full.names = TRUE,
                      pattern = ".narrowPeak") |> 
  lapply(read_narrowpeaks)

target <- c("CTCF", "H3K27ac", "RAD21", "YAP1")
condition <- c("control", "sorbitol")
names(cutntag) <- paste0(rep(target, each = 2), "_", condition)
lapply(cutntag, length)

## Load loops
loops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions() |> 
  as.data.frame() |> 
  mariner::as_ginteractions()

## separate loops into gained, lost, static
gainedLoops <- loops[loops$padj <= 0.05 & loops$log2FoldChange > 1]

lostLoops <- loops[loops$padj <= 0.05 & loops$log2FoldChange < -1]

# Helper Functions --------------------------------------------------------
## Function to convert GInteractions to GRanges
extractAnchors <- function(gi) {
  # Extract both anchors and combine them
  anchor1 <- anchors(gi, type="first")
  anchor2 <- anchors(gi, type="second")
  
  # Combine and remove duplicates
  c(anchor1, anchor2) |> 
    unique()
}

## Function to calculate bootstrap confidence intervals
calculateCI <- function(anchors, peaks, n_bootstrap = 1000, conf_level = 0.95) {
  # Total number of anchors
  n_anchors <- length(anchors)
  
  # Perform bootstrap resampling
  boot_props <- replicate(n_bootstrap, {
    # Sample anchors with replacement
    boot_idx <- sample(n_anchors, replace = TRUE)
    boot_anchors <- anchors[boot_idx]
    
    # Calculate proportion for this bootstrap sample
    hasOverlap <- countOverlaps(boot_anchors, peaks) > 0
    sum(hasOverlap) / n_anchors
  })
  
  # Calculate confidence intervals
  quantiles <- quantile(boot_props, probs = c((1-conf_level)/2, 1-(1-conf_level)/2))
  return(data.frame(
    lower = quantiles[1] * 100,  # Convert to percentage
    upper = quantiles[2] * 100
  ))
}

# Function to include confidence intervals during overlap
calculateOverlaps <- function(anchors, peaks_list, category) {
  # count overlaps (query = anchors, subject = peaks)
  overlaps <- lapply(peaks_list, function(peaks) {
    hasOverlap <- countOverlaps(anchors, peaks) > 0
    sum(hasOverlap) / length(anchors)
  })
  
  # Calculate confidence intervals for each comparison
  cis <- lapply(peaks_list, function(peaks) {
    calculateCI(anchors, peaks)
  })
  
  # Create df with CIs
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

## Calculate proportions for each category
results <- rbind(
  calculateOverlaps(gained_anchors, cutntag, "Gained"),
  calculateOverlaps(lost_anchors, cutntag, "Lost")
)

## Convert proportions to percentages
plot_data <- results |>
  mutate(
    Category = factor(Category, levels = c("Lost", "Gained")),
    Target = factor(Target, levels = c("CTCF", "RAD21", "H3K27ac", "YAP1"))
  )

# Create separate plots for each protein with flipped axes ----------------

# Function to create a plot for a specific protein with flipped axes
create_protein_plot <- function(data, protein_name) {
  protein_data <- data |> filter(Target == protein_name)
  
  ggplot(protein_data, aes(x = Category, y = Percentage, 
                           fill = Condition, group = Condition)) +
    geom_hline(yintercept = seq(0, 100, 25), 
               color = "gray90", linetype = "dashed") +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
                  position = position_dodge(width = 0.7),
                  width = 0.25) +
    # Add percentage labels on plot
    geom_text(aes(label = sprintf("%.1f%%", Percentage),
                  y = Percentage + 5),  
              position = position_dodge(width = 0.7),
              vjust = 0, size = 3) +
    scale_fill_manual(values = c("control" = "#619CFF", "sorbitol" = "#F8766D")) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),  
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(face = "bold"),
      axis.title.x = element_blank(),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "bottom"
    ) +
    labs(y = "Percentage of anchors bound (%)",
         title = protein_name) +
    scale_y_continuous(
      limits = c(0, 100),
      expand = expansion(mult = c(0, 0.15)),
      breaks = seq(0, 100, 25)
    )
}

# Create and save individual plots only
for (protein in levels(plot_data$Target)) {
  protein_plot <- create_protein_plot(plot_data, protein)
  ggsave(
    filename = paste0("plots/cutntag_HiCanchors_", protein, ".pdf"), 
    plot = protein_plot,
    width = 5, 
    height = 6
  )
}
