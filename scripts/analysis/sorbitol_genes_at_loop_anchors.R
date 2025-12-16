# Analysis: Overlap between differentially expressed genes and loop anchors
# Purpose: Verify claim that most sorbitol-responsive genes are NOT at gained/lost loop anchors

library(GenomicRanges)
library(InteractionSet)
library(DESeq2)
library(tidyverse)
library(mariner)

# Load data ---------------------------------------------------------------

dds <- readRDS("data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")

diff_loopCounts <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions()
diff_loopCounts <- keepStandardChromosomes(diff_loopCounts, pruning.mode = "coarse")

# Define loop categories --------------------------------------------------

gainedLoops <- diff_loopCounts[mcols(diff_loopCounts)$padj < 0.1 & 
                                 mcols(diff_loopCounts)$log2FoldChange > 0]
lostLoops   <- diff_loopCounts[mcols(diff_loopCounts)$padj < 0.1 & 
                                 mcols(diff_loopCounts)$log2FoldChange < 0]

# Extract loop anchors ----------------------------------------------------

gained_anchors <- c(anchors(gainedLoops, "first"), 
                    anchors(gainedLoops, "second")) |> 
  unique()

lost_anchors <- c(anchors(lostLoops, "first"), 
                  anchors(lostLoops, "second")) |> 
  unique()

# Get all loop anchors (for context)
all_anchors <- c(anchors(diff_loopCounts, "first"),
                 anchors(diff_loopCounts, "second")) |> 
  unique()

# Get DGE results as GRanges with gene info -------------------------------

dge_gr <- results(dds, format = "GRanges")
mcols(dge_gr) <- rowData(dds)

# Create promoter regions (TSS +/- some window, standard is +/- 1kb or 2kb)
# Using same approach as in your Panel C code
dge_gr_prom <- promoters(dge_gr, upstream = 2000, downstream = 2000) |> 
  keepStandardChromosomes(pruning.mode = "coarse")
seqlevelsStyle(dge_gr_prom) <- "UCSC"

# Define differentially expressed genes -----------------------------------

# Get all time coefficients
time_coefficients <- resultsNames(dds) |> 
  str_subset("Time")

# Extract DEGs at ANY timepoint (padj < 0.05, |LFC| > 2)
all_degs <- lapply(time_coefficients, function(coef) {
  results(dds, name = coef) |>
    as.data.frame() |>
    filter(padj < 0.05 & abs(log2FoldChange) > 2) |>
    rownames()
}) |>
  unlist() |>
  unique()

# Separate upregulated and downregulated genes
upreg_genes <- lapply(time_coefficients, function(coef) {
  results(dds, name = coef) |>
    as.data.frame() |>
    filter(padj < 0.05 & log2FoldChange > 2) |>
    rownames()
}) |>
  unlist() |>
  unique()

downreg_genes <- lapply(time_coefficients, function(coef) {
  results(dds, name = coef) |>
    as.data.frame() |>
    filter(padj < 0.05 & log2FoldChange < -2) |>
    rownames()
}) |>
  unlist() |>
  unique()

# Create GRanges of DEG promoters -----------------------------------------

all_deg_prom <- dge_gr_prom[names(dge_gr_prom) %in% all_degs]
upreg_prom <- dge_gr_prom[names(dge_gr_prom) %in% upreg_genes]
downreg_prom <- dge_gr_prom[names(dge_gr_prom) %in% downreg_genes]

# Calculate overlaps ------------------------------------------------------

# DEGs overlapping with gained loop anchors
deg_at_gained <- subsetByOverlaps(all_deg_prom, gained_anchors)
upreg_at_gained <- subsetByOverlaps(upreg_prom, gained_anchors)
downreg_at_gained <- subsetByOverlaps(downreg_prom, gained_anchors)

# DEGs overlapping with lost loop anchors
deg_at_lost <- subsetByOverlaps(all_deg_prom, lost_anchors)
upreg_at_lost <- subsetByOverlaps(upreg_prom, lost_anchors)
downreg_at_lost <- subsetByOverlaps(downreg_prom, lost_anchors)

# DEGs overlapping with ANY differential loop anchor (gained OR lost)
diff_anchors <- c(gained_anchors, lost_anchors) |> unique()
deg_at_any_diff <- subsetByOverlaps(all_deg_prom, diff_anchors)

# DEGs NOT at any differential loop anchor
# Use countOverlaps to find which DEGs have NO overlap with differential anchors
overlap_counts <- countOverlaps(all_deg_prom, diff_anchors)
deg_NOT_at_diff <- all_deg_prom[overlap_counts == 0]

# Summary statistics ------------------------------------------------------

results_summary <- tibble(
  category = c(
    "Total DEGs",
    "DEGs at gained anchors",
    "DEGs at lost anchors", 
    "DEGs at ANY differential anchor",
    "DEGs NOT at differential anchors",
    "",
    "Upregulated genes total",
    "Upregulated at gained anchors",
    "Upregulated at lost anchors",
    "",
    "Downregulated genes total", 
    "Downregulated at gained anchors",
    "Downregulated at lost anchors"
  ),
  count = c(
    length(all_degs),
    length(deg_at_gained),
    length(deg_at_lost),
    length(deg_at_any_diff),
    length(deg_NOT_at_diff),
    NA,
    length(upreg_genes),
    length(upreg_at_gained),
    length(upreg_at_lost),
    NA,
    length(downreg_genes),
    length(downreg_at_gained),
    length(downreg_at_lost)
  ),
  percentage = c(
    100,
    round(100 * length(deg_at_gained) / length(all_degs), 1),
    round(100 * length(deg_at_lost) / length(all_degs), 1),
    round(100 * length(deg_at_any_diff) / length(all_degs), 1),
    round(100 * length(deg_NOT_at_diff) / length(all_degs), 1),
    NA,
    100,
    round(100 * length(upreg_at_gained) / length(upreg_genes), 1),
    round(100 * length(upreg_at_lost) / length(upreg_genes), 1),
    NA,
    100,
    round(100 * length(downreg_at_gained) / length(downreg_genes), 1),
    round(100 * length(downreg_at_lost) / length(downreg_genes), 1)
  )
)


# Create visualization ----------------------------------------------------

# Create plots directory if it doesn't exist
if (!dir.exists("plots")) dir.create("plots")

# Color palette
color_at_anchor <- "#7570B3"
color_not_at_anchor <- "#E7298A"

# Prepare data for the key comparison
plot_data <- tibble(
  location = c("At anchors", "Not at anchors"),
  count = c(length(deg_at_any_diff), length(deg_NOT_at_diff)),
  percentage = c(
    100 * length(deg_at_any_diff) / length(all_degs),
    100 * length(deg_NOT_at_diff) / length(all_degs)
  )
) |>
  mutate(location = factor(location, levels = location))

# Create single focused plot
p <- ggplot(plot_data, aes(x = location, y = percentage, fill = location)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("n = %d\n%.1f%%", count, percentage)), 
            vjust = -0.3, size = 5) +
  scale_fill_manual(values = c(color_at_anchor, color_not_at_anchor)) +
  labs(
    title = "Sorbitol-responsive genes at differential loop anchors",
    x = NULL,
    y = "Percentage of all DEGs (%)"
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave("plots/sorbitol_genes_at_loop_anchors.pdf", p, width = 8, height = 6)