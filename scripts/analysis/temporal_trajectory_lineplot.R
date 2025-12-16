# Transcriptional Response at Differential Loop Anchors

library(DESeq2)
library(InteractionSet)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(tidyverse)
library(mariner)
library(rstatix)

# color palette -----------------------------------------------------------

color_upregulated <- "#F8766D"
color_downregulated <- "#619CFF"
color_static <- "#999999"

loop_colors <- c(
  static = color_static,
  gained = color_upregulated,
  lost = color_downregulated
)

# load data ---------------------------------------------------------------

dds <- readRDS("data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")

diff_loopCounts <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions()
diff_loopCounts <- keepStandardChromosomes(diff_loopCounts, pruning.mode = "coarse")

# prepare data ------------------------------------------------------------

dge_gr <- results(dds, format = "GRanges")
mcols(dge_gr) <- rowData(dds)
dge_gr_prom <- promoters(dge_gr) |> 
  keepStandardChromosomes(pruning.mode = "coarse")
seqlevelsStyle(dge_gr_prom) <- "UCSC"

as_gi <- function(x) x |> as.data.frame() |> as_ginteractions()

gained <- as_gi(diff_loopCounts) |> 
  subset(padj <= 0.1 & log2FoldChange > 0)
lost <- as_gi(diff_loopCounts) |> 
  subset(padj <= 0.1 & log2FoldChange < 0)
static <- as_gi(diff_loopCounts) |> 
  subset(padj > 0.1)

gained_df <- as.data.frame(subsetByOverlaps(dge_gr_prom, gained)) |> 
  mutate(type = "gained")
lost_df <- as.data.frame(subsetByOverlaps(dge_gr_prom, lost)) |> 
  mutate(type = "lost")
static_df <- as.data.frame(subsetByOverlaps(dge_gr_prom, static)) |> 
  mutate(type = "static")

combined_raw <- bind_rows(gained_df, lost_df, static_df)
combined_raw$padj_lrt <- p.adjust(combined_raw$LRTPvalue, method = "BH")

combined <- combined_raw |>
  filter(padj_lrt < 0.1) |>
  dplyr::select(starts_with("Time_"), type, symbol, gene_id) |>
  pivot_longer(
    cols = starts_with("Time"),
    values_to = "log2FoldChange",
    names_to = "timepoint"
  ) |>
  mutate(
    timepoint = recode(
      timepoint,
      "Time_1h_vs_0h" = "1",
      "Time_3h_vs_0h" = "3",
      "Time_6h_vs_0h" = "6",
      "Time_9h_vs_0h" = "9",
      "Time_12h_vs_0h" = "12",
      "Time_24h_vs_0h" = "24"
    ),
    timepoint_numeric = as.numeric(as.character(timepoint)),
    timepoint = factor(timepoint, levels = c("1", "3", "6", "9", "12", "24")),
    type = factor(type, levels = c("static", "gained", "lost"))
  )

# statistical analysis ----------------------------------------------------

stat_test <- combined |>
  filter(type == "gained") |>
  group_by(timepoint, timepoint_numeric) |>
  t_test(log2FoldChange ~ 1, mu = 0) |>
  adjust_pvalue(method = "BH") |>
  add_significance("p.adj") |>
  mutate(
    y_position = 0.40,  # Moved up from 0.25 to 0.40
    significance_star = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

stat_annotations <- stat_test |>
  filter(significance_star != "ns")

# prepare data for direct labels ------------------------------------------

label_data <- combined |>
  filter(timepoint_numeric == 24) |>
  group_by(type) |>
  summarise(
    x = 24,
    y = median(log2FoldChange),
    label = case_when(
      type == "static" ~ "Static",
      type == "gained" ~ "Gained",
      type == "lost" ~ "Lost"
    ),
    .groups = "drop"
  ) |>
  mutate(
    # Position labels to the right and slightly below each point
    x_adjusted = x + 0.5,
    y_adjusted = y - 0.03
  )

# create line plot --------------------------------------------------------

line_plot <- ggplot(combined, aes(x = timepoint_numeric, y = log2FoldChange, color = type, fill = type)) +
  stat_summary(
    fun = median,
    geom = "line",
    linewidth = 1.5,
    aes(group = type),
    position = position_dodge(width = 0.5)
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 4,
    shape = 21,
    fill = "white",
    stroke = 1.5,
    position = position_dodge(width = 0.5)
  ) +
  geom_hline(yintercept = 0, color = "gray30", linetype = "dashed", linewidth = 0.4) +
  {if(nrow(stat_annotations) > 0) {
    geom_text(
      data = stat_annotations,
      aes(x = timepoint_numeric, y = y_position, label = significance_star),
      color = color_upregulated,
      size = 6,
      inherit.aes = FALSE,
      show.legend = FALSE,
      fontface = "bold"
    )
  }} +
  geom_text(
    data = label_data,
    aes(x = x_adjusted, y = y_adjusted, label = label, color = type),
    hjust = 0,
    vjust = 1,
    size = 4.5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_color_manual(values = loop_colors) +
  scale_fill_manual(values = loop_colors) +
  scale_x_continuous(
    breaks = c(1, 3, 6, 9, 12, 24),
    labels = c("1", "3", "6", "9", "12", "24")
  ) +
  scale_y_continuous(
    breaks = seq(-0.5, 0.5, 0.25)
  ) +
  coord_cartesian(
    xlim = c(1, 28),
    ylim = c(-0.5, 0.5)
  ) +
  labs(
    x = "Hours after hyperosmotic stress",
    y = "RNA-seq log2FC(sorbitol/control)"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "none",
    plot.margin = margin(15, 15, 15, 15, "pt")
  )

# save outputs ------------------------------------------------------------

ggsave(
  "plots/temporal_trajectory_lineplot.pdf",
  plot = line_plot,
  width = 8,
  height = 6,
  units = "in"
)