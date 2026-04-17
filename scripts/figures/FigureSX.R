# FigureSX.R ------------------------------------------------
# Author:      JP Flores
# Date:        2026-04-17
# Project:     STRS
# Description: Single-row supplementary figure showing 3D chromatin domain
#              organization changes under hyperosmotic stress. Panel A: example
#              chr6 locus with PC1 eigenvector tracks + stacked Hi-C rectangles
#              (control on top, sorbitol below) + genes + genome label.
#              Panel B: insulation score profile. Panel C: differential saddle
#              plot (Sorbitol − Control).
# Input:       data/processed/hic/maps/*.hic
#              data/processed/hic/eigenvectors/*.cis.vecs.tsv
#              data/processed/hic/domain_rds/*.rds  (from domainAnalysis.R)
# Output:      figures/FigureSX.pdf
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

project_dir  <- "/work/users/j/p/jpflores/projects/STRS"
hic_dir      <- file.path(project_dir, "data/processed/hic/maps")
eigvec_dir   <- file.path(project_dir, "data/processed/hic/eigenvectors")
rds_dir      <- file.path(project_dir, "data/processed/hic/domain_rds")
output_dir   <- file.path(project_dir, "figures")

# .hic megaMap files
hic_control  <- file.path(hic_dir, "YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic")
hic_sorbitol <- file.path(hic_dir, "YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic")

# Example locus — chr6 1–6Mb (Amat et al. Fig 1D)
locus_chrom <- "chr6"
locus_start <- 1e6
locus_end   <- 6e6

resolution    <- 10e3   # Hi-C resolution (bp)
n_saddle_bins <- 50     # must match --n-bins in cooltools saddle
map_zrange    <- c(0, 100)

# Layout constants
map_width   <- 2.5      # width of each Hi-C rectangle panel
map_height  <- 0.75     # height of each Hi-C rectangle
eig_height  <- 0.35     # height of eigenvector signal track
gap         <- 0.05     # gap between stacked elements
x_locus     <- 0.5     # x position of locus panel

# Colors consistent with existing STRS figures
col_control  <- "#4C72B0"
col_sorbitol <- "#DD8452"


# Libraries ---------------------------------------------------------------

library(patchwork)
library(dplyr)
library(ggplot2)
library(plotgardener)
library(readr)
library(scales)
library(tidyr)


# Load data ---------------------------------------------------------------

## Pre-computed insulation profile and saddle data from domainAnalysis.R
avg_profile <- readRDS(file.path(rds_dir, "insulation_avg_profile.rds"))
saddle_tidy <- readRDS(file.path(rds_dir, "saddle_tidy.rds"))

## Eigenvector TSVs — filtered to locus for signal tracks
eig_ctrl_locus <- read_tsv(
  file.path(eigvec_dir, "control.cis.vecs.tsv"),
  show_col_types = FALSE
) |>
  filter(chrom == locus_chrom,
         start >= locus_start,
         end   <= locus_end)

eig_sorb_locus <- read_tsv(
  file.path(eigvec_dir, "sorbitol.cis.vecs.tsv"),
  show_col_types = FALSE
) |>
  filter(chrom == locus_chrom,
         start >= locus_start,
         end   <= locus_end)


# Wrangle data ------------------------------------------------------------

## Insulation profile line labels at bin -40 (-400kb)
line_labels <- avg_profile |>
  filter(bin == -45) |>
  mutate(label = case_when(condition == "control"  ~ "Control",
                           condition == "sorbitol" ~ "Sorbitol 1h"))

## Differential saddle: log2(sorbitol) - log2(control)
saddle_wide <- saddle_tidy |>
  pivot_wider(names_from = condition, values_from = log2_oe) |>
  mutate(diff = `Sorbitol 1h` - Control)

saddle_diff_limit <- min(max(abs(saddle_wide$diff), na.rm = TRUE), 1)

## Eigenvector range (finite values only)
eig_vals  <- c(eig_ctrl_locus$E1, eig_sorb_locus$E1)
eig_vals  <- eig_vals[is.finite(eig_vals)]
eig_range <- c(min(eig_vals), max(eig_vals))


# Visualization -----------------------------------------------------------

## --- ggplot: insulation score profile -----------------------------------
p_insulation <- ggplot(avg_profile,
                       aes(x = bin * 10, y = mean_insulation,
                           color = condition)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_text(
    data        = line_labels,
    aes(label   = label, color = condition),
    nudge_y     = ifelse(line_labels$condition == "control", -0.010, -0.006),
    hjust       = 0.5,
    size        = 1.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("control"  = col_control,
                                "sorbitol" = col_sorbitol)) +
  labs(x = "Distance from boundary (kb)",
       y = "Mean insulation score (log2)") +
  theme_bw() +
  theme(legend.position = "none",
        text            = element_text(size = 7))

## --- ggplots: saddle plots (control, sorbitol, differential) ------------

saddle_limit      <- min(max(abs(saddle_tidy$log2_oe), na.rm = TRUE), 2)
saddle_diff_limit <- min(max(abs(saddle_wide$diff),    na.rm = TRUE), 1)

## Shared saddle theme
saddle_theme <- theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid   = element_blank(),
        text         = element_text(size = 7),
        plot.title   = element_text(hjust = 0.5, size = 7))

## Control saddle
p_saddle_ctrl <- saddle_tidy |>
  filter(condition == "Control") |>
  ggplot(aes(x = col, y = row, fill = log2_oe)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = rev(RColorBrewer::brewer.pal(11, "RdBu")),
    limits = c(-saddle_limit, saddle_limit),
    oob    = scales::squish, name = "log2(O/E)",
    guide  = guide_colorbar(barwidth = 0.4, barheight = 4,
                            ticks.linewidth = 0.3,
                            frame.linewidth = 0.3)
  ) +
  scale_x_continuous(breaks = c(1, n_saddle_bins / 2, n_saddle_bins),
                     labels = c("B", "", "A")) +
  scale_y_reverse(breaks = c(1, n_saddle_bins / 2, n_saddle_bins),
                  labels = c("B", "", "A")) +
  labs(x = "PC1 quantile", y = "PC1 quantile", title = "Control") +
  saddle_theme +
  theme(legend.margin = margin(0, 0, 0, 2))

## Sorbitol saddle
p_saddle_sorb <- saddle_tidy |>
  filter(condition == "Sorbitol 1h") |>
  ggplot(aes(x = col, y = row, fill = log2_oe)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = rev(RColorBrewer::brewer.pal(11, "RdBu")),
    limits = c(-saddle_limit, saddle_limit),
    oob    = scales::squish, name = "log2(O/E)",
    guide  = guide_colorbar(barwidth = 0.4, barheight = 4,
                            ticks.linewidth = 0.3,
                            frame.linewidth = 0.3)
  ) +
  scale_x_continuous(breaks = c(1, n_saddle_bins / 2, n_saddle_bins),
                     labels = c("B", "", "A")) +
  scale_y_reverse(breaks = c(1, n_saddle_bins / 2, n_saddle_bins),
                  labels = c("B", "", "A")) +
  labs(x = "PC1 quantile", y = NULL, title = "Sorbitol 1h") +
  saddle_theme +
  theme(legend.margin = margin(0, 0, 0, 2))

## Differential saddle
p_saddle_diff <- saddle_wide |>
  select(row, col, log2_oe = diff) |>
  ggplot(aes(x = col, y = row, fill = log2_oe)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = rev(RColorBrewer::brewer.pal(11, "RdBu")),
    limits = c(-saddle_diff_limit, saddle_diff_limit),
    oob    = scales::squish,
    name = expression(Delta*log[2](O/E)),
    guide  = guide_colorbar(barwidth = 0.4, barheight = 4,
                            ticks.linewidth = 0.3,
                            frame.linewidth = 0.3)
  ) +
  scale_x_continuous(breaks = c(1, n_saddle_bins / 2, n_saddle_bins),
                     labels = c("B", "", "A")) +
  scale_y_reverse(breaks = c(1, n_saddle_bins / 2, n_saddle_bins),
                  labels = c("B", "", "A")) +
  labs(x = "PC1 quantile", y = "PC1 quantile",
       title = "Sorbitol − Control") +
  saddle_theme +
  theme(legend.margin = margin(0, 0, 0, 2))

## Compose: top row = control + sorbitol, bottom row = differential centered
p_saddle_composed <- (p_saddle_ctrl | p_saddle_sorb) /
  (plot_spacer() | p_saddle_diff | plot_spacer()) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


# Page layout -------------------------------------------------------------

## y positions — stacked from top
y_label    <- 0.15
y_eig_ctrl <- 0.35
y_map_ctrl <- y_eig_ctrl + eig_height + gap
y_eig_sorb <- y_map_ctrl + map_height + gap
y_map_sorb <- y_eig_sorb + eig_height + gap
y_genes    <- y_map_sorb + map_height + gap
y_genome   <- y_genes + 0.45 + gap

## x positions for ggplot panels — to the right of locus
x_ins      <- x_locus + map_width + 0.4
x_saddle   <- x_ins + 2.8

page_h     <- y_genome + 0.35
page_w     <- x_saddle + 4.2

## plotgardener locus parameters
params_locus <- pgParams(
  chrom      = locus_chrom,
  chromstart = as.integer(locus_start),
  chromend   = as.integer(locus_end),
  assembly   = "hg38",
  resolution = as.integer(resolution),
  norm       = "SCALE",
  zrange     = map_zrange,
  width      = map_width,
  height     = map_height,
  fontsize   = 6
)

pdf(file.path(output_dir, "FigureSX.pdf"),
    width = page_w, height = page_h)

pageCreate(width = page_w, height = page_h, default.units = "inches",
           showGuides = F)

## --- Panel A label ------------------------------------------------------
plotText("A", x = 0.15, y = y_label, fontsize = 11, fontface = "bold",
         default.units = "inches")

## --- Control eigenvector track ------------------------------------------
plotSignal(
  data      = eig_ctrl_locus |>
    select(chrom, start, end, score = E1) |>
    filter(!is.na(score)) |>
    mutate(end = end - 1L),
  params    = params_locus,
  x = x_locus, y = y_eig_ctrl,
  width = map_width, height = eig_height,
  range     = eig_range,
  fill      = col_control, linecolor = col_control,
  baseline  = TRUE, default.units = "inches"
)
plotText("untreated", x = x_locus, y = y_eig_ctrl,
         just = c("left", "top"), fontsize = 6, fontcolor = col_control,
         default.units = "inches")
plotText("Eigenvalue", x = x_locus - 0.1, y = y_eig_ctrl + eig_height / 2,
         rot = 90, just = "center", fontsize = 4,
         default.units = "inches")

## --- Control Hi-C rectangle ---------------------------------------------
hic_ctrl_plot <- plotHicRectangle(
  data   = hic_control,
  params = params_locus,
  x      = x_locus,
  y      = y_map_ctrl
)

annoHeatmapLegend(
  hic_ctrl_plot,
  x = x_locus + map_width + 0.05, y = y_map_ctrl,
  width = 0.05, height = map_height,
  orientation = "v", fontsize = 5, fontcolor = "black",
  just = c("left", "top"), default.units = "inches"
)

## --- Sorbitol eigenvector track -----------------------------------------
plotSignal(
  data      = eig_sorb_locus |>
    select(chrom, start, end, score = E1) |>
    filter(!is.na(score)) |>
    mutate(end = end - 1L),
  params    = params_locus,
  x = x_locus, y = y_eig_sorb,
  width = map_width, height = eig_height,
  range     = eig_range,
  fill      = col_sorbitol, linecolor = col_sorbitol,
  baseline  = TRUE, default.units = "inches"
)
plotText("200mM sorbitol", x = x_locus, y = y_eig_sorb,
         just = c("left", "top"), fontsize = 6, fontcolor = col_sorbitol,
         default.units = "inches")

## --- Sorbitol Hi-C rectangle --------------------------------------------
plotHicRectangle(
  data   = hic_sorbitol,
  params = params_locus,
  x      = x_locus,
  y      = y_map_sorb
)

## --- Genes --------------------------------------------------------------
plotGenes(
  params = params_locus,
  x      = x_locus,
  y      = y_genes,
  height = 0.4
)

## --- Genome label -------------------------------------------------------
plotGenomeLabel(
  params = params_locus,
  x      = x_locus,
  y      = y_genome,
  length = map_width
)

## --- Panel B: insulation profile ----------------------------------------
plotText("B", x = x_ins - 0.2, y = y_label, fontsize = 11,
         fontface = "bold", default.units = "inches")
plotGG(
  plot  = p_insulation,
  x     = x_ins, y = y_eig_ctrl,
  width = 2.6, height = page_h - y_eig_ctrl - 0.2,
  default.units = "inches"
)

## --- Panel C: composed saddle plots -------------------------------------
plotText("C", x = x_saddle - 0.2, y = y_label, fontsize = 11,
         fontface = "bold", default.units = "inches")
plotGG(
  plot  = p_saddle_composed,
  x     = x_saddle, y = y_eig_ctrl,
  width = 3.8, height = page_h - y_eig_ctrl - 0.2,
  default.units = "inches"
)

dev.off()

message("Saved: ", file.path(output_dir, "FigureSX.pdf"))


# Session info ------------------------------------------------------------

sessionInfo()