# FigureS2 ----------------------------------------------------------------
# Author:      JP Flores
# Date:        2026-04-20
# Project:     STRS
# Description: Supplementary figure showing 3D chromatin domain organization
#              changes under hyperosmotic stress. Panel A: example chr6 locus
#              with Hi-C rectangles (control on top, sorbitol below) + PC1
#              eigenvector tracks + genes + genome label. Panel B: saddle plots
#              (control, sorbitol, differential). Panel C: per-replicate
#              compartment strength boxplot faceted by genotype. Panel D: mean
#              insulation score profile across chromatin domain boundaries.
#              NOTE: "TAD" is reserved elsewhere in this manuscript for YAP1's
#              transcriptional activation domain (eGFP-YAP1dTAD), so insulation-
#              derived domains are called "chromatin domains" here, not TADs.
# Input:       - data/processed/hic/maps/*.hic
#              - data/processed/hic/eigenvectors/*.cis.vecs.tsv
#              - data/processed/hic/domain_rds/*.rds  (from domainAnalysis.R)
#              - data/processed/hic/compartment_strength/saddle_*.saddledata.tsv
# Output:      - figures/FigureS2.pdf
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

project_dir   <- "/work/users/j/p/jpflores/projects/STRS"
hic_dir       <- file.path(project_dir, "data/processed/hic/maps")
eigvec_dir    <- file.path(project_dir, "data/processed/hic/eigenvectors")
rds_dir       <- file.path(project_dir, "data/processed/hic/domain_rds")
saddle_dir    <- file.path(project_dir, "data/processed/hic/compartment_strength")
output_dir    <- file.path(project_dir, "figures")
output_pdf    <- file.path(output_dir, "FigureS2.pdf")

## .hic megaMap files
hic_control  <- file.path(hic_dir, "YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic")
hic_sorbitol <- file.path(hic_dir, "YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic")

## Example locus — chr6 1–6Mb
locus_chrom <- "chr6"
locus_start <- 1e6
locus_end   <- 6e6

resolution    <- 10e3   ## Hi-C resolution (bp)
n_saddle_bins <- 50L    ## must match --n-bins in cooltools saddle
map_zrange    <- c(0, 100)

## Compartment strength corner bins
n_corner <- 3L

## Colors consistent with existing STRS figures
col_control  <- "#4C72B0"
col_sorbitol <- "#DD8452"
col_ctrl_lt  <- "#89A8D0"   ## lighter hue for boxplot fills
col_sorb_lt  <- "#EDB48E"

## Condition / genotype factor levels
condition_levels <- c("control", "sorbitol")
genotype_levels  <- c("eGFP-YAP", "YAPdTAD")

## Layout constants (inches)
map_width  <- 2.5    ## width of Hi-C rectangle panel
map_height <- 0.75   ## height of each Hi-C rectangle
eig_height <- 0.35   ## height of eigenvector signal track
gap        <- 0.05   ## gap between stacked elements
x_locus    <- 0.5   ## x origin of panel A


# Libraries ---------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(plotgardener)
library(purrr)
library(RColorBrewer)
library(readr)
library(scales)
library(tibble)
library(tidyr)


# Load data ---------------------------------------------------------------

## Pre-computed insulation profile and saddle data from domainAnalysis.R
avg_profile <- readRDS(file.path(rds_dir, "insulation_avg_profile.rds"))
saddle_tidy <- readRDS(file.path(rds_dir, "saddle_tidy.rds"))

## Eigenvector TSVs — filtered to example locus for signal tracks
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

## Per-replicate saddle TSVs for compartment strength boxplot
sample_meta <- tribble(
  ~label,                ~genotype,   ~condition,
  "eGFP-YAP_ctrl_rep1", "eGFP-YAP",  "control",
  "eGFP-YAP_ctrl_rep2", "eGFP-YAP",  "control",
  "eGFP-YAP_ctrl_rep3", "eGFP-YAP",  "control",
  "eGFP-YAP_ctrl_rep4", "eGFP-YAP",  "control",
  "eGFP-YAP_sorb_rep1", "eGFP-YAP",  "sorbitol",
  "eGFP-YAP_sorb_rep2", "eGFP-YAP",  "sorbitol",
  "eGFP-YAP_sorb_rep3", "eGFP-YAP",  "sorbitol",
  "eGFP-YAP_sorb_rep4", "eGFP-YAP",  "sorbitol",
  "YAPdTAD_ctrl_rep1",  "YAPdTAD",   "control",
  "YAPdTAD_sorb_rep1",  "YAPdTAD",   "sorbitol"
)


# Wrangle data ------------------------------------------------------------

## --- Panel A: eigenvector range for shared y-axis -----------------------
eig_vals  <- c(eig_ctrl_locus$E1, eig_sorb_locus$E1)
eig_vals  <- eig_vals[is.finite(eig_vals)]
eig_range <- c(min(eig_vals), max(eig_vals))

## --- Panel B: saddle plots ----------------------------------------------
saddle_limit <- min(max(abs(saddle_tidy$log2_oe), na.rm = TRUE), 2)

saddle_wide <- saddle_tidy |>
  pivot_wider(names_from = condition, values_from = log2_oe) |>
  mutate(diff = `Sorbitol 1h` - Control)

saddle_diff_limit <- min(max(abs(saddle_wide$diff), na.rm = TRUE), 1)

## --- Panel B: compartment strength per replicate ------------------------
calc_compartment_strength <- function(mat, n_corner) {
  n   <- nrow(mat)
  idx <- seq_len(n_corner)
  BB  <- mean(mat[idx,                    idx                   ], na.rm = TRUE)
  AA  <- mean(mat[(n - n_corner + 1):n,   (n - n_corner + 1):n ], na.rm = TRUE)
  BA  <- mean(mat[idx,                    (n - n_corner + 1):n  ], na.rm = TRUE)
  AB  <- mean(mat[(n - n_corner + 1):n,   idx                   ], na.rm = TRUE)
  (AA + BB) / (AB + BA)
}

strength_df <- sample_meta |>
  mutate(
    tsv_path  = file.path(saddle_dir,
                          paste0("saddle_", label, ".saddledata.tsv")),
    strength  = map_dbl(tsv_path, \(path) {
      mat <- read_tsv(path, col_names = FALSE, show_col_types = FALSE) |>
        as.matrix()
      mat <- mat[2:(nrow(mat) - 1), 2:(ncol(mat) - 1)]
      calc_compartment_strength(mat, n_corner)
    }),
    condition = factor(condition, levels = condition_levels),
    genotype  = factor(genotype,  levels = genotype_levels)
  )

## --- Panel C: insulation score line labels ------------------------------
line_labels <- avg_profile |>
  filter(bin == -45) |>
  mutate(label = case_when(condition == "control"  ~ "Control",
                           condition == "sorbitol" ~ "Sorbitol 1h"))


# Visualization -----------------------------------------------------------

## Shared saddle theme
saddle_theme <- theme_bw() +
  theme(panel.grid    = element_blank(),
        text          = element_text(size = 5),
        plot.title    = element_text(hjust = 0.5, size = 5,
                                     margin = margin(0, 0, 1, 0)),
        plot.margin   = margin(1, 2, 1, 2),
        legend.margin = margin(0, 0, 0, -8),
        legend.text   = element_text(size = 4),
        legend.title  = element_text(size = 4))

## --- Panel B: control saddle --------------------------------------------
p_saddle_ctrl <- saddle_tidy |>
  filter(condition == "Control") |>
  ggplot(aes(x = col, y = row, fill = log2_oe)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = rev(brewer.pal(11, "RdBu")),
    limits = c(-saddle_limit, saddle_limit),
    oob    = scales::squish, name = "log2(O/E)",
    guide  = guide_colorbar(barwidth = 0.25, barheight = 1.5,
                            ticks.linewidth = 0.3, frame.linewidth = 0.3)
  ) +
  scale_x_continuous(breaks = c(1, n_saddle_bins / 2, n_saddle_bins),
                     labels = c("B", "", "A")) +
  scale_y_reverse(breaks = c(1, n_saddle_bins / 2, n_saddle_bins),
                  labels = c("B", "", "A")) +
  labs(x = "PC1 quantile", y = "PC1 quantile", title = "Control") +
  coord_fixed(ratio = 1) +
  saddle_theme

## --- Panel B: sorbitol saddle -------------------------------------------
p_saddle_sorb <- saddle_tidy |>
  filter(condition == "Sorbitol 1h") |>
  ggplot(aes(x = col, y = row, fill = log2_oe)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = rev(brewer.pal(11, "RdBu")),
    limits = c(-saddle_limit, saddle_limit),
    oob    = scales::squish, name = "log2(O/E)",
    guide  = guide_colorbar(barwidth = 0.25, barheight = 1.5,
                            ticks.linewidth = 0.3, frame.linewidth = 0.3)
  ) +
  scale_x_continuous(breaks = c(1, n_saddle_bins / 2, n_saddle_bins),
                     labels = c("B", "", "A")) +
  scale_y_reverse(breaks = c(1, n_saddle_bins / 2, n_saddle_bins),
                  labels = c("B", "", "A")) +
  labs(x = "PC1 quantile", y = "PC1 quantile", title = "Sorbitol 1h") +
  coord_fixed(ratio = 1) +
  saddle_theme

## --- Panel B: differential saddle ---------------------------------------
p_saddle_diff <- saddle_wide |>
  dplyr::select(row, col, log2_oe = diff) |>
  ggplot(aes(x = col, y = row, fill = log2_oe)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = rev(brewer.pal(11, "RdBu")),
    limits = c(-saddle_diff_limit, saddle_diff_limit),
    oob    = scales::squish,
    name   = "log2(O/E)",
    guide  = guide_colorbar(barwidth = 0.25, barheight = 1.5,
                            ticks.linewidth = 0.3, frame.linewidth = 0.3)
  ) +
  scale_x_continuous(breaks = c(1, n_saddle_bins / 2, n_saddle_bins),
                     labels = c("B", "", "A")) +
  scale_y_reverse(breaks = c(1, n_saddle_bins / 2, n_saddle_bins),
                  labels = c("B", "", "A")) +
  labs(x = "PC1 quantile", y = "PC1 quantile",
       title = "Sorbitol \u2212 Control") +
  coord_fixed(ratio = 1) +
  saddle_theme

## Compose saddle: top row ctrl + sorb, bottom row diff centered
## (kept for reference — individual placement used in page layout instead)
# p_saddles_composed <- (p_saddle_ctrl | p_saddle_sorb) /
#   (plot_spacer() | p_saddle_diff | plot_spacer()) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "right")

## --- Panel B: compartment strength boxplot ------------------------------
fill_vals <- c("control" = col_ctrl_lt, "sorbitol" = col_sorb_lt)

p_strength <- ggplot(strength_df,
                     aes(x = condition, y = strength, fill = condition)) +
  geom_boxplot(
    data          = \(d) filter(d, genotype == "eGFP-YAP"),
    width         = 0.5,
    outlier.shape = NA,
    color         = "black",
    linewidth     = 0.4
  ) +
  geom_point(
    aes(shape = genotype),
    size     = 2.5,
    stroke   = 0.5,
    color    = "black",
    position = position_jitter(width = 0.08, seed = 42L)
  ) +
  scale_fill_manual(values = fill_vals) +
  scale_shape_manual(values = c("eGFP-YAP" = 21, "YAPdTAD" = 23)) +
  facet_grid(. ~ genotype, scales = "free_x") +
  labs(x = NULL,
       y = "Compartment strength\n(AA + BB) / (AB + BA)",
       title = "Genome-wide compartment strength") +
  theme_bw() +
  theme(
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 7),
    axis.text.y        = element_text(size = 7),
    axis.title.y       = element_text(size = 8),
    plot.title         = element_text(size = 8, hjust = 0.5),
    strip.text         = element_text(size = 8, face = "italic"),
    strip.background   = element_rect(fill = "grey92", color = "grey70")
  )

## --- Panel C: insulation score profile ----------------------------------
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
  labs(x = "Distance from chromatin domain boundary (kb)",
       y = "Mean insulation score (log2)") +
  theme_bw() +
  theme(legend.position = "none",
        text            = element_text(size = 7))


# Page layout -------------------------------------------------------------

## Panel A y positions
y_label    <- 0.15
y_map_ctrl <- 0.35
y_eig_ctrl <- y_map_ctrl + map_height + gap
y_map_sorb <- y_eig_ctrl + eig_height + gap
y_eig_sorb <- y_map_sorb + map_height + gap
y_genes    <- y_eig_sorb + eig_height + gap
y_genome   <- y_genes + 0.45 + gap

## Panel A total content height — C and D bottom aligns with genome label
panel_a_height <- y_genome + 0.3 - y_map_ctrl

## Panel B: saddle plots stacked, starting just below the B label
saddle_size   <- 1.3
saddle_gap    <- -0.08
y_saddle_ctrl <- y_label + 0.05
y_saddle_sorb <- y_saddle_ctrl + saddle_size + saddle_gap
y_saddle_diff <- y_saddle_sorb + saddle_size + saddle_gap

## Page height: tallest of panel A bottom or bottom of third saddle
page_h <- max(y_genome + 0.35,
              y_saddle_diff + saddle_size + 0.2)

## Panel C: compartment strength boxplot
b_box_w  <- 2.4

## Panel D: insulation profile
c_width  <- 2.6

## Panel x origins
x_panel_b  <- x_locus + map_width + 0.5
x_panel_c  <- x_panel_b + saddle_size + 0.5
x_panel_d  <- x_panel_c + b_box_w + 0.3

page_w     <- x_panel_d + c_width + 0.3

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


# Save outputs ------------------------------------------------------------

pdf(output_pdf, width = page_w, height = page_h)

pageCreate(width = page_w, height = page_h, default.units = "inches",
           showGuides = FALSE)

## ---- Panel A label -----------------------------------------------------
plotText("A", x = 0.15, y = y_label, fontsize = 11, fontface = "bold",
         default.units = "inches")

## ---- Panel A: control Hi-C rectangle -----------------------------------
hic_ctrl_plot <- plotHicRectangle(
  data   = hic_control,
  params = params_locus,
  x      = x_locus,
  y      = y_map_ctrl
)
plotText("untreated", x = x_locus, y = y_map_ctrl,
         just = c("left", "top"), fontsize = 6, fontcolor = col_control,
         default.units = "inches")
annoHeatmapLegend(
  hic_ctrl_plot,
  x = x_locus + map_width + 0.05, y = y_map_ctrl,
  width = 0.05, height = map_height,
  orientation = "v", fontsize = 5, fontcolor = "black",
  just = c("left", "top"), default.units = "inches"
)

## ---- Panel A: control eigenvector track --------------------------------
plotSignal(
  data      = eig_ctrl_locus |>
    dplyr::select(chrom, start, end, score = E1) |>
    filter(!is.na(score)) |>
    mutate(end = end - 1L),
  params    = params_locus,
  x = x_locus, y = y_eig_ctrl,
  width = map_width, height = eig_height,
  range     = eig_range,
  fill      = col_control, linecolor = col_control,
  baseline  = TRUE, default.units = "inches"
)
plotText("Eigenvalue", x = x_locus - 0.05, y = y_eig_ctrl + eig_height / 2,
         rot = 90, just = "center", fontsize = 4, default.units = "inches")

## ---- Panel A: sorbitol Hi-C rectangle ----------------------------------
plotHicRectangle(
  data   = hic_sorbitol,
  params = params_locus,
  x      = x_locus,
  y      = y_map_sorb
)
plotText("200mM sorbitol", x = x_locus, y = y_map_sorb,
         just = c("left", "top"), fontsize = 6, fontcolor = col_sorbitol,
         default.units = "inches")

## ---- Panel A: sorbitol eigenvector track -------------------------------
plotSignal(
  data      = eig_sorb_locus |>
    dplyr::select(chrom, start, end, score = E1) |>
    filter(!is.na(score)) |>
    mutate(end = end - 1L),
  params    = params_locus,
  x = x_locus, y = y_eig_sorb,
  width = map_width, height = eig_height,
  range     = eig_range,
  fill      = col_sorbitol, linecolor = col_sorbitol,
  baseline  = TRUE, default.units = "inches"
)
plotText("Eigenvalue", x = x_locus - 0.05, y = y_eig_sorb + eig_height / 2,
         rot = 90, just = "center", fontsize = 4, default.units = "inches")

## ---- Panel A: genes + genome label -------------------------------------
plotGenes(
  params = params_locus,
  x      = x_locus,
  y      = y_genes,
  height = 0.4
)
plotGenomeLabel(
  params = params_locus,
  x      = x_locus,
  y      = y_genome,
  length = map_width
)

## ---- Panel B label -----------------------------------------------------
plotText("B", x = x_panel_b - 0.2, y = y_label, fontsize = 11,
         fontface = "bold", default.units = "inches")

## ---- Panel B: control saddle (aligned with control Hi-C in A) ----------
plotGG(
  plot   = p_saddle_ctrl,
  x      = x_panel_b,
  y      = y_saddle_ctrl,
  width  = saddle_size,
  height = saddle_size,
  default.units = "inches"
)

## ---- Panel B: sorbitol saddle (aligned with sorbitol Hi-C in A) --------
plotGG(
  plot   = p_saddle_sorb,
  x      = x_panel_b,
  y      = y_saddle_sorb,
  width  = saddle_size,
  height = saddle_size,
  default.units = "inches"
)

## ---- Panel B: differential saddle (below sorbitol block) ---------------
plotGG(
  plot   = p_saddle_diff,
  x      = x_panel_b,
  y      = y_saddle_diff,
  width  = saddle_size,
  height = saddle_size,
  default.units = "inches"
)

## ---- Panel C label -----------------------------------------------------
plotText("C", x = x_panel_c - 0.2, y = y_label, fontsize = 11,
         fontface = "bold", default.units = "inches")

## ---- Panel C: compartment strength boxplot -----------------------------
plotGG(
  plot   = p_strength,
  x      = x_panel_c,
  y      = y_map_ctrl,
  width  = b_box_w,
  height = panel_a_height,
  default.units = "inches"
)

## ---- Panel D label -----------------------------------------------------
plotText("D", x = x_panel_d - 0.2, y = y_label, fontsize = 11,
         fontface = "bold", default.units = "inches")

## ---- Panel D: insulation score profile ---------------------------------
plotGG(
  plot   = p_insulation,
  x      = x_panel_d,
  y      = y_map_ctrl,
  width  = c_width,
  height = panel_a_height,
  default.units = "inches"
)

dev.off()

message("Saved: ", output_pdf)


# Session info ------------------------------------------------------------

sessionInfo()