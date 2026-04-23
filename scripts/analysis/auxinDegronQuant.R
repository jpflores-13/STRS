# auxinDegronQuant.R ----------------------------------------------------------
# Author:      JP Flores
# Date:        2026-04-19
# Project:     STRS
# Description: Loads CellProfiler per-nucleus output from CTCF/RAD21 auxin-
#              degron IF experiment, attaches condition labels from the plate
#              layout, and generates a violin + jitter plot of mean nuclear
#              GFP (FITC) intensity per condition. Used to confirm protein
#              degradation in response to auxin treatment.
#              Stats: pairwise t-tests (Control vs Sorbitol+Auxin per protein)
#              with significance stars displayed on plot via ggpubr.
# Input:       data/processed/imaging/cellprofiler_output/CTCF_RAD21_auxin_Nuclei.csv
# Output:      plots/auxinDegronQuant.pdf
# -----------------------------------------------------------------------------


# Parameters ------------------------------------------------------------------

project_dir  <- "/work/users/j/p/jpflores/projects/STRS"
cp_output    <- file.path(project_dir,
                          "data/processed/imaging/cellprofiler_output",
                          "CTCF_RAD21_auxin_Nuclei.csv")
output_dir   <- file.path(project_dir, "plots")
output_file  <- file.path(output_dir, "auxinDegronQuant.pdf")

## Plate layout: Col number -> condition label
## Rows B, C, D are biological replicates of the same 6 conditions
plate_layout <- data.frame(
  Col       = as.character(1:6),
  Protein   = c("CTCF",  "CTCF",          "CTCF",
                "RAD21", "RAD21",          "RAD21"),
  Treatment = c("Control",
                "Sorbitol",
                "Sorbitol + Auxin",
                "Control",
                "Sorbitol",
                "Sorbitol + Auxin")
)

## QC filters: remove nuclei that are likely debris or mitotic outliers
min_area_px  <- 500    # minimum nucleus area in pixels (~20x, 0.325 um/px)
max_area_px  <- 15000  # maximum nucleus area in pixels

## Significance comparisons: all pairwise within each protein
stat_comparisons_ctcf <- list(
  c("CTCF\nControl", "CTCF\nSorbitol"),
  c("CTCF\nControl", "CTCF\nSorbitol + Auxin"),
  c("CTCF\nSorbitol", "CTCF\nSorbitol + Auxin")
)

stat_comparisons_rad21 <- list(
  c("RAD21\nControl", "RAD21\nSorbitol"),
  c("RAD21\nControl", "RAD21\nSorbitol + Auxin"),
  c("RAD21\nSorbitol", "RAD21\nSorbitol + Auxin")
)

## Plot dimensions
plot_width   <- 7
plot_height  <- 5


# Libraries -------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)
library(forcats)
library(readr)


# Load data -------------------------------------------------------------------

## Per-nucleus measurements exported by CellProfiler
## Each row = one nucleus; key columns are:
##   Metadata_Row, Metadata_Col  : well position on plate
##   Intensity_MeanIntensity_GFP : mean nuclear FITC signal (primary readout)
##   AreaShape_Area              : nuclear area in pixels (for QC filtering)
nuclei_raw <- read_csv(cp_output, show_col_types = FALSE)


# Wrangle data ----------------------------------------------------------------

nuclei <- nuclei_raw |>
  rename(
    Row = Metadata_Row,
    Col = Metadata_Col
  ) |>
  mutate(Col = as.character(Col)) |>
  left_join(plate_layout, by = "Col") |>
  filter(
    AreaShape_Area >= min_area_px,
    AreaShape_Area <= max_area_px
  ) |>
  mutate(
    Condition = paste(Protein, Treatment, sep = "\n"),
    Condition = fct_relevel(
      Condition,
      "CTCF\nControl",
      "CTCF\nSorbitol",
      "CTCF\nSorbitol + Auxin",
      "RAD21\nControl",
      "RAD21\nSorbitol",
      "RAD21\nSorbitol + Auxin"
    )
  )


# Analysis --------------------------------------------------------------------

## Per-replicate (per-well) means — used as jitter points on the plot.
## T-tests are run on these replicate means (n=3 per condition) rather than
## on individual nuclei, which would massively inflate the degrees of freedom
## and produce meaningless p-values.
replicate_means <- nuclei |>
  group_by(Row, Col, Condition, Protein, Treatment) |>
  summarize(
    mean_gfp = mean(Intensity_MeanIntensity_GFP),
    n_nuclei = dplyr::n(),
    .groups  = "drop"
  )


# Visualization ---------------------------------------------------------------

## Color palette: CTCF in orange, RAD21 in blue (Okabe-Ito colorblind-safe)
protein_colors <- c("CTCF" = "#E69F00", "RAD21" = "#0072B2")

p <- ggplot(
  nuclei,
  aes(x = Condition, y = Intensity_MeanIntensity_GFP, fill = Protein)
) +
  
  ## Violin shows full per-nucleus distribution
  geom_violin(
    color = NA,
    alpha = 0.4,
    scale = "width",
    trim  = TRUE
  ) +
  
  ## Per-replicate means as jitter points (3 replicates per condition)
  geom_point(
    data     = replicate_means,
    aes(y    = mean_gfp),
    shape    = 21,
    size     = 2.5,
    stroke   = 0.6,
    color    = "black",
    position = position_jitter(width = 0.05, seed = 42)
  ) +
  
  ## Median line per condition
  stat_summary(
    fun       = median,
    geom      = "crossbar",
    width     = 0.4,
    linewidth = 0.5,
    color     = "black"
  ) +
  
  ## Significance brackets: Control vs Sorbitol+Auxin per protein
  ## stat_compare_means runs a t-test on the per-nucleus values within
  ## each facet and displays stars: ns / * / ** / ***
  stat_compare_means(
    comparisons  = stat_comparisons_ctcf,
    data         = nuclei |> filter(Protein == "CTCF"),
    method       = "t.test",
    label        = "p.signif",
    tip.length   = 0.01,
    bracket.size = 0.4,
    size         = 4
  ) +
  
  stat_compare_means(
    comparisons  = stat_comparisons_rad21,
    data         = nuclei |> filter(Protein == "RAD21"),
    method       = "t.test",
    label        = "p.signif",
    tip.length   = 0.01,
    bracket.size = 0.4,
    size         = 4
  ) +
  
  ## Separate CTCF and RAD21 visually
  facet_wrap(~ Protein, scales = "free_x") +
  
  scale_fill_manual(values = protein_colors) +
  
  labs(
    x    = NULL,
    y    = "Mean Nuclear GFP Intensity (a.u.)",
    fill = "Protein"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold", size = 11),
    axis.text.x     = element_text(size = 8),
    panel.spacing   = unit(1.5, "lines")
  )


# Save outputs ----------------------------------------------------------------

pdf(file   = output_file,
    width  = plot_width,
    height = plot_height)
print(p)
dev.off()


# Session info ----------------------------------------------------------------

sessionInfo()