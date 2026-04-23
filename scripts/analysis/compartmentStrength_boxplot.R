# compartmentStrength_boxplot.R -------------------------------------------
# Author:      JP Flores
# Date:        2026-04-19
# Project:     STRS
# Description: Computes genome-wide compartment strength per replicate from
#              per-replicate cooltools saddle matrices and produces a boxplot
#              comparing eGFP-YAP (control and sorbitol) vs YAPdTAD (control
#              and sorbitol). Compartment strength = (AA + BB) / (AB + BA),
#              where AA/BB are mean corner values and AB/BA are off-diagonal
#              corner values of the saddle matrix. Faceted by genotype.
# Input:       data/processed/hic/compartment_strength/saddle_*.saddledata.tsv
# Output:      plots/compartmentStrength_boxplot.pdf
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

project_dir   <- "/work/users/j/p/jpflores/projects/STRS"
saddle_dir    <- file.path(project_dir, "data/processed/hic/compartment_strength")
output_dir    <- file.path(project_dir, "plots")
output_pdf    <- file.path(output_dir, "compartmentStrength_boxplot.pdf")

## Number of corner bins used to compute AA, BB, AB, BA averages.
## Using the outermost 3 bins on each side (top-left = BB, top-right = BA,
## bottom-left = AB, bottom-right = AA) matches the convention in
## Flyamer et al. 2017 and the figure style shown.
n_corner      <- 3L
n_saddle_bins <- 50L    # must match --n-bins used in cooltools saddle

## Colors -- consistent with existing STRS figures
col_yap_ctrl  <- "#89A8D0"
col_yap_sorb  <- "#EDB48E"
col_dtad_ctrl <- "#89A8D0"
col_dtad_sorb <- "#EDB48E"

## Condition factor levels (x-axis within each facet)
condition_levels <- c("control", "+ sorbitol")

## Genotype factor levels (facet order)
genotype_levels  <- c("eGFP-YAP", "YAPdTAD")


# Libraries ---------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
library(tibble)
library(tidyr)


# Load data ---------------------------------------------------------------

## Metadata table mapping file labels to group assignments
sample_meta <- tribble(
  ~label,               ~genotype,   ~condition,
  "eGFP-YAP_ctrl_rep1", "eGFP-YAP", "control",
  "eGFP-YAP_ctrl_rep2", "eGFP-YAP", "control",
  "eGFP-YAP_ctrl_rep3", "eGFP-YAP", "control",
  "eGFP-YAP_ctrl_rep4", "eGFP-YAP", "control",
  "eGFP-YAP_sorb_rep1", "eGFP-YAP", "+ sorbitol",
  "eGFP-YAP_sorb_rep2", "eGFP-YAP", "+ sorbitol",
  "eGFP-YAP_sorb_rep3", "eGFP-YAP", "+ sorbitol",
  "eGFP-YAP_sorb_rep4", "eGFP-YAP", "+ sorbitol",
  "YAPdTAD_ctrl_rep1",  "YAPdTAD",  "control",
  "YAPdTAD_sorb_rep1",  "YAPdTAD",  "+ sorbitol"
)

## Helper: compute compartment strength from a saddle matrix
## Saddle convention (cooltools): rows/cols ordered B -> A (low -> high PC1)
## So: top-left     = BB (low x low)
##     bottom-right = AA (high x high)
##     top-right    = BA
##     bottom-left  = AB
calc_compartment_strength <- function(mat, n_corner) {
  n   <- nrow(mat)
  idx <- seq_len(n_corner)
  
  BB  <- mean(mat[idx,          idx                      ], na.rm = TRUE)
  AA  <- mean(mat[(n - n_corner + 1):n,
                  (n - n_corner + 1):n                   ], na.rm = TRUE)
  BA  <- mean(mat[idx,          (n - n_corner + 1):n     ], na.rm = TRUE)
  AB  <- mean(mat[(n - n_corner + 1):n, idx              ], na.rm = TRUE)
  
  (AA + BB) / (AB + BA)
}

## Read all per-replicate saddle TSVs and compute compartment strength
strength_df <- sample_meta |>
  mutate(
    tsv_path  = file.path(saddle_dir,
                          paste0("saddle_", label, ".saddledata.tsv")),
    strength  = map_dbl(tsv_path, \(path) {
      mat <- read_tsv(path, col_names = FALSE,
                      show_col_types = FALSE) |>
        as.matrix()
      ## Trim cooltools outlier bins (first and last row/col)
      mat <- mat[2:(nrow(mat) - 1), 2:(ncol(mat) - 1)]
      calc_compartment_strength(mat, n_corner)
    }),
    condition = factor(condition, levels = condition_levels),
    genotype  = factor(genotype,  levels = genotype_levels)
  )


# Visualization -----------------------------------------------------------

fill_vals <- c(
  "control"    = col_yap_ctrl,
  "+ sorbitol" = col_yap_sorb
)

p_strength <- ggplot(strength_df,
                     aes(x = condition, y = strength, fill = condition)) +
  ## Boxplot only for groups with >1 rep (eGFP-YAP); outliers hidden
  geom_boxplot(
    data      = \(d) filter(d, genotype == "eGFP-YAP"),
    width     = 0.5,
    outlier.shape = NA,
    color     = "black",
    linewidth = 0.4
  ) +
  ## Individual points for all samples
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
  labs(
    x     = NULL,
    y     = "Compartment strength\n(AA + BB) / (AB + BA)",
    title = "Genome-wide compartment strength"
  ) +
  theme_bw() +
  theme(
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 8),
    axis.text.y        = element_text(size = 8),
    axis.title.y       = element_text(size = 9),
    plot.title         = element_text(size = 10, hjust = 0.5),
    strip.text         = element_text(size = 9, face = "italic"),
    strip.background   = element_rect(fill = "grey92", color = "grey70")
  )


# Save outputs ------------------------------------------------------------

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

pdf(output_pdf, width = 5, height = 4.5)
print(p_strength)
dev.off()

saveRDS(strength_df,
        file.path(project_dir,
                  "data/processed/hic/compartment_strength",
                  "compartmentStrength_perRep.rds"))


# Session info ------------------------------------------------------------

sessionInfo()