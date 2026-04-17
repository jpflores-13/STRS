# Domain Analysis: Insulation Scores and Compartment Eigenvectors -----------
# Author:      JP Flores
# Date:        2026-04-16
# Project:     STRS
# Description: Visualizes changes in 3D chromatin domain organization between
#              untreated and 1h sorbitol-treated HEK293T eGFP-YAP1 cells.
#              All quantitative outputs (insulation scores, eigenvectors, saddle
#              matrices) are imported directly from cooltools, following the
#              analytical framework of Amat et al. 2019.
# Input:       data/processed/hic/insulation/*.tsv  (cooltools insulation)
#              data/processed/hic/eigenvectors/*.tsv (cooltools eigs-cis)
#              data/processed/hic/eigenvectors/*.saddledata.tsv (cooltools saddle)
# Output:      plots/*.pdf
#              data/processed/hic/domain_rds/*.rds
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

project_dir    <- "/work/users/j/p/jpflores/projects/STRS"
insulation_dir <- file.path(project_dir, "data/processed/hic/insulation")
eigvec_dir     <- file.path(project_dir, "data/processed/hic/eigenvectors")
output_dir     <- file.path(project_dir, "plots")
rds_dir        <- file.path(project_dir, "data/processed/hic/domain_rds")

# Autosomes only — must match view file order used in SLURM script
chroms <- c(paste0("chr", 1:7), "chr8", "chr9", "chr11", "chr10",
            paste0("chr", 12:18), "chr20", "chr19", "chr22", "chr21")

resolution    <- 10e3       # 10kb — resolution of insulation scores
window        <- 250e3      # insulation window (bp)
n_saddle_bins <- 50         # must match --n-bins in cooltools saddle
flank_bins    <- 50         # bins flanking boundary for avg profile

boundary_col  <- paste0("is_boundary_", as.integer(window))
score_col     <- paste0("log2_insulation_score_", as.integer(window))

# Colors consistent with existing STRS figures
col_control  <- "#4C72B0"
col_sorbitol <- "#DD8452"
col_AtoB     <- "#C44E52"
col_BtoA     <- "#4C72B0"
col_stable   <- "grey70"


# Libraries ---------------------------------------------------------------

library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(glue)
library(purrr)
library(readr)
library(stringr)
library(tidyr)


# Utility scripts ---------------------------------------------------------

source(file.path(project_dir, "scripts/utils/ggplot2_pgTheme.R"))


# Load data ---------------------------------------------------------------

## Insulation score TSVs — one row per 10kb genomic bin
ins_control  <- read_tsv(
  file.path(insulation_dir,
            glue("control_insulation_{as.integer(window)}bp.tsv")),
  show_col_types = FALSE
)

ins_sorbitol <- read_tsv(
  file.path(insulation_dir,
            glue("sorbitol_insulation_{as.integer(window)}bp.tsv")),
  show_col_types = FALSE
)

## Eigenvector TSVs from cooltools eigs-cis
## Columns: chrom, start, end, E1, E2, E3
eig_control  <- read_tsv(
  file.path(eigvec_dir, "control.cis.vecs.tsv"),
  show_col_types = FALSE
)

eig_sorbitol <- read_tsv(
  file.path(eigvec_dir, "sorbitol.cis.vecs.tsv"),
  show_col_types = FALSE
)

## Saddle data — extracted from cooltools .npz output
## saddledata is a (n_bins+2) x (n_bins+2) matrix including outlier bins;
## we trim the first and last row/col (outliers) to get the n_bins x n_bins core
saddle_control  <- read_tsv(
  file.path(eigvec_dir, "saddle_control.saddledata.tsv"),
  col_names = FALSE,
  show_col_types = FALSE
) |> as.matrix()

saddle_sorbitol <- read_tsv(
  file.path(eigvec_dir, "saddle_sorbitol.saddledata.tsv"),
  col_names = FALSE,
  show_col_types = FALSE
) |> as.matrix()

## Trim outlier bins (first and last row/col added by cooltools)
saddle_control  <- saddle_control[2:(nrow(saddle_control) - 1),
                                  2:(ncol(saddle_control) - 1)]
saddle_sorbitol <- saddle_sorbitol[2:(nrow(saddle_sorbitol) - 1),
                                   2:(ncol(saddle_sorbitol) - 1)]


# Wrangle data ------------------------------------------------------------

## --- Insulation: join conditions on bin coordinates ----------------------

ins_joined <- ins_control |>
  select(chrom, start, end,
         score_ctrl    = all_of(score_col),
         boundary_ctrl = all_of(boundary_col)) |>
  inner_join(
    ins_sorbitol |>
      select(chrom, start, end,
             score_sorb    = all_of(score_col),
             boundary_sorb = all_of(boundary_col)),
    by = c("chrom", "start", "end")
  ) |>
  filter(chrom %in% chroms,
         !is.na(score_ctrl),
         !is.na(score_sorb))

## Boundary counts per condition
boundary_counts <- tibble(
  condition    = c("Control", "Sorbitol 1h"),
  n_boundaries = c(
    sum(ins_control[[boundary_col]],  na.rm = TRUE),
    sum(ins_sorbitol[[boundary_col]], na.rm = TRUE)
  )
)

## --- Average insulation profile around control-defined boundaries --------
# For each control boundary, extract insulation scores ±flank_bins bins
# from both conditions, then average — replicates Amat et al. Figure 2B.

boundary_idx <- which(ins_joined$boundary_ctrl)

extract_window <- function(scores, idx, flank) {
  map(idx, \(i) {
    left  <- i - flank
    right <- i + flank
    if (left < 1 || right > length(scores)) return(NULL)
    scores[left:right]
  }) |>
    keep(\(x) !is.null(x)) |>
    do.call(rbind, args = _)
}

win_ctrl <- extract_window(ins_joined$score_ctrl, boundary_idx, flank_bins)
win_sorb <- extract_window(ins_joined$score_sorb, boundary_idx, flank_bins)

avg_profile <- tibble(
  bin      = seq(-flank_bins, flank_bins),
  control  = colMeans(win_ctrl, na.rm = TRUE),
  sorbitol = colMeans(win_sorb, na.rm = TRUE)
) |>
  pivot_longer(cols      = c(control, sorbitol),
               names_to  = "condition",
               values_to = "mean_insulation")

## --- Eigenvectors: join and classify A/B compartment transitions ---------

eig_joined <- eig_control |>
  select(chrom, start, end, E1_ctrl = E1) |>
  inner_join(
    eig_sorbitol |> select(chrom, start, end, E1_sorb = E1),
    by = c("chrom", "start", "end")
  ) |>
  filter(chrom %in% chroms,
         !is.na(E1_ctrl),
         !is.na(E1_sorb)) |>
  mutate(
    compartment_class = case_when(
      E1_ctrl > 0 & E1_sorb > 0  ~ "A → A (stable)",
      E1_ctrl < 0 & E1_sorb < 0  ~ "B → B (stable)",
      E1_ctrl > 0 & E1_sorb < 0  ~ "A → B",
      E1_ctrl < 0 & E1_sorb > 0  ~ "B → A",
      TRUE                        ~ NA_character_
    )
  )

compartment_summary <- eig_joined |>
  filter(!is.na(compartment_class)) |>
  count(compartment_class) |>
  mutate(pct = n / sum(n) * 100)

## --- Saddle matrices: log2 O/E ------------------------------------------
# cooltools saddle already outputs mean O/E per quantile bin pair.
# Take log2 for symmetric visualization around 0.

saddle_ctrl_log2 <- log2(saddle_control)
saddle_sorb_log2 <- log2(saddle_sorbitol)

tidy_saddle <- function(mat, condition_label) {
  n <- nrow(mat)
  expand_grid(row = seq_len(n), col = seq_len(n)) |>
    mutate(log2_oe   = as.vector(mat),
           condition = condition_label)
}

saddle_tidy <- bind_rows(
  tidy_saddle(saddle_ctrl_log2, "Control"),
  tidy_saddle(saddle_sorb_log2, "Sorbitol 1h")
) |>
  mutate(condition = factor(condition, levels = c("Control", "Sorbitol 1h")))


# Visualization -----------------------------------------------------------

## --- Plot 1: Boundary counts bar chart ----------------------------------
p_boundary_counts <- ggplot(boundary_counts,
                            aes(x = condition, y = n_boundaries,
                                fill = condition)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = n_boundaries), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("Control"     = col_control,
                                "Sorbitol 1h" = col_sorbitol)) +
  labs(x     = NULL,
       y     = "Number of TAD boundaries",
       title = glue("TAD boundary calls ({window / 1e3}kb insulation window)")) +
  pgTheme() +
  theme(legend.position = "none")

## --- Plot 2: Average insulation profile around control boundaries --------
# Each bin = 10kb; multiply by 10 to get kb on x-axis
p_insulation_profile <- ggplot(avg_profile,
                               aes(x = bin * 10, y = mean_insulation,
                                   color = condition)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  scale_color_manual(values = c("control"  = col_control,
                                 "sorbitol" = col_sorbitol),
                     labels  = c("control"  = "Control",
                                 "sorbitol" = "Sorbitol 1h")) +
  labs(x     = "Distance from boundary (kb)",
       y     = "Mean insulation score (log2)",
       color = NULL,
       title = "Average insulation at control-defined TAD boundaries") +
  pgTheme() +
  theme(legend.position = c(0.8, 0.85))

## --- Plot 3: Compartment eigenvector scatter ----------------------------
p_compartment_scatter <- eig_joined |>
  filter(!is.na(compartment_class)) |>
  ggplot(aes(x = E1_ctrl, y = E1_sorb, color = compartment_class)) +
  geom_point(alpha = 0.25, size = 0.6, stroke = 0) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "black") +
  scale_color_manual(
    values = c("A → A (stable)" = col_stable,
               "B → B (stable)" = col_stable,
               "A → B"          = col_AtoB,
               "B → A"          = col_BtoA)
  ) +
  labs(x     = "PC1 — Control",
       y     = "PC1 — Sorbitol 1h",
       color = "Compartment transition",
       title = "A/B compartment eigenvector comparison") +
  pgTheme() +
  theme(legend.position = "right")

## --- Plot 4: Saddle plots -----------------------------------------------
saddle_limit <- min(max(abs(saddle_tidy$log2_oe), na.rm = TRUE), 2)

p_saddle <- ggplot(saddle_tidy,
                   aes(x = col, y = row, fill = log2_oe)) +
  geom_tile() +
  facet_wrap(~ condition, nrow = 1) +
  scale_fill_gradientn(
    colors = rev(RColorBrewer::brewer.pal(11, "RdBu")),
    limits = c(-saddle_limit, saddle_limit),
    oob    = scales::squish,
    name   = "log2(O/E)"
  ) +
  scale_x_continuous(breaks = c(1, n_saddle_bins / 2, n_saddle_bins),
                     labels = c("B", "", "A")) +
  scale_y_reverse(breaks = c(1, n_saddle_bins / 2, n_saddle_bins),
                  labels = c("B", "", "A")) +
  labs(x     = "PC1 quantile",
       y     = "PC1 quantile",
       title = "Compartment strength (saddle plot)") +
  pgTheme() +
  theme(aspect.ratio = 1,
        panel.grid   = element_blank(),
        strip.text   = element_text(size = 10))


# Save outputs ------------------------------------------------------------

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir,    showWarnings = FALSE, recursive = TRUE)

## RDS objects
saveRDS(ins_joined,          file.path(rds_dir, "insulation_joined.rds"))
saveRDS(avg_profile,         file.path(rds_dir, "insulation_avg_profile.rds"))
saveRDS(eig_joined,          file.path(rds_dir, "eigenvectors_joined.rds"))
saveRDS(compartment_summary, file.path(rds_dir, "compartment_summary.rds"))
saveRDS(saddle_tidy,         file.path(rds_dir, "saddle_tidy.rds"))

## Plots
pdf(file.path(output_dir, "boundary_counts.pdf"),     width = 3.5, height = 4)
print(p_boundary_counts)
dev.off()

pdf(file.path(output_dir, "insulation_profile.pdf"),  width = 5,   height = 4)
print(p_insulation_profile)
dev.off()

pdf(file.path(output_dir, "compartment_scatter.pdf"), width = 5,   height = 5)
print(p_compartment_scatter)
dev.off()

pdf(file.path(output_dir, "saddle_plots.pdf"),        width = 8,   height = 4)
print(p_saddle)
dev.off()


# Session info ------------------------------------------------------------

sessionInfo()