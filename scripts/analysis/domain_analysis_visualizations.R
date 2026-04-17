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

library(GenomicRanges)
library(IRanges)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(glue)
library(purrr)
library(readr)
library(scales)
library(stringr)
library(tidyr)



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

## Saddle data TSVs extracted from cooltools .npz output
## saddledata is (n_bins+2) x (n_bins+2) — includes outlier bins on each edge
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

## Trim outlier bins (first and last row/col) to get n_saddle_bins x n_saddle_bins core
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
# cooltools saddle outputs mean O/E per quantile bin pair.
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


# Analysis: boundary overlap using GenomicRanges --------------------------

## Build GRanges from control and sorbitol boundary calls.
## findOverlaps with ±1 bin (±10kb) slop classifies boundaries as:
##   shared    — present in both conditions
##   lost      — control only
##   gained    — sorbitol only

ctrl_boundaries <- ins_control |>
  filter(.data[[boundary_col]] == TRUE, !is.na(.data[[score_col]])) |>
  filter(chrom %in% chroms)

sorb_boundaries <- ins_sorbitol |>
  filter(.data[[boundary_col]] == TRUE, !is.na(.data[[score_col]])) |>
  filter(chrom %in% chroms)

gr_ctrl <- GRanges(
  seqnames = ctrl_boundaries$chrom,
  ranges   = IRanges(start = ctrl_boundaries$start,
                     end   = ctrl_boundaries$end)
)

gr_sorb <- GRanges(
  seqnames = sorb_boundaries$chrom,
  ranges   = IRanges(start = sorb_boundaries$start,
                     end   = sorb_boundaries$end)
)

slop_bp    <- as.integer(resolution)
shared_hits <- findOverlaps(gr_ctrl + slop_bp, gr_sorb + slop_bp,
                            ignore.strand = TRUE)

n_shared <- length(unique(queryHits(shared_hits)))
n_lost   <- length(gr_ctrl) - n_shared
n_gained <- length(gr_sorb) - length(unique(subjectHits(shared_hits)))

boundary_overlap <- tibble(
  category = factor(
    c("Lost\n(control only)", "Shared", "Gained\n(sorbitol only)"),
    levels = c("Lost\n(control only)", "Shared", "Gained\n(sorbitol only)")
  ),
  n        = c(n_lost, n_shared, n_gained),
  fill_col = c(col_sorbitol, "grey60", col_control)
)


# Visualization -----------------------------------------------------------

## --- Plot 1: Average insulation profile around control boundaries --------
# Labels placed at x = -400kb where the lines are clearly vertically separated
line_labels <- avg_profile |>
  filter(bin == -40) |>
  mutate(label = case_when(condition == "control"  ~ "Control",
                           condition == "sorbitol" ~ "Sorbitol 1h"))

p_insulation_profile <- ggplot(avg_profile,
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
    size        = 3.5,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("control"  = col_control,
                                "sorbitol" = col_sorbitol)) +
  labs(x     = "Distance from boundary (kb)",
       y     = "Mean insulation score (log2)") +
  theme_bw() +
  theme(legend.position = "none")

## --- Plot 2: TAD boundary overlap bar chart ------------------------------
p_boundary_overlap <- ggplot(boundary_overlap,
                             aes(x = category, y = n, fill = fill_col)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = scales::comma(n)), vjust = -0.4, size = 3.5) +
  scale_fill_identity() +
  labs(x = NULL,
       y = "Number of TAD boundaries") +
  
  theme(panel.grid       = element_blank(),
        panel.background = element_blank(),
        plot.background  = element_blank())

## --- Plot 3: Compartment eigenvector scatter ----------------------------
compartment_levels <- c("B → B (stable)", "A → A (stable)", "B → A", "A → B")

p_compartment_scatter <- eig_joined |>
  filter(!is.na(compartment_class)) |>
  mutate(compartment_class = factor(compartment_class,
                                    levels = compartment_levels)) |>
  arrange(compartment_class) |>
  ggplot(aes(x = E1_ctrl, y = E1_sorb, color = compartment_class)) +
  geom_point(alpha = 0.25, size = 0.6, stroke = 0) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "black") +
  scale_color_manual(
    values = c("A → A (stable)" = col_stable,
               "B → B (stable)" = col_stable,
               "A → B"          = col_AtoB,
               "B → A"          = col_BtoA),
    labels = c("A → A (stable)" = "A → A (stable)",
               "B → B (stable)" = "B → B (stable)",
               "A → B"          = "A → B (switching)",
               "B → A"          = "B → A (switching)")
  ) +
  guides(color = guide_legend(
    override.aes = list(alpha = 1, size = 3),
    title        = "Compartment"
  )) +
  labs(x = "PC1 (Control)",
       y = "PC1 (Sorbitol 1h)") +
  
  theme(legend.position  = "right",
        panel.grid       = element_blank(),
        panel.background = element_blank(),
        plot.background  = element_blank())

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
  scale_x_continuous(
    breaks = c(1, n_saddle_bins / 2, n_saddle_bins),
    labels = c("B", "", "A")
  ) +
  scale_y_reverse(
    breaks = c(1, n_saddle_bins / 2, n_saddle_bins),
    labels = c("B", "", "A")
  ) +
  labs(x     = "Genomic bins ranked by PC1",
       y     = "Genomic bins ranked by PC1",
       title = "Compartment strength") +
  
  theme(aspect.ratio     = 1,
        panel.grid       = element_blank(),
        panel.background = element_blank(),
        plot.background  = element_blank(),
        strip.text       = element_text(size = 10))


# Save outputs ------------------------------------------------------------

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir,    showWarnings = FALSE, recursive = TRUE)

## RDS objects
saveRDS(ins_joined,          file.path(rds_dir, "insulation_joined.rds"))
saveRDS(avg_profile,         file.path(rds_dir, "insulation_avg_profile.rds"))
saveRDS(boundary_overlap,    file.path(rds_dir, "boundary_overlap.rds"))
saveRDS(eig_joined,          file.path(rds_dir, "eigenvectors_joined.rds"))
saveRDS(compartment_summary, file.path(rds_dir, "compartment_summary.rds"))
saveRDS(saddle_tidy,         file.path(rds_dir, "saddle_tidy.rds"))

## Plots
pdf(file.path(output_dir, "insulation_profile.pdf"),  width = 5,   height = 4)
print(p_insulation_profile)
dev.off()

pdf(file.path(output_dir, "boundary_overlap.pdf"),    width = 4,   height = 4)
print(p_boundary_overlap)
dev.off()

pdf(file.path(output_dir, "compartment_scatter.pdf"), width = 5,   height = 5)
print(p_compartment_scatter)
dev.off()

pdf(file.path(output_dir, "saddle_plots.pdf"),        width = 8,   height = 4)
print(p_saddle)
dev.off()


# Session info ------------------------------------------------------------

sessionInfo()