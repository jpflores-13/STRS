# ##############################################################################
# filename:    eisa_auxinRAD21_barplot.R
# author:      JP Flores
# date:        2026-05-29
# project:     STRS
# description: Grouped barplot showing mean delta-intron (nascent RNA log2FC
#              from eisaR +/- SEM) at genes whose promoters overlap sorbitol-
#              induced (gained) vs pre-existing (lost) loop anchors, comparing
#              sorbitol 6h vs untreated and sorbitol+auxin 6h vs untreated.
#              Mirrors auxinRAD21_gainedLoops_barplot.R but uses delta-intron
#              (nascent RNA from eisaR) rather than total RNA log2FC. Strip
#              labels below x-axis. Paired Wilcoxon signed-rank tests per loop
#              group.
# input:       data/processed/rna/AuxinRAD21/output/EISA/eisa_auxinRAD21_results.rds
# output:      plots/eisa_auxinRAD21_barplot.pdf
#              data/processed/rna/AuxinRAD21/output/EISA/gg_eisa_auxinRAD21_barplot.rds
# ##############################################################################


# Libraries -------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)


# Utility scripts -------------------------------------------------------------

source("scripts/utils/ggplot2_pgTheme.R")


# Parameters ------------------------------------------------------------------

project_dir <- "/work/users/j/p/jpflores/projects/STRS"
plots_dir   <- file.path(project_dir, "plots")
eisa_dir    <- file.path(project_dir,
                         "data/processed/rna/AuxinRAD21/output/EISA/per_replicate")

eisa_rds   <- file.path(eisa_dir, "eisa_auxinRAD21_results.rds")
gg_rds_out <- file.path(eisa_dir, "gg_eisa_auxinRAD21_barplot.rds")

color_sorbitol <- "#F8766D"
color_auxin    <- "#009E73"   # Okabe-Ito green

plot_width  <- 6
plot_height <- 4

output_pdf <- file.path(plots_dir, "eisa_auxinRAD21_barplot.pdf")


# Load data -------------------------------------------------------------------

eisa_df <- readRDS(eisa_rds)

message("EISA AuxinRAD21 results: ", nrow(eisa_df), " genes")


# Wrangle data ----------------------------------------------------------------

gained_df <- eisa_df |> dplyr::filter(loop_category == "gained")
lost_df   <- eisa_df |> dplyr::filter(loop_category == "lost")

bar_raw_df <- dplyr::bind_rows(
  data.frame(
    gene_id    = gained_df$gene_id,
    log2FC     = gained_df$lfc_intron_sorb,
    condition  = "Sorbitol",
    loop_group = "Sorbitol-Induced Loops"
  ),
  data.frame(
    gene_id    = gained_df$gene_id,
    log2FC     = gained_df$lfc_intron_auxin,
    condition  = "Sorbitol + Auxin",
    loop_group = "Sorbitol-Induced Loops"
  ),
  data.frame(
    gene_id    = lost_df$gene_id,
    log2FC     = lost_df$lfc_intron_sorb,
    condition  = "Sorbitol",
    loop_group = "Pre-Existing Loops"
  ),
  data.frame(
    gene_id    = lost_df$gene_id,
    log2FC     = lost_df$lfc_intron_auxin,
    condition  = "Sorbitol + Auxin",
    loop_group = "Pre-Existing Loops"
  )
) |>
  dplyr::filter(!is.na(log2FC))

bar_raw_df$condition  <- factor(bar_raw_df$condition,
                                levels = c("Sorbitol", "Sorbitol + Auxin"))
bar_raw_df$loop_group <- factor(bar_raw_df$loop_group,
                                levels = c("Sorbitol-Induced Loops",
                                           "Pre-Existing Loops"))


# Statistics ------------------------------------------------------------------

format_pval <- function(p) {
  if (p < 0.001) "p < 0.001" else paste0("p = ", round(p, 3))
}

run_paired_wilcox <- function(df, loop_grp) {
  sub    <- df |> dplyr::filter(loop_group == loop_grp)
  sorb   <- sub |> dplyr::filter(condition == "Sorbitol")
  auxin  <- sub |> dplyr::filter(condition == "Sorbitol + Auxin")
  common <- intersect(sorb$gene_id, auxin$gene_id)
  lfc_s  <- sorb$log2FC[match(common,  sorb$gene_id)]
  lfc_a  <- auxin$log2FC[match(common, auxin$gene_id)]
  valid  <- !is.na(lfc_s) & !is.na(lfc_a)
  wt     <- wilcox.test(lfc_s[valid], lfc_a[valid],
                        paired = TRUE, alternative = "two.sided")
  list(wt = wt, n = sum(valid))
}

wt_gained <- run_paired_wilcox(bar_raw_df, "Sorbitol-Induced Loops")
wt_lost   <- run_paired_wilcox(bar_raw_df, "Pre-Existing Loops")

message("Sorbitol-Induced Loops - paired Wilcoxon (delta-intron):")
message("  n = ", wt_gained$n, " genes | p = ", wt_gained$wt$p.value)
message("Pre-Existing Loops - paired Wilcoxon (delta-intron):")
message("  n = ", wt_lost$n,   " genes | p = ", wt_lost$wt$p.value)


# Summary statistics ----------------------------------------------------------

bar_summary <- bar_raw_df |>
  dplyr::group_by(loop_group, condition) |>
  dplyr::summarise(
    mean_lfc = mean(log2FC, na.rm = TRUE),
    se_lfc   = sd(log2FC,   na.rm = TRUE) / sqrt(dplyr::n()),
    n        = dplyr::n(),
    .groups  = "drop"
  )

bracket_offset <- diff(range(c(bar_summary$mean_lfc - bar_summary$se_lfc,
                               bar_summary$mean_lfc + bar_summary$se_lfc),
                             na.rm = TRUE)) * 0.06

sig_df <- bar_summary |>
  dplyr::group_by(loop_group) |>
  dplyr::summarise(
    ## Always place bracket above y = 0 — use whichever is higher:
    ## the tallest bar top or zero itself, then add offset
    y_bracket = max(max(mean_lfc + se_lfc, na.rm = TRUE), 0) + bracket_offset,
    .groups   = "drop"
  ) |>
  dplyr::mutate(
    pval_label = c(format_pval(wt_gained$wt$p.value),
                   format_pval(wt_lost$wt$p.value)),
    n_label    = c(paste0("n = ", wt_gained$n, " genes"),
                   paste0("n = ", wt_lost$n,   " genes")),
    y_text     = y_bracket + bracket_offset * 0.8,
    y_n        = y_bracket + bracket_offset * 2.2
  )


# Visualization ---------------------------------------------------------------

p_bar <- ggplot(bar_summary,
                aes(x = condition, y = mean_lfc, fill = condition)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray40", linewidth = 0.3) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = mean_lfc - se_lfc,
        ymax = mean_lfc + se_lfc),
    width     = 0.22,
    linewidth = 0.4
  ) +
  geom_segment(data = sig_df,
               aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket),
               inherit.aes = FALSE, linewidth = 0.3, color = "black") +
  geom_segment(data = sig_df,
               aes(x = 1, xend = 1,
                   y = y_bracket, yend = y_bracket - bracket_offset * 0.4),
               inherit.aes = FALSE, linewidth = 0.3, color = "black") +
  geom_segment(data = sig_df,
               aes(x = 2, xend = 2,
                   y = y_bracket, yend = y_bracket - bracket_offset * 0.4),
               inherit.aes = FALSE, linewidth = 0.3, color = "black") +
  geom_text(data = sig_df,
            aes(x = 1.5, y = y_text, label = pval_label),
            inherit.aes = FALSE, size = 6 / .pt, color = "black") +
  geom_text(data = sig_df,
            aes(x = 1.5, y = y_n, label = n_label),
            inherit.aes = FALSE, size = 5 / .pt, color = "gray40") +
  scale_fill_manual(values = c(
    "Sorbitol"         = color_sorbitol,
    "Sorbitol + Auxin" = color_auxin
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.30))) +
  facet_wrap(~loop_group, strip.position = "bottom") +
  labs(
    x = NULL,
    y = expression(Mean~Delta*"intron " * log[2]~FC~(treated / untreated))
  ) +
  theme_classic() +
  theme(
    strip.placement    = "outside",
    strip.background   = element_rect(fill = NA, color = "black",
                                      linewidth = 0.3),
    strip.text.x       = element_text(size = 7, face = "bold"),
    legend.position    = "none",
    axis.text          = element_text(size = 7, color = "black"),
    axis.title         = element_text(size = 8.5),
    axis.line          = element_line(linewidth = 0.3),
    axis.ticks         = element_line(linewidth = 0.3),
    plot.margin        = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    panel.spacing      = unit(1.5, "lines")
  )


# Save outputs ----------------------------------------------------------------

ggsave(output_pdf, p_bar,
       width  = plot_width,
       height = plot_height,
       units  = "in")

saveRDS(p_bar, file = gg_rds_out)

message("Barplot saved: ", output_pdf)
message("ggplot RDS saved: ", gg_rds_out)


# Session info ----------------------------------------------------------------

sessionInfo()