# ##############################################################################
# filename:    eisa_auxinRAD21_histogram.R
# author:      JP Flores
# date:        2026-05-29
# project:     STRS
# description: Overlapping density plot showing the distribution of delta-intron
#              (nascent RNA log2FC from eisaR) at genes whose promoters overlap
#              sorbitol-induced (gained) loop anchors. Two distributions are
#              overlaid: sorbitol 6h vs untreated and sorbitol+auxin 6h vs
#              untreated. Mirrors auxinRAD21_gainedLoops_histogram.R but uses
#              nascent RNA (lfc_intron from eisaR) rather than total RNA log2FC.
#              Centered on 0 with a dashed reference line. Paired Wilcoxon
#              signed-rank test.
# input:       data/processed/rna/AuxinRAD21/output/EISA/eisa_auxinRAD21_results.rds
# output:      plots/eisa_auxinRAD21_histogram.pdf
#              data/processed/rna/AuxinRAD21/output/EISA/gg_eisa_auxinRAD21_histogram.rds
# ##############################################################################


# Libraries -------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(stringr)


# Utility scripts -------------------------------------------------------------

source("scripts/utils/ggplot2_pgTheme.R")


# Parameters ------------------------------------------------------------------

project_dir <- "/work/users/j/p/jpflores/projects/STRS"
plots_dir   <- file.path(project_dir, "plots")
eisa_dir    <- file.path(project_dir,
                         "data/processed/rna/AuxinRAD21/output/EISA/per_replicate")

eisa_rds   <- file.path(eisa_dir, "eisa_auxinRAD21_results.rds")
gg_rds_out <- file.path(eisa_dir, "gg_eisa_auxinRAD21_histogram.rds")

color_sorbitol <- "#F8766D"
color_auxin    <- "#009E73"   # Okabe-Ito green

plot_width  <- 5
plot_height <- 4

output_pdf <- file.path(plots_dir, "eisa_auxinRAD21_histogram.pdf")


# Load data -------------------------------------------------------------------

eisa_df <- readRDS(eisa_rds)

message("EISA AuxinRAD21 results: ", nrow(eisa_df), " genes")
message("  Gained anchor genes: ", sum(eisa_df$loop_category == "gained"))


# Wrangle data ----------------------------------------------------------------

gained_df <- eisa_df |> dplyr::filter(loop_category == "gained")

hist_df <- dplyr::bind_rows(
  data.frame(
    gene_id   = gained_df$gene_id,
    log2FC    = gained_df$lfc_intron_sorb,
    condition = "Sorbitol"
  ),
  data.frame(
    gene_id   = gained_df$gene_id,
    log2FC    = gained_df$lfc_intron_auxin,
    condition = "Sorbitol + Auxin"
  )
) |>
  dplyr::filter(!is.na(log2FC))

hist_df$condition <- factor(hist_df$condition,
                            levels = c("Sorbitol", "Sorbitol + Auxin"))


# Statistics ------------------------------------------------------------------

## Paired Wilcoxon signed-rank test — same gene measured in both conditions
common_genes     <- intersect(
  gained_df$gene_id[!is.na(gained_df$lfc_intron_sorb)],
  gained_df$gene_id[!is.na(gained_df$lfc_intron_auxin)]
)
lfc_sorb_paired  <- gained_df$lfc_intron_sorb[
  match(common_genes, gained_df$gene_id)]
lfc_auxin_paired <- gained_df$lfc_intron_auxin[
  match(common_genes, gained_df$gene_id)]

valid_pairs      <- !is.na(lfc_sorb_paired) & !is.na(lfc_auxin_paired)
lfc_sorb_paired  <- lfc_sorb_paired[valid_pairs]
lfc_auxin_paired <- lfc_auxin_paired[valid_pairs]
n_genes          <- sum(valid_pairs)

wt <- wilcox.test(lfc_sorb_paired, lfc_auxin_paired,
                  paired      = TRUE,
                  alternative = "two.sided")

message("Paired Wilcoxon signed-rank test (delta-intron: sorbitol vs auxin):")
message("  n = ", n_genes, " genes")
message("  W = ", wt$statistic)
message("  p = ", wt$p.value)

p_label <- if (wt$p.value < 0.001) {
  "italic(p) < 0.001"
} else {
  paste0("italic(p) == ", round(wt$p.value, 3))
}
n_label <- paste0("n = ", n_genes, " genes")


# Visualization ---------------------------------------------------------------

dens_sorb  <- density(lfc_sorb_paired,  na.rm = TRUE)
dens_auxin <- density(lfc_auxin_paired, na.rm = TRUE)

sorb_label_x  <- dens_sorb$x[which.max(dens_sorb$y)]   + 0.6
sorb_label_y  <- max(dens_sorb$y) * 0.85
auxin_label_x <- dens_auxin$x[which.max(dens_auxin$y)] - 2.5
auxin_label_y <- max(dens_auxin$y) * 0.65

p_hist <- ggplot(hist_df, aes(x = log2FC, fill = condition)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray75", linewidth = 0.3) +
  geom_density(alpha = 0.4, color = NA) +
  annotate("text",
           x        = sorb_label_x,
           y        = sorb_label_y,
           label    = "Sorbitol",
           color    = color_sorbitol,
           fontface = "bold",
           size     = 7 / .pt,
           hjust    = 0) +
  annotate("text",
           x          = auxin_label_x,
           y          = auxin_label_y,
           label      = "Sorbitol +\nAuxin",
           color      = color_auxin,
           fontface   = "bold",
           size       = 7 / .pt,
           hjust      = 0,
           lineheight = 0.8) +
  annotate("text",
           x     = 4.4,
           y     = max(dens_sorb$y) * 0.95,
           label = p_label,
           parse = TRUE,
           size  = 7 / .pt,
           hjust = 1,
           color = "black") +
  annotate("text",
           x     = 4.4,
           y     = max(dens_sorb$y) * 0.83,
           label = n_label,
           size  = 6 / .pt,
           hjust = 1,
           color = "gray40") +
  scale_fill_manual(values = c(
    "Sorbitol"         = color_sorbitol,
    "Sorbitol + Auxin" = color_auxin
  )) +
  coord_cartesian(xlim = c(-4.5, 4.5)) +
  scale_x_continuous(limits = c(-4.5, 4.5)) +
  labs(
    x     = expression(Delta*"intron " * log[2]~FC~(treated / untreated)),
    y     = "Density",
    title = "Nascent RNA at Sorbitol-Induced Loop Anchors"
  ) +
  theme_classic() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 9),
    legend.position = "none",
    axis.text       = element_text(size = 7, color = "black"),
    axis.title      = element_text(size = 8.5),
    axis.line       = element_line(linewidth = 0.3),
    axis.ticks      = element_line(linewidth = 0.3),
    plot.margin     = margin(5.5, 5.5, 5.5, 5.5, "pt")
  )


# Save outputs ----------------------------------------------------------------

ggsave(output_pdf, p_hist,
       width  = plot_width,
       height = plot_height,
       units  = "in")

saveRDS(p_hist, file = gg_rds_out)

message("Histogram saved: ", output_pdf)
message("ggplot RDS saved: ", gg_rds_out)


# Session info ----------------------------------------------------------------

sessionInfo()