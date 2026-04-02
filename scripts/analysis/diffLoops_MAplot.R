# ##############################################################################
# filename:    diffLoops_MAplot.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: MA plot and density companion for differential Hi-C loops;
#              colors gained loops red and lost loops blue
# ##############################################################################

# Libraries ----
library(DESeq2)
library(InteractionSet)
library(ggplot2)
library(dplyr)
library(glue)
library(stringr)
library(purrr)
library(apeglm)
library(RColorBrewer)
library(plotgardener)

# Parameters ----
diff_loops_rds     <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
output_pdf         <- "plots/diffLoops_MAplot.pdf"
output_density_pdf <- "plots/diffLoops_MAplot_density.pdf"

# Data import ----
res <- readRDS(diff_loops_rds) |> interactions()

# Analysis ----
res_df <- res |>
  as.data.frame() |>
  dplyr::select(baseMean, log2FoldChange, padj) |>
  mutate(isDE = case_when(
    log2FoldChange > 0 & padj < 0.1 ~ "TRUE - upreg",
    log2FoldChange < 0 & padj < 0.1 ~ "TRUE - downreg",
    padj > 0.1           ~ "FALSE",
    is.character("NA")   ~ "FALSE"
  )) |>
  arrange(isDE)

# Visualization ----

## MA plot ----
res_gg <- res_df |>
  ggplot(aes(x = baseMean, y = log2FoldChange, color = isDE)) +
  geom_point(alpha = 1) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  scale_color_manual(values = c("TRUE - upreg"   = "#F8766D",
                                "TRUE - downreg" = "#619CFF",
                                "FALSE"          = "grey80")) +
  ylim(c(-4, 4)) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100)) +
  labs(y = "Hi-C Contact log2FoldChange (sorbitol/untreated)",
       x = "mean of normalized counts") +
  theme_classic() +
  theme(legend.position = "NONE") +
  annotate(geom = "text", label = "Gained", x = 100, y = 3.5,  color = "#F8766D") +
  annotate(geom = "text",
           label = paste0("n = ", nrow(subset(res_df, isDE == "TRUE - upreg"))),
           x = 100, y = 3.1, color = "#F8766D") +
  annotate(geom = "text", label = "Lost",   x = 100, y = -3.5, color = "#619CFF") +
  annotate(geom = "text",
           label = paste0("n = ", nrow(subset(res_df, isDE == "TRUE - downreg"))),
           x = 100, y = -3.9, color = "#619CFF")

## Density companion ----
densityMA <- ggplot(res_df, aes(y = log2FoldChange)) +
  geom_density(color = "#619CFF", fill = 4, alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  ylim(c(-4, 4)) +
  xlim(c(0, 1.5)) +
  theme_classic() +
  theme(legend.position = "NONE",
        axis.text   = element_blank(),
        axis.title  = element_blank(),
        axis.ticks  = element_blank(),
        axis.line.x = element_blank())

# Save outputs ----
ggsave(output_pdf, plot = res_gg, device = "pdf", width = 8, height = 5)
ggsave(output_density_pdf, plot = densityMA, device = "pdf", width = 3, height = 5)

sessionInfo()
