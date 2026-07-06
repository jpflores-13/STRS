# ##############################################################################
# filename:    FigureS4.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Supplementary Figure S4 — matchRanges diagnostic plots showing
#              covariate matching quality for aggregate contacts and loop size;
#              three-panel plotgardener figure
# ##############################################################################

# Libraries ----
library(InteractionSet)
library(strawr)
library(tidyverse)
library(glue)
library(nullranges)
library(data.table)
library(mariner)
library(plotgardener)

# Parameters ----
hic_dir        <- "data/processed/hic/hg38/220715_dietJuicerCore/output"
diff_loops_rds <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
output_pdf     <- "figures/FigureS4.pdf"
page_width     <- 11
page_height    <- 4

# Data import ----
hicFiles <- list.files(hic_dir,
                       full.names = TRUE,
                       recursive  = TRUE,
                       pattern    = "inter_30.hic")

noDroso_loops <- readRDS(diff_loops_rds) |>
  interactions() |>
  pullHicPixels(binSize = 10e3,
                files   = hicFiles,
                half    = "both",
                norm    = "VC_SQRT",
                matrix  = "observed")

# Analysis ----
mcols(noDroso_loops)$loop_size <- pairdist(noDroso_loops)

mcols(noDroso_loops)$loop_type <- case_when(
  mcols(noDroso_loops)$padj < 0.05 & mcols(noDroso_loops)$log2FoldChange > 1 &
    mcols(noDroso_loops)$loop_size >= 150000 ~ "truegained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange > 0 &
    mcols(noDroso_loops)$loop_size >= 150000 ~ "gained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange < 0 ~ "lost",
  mcols(noDroso_loops)$padj > 0.1 ~ "static",
  TRUE ~ "other")

mcols(noDroso_loops)$sorb_contacts <-
  counts(noDroso_loops)[, "YAPP_HEK_sorbitol_4_2_inter_30.hic"] +
  counts(noDroso_loops)[, "YAPP_HEK_sorbitol_5_2_inter_30.hic"] +
  counts(noDroso_loops)[, "YAPP_HEK_sorbitol_6_2_inter_30.hic"]

mcols(noDroso_loops)$cont_contacts <-
  counts(noDroso_loops)[, "YAPP_HEK_control_1_2_inter_30.hic"] +
  counts(noDroso_loops)[, "YAPP_HEK_control_2_2_inter_30.hic"] +
  counts(noDroso_loops)[, "YAPP_HEK_control_3_2_inter_30.hic"]

mcols(noDroso_loops)$agg_contacts <-
  mcols(noDroso_loops)$sorb_contacts + mcols(noDroso_loops)$cont_contacts

mcols(noDroso_loops)$loop_size <- log(mcols(noDroso_loops)$loop_size)

mcols(noDroso_loops)$agg_contacts[mcols(noDroso_loops)$agg_contacts == 0] <- NA
noDroso_loops <- interactions(noDroso_loops) |>
  as.data.frame() |>
  na.omit() |>
  as_ginteractions()

mcols(noDroso_loops)$agg_contacts <- log((mcols(noDroso_loops)$agg_contacts + 1))

focal   <- noDroso_loops[!noDroso_loops$loop_type %in% c("static", "lost", "other")]
pool    <- noDroso_loops[noDroso_loops$loop_type  %in% c("static", "lost", "other")]

nullSet <- matchRanges(focal  = focal,
                       pool   = pool,
                       covar  = ~ agg_contacts + loop_size,
                       method = "stratified",
                       replace = FALSE)

# Visualization ----
pdf(output_pdf, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = FALSE)

plotText(label = "A", fontsize = 14, fontface = "bold",
         x = 0.1, y = 0.1, just = c("left", "top"))

plotA_clean <- plotCovariate(nullSet) +
  labs(x = "Aggregate Contacts", y = "Density") +
  scale_fill_discrete(name = NULL,
                      labels = c("f" = "Focal", "m" = "Matched", "p" = "Pool"))
plotGG(plot = plotA_clean, x = 0.25, y = 0.25,
       width = 3.25, height = 3.25, just = c("left", "top"))

plotText(label = "B", fontsize = 14, fontface = "bold",
         x = 3.85, y = 0.1, just = c("left", "top"))

plotB_clean <- plotCovariate(nullSet, covar = "loop_size") +
  labs(x = "Loop Size", y = "Density") +
  scale_fill_discrete(name = NULL,
                      labels = c("f" = "Focal", "m" = "Matched", "p" = "Pool"))
plotGG(plot = plotB_clean, x = 4.0, y = 0.25,
       width = 3.25, height = 3.25, just = c("left", "top"))

plotText(label = "C", fontsize = 14, fontface = "bold",
         x = 7.6, y = 0.1, just = c("left", "top"))

propensity_plot <- plotPropensity(nullSet, sets = c("f", "p", "m"), log = "x") +
  scale_x_log10(breaks = scales::breaks_log(n = 5),
                labels = scales::label_log()) +
  labs(title = "~ Aggregate Contacts + Loop Size",
       x = "log(Propensity Score)",
       y = "Density") +
  scale_fill_discrete(name = NULL,
                      labels = c("f" = "Focal", "m" = "Matched", "p" = "Pool")) +
  theme(plot.title = element_text(hjust = 0.5))

plotGG(plot = propensity_plot, x = 7.75, y = 0.25,
       width = 3.25, height = 3.25, just = c("left", "top"))

# Save outputs ----
dev.off()

sessionInfo()
