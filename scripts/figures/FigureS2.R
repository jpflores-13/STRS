## Create Figure S2: matchRanges diagnostic plots
## Three-panel figure showing covariate matching quality

# Load required libraries -------------------------------------------------
library(InteractionSet)
library(strawr)
library(tidyverse)
library(glue)
library(nullranges)
library(data.table)
library(mariner)
library(plotgardener)

#------------------------- Read in and filter data -----------------------#

## Load in each HiC replicate 
hicFiles <- list.files("/proj/phanstiel_lab/Data/processed/YAPP/hic/hg38/220715_dietJuicerCore/output/",
                       full.names = TRUE,
                       recursive = TRUE,
                       pattern = "inter_30.hic")

## Load in loops
noDroso_loops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions() |> 
  pullHicPixels(binSize = 10e3,
                files = hicFiles,
                half = "both",
                norm = "VC_SQRT",
                matrix = "observed")

## create an mcol for loop size
mcols(noDroso_loops)$loop_size <- pairdist(noDroso_loops)

## create an mcol for loop type
mcols(noDroso_loops)$loop_type <- case_when(
  mcols(noDroso_loops)$padj < 0.05 & mcols(noDroso_loops)$log2FoldChange > 1 &
    mcols(noDroso_loops)$loop_size >= 150000 ~ "truegained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange > 0 &
    mcols(noDroso_loops)$loop_size >= 150000 ~ "gained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange < 0 ~ "lost",
  mcols(noDroso_loops)$padj > 0.1 ~ "static",
  is.character("NA") ~ "other")

## create aggregated sorb counts mcol
mcols(noDroso_loops)$sorb_contacts <- counts(noDroso_loops)[,"YAPP_HEK_sorbitol_4_2_inter_30.hic"] +
  counts(noDroso_loops)[,"YAPP_HEK_sorbitol_5_2_inter_30.hic"] + 
  counts(noDroso_loops)[,"YAPP_HEK_sorbitol_6_2_inter_30.hic"]

## create aggregated cont counts mcol
mcols(noDroso_loops)$cont_contacts <- counts(noDroso_loops)[,"YAPP_HEK_control_1_2_inter_30.hic"] +
  counts(noDroso_loops)[,"YAPP_HEK_control_2_2_inter_30.hic"] + 
  counts(noDroso_loops)[,"YAPP_HEK_control_3_2_inter_30.hic"]

## create aggregated sorb+cont counts mcol
mcols(noDroso_loops)$agg_contacts <- mcols(noDroso_loops)$sorb_contacts + mcols(noDroso_loops)$cont_contacts

## log transform loop_size to get a normal distribution for matchRanges
mcols(noDroso_loops)$loop_size <- log(mcols(noDroso_loops)$loop_size)

## convert all `0` agg_contact values to NAs and log transform for matchRanges
mcols(noDroso_loops)$agg_contacts[mcols(noDroso_loops)$agg_contacts == 0] <- NA
noDroso_loops <- interactions(noDroso_loops) |>
  as.data.frame() |>
  na.omit() |>
  as_ginteractions()

mcols(noDroso_loops)$agg_contacts <- log((mcols(noDroso_loops)$agg_contacts + 1))

# Matched set for gained YAPP loops ---------------------------------------
focal <- noDroso_loops[!noDroso_loops$loop_type %in% c("static", "lost","other")] 
pool <- noDroso_loops[noDroso_loops$loop_type %in% c("static", "lost","other")] 

nullSet <- matchRanges(focal = focal,
                       pool = pool,
                       covar = ~ agg_contacts + loop_size, 
                       method = 'stratified',
                       replace = FALSE)

# Create output directory -------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# Initialize PDF ----------------------------------------------------------
pdf("figures/FigureS2.pdf", width = 11, height = 4)

# Create page -------------------------------------------------------------
pageCreate(width = 11, height = 4, showGuides = FALSE)

# Panel A: Aggregate contacts covariate -----------------------------------
plotText(
  label = "A",
  fontsize = 14,
  fontface = "bold",
  x = 0.1, y = 0.1,
  just = c("left", "top")
)

## Capture the plotCovariate output and place it
plotA <- plotGG(
  plot = plotCovariate(nullSet),
  x = 0.25, y = 0.25,
  width = 3.25, height = 3.25,
  just = c("left", "top")
)

# Panel B: Loop size covariate --------------------------------------------
plotText(
  label = "B",
  fontsize = 14,
  fontface = "bold",
  x = 3.85, y = 0.1,
  just = c("left", "top")
)

plotB <- plotGG(
  plot = plotCovariate(nullSet, covar = "loop_size"),
  x = 4.0, y = 0.25,
  width = 3.25, height = 3.25,
  just = c("left", "top")
)

# Panel C: Propensity scores ----------------------------------------------
plotText(
  label = "C",
  fontsize = 14,
  fontface = "bold",
  x = 7.6, y = 0.1,
  just = c("left", "top")
)

## Create the propensity plot with custom x-axis breaks
propensity_plot <- plotPropensity(nullSet, sets = c('f', 'p', 'm'), log = 'x') +
  scale_x_log10(
    breaks = scales::breaks_log(n = 5),
    labels = scales::label_log()
  )

plotC <- plotGG(
  plot = propensity_plot,
  x = 7.75, y = 0.25,
  width = 3.25, height = 3.25,
  just = c("left", "top")
)

# Close device ------------------------------------------------------------
dev.off()

cat("\n===========================================\n")
cat("FIGURE CREATED: figures/FigureS2.pdf\n")
cat("===========================================\n")
cat("Panel A: Aggregate contacts covariate matching\n")
cat("Panel B: Loop size covariate matching\n")
cat("Panel C: Propensity score distributions\n")
cat("===========================================\n\n")