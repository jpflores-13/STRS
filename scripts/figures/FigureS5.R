## Figure S5
pdf("figures/FigureS5.pdf",
    width = 10,
    height = 4)

library(DESeq2)
library(InteractionSet)
library(glue)
library(apeglm)
library(RColorBrewer)
library(plotgardener)
library(tidyverse)
library(dbscan)
library(data.table)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(mariner)
library(plyranges)
source("scripts/utils/ggplot2_pgTheme.R")

# Load data ---------------------------------------------------------------

## Load in differential loops
diffLoops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds")
normalizedAPAs <- list.files("data/processed/hic/normalizedAPA",
                             full.names = TRUE) |> 
  lapply(readRDS)

names(normalizedAPAs) <- list.files("data/processed/hic/normalizedAPA/",
                                    recursive = TRUE)

# Load read depth data ----------------------------------------------------

## Load the updated condition averages for read depth annotations
depth_data <- read_csv("data/raw/hic/sequencing_depth_results/updated_condition_averages.csv",
                       show_col_types = FALSE)

## Create a lookup function to get read depth in millions
get_depth_millions <- function(cell_type, condition) {
  condition_standard <- case_when(
    condition %in% c("Control", "Untreated", "untreated") ~ "Control",
    condition %in% c("Sorbitol", "NaCl_Treated", "Treated", "+ sorbitol", "+ NaCl") ~ "Treated",
    TRUE ~ condition
  )
  
  depth_row <- depth_data |>
    filter(Cell_Type == cell_type, 
           condition_standard == condition_standard)
  
  if (nrow(depth_row) > 0) {
    depth_millions <- depth_row$average_depth_per_replicate
    return(sprintf("%.0fM", depth_millions))
  } else {
    warning(glue("Could not find depth for {cell_type} - {condition_standard}"))
    return("N/A")
  }
}

## Pre-calculate depths for the three conditions
depths <- list(
  hek_egfp = get_depth_millions("HEK293_eGFP_YAP", "Control"),
  hek_wt = get_depth_millions("HEK293_WT", "Control"),
  hek_dtad = "502M"  # Average of 513.76 and 489.24
)

## Print depths to verify
cat("Calculated depths:\n")
cat("eGFP-YAP:", depths$hek_egfp, "\n")
cat("WT:", depths$hek_wt, "\n")
cat("dTAD:", depths$hek_dtad, "\n")

# Create gained & lost loop lists -----------------------------------------

## gained loops with a p-adj. value of < 0.1 and a (+) log2FC
gained_adj <- 
  diffLoops[which(rowData(diffLoops)$padj < 0.1 &
                    rowData(diffLoops)$log2FoldChange > 0)]

## lost loops with a p-adj. value of < 0.1 and a (-) log2FC
lost_adj <- 
  diffLoops[which(rowData(diffLoops)$padj < 0.1 &
                    rowData(diffLoops)$log2FoldChange < 0)]

# Define APA parameters ---------------------------------------------------

apaParams <- pgParams(assembly = "hg38",
                      height = 0.75,
                      width = 0.75,
                      fontsize = 5,
                      palette = colorRampPalette(brewer.pal(n = 9, "YlGnBu")))

## Make page
pageCreate(width = 10,
           height = 4,
           showGuides = FALSE)

# Add panel label ---------------------------------------------------------

plotText(label = "A", 
         x = 0.25, 
         y = 0.25, 
         just = c("left", "top"),
         fontface = "bold", 
         fontsize = 14)

# Visualize normalized gained APAs ----------------------------------------

## Combine APA matrices to pull out the max value for zrange max
wt_gainedLoops_combined <- c(normalizedAPAs$wt_gainedLoops_contAPA_normalized.rds,
                             normalizedAPAs$wt_gainedLoops_sorbAPA_normalized.rds)

oe_gainedLoops_combined <- c(normalizedAPAs$oe_gainedLoops_contAPA_normalized.rds,
                             normalizedAPAs$oe_gainedLoops_sorbAPA_normalized.rds)

dtad_gainedLoops_combined <- c(normalizedAPAs$dtad_gainedLoops_contAPA_normalized.rds,
                               normalizedAPAs$dtad_gainedLoops_sorbAPA_normalized.rds)

# Plot eGFP-YAP gained loops ----------------------------------------------

plotMatrix(normalizedAPAs$oe_gainedLoops_contAPA_normalized.rds,
           x = 0.5,
           y = 0.5,
           zrange = c(0, max(oe_gainedLoops_combined)),
           params = apaParams) |> 
  annoHeatmapLegend(x = 1.3,
                    y = 0.5,
                    width = 0.05,
                    height = 0.75,
                    fontcolor = "black",
                    fontsize = 5)

plotMatrix(normalizedAPAs$oe_gainedLoops_sorbAPA_normalized.rds,
           x = 0.5,
           y = 1.3,
           zrange = c(0, max(oe_gainedLoops_combined)),
           params = apaParams) 

plotText(label = depths$hek_egfp,
         fontcolor = "black",
         fontsize = 4,
         x = 1.25,
         y = 2.05,
         just = c("right", "bottom"))

# Plot WT gained loops ----------------------------------------------------

plotMatrix(normalizedAPAs$wt_gainedLoops_contAPA_normalized.rds,
           x = 1.45,
           y = 0.5,
           zrange = c(0, max(wt_gainedLoops_combined)),
           params = apaParams) |> 
  annoHeatmapLegend(x = 2.25,
                    y = 0.5,
                    width = 0.05,
                    height = 0.75,
                    fontcolor = "black",
                    fontsize = 5)

plotMatrix(normalizedAPAs$wt_gainedLoops_sorbAPA_normalized.rds,
           x = 1.45,
           y = 1.3,
           zrange = c(0, max(wt_gainedLoops_combined)),
           params = apaParams) 

plotText(label = depths$hek_wt,
         fontcolor = "black",
         fontsize = 4,
         x = 2.20,
         y = 2.05,
         just = c("right", "bottom"))

# Plot dTAD gained loops --------------------------------------------------

plotMatrix(normalizedAPAs$dtad_gainedLoops_contAPA_normalized.rds,
           x = 2.4,
           y = 0.5,
           zrange = c(0, max(dtad_gainedLoops_combined)),
           params = apaParams) |> 
  annoHeatmapLegend(x = 3.2,
                    y = 0.5,
                    width = 0.05,
                    height = 0.75,
                    fontcolor = "black",
                    fontsize = 5)

plotMatrix(normalizedAPAs$dtad_gainedLoops_sorbAPA_normalized.rds,
           x = 2.4,
           y = 1.3,
           zrange = c(0, max(dtad_gainedLoops_combined)),
           params = apaParams)

plotText(label = depths$hek_dtad,
         fontcolor = "black",
         fontsize = 4,
         x = 3.15,
         y = 2.05,
         just = c("right", "bottom"))

# Add labels for gained loops ---------------------------------------------

plotText(label = "HEK293 eGFP-YAP", fontcolor = "black",
         x = 0.875,
         y = 2.15,
         params = apaParams,
         just = c("center", "top"))

plotText(label = "HEK293 WT", fontcolor = "black",
         x = 1.825,
         y = 2.15,
         params = apaParams,
         just = c("center", "top"))

plotText(label = "HEK293 eGFP-YAPdTAD", fontcolor = "black",
         x = 2.775,
         y = 2.15,
         params = apaParams,
         just = c("center", "top"))

plotText(label = "untreated", fontcolor = "black",
         x = 0.45,
         y = 0.875,
         params = apaParams,
         rot = 90,
         just = c("center", "bottom"))

plotText(label = "+ sorbitol", fontcolor = "black",
         x = 0.45,
         y = 1.675,
         params = apaParams,
         rot = 90,
         just = c("center", "bottom"))

# Visualize normalized lost APAs ------------------------------------------

## Combine APA matrices to pull out the max value for zrange max
wt_lostLoops_combined <- c(normalizedAPAs$wt_lostLoops_contAPA_normalized.rds,
                           normalizedAPAs$wt_lostLoops_sorbAPA_normalized.rds)

oe_lostLoops_combined <- c(normalizedAPAs$oe_lostLoops_contAPA_normalized.rds,
                           normalizedAPAs$oe_lostLoops_sorbAPA_normalized.rds)

dtad_lostLoops_combined <- c(normalizedAPAs$dtad_lostLoops_contAPA_normalized.rds,
                             normalizedAPAs$dtad_lostLoops_sorbAPA_normalized.rds)

# Plot eGFP-YAP lost loops ------------------------------------------------

plotMatrix(normalizedAPAs$oe_lostLoops_contAPA_normalized.rds,
           x = 3.5,
           y = 0.5,
           zrange = c(0, max(oe_lostLoops_combined)),
           params = apaParams) |>
  annoHeatmapLegend(x = 4.3,
                    y = 0.5,
                    width = 0.05,
                    height = 0.75,
                    fontcolor = "black",
                    fontsize = 5)

plotMatrix(normalizedAPAs$oe_lostLoops_sorbAPA_normalized.rds,
           x = 3.5,
           y = 1.3,
           zrange = c(0, max(oe_lostLoops_combined)),
           params = apaParams)

plotText(label = depths$hek_egfp,
         fontcolor = "black",
         fontsize = 4,
         x = 4.25,
         y = 2.05,
         just = c("right", "bottom"))

# Plot WT lost loops ------------------------------------------------------

plotMatrix(normalizedAPAs$wt_lostLoops_contAPA_normalized.rds,
           x = 4.45,
           y = 0.5,
           zrange = c(0, max(wt_lostLoops_combined)),
           params = apaParams) |>
  annoHeatmapLegend(x = 5.25,
                    y = 0.5,
                    width = 0.05,
                    height = 0.75,
                    fontcolor = "black",
                    fontsize = 5)

plotMatrix(normalizedAPAs$wt_lostLoops_sorbAPA_normalized.rds,
           x = 4.45,
           y = 1.3,
           zrange = c(0, max(wt_lostLoops_combined)),
           params = apaParams)

plotText(label = depths$hek_wt,
         fontcolor = "black",
         fontsize = 4,
         x = 5.20,
         y = 2.05,
         just = c("right", "bottom"))

# Plot dTAD lost loops ----------------------------------------------------

plotMatrix(normalizedAPAs$dtad_lostLoops_contAPA_normalized.rds,
           x = 5.4,
           y = 0.5,
           zrange = c(0, max(dtad_lostLoops_combined)),
           params = apaParams) |>
  annoHeatmapLegend(x = 6.2,
                    y = 0.5,
                    width = 0.05,
                    height = 0.75,
                    fontcolor = "black",
                    fontsize = 5)

plotMatrix(normalizedAPAs$dtad_lostLoops_sorbAPA_normalized.rds,
           x = 5.4,
           y = 1.3,
           zrange = c(0, max(dtad_lostLoops_combined)),
           params = apaParams) 

plotText(label = depths$hek_dtad,
         fontcolor = "black",
         fontsize = 4,
         x = 6.15,
         y = 2.05,
         just = c("right", "bottom"))

# Add labels for lost loops -----------------------------------------------

plotText(label = "HEK293 eGFP-YAP", fontcolor = "black",
         x = 3.875,
         y = 2.15,
         params = apaParams,
         just = c("center", "top"))

plotText(label = "HEK293 WT", fontcolor = "black",
         x = 4.825,
         y = 2.15,
         params = apaParams,
         just = c("center", "top"))

plotText(label = "HEK293 eGFP-YAPdTAD", fontcolor = "black",
         x = 5.775,
         y = 2.15,
         params = apaParams,
         just = c("center", "top"))

# Add dividing lines and section labels -----------------------------------

plotSegments(
  x0 = 0.5,
  y0 = 0.4,
  x1 = 3.3,
  y1 = 0.4,
  default.units = "inches",
  linecolor = "#F8766D",
  lwd = 1,
  lty = 1)

plotSegments(
  x0 = 3.5,
  y0 = 0.4,
  x1 = 6.3,
  y1 = 0.4,
  default.units = "inches",
  linecolor = "#619CFF",
  lwd = 1,
  lty = 1)

plotText(label = "gained loops",
         x = 1.9,
         y = 0.3,
         fontcolor = "#F8766D",
         fontsize = 10)

plotText(label = "lost loops",
         x = 4.9,
         y = 0.3,
         fontcolor = "#619CFF",
         fontsize = 10)

dev.off()