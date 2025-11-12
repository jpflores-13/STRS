## Figure 1 

pdf("figures/Figure1.pdf",
    width = 8.5,
    height = 11)

## Perform differential loops analysis with DESeq2

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
                             full.names = T) |> 
  lapply(readRDS)

names(normalizedAPAs) <- list.files("data/processed/hic/normalizedAPA/",
                                    recursive = T)

## Calculate loop size
rowData(diffLoops)$loop_size <- pairdist(diffLoops)

# Load read depth data ----------------------------------------------------

## Load the eight condition averages for read depth annotations
depth_data <- read_csv("data/raw/hic/sequencing_depth_results/eight_condition_averages.csv",
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

## Pre-calculate depths for all conditions
depths <- list(
  hek_egfp = get_depth_millions("HEK293_eGFP_YAP", "Control"),
  hek_wt = get_depth_millions("HEK293_WT", "Control"),
  hct116 = get_depth_millions("HCT116", "Control"),
  t47d = get_depth_millions("T47D", "Control")
)

# Create gained & lost loop lists -----------------------------------------

## gained loops with a p-adj. value of < 0.1 and a (+) log2FC
gained_adj <- 
  diffLoops[which(rowData(diffLoops)$padj < 0.1 &
                    rowData(diffLoops)$log2FoldChange >0)]

## lost loops with a p-adj. value of < 0.1 and a (-) log2FC
lost_adj <- 
  diffLoops[which(rowData(diffLoops)$padj < 0.1 &
                    rowData(diffLoops)$log2FoldChange < 0)]

## filter for the 100 best gained & lost loops
bestGained <- head(gained_adj[order(rowData(gained_adj)$padj, decreasing = F)], 100)
bestLost <- head(lost_adj[order(rowData(lost_adj)$padj, decreasing = F)], 100)

# top 100 gained & lost loops as GRanges
loopRegions_gained <-
  GRanges(seqnames = as.character(seqnames(anchors(x = bestGained, "first"))),
          ranges = IRanges(start = start(anchors(bestGained, "first")),
                           end = end(anchors(bestGained, "second"))),
          mcols = mcols(bestGained))

loopRegions_lost <-
  GRanges(seqnames = as.character(seqnames(anchors(x = bestLost, "first"))),
          ranges = IRanges(start = start(anchors(bestLost, "first")),
                           end = end(anchors(bestLost, "second"))),
          mcols = mcols(bestLost))

## Expand regions by buffer
buffer <- 250e3
loopRegions_gained_buffed <- loopRegions_gained + buffer
loopRegions_lost_buffed <- loopRegions_lost + buffer

# load loop lists --------------------------------------------------------
lostLoops <- diffLoops[which(rowData(diffLoops)$padj < 0.1 &
                               rowData(diffLoops)$log2FoldChange < 0)] 

gainedLoops <- diffLoops[which(rowData(diffLoops)$padj < 0.1 &
                                 rowData(diffLoops)$log2FoldChange > 0)]

## Define parameters
pix <- 0.25
gainedLostParams <- pgParams(assembly = "hg38",
                             resolution = 10e3,
                             chrom = as.character(seqnames(loopRegions_gained_buffed))[36],
                             chromstart = start(loopRegions_gained_buffed)[36],
                             chromend = end(loopRegions_gained_buffed)[36],
                             zrange = c(0,100),
                             norm = "SCALE",
                             x = 0,
                             width = 1.5,
                             height = 0.75,
                             fontsize = 5)

gainedParams <- pgParams(assembly = "hg38",
                         resolution = 10e3,
                         chrom = as.character(seqnames(loopRegions_gained_buffed))[38],
                         chromstart = start(loopRegions_gained_buffed)[38],
                         chromend = end(loopRegions_gained_buffed)[38],
                         zrange = c(0,100),
                         norm = "SCALE",
                         x = 0,
                         width = 1.5,
                         height = 0.75,
                         fontsize = 5)

lostParams <- pgParams(assembly = "hg38",
                       resolution = 10e3,
                       chrom = as.character(seqnames(loopRegions_lost_buffed))[20],
                       chromstart = start(loopRegions_lost_buffed)[20],
                       chromend = end(loopRegions_lost_buffed)[20],
                       zrange = c(0,100),
                       norm = "SCALE",
                       x = 0,
                       width = 1.5,
                       height = 0.75,
                       fontsize = 5)

## Make page
pageCreate(width = 8, height = 5.75,
           showGuides = F)

## Plot gained + lost loop example region

gainedLost_control <- plotHicRectangle(data = "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic",
                                       params = gainedLostParams,
                                       y = pix)

annoHeatmapLegend(gainedLost_control,
                  params = gainedLostParams,
                  orientation = "v",
                  fontcolor = "black",
                  digits = 2,
                  x = gainedLostParams$width + 0.05,
                  y = pix,
                  width = 0.05,
                  just = c("left", "top"),
                  default.units = "inches")

annoPixels(gainedLost_control,
           data = interactions(lostLoops),
           type = "arrow",
           col = "#619CFF",
           lwd = 0.001,
           lty = 0.001,
           shift = 10)

annoPixels(gainedLost_control,
           data = interactions(gainedLoops),
           type = "arrow",
           col = "#F8766D",
           lwd = 0.001,
           lty = 0.001,
           shift = 10)

## Plot bottom Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls

gainedLost_sorb <- 
  plotHicRectangle(data = "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic", 
                   params = gainedLostParams,
                   y = gainedLostParams$height + 0.3)

annoPixels(gainedLost_sorb,
           data = interactions(lostLoops),
           type = "arrow",
           col = "#619CFF",
           lwd = 0.001,
           lty = 0.001,
           shift = 10)

annoPixels(gainedLost_sorb,
           data = interactions(gainedLoops),
           type = "arrow",
           col = "#F8766D",
           lwd = 0.001,
           lty = 0.001,
           shift = 10)

## Plot genes
plotGenes(param = gainedLostParams,
          chrom = gainedLostParams$chrom,
          y = as.numeric(str_remove(gainedLost_sorb$y, "inches")) + 0.75,
          height = 0.5)

## Plot genome label
plotGenomeLabel(params = gainedLostParams,
                length = gainedLostParams$width,
                y = 2.3)

plotText(label = "untreated",
         params = gainedLostParams,
         x = 0,
         y = pix,
         just = c("top", "left"))

plotText(label = "200mM sorbitol",
         params = gainedLostParams,
         x = 0,
         y = pix + 0.8,
         just = c("top", "left"))

## Plot gained loop example region 

gained_control <- plotHicRectangle(data = "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic",
                                   params = gainedParams,
                                   x = gainedLostParams$width + 0.15,
                                   y = pix)

annoHeatmapLegend(gained_control,
                  params = gainedParams,
                  orientation = "v",
                  fontcolor = "black",
                  digits = 2,
                  x = as.numeric(str_remove(gained_control$x, "inches")) + 
                    gainedLostParams$width + 0.05,
                  y = pix,
                  width = 0.05,
                  just = c("left", "top"),
                  default.units = "inches")

annoPixels(gained_control,
           data = interactions(lostLoops),
           type = "arrow",
           col = "#619CFF",
           lwd = 0.001,
           lty = 0.001,
           shift = 10)

annoPixels(gained_control,
           data = interactions(gainedLoops),
           type = "arrow",
           col = "#F8766D",
           lwd = 0.001,
           lty = 0.001,
           shift = 10)

## Plot bottom Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls

gained_sorb <- 
  plotHicRectangle(data = "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic", 
                   params = gainedParams,
                   x = gainedLostParams$width + 0.15,
                   y = gainedLostParams$height + 0.3)

annoPixels(gained_sorb,
           data = interactions(lostLoops),
           type = "arrow",
           col = "#619CFF",
           lwd = 0.001,
           lty = 0.001,
           shift = 10)

annoPixels(gained_sorb,
           data = interactions(gainedLoops),
           type = "arrow",
           col = "#F8766D",
           lwd = 0.001,
           lty = 0.001,
           shift = 10)

## Plot genes
plotGenes(param = gainedParams,
          chrom = gainedParams$chrom,
          y = as.numeric(str_remove(gained_sorb$y, "inches")) + 0.75,
          x = as.numeric(str_remove(gained_control$x, "inches")),
          height = 0.5)

## Plot genome label
plotGenomeLabel(params = gainedParams,
                length = gainedParams$width,
                x = as.numeric(str_remove(gained_control$x, "inches")),
                y = 2.3)

plotText(label = "untreated",
         params = gainedParams,
         x = as.numeric(str_remove(gained_control$x, "inches")),
         y = pix,
         just = c("top", "left"))

plotText(label = "200mM sorbitol",
         params = gainedParams,
         x = as.numeric(str_remove(gained_control$x, "inches")),
         y = pix + 0.8,
         just = c("top", "left"))

## Plot lost loop example region

lost_control <- plotHicRectangle(data = "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic",
                                 params = lostParams,
                                 x = gainedLostParams$width + 
                                   as.numeric(str_remove(gained_control$width, "inches")) + 0.3,
                                 y = pix)

annoHeatmapLegend(lost_control,
                  params = lostParams,
                  orientation = "v",
                  fontcolor = "black",
                  digits = 2,
                  x = as.numeric(str_remove(lost_control$x, "inches")) + 
                    gainedLostParams$width + 0.05,
                  y = pix,
                  width = 0.05,
                  just = c("left", "top"),
                  default.units = "inches")

annoPixels(lost_control,
           data = interactions(lostLoops),
           type = "arrow",
           col = "#619CFF",
           lwd = 0.001,
           lty = 0.001,
           shift = 10)

annoPixels(lost_control,
           data = interactions(gainedLoops),
           type = "arrow",
           col = "#F8766D",
           lwd = 0.001,
           lty = 0.001,
           shift = 10)

## Plot bottom Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls

lost_sorb <- 
  plotHicRectangle(data = "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic", 
                   params = lostParams,
                   x = gainedLostParams$width + 
                     as.numeric(str_remove(gained_control$width, "inches")) + 0.3,
                   y = gainedLostParams$height + 0.3)

annoPixels(lost_sorb,
           data = interactions(lostLoops),
           type = "arrow",
           col = "#619CFF",
           lwd = 0.001,
           lty = 0.001,
           shift = 10)

annoPixels(lost_sorb,
           data = interactions(gainedLoops),
           type = "arrow",
           col = "#F8766D",
           lwd = 0.001,
           lty = 0.001,
           shift = 10)

## Plot genes
plotGenes(param = lostParams,
          chrom = lostParams$chrom,
          y = as.numeric(str_remove(lost_sorb$y, "inches")) + 0.75,
          x = as.numeric(str_remove(lost_control$x, "inches")),
          height = 0.5)

## Plot genome label
plotGenomeLabel(params = lostParams,
                length = lostParams$width,
                x = as.numeric(str_remove(lost_control$x, "inches")),
                y = 2.3)

plotText(label = "untreated",
         params = lostParams,
         x = as.numeric(str_remove(lost_control$x, "inches")),
         y = pix,
         just = c("top", "left"))

plotText(label = "200mM sorbitol",
         params = lostParams,
         x = as.numeric(str_remove(lost_control$x, "inches")),
         y = pix + 0.8,
         just = c("top", "left"))

# manipulate dataset for custom MA plot -----------------------------------
res <- interactions(diffLoops)
res_df <- res |>
  as.data.frame() |>
  dplyr::select(baseMean, log2FoldChange, padj) |> 
  mutate(isDE = case_when(
    log2FoldChange > 0 & padj < 0.1 ~ "TRUE - upreg",
    log2FoldChange < 0 & padj < 0.1 ~ "TRUE - downreg",
    padj > 0.1  ~ "FALSE",
    is.character("NA") ~ "FALSE"
  )) |> 
  arrange(isDE)

# custom MA plot ----------------------------------------------------------

res_gg <- res_df |> 
  ggplot(aes(x = baseMean, y = log2FoldChange, color = isDE)) +
  geom_point(alpha = 1,
             size = 0.5) +
  geom_hline(yintercept = 0,
             linetype = 2,
             color = "grey40") +
  scale_color_manual(values = c("TRUE - upreg" = "#F8766D",
                                "TRUE - downreg" = "#619CFF",
                                "FALSE" = "grey80")) +
  ylim(c(-4,4)) +
  scale_x_log10(breaks=c(1,2,5,10,20,50,100)) +
  labs(y = "Hi-C Contact log2FoldChange (sorbitol/untreated)",
       x = "mean of normalized counts") +
  annotate(geom = "text",
           label = "Gained",
           x = 100,
           y = 3.5,
           color = "#F8766D",
           size  = 5/.pt) +
  annotate(geom = "text",
           label = paste0("n = ",
                          nrow(subset(res_df,
                                      isDE == "TRUE - upreg"))),
           x = 100,
           y = 3.1,
           color = "#F8766D",
           size  = 5/.pt) +
  annotate(geom = "text",
           label = "Lost",
           x = 100,
           y = -3.5,
           color = "#619CFF",
           size  = 5/.pt) +
  annotate(geom = "text",
           label = paste0("n = ",
                          nrow(subset(res_df,
                                      isDE == "TRUE - downreg"))),
           x = 100,
           y = -3.9,
           color = "#619CFF",
           size  = 5/.pt) +
  theme_pg +
  theme(
    legend.position = "NONE",
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.ticks.length = unit(-.1, "cm"),
    aspect.ratio = 1)

# Plot density of log2FoldChange ------------------------------------------

densityMA <- ggplot(res_df, aes(y = log2FoldChange)) +
  geom_density(
    color = "#9370DB",
    fill = "#9370DB",
    alpha = 0.25) +
  geom_hline(yintercept = 0,
             linetype = 2,
             color = "grey40") +
  ylim(c(-4,4)) +
  xlim(c(0, 2)) +
  theme_classic() +
  theme_pg +
  theme(
    legend.position = "NONE",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(-.1, "cm"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank())

plotGG(res_gg,
       params = gainedLostParams,
       x = gainedLostParams$width + 
         as.numeric(str_remove(gained_control$width, "inches")) + 2.25,
       y = pix,
       height = 2.3,
       width = 2.3,
       just = c("top", "left"))

plotGG(densityMA,
       x = gainedLostParams$width + 
         as.numeric(str_remove(gained_control$width, "inches")) + 4.6,
       y = pix,
       height = 2.05,
       width = 0.5,
       just = c("top", "left"))

## Visualize eGFP-YAP gained & lost loop APA plots comparing 
## HEK293 eGFP-YAP/WT, HCT116, & T47D genotypes

apaParams <- pgParams(assembly = "hg38",
                      height = 0.75,
                      width = 0.75,
                      fontsize = 5,
                      palette = colorRampPalette(brewer.pal(n = 9, "YlGnBu")))

# Visualize normalized gained APAs ----------------------------------------

## combine APA matrices to pull out the max value for zrange max
wt_gainedLoops_combined <- c(normalizedAPAs$wt_gainedLoops_contAPA_normalized.rds,
                             normalizedAPAs$wt_gainedLoops_sorbAPA_normalized.rds)

oe_gainedLoops_combined <- c(normalizedAPAs$oe_gainedLoops_contAPA_normalized.rds,
                             normalizedAPAs$oe_gainedLoops_sorbAPA_normalized.rds)

hct_gainedLoops_combined <- c(normalizedAPAs$hct_gainedLoops_contAPA_normalized.rds,
                              normalizedAPAs$hct_gainedLoops_sorbAPA_normalized.rds)

amat_gainedLoops_combined <- c(normalizedAPAs$amat_gainedLoops_contAPA_normalized.rds,
                               normalizedAPAs$amat_gainedLoops_sorbAPA_normalized.rds)

plotMatrix(normalizedAPAs$oe_gainedLoops_contAPA_normalized.rds,
           x = 0,
           y = 3.1,
           zrange = c(0, max(oe_gainedLoops_combined)),
           params = apaParams) |> 
  annoHeatmapLegend(x = 0.8,
                    y = 3.1,
                    width = 0.05,
                    height = 0.75,
                    fontcolor = "black",
                    fontsize = 5)

plotMatrix(normalizedAPAs$oe_gainedLoops_sorbAPA_normalized.rds,
           x = 0,
           y = 3.9,
           zrange = c(0, max(oe_gainedLoops_combined)),
           params = apaParams) 

plotText(label = depths$hek_egfp,
         fontcolor = "black",
         fontsize = 4,
         x = 0.75 - 0.05,
         y = 4.65 - 0.05,
         just = c("right", "bottom"))

plotMatrix(normalizedAPAs$wt_gainedLoops_contAPA_normalized.rds,
           x = 0.95,
           y = 3.1,
           zrange = c(0, max(wt_gainedLoops_combined)),
           params = apaParams)  |> 
  annoHeatmapLegend(x = 1.75,
                    y = 3.1,
                    width = 0.05,
                    height = 0.75,
                    fontcolor = "black",
                    fontsize = 5)

plotMatrix(normalizedAPAs$wt_gainedLoops_sorbAPA_normalized.rds,
           x = 0.95,
           y = 3.9,
           zrange = c(0, max(wt_gainedLoops_combined)),
           params = apaParams) 

plotText(label = depths$hek_wt,
         fontcolor = "black",
         fontsize = 4,
         x = 1.70 - 0.05,
         y = 4.65 - 0.05,
         just = c("right", "bottom"))


# HCT section below -------------------------------------------------------

plotMatrix(normalizedAPAs$hct_gainedLoops_contAPA_normalized.rds,
           x = 1.9,
           y = 3.1,
           zrange = c(0, max(hct_gainedLoops_combined)),
           params = apaParams) |>
  annoHeatmapLegend(x = 2.7,
                    y = 3.1,
                    width = 0.05,
                    height = 0.75,
                    fontcolor = "black", fontsize = 5)

plotMatrix(normalizedAPAs$hct_gainedLoops_sorbAPA_normalized.rds,
           x = 1.9,
           y = 3.9,
           zrange = c(0, max(hct_gainedLoops_combined)),
           params = apaParams)

plotText(label = depths$hct116,
         fontcolor = "black",
         fontsize = 4,
         x = 2.65 - 0.05,
         y = 4.65 - 0.05,
         just = c("right", "bottom"))


# HCT Section above -------------------------------------------------------

plotMatrix(normalizedAPAs$amat_gainedLoops_contAPA_normalized.rds,
           x = 2.85,
           y = 3.1,
           zrange = c(0, max(amat_gainedLoops_combined)),
           params = apaParams) |> 
  annoHeatmapLegend(x = 3.65,
                    y = 3.1,
                    width = 0.05,
                    height = 0.75,
                    fontcolor = "black",
                    fontsize = 5)

plotMatrix(normalizedAPAs$amat_gainedLoops_sorbAPA_normalized.rds,
           x = 2.85,
           y = 3.9,
           zrange = c(0, max(amat_gainedLoops_combined)),
           params = apaParams)

plotText(label = depths$t47d,
         fontcolor = "black",
         fontsize = 4,
         x = 3.60 - 0.05,
         y = 4.65 - 0.05,
         just = c("right", "bottom"))

plotText(label = "HEK293 eGFP-YAP", fontcolor = "black",
         x = 0.05,
         y = 4.7,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "HEK293 WT", fontcolor = "black",
         x = 1.125,
         y = 4.7,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "HCT116", fontcolor = "black",
         x = 2.15,
         y = 4.7,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "T47D (Amat et al 2019)", fontcolor = "black",
         x = 2.85,
         y = 4.7,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "untreated", fontcolor = "black",
         x = 0.025,
         y = 3.125,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "+ sorbitol", fontcolor = "black",
         x = 0.025,
         y = 3.925,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "untreated", fontcolor = "black",
         x = 0.95 + 0.025,
         y = 3.125,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "+ sorbitol", fontcolor = "black",
         x = 0.95 + 0.025,
         y = 3.925,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "untreated", fontcolor = "black",
         x = 1.925,
         y = 3.125,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "+ sorbitol", fontcolor = "black",
         x = 1.925,
         y = 3.925,
         params = apaParams,
         just = c ("top", "left"))

plotText(label = "untreated", fontcolor = "black",
         x = 2.875,
         y = 3.125,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "+ NaCl", fontcolor = "black",
         x = 2.875,
         y = 3.925,
         params = apaParams,
         just = c("top", "left"))

# Visualize normalized lost APAs ------------------------------------------

## combine APA matrices to pull out the max value for zrange max
wt_lostLoops_combined <- c(normalizedAPAs$wt_lostLoops_contAPA_normalized.rds,
                           normalizedAPAs$wt_lostLoops_sorbAPA_normalized.rds)

oe_lostLoops_combined <- c(normalizedAPAs$oe_lostLoops_contAPA_normalized.rds,
                           normalizedAPAs$oe_lostLoops_sorbAPA_normalized.rds)

hct_lostLoops_combined <- c(normalizedAPAs$hct_lostLoops_contAPA_normalized.rds,
                            normalizedAPAs$hct_lostLoops_sorbAPA_normalized.rds)

amat_lostLoops_combined <- c(normalizedAPAs$amat_lostLoops_contAPA_normalized.rds,
                             normalizedAPAs$amat_lostLoops_sorbAPA_normalized.rds)

# Visualize APAs ----------------------------------------------------------

plotMatrix(normalizedAPAs$oe_lostLoops_contAPA_normalized.rds,
           x = 3.8,
           y = 3.1,
           zrange = c(0, max(oe_lostLoops_combined)),
           params = apaParams) |>
  annoHeatmapLegend(x = 4.6,
                    y = 3.1,
                    width = 0.05,
                    height = 0.75,
                    fontcolor = "black",
                    fontsize = 5)

plotMatrix(normalizedAPAs$oe_lostLoops_sorbAPA_normalized.rds,
           x = 3.8,
           y = 3.9,
           zrange = c(0, max(oe_lostLoops_combined)),
           params = apaParams)

plotText(label = depths$hek_egfp,
         fontcolor = "black",
         fontsize = 4,
         x = 4.55 - 0.05,
         y = 4.65 - 0.05,
         just = c("right", "bottom"))

plotMatrix(normalizedAPAs$wt_lostLoops_contAPA_normalized.rds,
           x = 4.75,
           y = 3.1,
           zrange = c(0, max(wt_lostLoops_combined)),
           params = apaParams)  |>
  annoHeatmapLegend(x = 5.55,
                    y = 3.1,
                    width = 0.05,
                    height = 0.75,
                    fontcolor = "black",
                    fontsize = 5)

plotMatrix(normalizedAPAs$wt_lostLoops_sorbAPA_normalized.rds,
           x = 4.75,
           y = 3.9,
           zrange = c(0, max(wt_lostLoops_combined)),
           params = apaParams)

plotText(label = depths$hek_wt,
         fontcolor = "black",
         fontsize = 4,
         x = 5.50 - 0.05,
         y = 4.65 - 0.05,
         just = c("right", "bottom"))


# HCT section below -------------------------------------------------------

plotMatrix(normalizedAPAs$hct_lostLoops_contAPA_normalized.rds,
           x = 5.7,
           y = 3.1,
           zrange = c(0, max(hct_lostLoops_combined)),
           params = apaParams) |>
  annoHeatmapLegend(x = 6.5,
                    y = 3.1,
                    width = 0.05,
                    height = 0.75,
                    fontcolor = "black",
                    fontsize = 5,
                    params = apaParams)

plotMatrix(normalizedAPAs$hct_lostLoops_sorbAPA_normalized.rds,
           x = 5.7,
           y = 3.9,
           zrange = c(0, max(hct_lostLoops_combined)),
           params = apaParams)

plotText(label = depths$hct116,
         fontcolor = "black",
         fontsize = 4,
         x = 6.45 - 0.05,
         y = 4.65 - 0.05,
         just = c("right", "bottom"))


# HCT section above -------------------------------------------------------

plotMatrix(normalizedAPAs$amat_lostLoops_contAPA_normalized.rds,
           x = 6.65,
           y = 3.1,
           zrange = c(0, max(amat_lostLoops_combined)),
           params = apaParams) |>
  annoHeatmapLegend(x = 7.45,
                    y = 3.1,
                    width = 0.05,
                    height = 0.75,
                    fontcolor = "black",
                    fontsize = 5,
                    params = apaParams)

plotMatrix(normalizedAPAs$amat_lostLoops_sorbAPA_normalized.rds,
           x = 6.65,
           y = 3.9,
           zrange = c(0, max(amat_lostLoops_combined)),
           params = apaParams) 

plotText(label = depths$t47d,
         fontcolor = "black",
         fontsize = 4,
         x = 7.40 - 0.05,
         y = 4.65 - 0.05,
         just = c("right", "bottom"))

plotText(label = "HEK293 eGFP-YAP", fontcolor = "black",
         x = 3.85,
         y = 4.7,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "HEK293 WT", fontcolor = "black",
         x = 4.9,
         y = 4.7,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "HCT116", fontcolor = "black",
         x = 5.9,
         y = 4.7,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "T47D (Amat et al 2019)", fontcolor = "black",
         x = 6.65,
         y = 4.7,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "untreated", fontcolor = "black",
         x = 3.825,
         y = 3.125,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "+ sorbitol", fontcolor = "black",
         x = 3.825,
         y = 3.925,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "untreated", fontcolor = "black",
         x = 4.75 + 0.025,
         y = 3.125,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "+ sorbitol", fontcolor = "black",
         x = 4.75 + 0.025,
         y = 3.925,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "untreated", fontcolor = "black",
         x = 5.725,
         y = 3.125,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "+ sorbitol", fontcolor = "black",
         x = 5.725,
         y = 3.925,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "untreated", fontcolor = "black",
         x = 6.675,
         y = 3.125,
         params = apaParams,
         just = c("top", "left"))

plotText(label = "+ NaCl", fontcolor = "black",
         x = 6.675,
         y = 3.925,
         params = apaParams,
         just = c("top", "left"))

plotSegments(
  x0 = 0,
  y0 = 3,
  x1 = 3.75,
  y1 = 3,
  default.units = "inches",
  linecolor = "#F8766D",
  lwd = 1,
  lty = 1)

plotSegments(
  x0 = 3.8,
  y0 = 3,
  x1 = 7.6,
  y1 = 3,
  default.units = "inches",
  linecolor = "#619CFF",
  lwd = 1,
  lty = 1)

plotText(label = "gained loops",
         x = 1.8,
         y = 2.9,
         fontcolor = "#F8766D",
         fontsize = 10)

plotText(label = "lost loops",
         x = 5.55,
         y = 2.9,
         fontcolor = "#619CFF",
         fontsize = 10)

plotText(label = "A",
         x = 0,
         y = 0.125,
         fontface = "bold",
         fontsize = 12)

plotText(label = "B",
         x = 5,
         y = 0.125,
         fontface = "bold",
         fontsize = 12)

plotText(label = "C",
         x = 0,
         y = 2.8,
         fontface = "bold",
         fontsize = 12)

dev.off()