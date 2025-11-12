## Lost differential loop visualizations (HiC rectangles)

# Load packages / data ----------------------------------------------------

library(plotgardener)
library(InteractionSet)
library(tidyverse)
library(glue)
library(dbscan)
library(data.table)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(mariner)
library(plyranges)
library(RColorBrewer)

diff_loopCounts <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions() 

## calculate loop size
mcols(diff_loopCounts)$loop_size <- pairdist(diff_loopCounts)

# Setting static and lost loops -------------------------------------------

## lost loops with a p-adj. value of < 0.1 and a (-) log2FC
lost_adj <- 
  diff_loopCounts[which(diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange < 0)]

## filter for the 100 best lost loops

## based off lowest padj value
bestLost <- head(lost_adj[order(lost_adj$padj, decreasing = F)], 100)

# top 100 lost loops
loopRegions_lost <-
  GRanges(seqnames = as.character(seqnames(anchors(x = bestLost, "first"))),
          ranges = IRanges(start = start(anchors(bestLost, "first")),
                           end = end(anchors(bestLost, "second"))),
          mcols = mcols(bestLost))

## Expand regions by buffer
buffer <- 200e3
loopRegions_lost_buffed <- loopRegions_lost + buffer

# load loop lists --------------------------------------------------------

lostLoops <- diff_loopCounts[which(diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange < 0)] 

gainedLoops <- diff_loopCounts[which(diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange >0)]

# Create Survey Plots -----------------------------------------------------
RColorBrewer::brewer.pal(n = 5,name = "YlGnBu")
# "#A1DAB4" "#41B6C4" "#2C7FB8" "#253494"

##make pdf
pdf(file = "plots/surveyPlot_lostLoops_10kb.pdf",
    width = 5.75,
    height = 5.5)

## Loop through each region
for(i in seq_along(loopRegions_lost_buffed)){
  
  ## Define parameters
  p <- pgParams(assembly = "hg38",
                resolution = 10e3,
                chrom = as.character(seqnames(loopRegions_lost_buffed))[i],
                chromstart = start(loopRegions_lost_buffed)[i],
                chromend = end(loopRegions_lost_buffed)[i],
                zrange = c(0,100),
                norm = "SCALE",
                x = 0.25,
                width = 5,
                length = 5,
                height = 2)
  
  
  # Begin Visualization -----------------------------------------------------
  ## Make page
  pageCreate(width = 5.75, height = 5.5,
             xgrid = 0, ygrid = 0, showGuides = F)
  
  ## Plot top Hi-C rectangle + lost loops
  
  control <- plotHicRectangle(data = "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic",
                              params = p,
                              y = 0.5)
  
  annoHeatmapLegend(control, orientation = "v",
                    fontsize = 8,
                    fontcolor = "black",
                    digits = 2,
                    x = 5.5,
                    y = 0.5,
                    width = 0.1,
                    height = 1.5,
                    just = c("left", "top"),
                    default.units = "inches")
  
  annoPixels(control,
             data = lostLoops,
             shift = 0.5,
             type = "arrow",
             col = "#005AB5")
  
  annoPixels(control,
             data = gainedLoops,
             shift = 0.5,
             type = "arrow",
             col = "#DC3220")
  
  ## Plot bottom Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls
  
  sorb <- 
    plotHicRectangle(data = "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic", 
                     params = p,
                     y = 2.6)
  
  annoHeatmapLegend(sorb, orientation = "v",
                    fontsize = 8,
                    fontcolor = "black",
                    digits = 2,
                    x = 5.5,
                    y = 2.6,
                    width = 0.1,
                    height = 1.5,
                    just = c("left", "top"),
                    default.units = "inches")
  
  annoPixels(sorb,
             data = lostLoops,
             shift = 0.5,
             type = "arrow",
             col = "#005AB5")
  
  annoPixels(sorb,
             data = gainedLoops,
             shift = 0.5,
             type = "arrow",
             col = "#DC3220")
  
  ## Plot genes
  plotGenes(param = p,
            chrom = p$chrom,
            x = 0.25,
            y = 4.7,
            height = 0.5)
  
  
  ## Plot genome label
  plotGenomeLabel(params = p,
                  x = 0.25,
                  y = 5.3)
  
  # Annotate Hi-C rectangles by treatment ------------------------------------
  
  plotText(label = "untreated",
           x = 0.25,
           y = 0.5,
           just = c("top", "left"))
  
  plotText(label = "+ sorbitol",
           x = 0.25,
           y = 2.6,
           just = c("top", "left"))
  
  annoHighlight(
    plot = control,
    fill = "lightgrey",
    chrom = as.character(seqnames(loopRegions_lost))[i],
    chromstart = start(loopRegions_lost)[i],
    chromend = start(loopRegions_lost)[i] + 10e3,
    default.units = "inches",
    y = 0.5, x = 0.25, height = 5.3, just = c("left", "top")
  )
  
  
  annoHighlight(
    plot = sorb,
    y = 0.5, x = 0.25, height = 5.3, just = c("left", "top"),
    fill = "lightgrey",
    chrom = as.character(seqnames(loopRegions_lost))[i],
    chromstart = end(loopRegions_lost)[i],
    chromend = end(loopRegions_lost)[i] - 10e3,
    default.units = "inches"
  )
  
}
dev.off()