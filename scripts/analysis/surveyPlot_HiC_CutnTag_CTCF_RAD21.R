## Gained differential loop visualizations (HiC rectangles) with Cut&Tag peaks

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
library(GenomicRanges)

# Data --------------------------------------------------------------------
diff_loopCounts <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions() 

## calculate loop size
mcols(diff_loopCounts)$loop_size <- pairdist(diff_loopCounts)

# Setting static and lost loops -------------------------------------------

## gained loops with a p-adj. value of < 0.1 and a (+) log2FC
gained_adj <- 
  diff_loopCounts[which(diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange > 0)]

## filter for the 100 best gained loops

## based off lowest padj value
bestGained <- head(gained_adj[order(gained_adj$padj, decreasing = F)], 100)

# top 100 gained loops
loopRegions_gained <-
  GRanges(seqnames = as.character(seqnames(anchors(x = bestGained, "first"))),
          ranges = IRanges(start = start(anchors(bestGained, "first")),
                           end = end(anchors(bestGained, "second"))),
          mcols = mcols(bestGained))

## Expand regions by buffer
buffer <- 200e3
loopRegions_gained_buffed <- loopRegions_gained + buffer

## Specify loop type

lostLoops <- diff_loopCounts[which(diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange < 0)] 

gainedLoops <- diff_loopCounts[which(diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange > 0)]

## Load in cutntag files to calculate signal range
cutntag <- list.files(path = "data/processed/cutntag/output/mergeSignal/",
                      full.names = T)

##make pdf
pdf(file = "plots/surveyPlot_HiC_CutnTag_CTCF_RAD21.pdf",
    width = 5.75,
    height = 7)

## Loop through each region
for(i in seq_along(loopRegions_gained_buffed)){
  
  ## Define parameters
  p <- pgParams(assembly = "hg38",
                resolution = 10e3,
                chrom = as.character(seqnames(loopRegions_gained_buffed))[i],
                chromstart = start(loopRegions_gained_buffed)[i],
                chromend = end(loopRegions_gained_buffed)[i],
                zrange = c(0,100),
                norm = "SCALE",
                x = 0.25,
                width = 5,
                length = 5,
                height = 2)
  
  
  # Begin Visualization -----------------------------------------------------
  ## Make page
  pageCreate(width = 5.75, height = 7,
             showGuides = F)
  
  ## Plot top Hi-C rectangle + lost loops
  
  control <- plotHicRectangle(data = "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic",
                              params = p,
                              y = 0.25)
  
  plotText(label = "untreated",
           x = 0.25,
           y = 0.25,
           just = c("top", "left"))
  
  annoHeatmapLegend(control, orientation = "v",
                    fontsize = 8,
                    fontcolor = "black",
                    digits = 2,
                    x = 5.5,
                    y = 0.25,
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
                     y = 2.35)
  
  plotText(label = "+ sorbitol",
           x = 0.25,
           y = 2.35,
           just = c("top", "left"))
  
  annoHeatmapLegend(sorb, orientation = "v",
                    fontsize = 8,
                    fontcolor = "black",
                    digits = 2,
                    x = 5.5,
                    y = 2.35,
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
  
  ## Calculate Signal Range for CTCF Signal tracks
  ctcf_signalRange <- calcSignalRange(data = str_subset(cutntag, "CTCF"),
                                      chrom = p$chrom,
                                      chromstart = p$chromstart,
                                      chromend = p$chromend,
                                      assembly = "hg38",
                                      negData = F)
  
  ## Plot CTCF control signal
  plotSignal("data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_CTCF_cont_0h.bw",
             params = p,
             x = 0.25,
             y = 4.45,
             height = 0.25,
             linecolor = "#7FCDBB",
             fill = "#7FCDBB",
             scale = T,
             range = ctcf_signalRange
  )
  
  plotText("CTCF",
           fontcolor = "#7FCDBB",
           rot = 90,
           y = 4.85,
           x = 0.1)
  
  ## Plot CTCF sorbitol signal
  plotSignal("data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h.bw",
             params = p,
             x = 0.25,
             y = 4.95,
             height = 0.25,
             linecolor = "#7FCDBB",
             fill = "#7FCDBB",
             scale = T,
             range = ctcf_signalRange
  )
  
  ## Calculate Signal Range for RAD21 Signal tracks
  rad21_signalRange <- calcSignalRange(data = str_subset(cutntag, "RAD21"),
                                       chrom = p$chrom,
                                       chromstart = p$chromstart,
                                       chromend = p$chromend,
                                       assembly = "hg38",
                                       negData = F)
  
  ## Plot RAD21 control signal
  plotSignal("data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_RAD21_cont_0h.bw",
             params = p,
             x = 0.25,
             y = 5.35,
             height = 0.25,
             linecolor = "#41B6C4",
             fill = "#41B6C4",
             scale = T,
             range = rad21_signalRange
  )
  
  plotText("RAD21",
           fontcolor = "#41B6C4",
           rot = 90,
           y = 5.75,
           x = 0.1)
  
  ## Plot RAD21 sorbitol signal
  plotSignal("data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_RAD21_sorbitol_1h.bw",
             params = p,
             x = 0.25,
             y = 5.85,
             height = 0.25,
             linecolor = "#41B6C4",
             fill = "#41B6C4",
             scale = T,
             range = rad21_signalRange
  )
  
  ## Plot genes
  plotGenes(param = p,
            chrom = p$chrom,
            x = 0.25,
            y = 6.15,
            height = 0.5)
  
  
  ## Plot genome label
  plotGenomeLabel(params = p,
                  x = 0.25,
                  y = 6.75)
  
  annoHighlight(
    plot = control,
    fill = "lightgrey",
    chrom = as.character(seqnames(loopRegions_gained))[i],
    chromstart = start(loopRegions_gained)[i],
    chromend = start(loopRegions_gained)[i] + 10e3,
    default.units = "inches",
    y = 0.25, x = 0.25, height = 6.5, just = c("left", "top")
  )
  
  
  annoHighlight(
    plot = sorb,
    y = 0.25, x = 0.25, height = 6.5, just = c("left", "top"),
    fill = "lightgrey",
    chrom = as.character(seqnames(loopRegions_gained))[i],
    chromstart = end(loopRegions_gained)[i],
    chromend = end(loopRegions_gained)[i] - 10e3,
    default.units = "inches"
  )
  
  
}
dev.off()
