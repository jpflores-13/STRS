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

# Function to read narrowPeak files ---------------------------------------
read_narrowpeaks <- function(file) {
  peaks <- read.table(file, sep = "\t", header = FALSE,
                      col.names = c("seqnames", "start", "end", "name", 
                                    "score", "strand", "signalValue", 
                                    "pValue", "qValue", "peak"))
  
  # Convert to GRanges
  GRanges(seqnames = peaks$seqnames,
          ranges = IRanges(start = peaks$start + 1, end = peaks$end), # Convert to 1-based
          strand = "*",
          score = peaks$score,
          signalValue = peaks$signalValue,
          pValue = peaks$pValue,
          qValue = peaks$qValue,
          peak = peaks$peak)
}

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

## Load cutntag narrowPeak data (no filtering initially)
cutntag_np <- list.files("data/processed/cutntag/output/peaks/",
                             full.names = T,
                             pattern = ".narrowPeak") |> 
  lapply(read_narrowpeaks)

target <- c("CTCF", "H3K27ac", "RAD21", "YAP1")
condition <- c("control", "sorbitol")
names(cutntag_np) <- paste0(rep(target, each = 2), "_", condition)

##make pdf
pdf(file = "plots/surveyPlot_HiC_CutnTag.pdf",
    width = 5.75,
    height = 9.75)

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
  pageCreate(width = 5.75, height = 9.75,
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
  
  ## Calculate Signal Range for H3K27ac Signal tracks
  h3k27ac_signalRange <- calcSignalRange(data = str_subset(cutntag, "H3K27ac"),
                                         chrom = p$chrom,
                                         chromstart = p$chromstart,
                                         chromend = p$chromend,
                                         assembly = "hg38",
                                         negData = F)
  
  ## Plot H3K27ac control peak ranges ABOVE signal
  plotRanges(
    data = cutntag_np$H3K27ac_control, 
    params = p,
    x = 0.25, y = 4.35, width = 5, height = 0.08,
    fill = "#C7E9B4",
    colorby = "signalValue",
    palette = colorRampPalette(c("white", "#C7E9B4"))(100),
    stroke = 0.05,
    strokecolor = "black",
    collapse = TRUE,
    just = c("left", "top"), default.units = "inches"
  )
  
  ## Plot H3K27ac control signal
  plotSignal("data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_H3K27ac_cont_0h.bw",
             params = p,
             x = 0.25,
             y = 4.4,
             height = 0.25,
             linecolor = "#C7E9B4",
             fill = "#C7E9B4",
             scale = T,
             range = h3k27ac_signalRange
  )
  
  plotText("H3K27ac",
           fontcolor = "#C7E9B4",
           rot = 90,
           y = 4.9,
           x = 0.1)
  
  ## Plot H3K27ac sorbitol peak ranges ABOVE signal
  plotRanges(
    data = cutntag_np$H3K27ac_sorbitol, 
    params = p,
    x = 0.25, y = 4.95, width = 5, height = 0.08,
    fill = "#C7E9B4",
    colorby = "signalValue",
    palette = colorRampPalette(c("white", "#C7E9B4"))(100),
    stroke = 0.05,
    strokecolor = "black",
    alpha = 0.7,
    collapse = TRUE,
    just = c("left", "top"), default.units = "inches"
  )
  
  ## Plot H3K27ac sorbitol signal
  plotSignal("data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_H3K27ac_sorbitol_1h.bw",
             params = p,
             x = 0.25,
             y = 5,
             height = 0.25,
             linecolor = "#C7E9B4",
             fill = "#C7E9B4",
             scale = T,
             range = h3k27ac_signalRange
  )
  
  ## Calculate Signal Range for CTCF Signal tracks
  ctcf_signalRange <- calcSignalRange(data = str_subset(cutntag, "CTCF"),
                                      chrom = p$chrom,
                                      chromstart = p$chromstart,
                                      chromend = p$chromend,
                                      assembly = "hg38",
                                      negData = F)
  
  ## Plot CTCF control peak ranges ABOVE signal
  plotRanges(
    data = cutntag_np$CTCF_control, 
    params = p,
    x = 0.25, y = 5.55, width = 5, height = 0.08,
    fill = "#7FCDBB",
    colorby = "signalValue",
    palette = colorRampPalette(c("white", "#7FCDBB"))(100),
    stroke = 0.05,
    strokecolor = "black",
    collapse = TRUE,
    just = c("left", "top"), default.units = "inches"
  )
  
  ## Plot CTCF control signal
  plotSignal("data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_CTCF_cont_0h.bw",
             params = p,
             x = 0.25,
             y = 5.6,
             height = 0.25,
             linecolor = "#7FCDBB",
             fill = "#7FCDBB",
             scale = T,
             range = ctcf_signalRange
  )
  
  plotText("CTCF",
           fontcolor = "#7FCDBB",
           rot = 90,
           y = 6,
           x = 0.1)
  
  ## Plot CTCF sorbitol peak ranges ABOVE signal
  plotRanges(
    data = cutntag_np$CTCF_sorbitol, 
    params = p,
    x = 0.25, y = 6.15, width = 5, height = 0.08,
    fill = "#7FCDBB",
    colorby = "signalValue",
    palette = colorRampPalette(c("white", "#7FCDBB"))(100),
    stroke = 0.05,
    strokecolor = "black",
    alpha = 0.7,
    collapse = TRUE,
    just = c("left", "top"), default.units = "inches"
  )
  
  ## Plot CTCF sorbitol signal
  plotSignal("data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h.bw",
             params = p,
             x = 0.25,
             y = 6.2,
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
  
  ## Plot RAD21 control peak ranges ABOVE signal
  plotRanges(
    data = cutntag_np$RAD21_control, 
    params = p,
    x = 0.25, y = 6.75, width = 5, height = 0.08,
    fill = "#41B6C4",
    colorby = "signalValue",
    palette = colorRampPalette(c("white", "#41B6C4"))(100),
    stroke = 0.05,
    strokecolor = "black",
    collapse = TRUE,
    just = c("left", "top"), default.units = "inches"
  )
  
  ## Plot RAD21 control signal
  plotSignal("data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_RAD21_cont_0h.bw",
             params = p,
             x = 0.25,
             y = 6.8,
             height = 0.25,
             linecolor = "#41B6C4",
             fill = "#41B6C4",
             scale = T,
             range = rad21_signalRange
  )
  
  plotText("RAD21",
           fontcolor = "#41B6C4",
           rot = 90,
           y = 7.2,
           x = 0.1)
  
  ## Plot RAD21 sorbitol peak ranges ABOVE signal
  plotRanges(
    data = cutntag_np$RAD21_sorbitol, 
    params = p,
    x = 0.25, y = 7.35, width = 5, height = 0.08,
    fill = "#41B6C4",
    colorby = "signalValue",
    palette = colorRampPalette(c("white", "#41B6C4"))(100),
    stroke = 0.05,
    strokecolor = "black",
    alpha = 0.7,
    collapse = TRUE,
    just = c("left", "top"), default.units = "inches"
  )
  
  ## Plot RAD21 sorbitol signal
  plotSignal("data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_RAD21_sorbitol_1h.bw",
             params = p,
             x = 0.25,
             y = 7.4,
             height = 0.25,
             linecolor = "#41B6C4",
             fill = "#41B6C4",
             scale = T,
             range = rad21_signalRange
  )
  
  ## Calculate Signal Range for YAP1 Signal tracks
  yap1_signalRange <- calcSignalRange(data = str_subset(cutntag, "YAP1"),
                                      chrom = p$chrom,
                                      chromstart = p$chromstart,
                                      chromend = p$chromend,
                                      assembly = "hg38",
                                      negData = F)
  
  ## Plot YAP1 control peak ranges ABOVE signal
  plotRanges(
    data = cutntag_np$YAP1_control, 
    params = p,
    x = 0.25, y = 7.95, width = 5, height = 0.08,
    fill = "#1D91C0",
    colorby = "signalValue",
    palette = colorRampPalette(c("white", "#1D91C0"))(100),
    stroke = 0.05,
    strokecolor = "black",
    collapse = TRUE,
    just = c("left", "top"), default.units = "inches"
  )
  
  ## Plot YAP1 control signal
  plotSignal("data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_YAP1_cont_0h.bw",
             params = p,
             x = 0.25,
             y = 8,
             height = 0.25,
             linecolor = "#1D91C0",
             fill = "#1D91C0",
             scale = T,
             range = yap1_signalRange
  )
  
  plotText("YAP1",
           fontcolor = "#1D91C0",
           rot = 90,
           y = 8.5,
           x = 0.1)
  
  ## Plot YAP1 sorbitol peak ranges ABOVE signal
  plotRanges(
    data = cutntag_np$YAP1_sorbitol, 
    params = p,
    x = 0.25, y = 8.55, width = 5, height = 0.08,
    fill = "#1D91C0",
    colorby = "signalValue",
    palette = colorRampPalette(c("white", "#1D91C0"))(100),
    stroke = 0.05,
    strokecolor = "black",
    alpha = 0.7,
    collapse = TRUE,
    just = c("left", "top"), default.units = "inches"
  )
  
  ## Plot YAP1 sorbitol signal
  plotSignal("data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_YAP1_sorbitol_1h.bw",
             params = p,
             x = 0.25,
             y = 8.6,
             height = 0.25,
             linecolor = "#1D91C0",
             fill = "#1D91C0",
             scale = T,
             range = yap1_signalRange
  )
  
  ## Plot genes
  plotGenes(param = p,
            chrom = p$chrom,
            x = 0.25,
            y = 8.9,
            height = 0.5)
  
  
  ## Plot genome label
  plotGenomeLabel(params = p,
                  x = 0.25,
                  y = 9.5)
  
  annoHighlight(
    plot = control,
    fill = "lightgrey",
    chrom = as.character(seqnames(loopRegions_gained))[i],
    chromstart = start(loopRegions_gained)[i],
    chromend = start(loopRegions_gained)[i] + 10e3,
    default.units = "inches",
    y = 0.25, x = 0.25, height = 9.25, just = c("left", "top")
  )
  
  
  annoHighlight(
    plot = sorb,
    y = 0.25, x = 0.25, height = 9.25, just = c("left", "top"),
    fill = "lightgrey",
    chrom = as.character(seqnames(loopRegions_gained))[i],
    chromstart = end(loopRegions_gained)[i],
    chromend = end(loopRegions_gained)[i] - 10e3,
    default.units = "inches"
  )
  
  
}
dev.off()
