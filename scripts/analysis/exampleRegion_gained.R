# ##############################################################################
# filename:    exampleRegion_gained.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Single-region Hi-C rectangle plot centered on a top gained loop;
#              control and sorbitol maps shown with gained/lost loop annotations
# ##############################################################################

# Libraries ----
library(plotgardener)
library(InteractionSet)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(mariner)
library(plyranges)
library(RColorBrewer)

# Parameters ----
diff_loops_rds  <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
hic_control     <- "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic"
hic_sorbitol    <- "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic"
output_pdf      <- "plots/exampleRegion_gained.pdf"
region_index    <- 38
buffer          <- 200e3
page_width      <- 5.75
page_height     <- 5.75

# Data import ----
diff_loopCounts <- readRDS(diff_loops_rds) |>
  interactions()

mcols(diff_loopCounts)$loop_size <- pairdist(diff_loopCounts)

# Analysis ----
gained_adj <- diff_loopCounts[
  diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange > 0
]

bestGained <- head(gained_adj[order(gained_adj$padj)], 100)

loopRegions_gained <- GRanges(
  seqnames = as.character(seqnames(anchors(x = bestGained, "first"))),
  ranges   = IRanges(start = start(anchors(bestGained, "first")),
                     end   = end(anchors(bestGained, "second"))),
  mcols    = mcols(bestGained)
)

loopRegions_gained_buffed <- loopRegions_gained + buffer

lostLoops   <- diff_loopCounts[diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange < 0]
gainedLoops <- diff_loopCounts[diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange > 0]

# Visualization ----
p <- pgParams(
  assembly   = "hg38",
  resolution = 10e3,
  chrom      = as.character(seqnames(loopRegions_gained_buffed))[region_index],
  chromstart = start(loopRegions_gained_buffed)[region_index],
  chromend   = end(loopRegions_gained_buffed)[region_index],
  zrange     = c(0, 100),
  norm       = "SCALE",
  x          = 0.25,
  width      = 5,
  length     = 5,
  height     = 2
)

# Save outputs ----
pdf(file = output_pdf, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height,
           xgrid = 0, ygrid = 0, showGuides = FALSE)

control <- plotHicRectangle(data = hic_control, params = p, y = 0.5)

annoHeatmapLegend(control, orientation = "v",
                  fontsize = 8, fontcolor = "black", digits = 2,
                  x = 5.5, y = 0.5, width = 0.1, height = 1.5,
                  just = c("left", "top"), default.units = "inches")

annoPixels(control, data = lostLoops,   shift = 0.5, type = "arrow", col = "#005AB5")
annoPixels(control, data = gainedLoops, shift = 0.5, type = "arrow", col = "#DC3220")

sorb <- plotHicRectangle(data = hic_sorbitol, params = p, y = 2.6)

annoPixels(sorb, data = lostLoops,   shift = 0.5, type = "arrow", col = "#005AB5")
annoPixels(sorb, data = gainedLoops, shift = 0.5, type = "arrow", col = "#DC3220")

plotGenes(param = p, chrom = p$chrom, x = 0.25, y = 4.7, height = 0.5)
plotGenomeLabel(params = p, x = 0.25, y = 5.3)

plotText(label = "untreated",      x = 0.25, y = 0.5, just = c("top", "left"))
plotText(label = "200mM sorbitol", x = 0.25, y = 2.6, just = c("top", "left"))

annoHighlight(
  plot = control, fill = "lightgrey",
  chrom      = as.character(seqnames(loopRegions_gained))[region_index],
  chromstart = start(loopRegions_gained)[region_index],
  chromend   = start(loopRegions_gained)[region_index] + 10e3,
  default.units = "inches",
  y = 0.5, x = 0.25, height = 4.8, just = c("left", "top")
)

annoHighlight(
  plot = sorb, fill = "lightgrey",
  chrom      = as.character(seqnames(loopRegions_gained))[region_index],
  chromstart = end(loopRegions_gained)[region_index],
  chromend   = end(loopRegions_gained)[region_index] - 10e3,
  default.units = "inches",
  y = 0.5, x = 0.25, height = 4.8, just = c("left", "top")
)

dev.off()

sessionInfo()
