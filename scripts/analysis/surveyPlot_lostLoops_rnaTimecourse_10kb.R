# ##############################################################################
# filename:    surveyPlot_lostLoops_rnaTimecourse_10kb.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Survey plots for top 100 lost loops with RNA-seq timecourse
#              signal tracks (0h–24h, YlGnBu palette); two Hi-C rectangles
#              plus seven stranded signal tracks and anchor highlights
# ##############################################################################

# Libraries ----
library(plotgardener)
library(InteractionSet)
library(mariner)

# Parameters ----
diff_loops_rds <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
hic_control    <- "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic"
hic_sorbitol   <- "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic"
rna_signal_dir <- "data/processed/rna/timecourse/output/mergeSignal/stranded/"
output_pdf     <- "plots/surveyPlot_lostLoops_rnaTimecourse_10kb.pdf"
n_loops        <- 100L
buffer         <- 200e3
page_width     <- 5.75
page_height    <- 9.1

tp_colors <- c("0h"  = "#C7E9B4", "1h"  = "#7FCDBB", "3h"  = "#41B6C4",
               "6h"  = "#1D91C0", "9h"  = "#225EA8", "12h" = "#253494",
               "24h" = "#081D58")

# Data import ----
diff_loopCounts <- readRDS(diff_loops_rds) |> interactions()
mcols(diff_loopCounts)$loop_size <- pairdist(diff_loopCounts)

rna <- list.files(rna_signal_dir, full.names = TRUE)

# Analysis ----
lost_adj  <- diff_loopCounts[diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange < 0]
bestLost  <- head(lost_adj[order(lost_adj$padj)], n_loops)

loopRegions_lost <- GRanges(
  seqnames = as.character(seqnames(anchors(bestLost, "first"))),
  ranges   = IRanges(start = start(anchors(bestLost, "first")),
                     end   = end(anchors(bestLost, "second"))),
  mcols    = mcols(bestLost)
)
loopRegions_lost_buffed <- loopRegions_lost + buffer

lostLoops   <- diff_loopCounts[diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange < 0]
gainedLoops <- diff_loopCounts[diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange > 0]

# Visualization ----
pdf(file = output_pdf, width = page_width, height = page_height)

for (i in seq_along(loopRegions_lost_buffed)) {

  p <- pgParams(assembly   = "hg38",
                resolution = 10e3,
                chrom      = as.character(seqnames(loopRegions_lost_buffed))[i],
                chromstart = start(loopRegions_lost_buffed)[i],
                chromend   = end(loopRegions_lost_buffed)[i],
                zrange     = c(0, 100),
                norm       = "SCALE",
                x          = 0.25,
                width      = 5,
                length     = 5,
                height     = 2)

  pageCreate(width = page_width, height = page_height, showGuides = FALSE)

  ## Hi-C panels ----
  control <- plotHicRectangle(data = hic_control, params = p, y = 0.25)
  plotText(label = "untreated", x = 0.25, y = 0.25, just = c("top", "left"))
  annoHeatmapLegend(control, orientation = "v",
                    fontsize = 8, fontcolor = "black", digits = 2,
                    x = 5.5, y = 0.25, width = 0.1, height = 1.5,
                    just = c("left", "top"), default.units = "inches")
  annoPixels(control, data = lostLoops,   shift = 0.5, type = "arrow", col = "#005AB5")
  annoPixels(control, data = gainedLoops, shift = 0.5, type = "arrow", col = "#DC3220")

  sorb <- plotHicRectangle(data = hic_sorbitol, params = p, y = 2.35)
  plotText(label = "+ sorbitol", x = 0.25, y = 2.35, just = c("top", "left"))
  annoHeatmapLegend(sorb, orientation = "v",
                    fontsize = 8, fontcolor = "black", digits = 2,
                    x = 5.5, y = 2.35, width = 0.1, height = 1.5,
                    just = c("left", "top"), default.units = "inches")
  annoPixels(sorb, data = lostLoops,   shift = 0.5, type = "arrow", col = "#005AB5")
  annoPixels(sorb, data = gainedLoops, shift = 0.5, type = "arrow", col = "#DC3220")

  ## RNA-seq signal tracks ----
  signalRange <- calcSignalRange(data = rna,
                                 chrom      = as.character(seqnames(loopRegions_lost_buffed))[i],
                                 chromstart = start(loopRegions_lost_buffed)[i],
                                 chromend   = end(loopRegions_lost_buffed)[i],
                                 assembly   = "hg38",
                                 negData    = FALSE)

  tp_y      <- c(4.4, 5.0, 5.6, 6.2, 6.8, 7.4, 8.0)
  tp_labels <- c("0h", "1h", "3h", "6h", "9h", "12h", "24h")
  tp_fwd    <- c(paste0(rna_signal_dir, "STRS_HEK293_WT_cont_0h_fwd.bw"),
                 paste0(rna_signal_dir, "STRS_HEK293_WT_sorb_1h_fwd.bw"),
                 paste0(rna_signal_dir, "STRS_HEK293_WT_sorb_3h_fwd.bw"),
                 paste0(rna_signal_dir, "STRS_HEK293_WT_sorb_6h_fwd.bw"),
                 paste0(rna_signal_dir, "STRS_HEK293_WT_sorb_9h_fwd.bw"),
                 paste0(rna_signal_dir, "STRS_HEK293_WT_sorb_12h_fwd.bw"),
                 paste0(rna_signal_dir, "STRS_HEK293_WT_sorb_24h_fwd.bw"))
  tp_rev    <- sub("_fwd.bw", "_rev.bw", tp_fwd)

  for (j in seq_along(tp_labels)) {
    col <- tp_colors[tp_labels[j]]
    plotSignal(tp_fwd[j], params = p, x = 0.25, y = tp_y[j], height = 0.25,
               linecolor = col, scale = TRUE, range = signalRange)
    plotSignal(tp_rev[j], params = p, x = 0.25, y = tp_y[j], height = 0.25,
               linecolor = col, scale = TRUE, range = signalRange)
    plotText(tp_labels[j], fontcolor = col, rot = 90,
             y = tp_y[j] + 0.1, x = 0.1)
  }

  ## Genes and genome label ----
  plotGenes(param = p, chrom = p$chrom, x = 0.25, y = 8.25, height = 0.5)
  plotGenomeLabel(params = p, x = 0.25, y = 8.85)

  ## Anchor highlights ----
  annoHighlight(plot = control, fill = "lightgrey",
                chrom      = as.character(seqnames(loopRegions_lost))[i],
                chromstart = start(loopRegions_lost)[i],
                chromend   = start(loopRegions_lost)[i] + 10e3,
                default.units = "inches",
                y = 0.25, x = 0.25, height = 9.1, just = c("left", "top"))

  annoHighlight(plot = sorb, fill = "lightgrey",
                chrom      = as.character(seqnames(loopRegions_lost))[i],
                chromstart = end(loopRegions_lost)[i],
                chromend   = end(loopRegions_lost)[i] - 10e3,
                default.units = "inches",
                y = 0.25, x = 0.25, height = 9.1, just = c("left", "top"))
}

dev.off()

sessionInfo()
