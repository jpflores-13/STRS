# ##############################################################################
# filename:    surveyPlot_HiC_CutnTag_CTCF_RAD21.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Survey plots for top 100 gained loops with CTCF and RAD21
#              CUT&Tag signal tracks (control and sorbitol); anchor highlights
# ##############################################################################

# Libraries ----
library(plotgardener)
library(InteractionSet)
library(mariner)
library(stringr)

# Parameters ----
diff_loops_rds     <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
hic_control        <- "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic"
hic_sorbitol       <- "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic"
cutntag_signal_dir <- "data/processed/cutntag/output/mergeSignal/"
output_pdf         <- "plots/surveyPlot_HiC_CutnTag_CTCF_RAD21.pdf"
n_loops            <- 100L
buffer             <- 200e3
page_width         <- 5.75
page_height        <- 7.0

# Data import ----
diff_loopCounts <- readRDS(diff_loops_rds) |> interactions()
mcols(diff_loopCounts)$loop_size <- pairdist(diff_loopCounts)

cutntag <- list.files(cutntag_signal_dir, full.names = TRUE)

# Analysis ----
gained_adj <- diff_loopCounts[diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange > 0]
bestGained <- head(gained_adj[order(gained_adj$padj)], n_loops)

loopRegions_gained <- GRanges(
  seqnames = as.character(seqnames(anchors(bestGained, "first"))),
  ranges   = IRanges(start = start(anchors(bestGained, "first")),
                     end   = end(anchors(bestGained, "second"))),
  mcols    = mcols(bestGained)
)
loopRegions_gained_buffed <- loopRegions_gained + buffer

lostLoops   <- diff_loopCounts[diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange < 0]
gainedLoops <- diff_loopCounts[diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange > 0]

# Visualization ----
pdf(file = output_pdf, width = page_width, height = page_height)

for (i in seq_along(loopRegions_gained_buffed)) {

  p <- pgParams(assembly   = "hg38",
                resolution = 10e3,
                chrom      = as.character(seqnames(loopRegions_gained_buffed))[i],
                chromstart = start(loopRegions_gained_buffed)[i],
                chromend   = end(loopRegions_gained_buffed)[i],
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

  ## CTCF tracks ----
  ctcf_range <- calcSignalRange(data = str_subset(cutntag, "CTCF"),
                                chrom = p$chrom, chromstart = p$chromstart,
                                chromend = p$chromend, assembly = "hg38",
                                negData = FALSE)

  plotSignal(paste0(cutntag_signal_dir, "STRS_HEK293_eGFP-YAP_CTCF_cont_0h.bw"),
             params = p, x = 0.25, y = 4.45, height = 0.25,
             linecolor = "#7FCDBB", fill = "#7FCDBB",
             scale = TRUE, range = ctcf_range)
  plotText("CTCF", fontcolor = "#7FCDBB", rot = 90, y = 4.85, x = 0.1)

  plotSignal(paste0(cutntag_signal_dir, "STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h.bw"),
             params = p, x = 0.25, y = 4.95, height = 0.25,
             linecolor = "#7FCDBB", fill = "#7FCDBB",
             scale = TRUE, range = ctcf_range)

  ## RAD21 tracks ----
  rad21_range <- calcSignalRange(data = str_subset(cutntag, "RAD21"),
                                 chrom = p$chrom, chromstart = p$chromstart,
                                 chromend = p$chromend, assembly = "hg38",
                                 negData = FALSE)

  plotSignal(paste0(cutntag_signal_dir, "STRS_HEK293_eGFP-YAP_RAD21_cont_0h.bw"),
             params = p, x = 0.25, y = 5.35, height = 0.25,
             linecolor = "#41B6C4", fill = "#41B6C4",
             scale = TRUE, range = rad21_range)
  plotText("RAD21", fontcolor = "#41B6C4", rot = 90, y = 5.75, x = 0.1)

  plotSignal(paste0(cutntag_signal_dir, "STRS_HEK293_eGFP-YAP_RAD21_sorbitol_1h.bw"),
             params = p, x = 0.25, y = 5.85, height = 0.25,
             linecolor = "#41B6C4", fill = "#41B6C4",
             scale = TRUE, range = rad21_range)

  ## Genes and genome label ----
  plotGenes(param = p, chrom = p$chrom, x = 0.25, y = 6.15, height = 0.5)
  plotGenomeLabel(params = p, x = 0.25, y = 6.75)

  ## Anchor highlights ----
  annoHighlight(plot = control, fill = "lightgrey",
                chrom      = as.character(seqnames(loopRegions_gained))[i],
                chromstart = start(loopRegions_gained)[i],
                chromend   = start(loopRegions_gained)[i] + 10e3,
                default.units = "inches",
                y = 0.25, x = 0.25, height = 6.5, just = c("left", "top"))

  annoHighlight(plot = sorb, fill = "lightgrey",
                chrom      = as.character(seqnames(loopRegions_gained))[i],
                chromstart = end(loopRegions_gained)[i],
                chromend   = end(loopRegions_gained)[i] - 10e3,
                default.units = "inches",
                y = 0.25, x = 0.25, height = 6.5, just = c("left", "top"))
}

dev.off()

sessionInfo()
