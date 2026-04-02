# ##############################################################################
# filename:    surveyPlot_HiC_CutnTag.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Survey plots for top 100 gained loops with CUT&Tag signal and
#              peak tracks for H3K27ac, CTCF, RAD21, and YAP1
#              (control and sorbitol); anchor highlights included
# ##############################################################################

# Libraries ----
library(plotgardener)
library(InteractionSet)
library(GenomicRanges)
library(mariner)
library(stringr)

# Parameters ----
diff_loops_rds    <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
hic_control       <- "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic"
hic_sorbitol      <- "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic"
cutntag_signal_dir <- "data/processed/cutntag/output/mergeSignal/"
cutntag_peaks_dir  <- "data/processed/cutntag/output/peaks/"
output_pdf        <- "plots/surveyPlot_HiC_CutnTag.pdf"
n_loops           <- 100L
buffer            <- 200e3
page_width        <- 5.75
page_height       <- 9.75

# Data import ----
diff_loopCounts <- readRDS(diff_loops_rds) |> interactions()
mcols(diff_loopCounts)$loop_size <- pairdist(diff_loopCounts)

cutntag <- list.files(cutntag_signal_dir, full.names = TRUE)

read_narrowpeaks <- function(file) {
  peaks <- read.table(file, sep = "\t", header = FALSE,
                      col.names = c("seqnames", "start", "end", "name",
                                    "score", "strand", "signalValue",
                                    "pValue", "qValue", "peak"))
  GRanges(seqnames    = peaks$seqnames,
          ranges      = IRanges(start = peaks$start + 1L, end = peaks$end),
          strand      = "*",
          score       = peaks$score,
          signalValue = peaks$signalValue,
          pValue      = peaks$pValue,
          qValue      = peaks$qValue,
          peak        = peaks$peak)
}

cutntag_np <- list.files(cutntag_peaks_dir, full.names = TRUE,
                          pattern = "\\.narrowPeak$") |>
  lapply(read_narrowpeaks)
names(cutntag_np) <- paste0(rep(c("CTCF", "H3K27ac", "RAD21", "YAP1"), each = 2),
                             "_", c("control", "sorbitol"))

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

  sorb <- plotHicRectangle(data = hic_sorbitol, params = p, y = 2.35)
  plotText(label = "+ sorbitol", x = 0.25, y = 2.35, just = c("top", "left"))
  annoHeatmapLegend(sorb, orientation = "v",
                    fontsize = 8, fontcolor = "black", digits = 2,
                    x = 5.5, y = 2.35, width = 0.1, height = 1.5,
                    just = c("left", "top"), default.units = "inches")

  ## H3K27ac tracks ----
  h3k27ac_range <- calcSignalRange(data = str_subset(cutntag, "H3K27ac"),
                                   chrom = p$chrom, chromstart = p$chromstart,
                                   chromend = p$chromend, assembly = "hg38",
                                   negData = FALSE)

  plotRanges(data = cutntag_np$H3K27ac_control, params = p,
             x = 0.25, y = 4.35, width = 5, height = 0.08,
             fill = "#C7E9B4", colorby = "signalValue",
             palette = colorRampPalette(c("white", "#C7E9B4"))(100),
             stroke = 0.05, strokecolor = "black", collapse = TRUE,
             just = c("left", "top"), default.units = "inches")
  plotSignal(paste0(cutntag_signal_dir, "STRS_HEK293_eGFP-YAP_H3K27ac_cont_0h.bw"),
             params = p, x = 0.25, y = 4.4, height = 0.25,
             linecolor = "#C7E9B4", fill = "#C7E9B4",
             scale = TRUE, range = h3k27ac_range)
  plotText("H3K27ac", fontcolor = "#C7E9B4", rot = 90, y = 4.9, x = 0.1)

  plotRanges(data = cutntag_np$H3K27ac_sorbitol, params = p,
             x = 0.25, y = 4.95, width = 5, height = 0.08,
             fill = "#C7E9B4", colorby = "signalValue",
             palette = colorRampPalette(c("white", "#C7E9B4"))(100),
             stroke = 0.05, strokecolor = "black", alpha = 0.7, collapse = TRUE,
             just = c("left", "top"), default.units = "inches")
  plotSignal(paste0(cutntag_signal_dir, "STRS_HEK293_eGFP-YAP_H3K27ac_sorbitol_1h.bw"),
             params = p, x = 0.25, y = 5.0, height = 0.25,
             linecolor = "#C7E9B4", fill = "#C7E9B4",
             scale = TRUE, range = h3k27ac_range)

  ## CTCF tracks ----
  ctcf_range <- calcSignalRange(data = str_subset(cutntag, "CTCF"),
                                chrom = p$chrom, chromstart = p$chromstart,
                                chromend = p$chromend, assembly = "hg38",
                                negData = FALSE)

  plotRanges(data = cutntag_np$CTCF_control, params = p,
             x = 0.25, y = 5.55, width = 5, height = 0.08,
             fill = "#7FCDBB", colorby = "signalValue",
             palette = colorRampPalette(c("white", "#7FCDBB"))(100),
             stroke = 0.05, strokecolor = "black", collapse = TRUE,
             just = c("left", "top"), default.units = "inches")
  plotSignal(paste0(cutntag_signal_dir, "STRS_HEK293_eGFP-YAP_CTCF_cont_0h.bw"),
             params = p, x = 0.25, y = 5.6, height = 0.25,
             linecolor = "#7FCDBB", fill = "#7FCDBB",
             scale = TRUE, range = ctcf_range)
  plotText("CTCF", fontcolor = "#7FCDBB", rot = 90, y = 6.0, x = 0.1)

  plotRanges(data = cutntag_np$CTCF_sorbitol, params = p,
             x = 0.25, y = 6.15, width = 5, height = 0.08,
             fill = "#7FCDBB", colorby = "signalValue",
             palette = colorRampPalette(c("white", "#7FCDBB"))(100),
             stroke = 0.05, strokecolor = "black", alpha = 0.7, collapse = TRUE,
             just = c("left", "top"), default.units = "inches")
  plotSignal(paste0(cutntag_signal_dir, "STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h.bw"),
             params = p, x = 0.25, y = 6.2, height = 0.25,
             linecolor = "#7FCDBB", fill = "#7FCDBB",
             scale = TRUE, range = ctcf_range)

  ## RAD21 tracks ----
  rad21_range <- calcSignalRange(data = str_subset(cutntag, "RAD21"),
                                 chrom = p$chrom, chromstart = p$chromstart,
                                 chromend = p$chromend, assembly = "hg38",
                                 negData = FALSE)

  plotRanges(data = cutntag_np$RAD21_control, params = p,
             x = 0.25, y = 6.75, width = 5, height = 0.08,
             fill = "#41B6C4", colorby = "signalValue",
             palette = colorRampPalette(c("white", "#41B6C4"))(100),
             stroke = 0.05, strokecolor = "black", collapse = TRUE,
             just = c("left", "top"), default.units = "inches")
  plotSignal(paste0(cutntag_signal_dir, "STRS_HEK293_eGFP-YAP_RAD21_cont_0h.bw"),
             params = p, x = 0.25, y = 6.8, height = 0.25,
             linecolor = "#41B6C4", fill = "#41B6C4",
             scale = TRUE, range = rad21_range)
  plotText("RAD21", fontcolor = "#41B6C4", rot = 90, y = 7.2, x = 0.1)

  plotRanges(data = cutntag_np$RAD21_sorbitol, params = p,
             x = 0.25, y = 7.35, width = 5, height = 0.08,
             fill = "#41B6C4", colorby = "signalValue",
             palette = colorRampPalette(c("white", "#41B6C4"))(100),
             stroke = 0.05, strokecolor = "black", alpha = 0.7, collapse = TRUE,
             just = c("left", "top"), default.units = "inches")
  plotSignal(paste0(cutntag_signal_dir, "STRS_HEK293_eGFP-YAP_RAD21_sorbitol_1h.bw"),
             params = p, x = 0.25, y = 7.4, height = 0.25,
             linecolor = "#41B6C4", fill = "#41B6C4",
             scale = TRUE, range = rad21_range)

  ## YAP1 tracks ----
  yap1_range <- calcSignalRange(data = str_subset(cutntag, "YAP1"),
                                chrom = p$chrom, chromstart = p$chromstart,
                                chromend = p$chromend, assembly = "hg38",
                                negData = FALSE)

  plotRanges(data = cutntag_np$YAP1_control, params = p,
             x = 0.25, y = 7.95, width = 5, height = 0.08,
             fill = "#1D91C0", colorby = "signalValue",
             palette = colorRampPalette(c("white", "#1D91C0"))(100),
             stroke = 0.05, strokecolor = "black", collapse = TRUE,
             just = c("left", "top"), default.units = "inches")
  plotSignal(paste0(cutntag_signal_dir, "STRS_HEK293_eGFP-YAP_YAP1_cont_0h.bw"),
             params = p, x = 0.25, y = 8.0, height = 0.25,
             linecolor = "#1D91C0", fill = "#1D91C0",
             scale = TRUE, range = yap1_range)
  plotText("YAP1", fontcolor = "#1D91C0", rot = 90, y = 8.5, x = 0.1)

  plotRanges(data = cutntag_np$YAP1_sorbitol, params = p,
             x = 0.25, y = 8.55, width = 5, height = 0.08,
             fill = "#1D91C0", colorby = "signalValue",
             palette = colorRampPalette(c("white", "#1D91C0"))(100),
             stroke = 0.05, strokecolor = "black", alpha = 0.7, collapse = TRUE,
             just = c("left", "top"), default.units = "inches")
  plotSignal(paste0(cutntag_signal_dir, "STRS_HEK293_eGFP-YAP_YAP1_sorbitol_1h.bw"),
             params = p, x = 0.25, y = 8.6, height = 0.25,
             linecolor = "#1D91C0", fill = "#1D91C0",
             scale = TRUE, range = yap1_range)

  ## Genes and genome label ----
  plotGenes(param = p, chrom = p$chrom, x = 0.25, y = 8.9, height = 0.5)
  plotGenomeLabel(params = p, x = 0.25, y = 9.5)

  ## Anchor highlights ----
  annoHighlight(plot = control, fill = "lightgrey",
                chrom      = as.character(seqnames(loopRegions_gained))[i],
                chromstart = start(loopRegions_gained)[i],
                chromend   = start(loopRegions_gained)[i] + 10e3,
                default.units = "inches",
                y = 0.25, x = 0.25, height = 9.25, just = c("left", "top"))

  annoHighlight(plot = sorb, fill = "lightgrey",
                chrom      = as.character(seqnames(loopRegions_gained))[i],
                chromstart = end(loopRegions_gained)[i],
                chromend   = end(loopRegions_gained)[i] - 10e3,
                default.units = "inches",
                y = 0.25, x = 0.25, height = 9.25, just = c("left", "top"))
}

dev.off()

sessionInfo()
