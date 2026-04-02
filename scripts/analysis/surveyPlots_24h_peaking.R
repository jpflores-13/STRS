# ##############################################################################
# filename:    surveyPlots_24h_peaking.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Survey plots for genes at gained loop anchors whose expression
#              peaks at 24h (continuously rising, significant at 24h); two Hi-C
#              panels with RNA-seq timecourse tracks; ordered by 24h log2FC
# ##############################################################################

# Libraries ----
library(plotgardener)
library(InteractionSet)
library(DESeq2)
library(GenomicRanges)
library(mariner)
library(glue)

# Parameters ----
dds_rds        <- "data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds"
diff_loops_rds <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
hic_control    <- "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic"
hic_sorbitol   <- "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic"
rna_signal_dir <- "data/processed/rna/timecourse/output/mergeSignal/stranded/"
output_pdf     <- "plots/surveyPlots_24h_peaking.pdf"
padj_24h_cutoff <- 0.1
max_genes       <- 150L
buffer          <- 400e3
page_width      <- 5.75
page_height     <- 10.5

tp_colors <- c("0h"  = "#C7E9B4", "1h"  = "#7FCDBB", "3h"  = "#41B6C4",
               "6h"  = "#1D91C0", "9h"  = "#225EA8", "12h" = "#253494",
               "24h" = "#081D58")

# Data import ----
dds <- readRDS(dds_rds)

diff_loopCounts <- readRDS(diff_loops_rds) |>
  interactions() |>
  keepStandardChromosomes(pruning.mode = "coarse")

rna_files <- list.files(rna_signal_dir, pattern = "\\.bw$", full.names = TRUE)

# Analysis ----
gainedLoops <- diff_loopCounts[diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange > 0]
lostLoops   <- diff_loopCounts[diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange < 0]

as_gi <- function(x) x |> as.data.frame() |> as_ginteractions()
gained_gi <- as_gi(gainedLoops)

dge_gr       <- results(dds, format = "GRanges")
mcols(dge_gr) <- rowData(dds)
dge_gr_prom  <- promoters(dge_gr) |>
  keepStandardChromosomes(pruning.mode = "coarse")
seqlevelsStyle(dge_gr_prom) <- "UCSC"

genes_at_gained <- subsetByOverlaps(dge_gr_prom, gained_gi)
message("Found ", length(genes_at_gained), " genes at gained loop anchors")

## Wald test results per timepoint ----
wald <- lapply(c("Time_1h_vs_0h", "Time_3h_vs_0h", "Time_6h_vs_0h",
                 "Time_9h_vs_0h", "Time_12h_vs_0h", "Time_24h_vs_0h"),
               function(nm) results(dds, name = nm))
names(wald) <- c("1h", "3h", "6h", "9h", "12h", "24h")

mcols(genes_at_gained)$wald_1h_log2FC  <- wald[["1h"]] [names(genes_at_gained), "log2FoldChange"]
mcols(genes_at_gained)$wald_3h_log2FC  <- wald[["3h"]] [names(genes_at_gained), "log2FoldChange"]
mcols(genes_at_gained)$wald_6h_log2FC  <- wald[["6h"]] [names(genes_at_gained), "log2FoldChange"]
mcols(genes_at_gained)$wald_9h_log2FC  <- wald[["9h"]] [names(genes_at_gained), "log2FoldChange"]
mcols(genes_at_gained)$wald_12h_log2FC <- wald[["12h"]][names(genes_at_gained), "log2FoldChange"]
mcols(genes_at_gained)$wald_24h_log2FC <- wald[["24h"]][names(genes_at_gained), "log2FoldChange"]
mcols(genes_at_gained)$wald_24h_padj   <- wald[["24h"]][names(genes_at_gained), "padj"]

## Filter for 24h-peaking genes ----
mc <- mcols(genes_at_gained)
genes_24h_peak <- genes_at_gained[
  !is.na(mc$wald_24h_padj)   & !is.na(mc$wald_24h_log2FC) &
  !is.na(mc$wald_1h_log2FC)  & !is.na(mc$wald_3h_log2FC)  &
  !is.na(mc$wald_6h_log2FC)  & !is.na(mc$wald_9h_log2FC)  &
  !is.na(mc$wald_12h_log2FC) &
  mc$wald_24h_padj   < padj_24h_cutoff &
  mc$wald_24h_log2FC > 0                &
  mc$wald_24h_log2FC > mc$wald_1h_log2FC  &
  mc$wald_24h_log2FC > mc$wald_3h_log2FC  &
  mc$wald_24h_log2FC > mc$wald_6h_log2FC  &
  mc$wald_24h_log2FC > mc$wald_9h_log2FC  &
  mc$wald_24h_log2FC > mc$wald_12h_log2FC
]
genes_24h_peak <- head(genes_24h_peak[order(-genes_24h_peak$wald_24h_log2FC)], max_genes)

message("Genes peaking at 24h: ", length(genes_24h_peak))
if (length(genes_24h_peak) > 0)
  message("Top genes: ", paste(head(genes_24h_peak$symbol, 10), collapse = ", "))

## Anchor / region helpers ----
find_overlapping_anchors <- function(gene_gr, loops_gi) {
  ov1 <- findOverlaps(gene_gr, anchors(loops_gi, "first"))
  ov2 <- findOverlaps(gene_gr, anchors(loops_gi, "second"))
  list(anchor1_loops = subjectHits(ov1),
       anchor2_loops = subjectHits(ov2),
       all_loops     = unique(c(subjectHits(ov1), subjectHits(ov2))))
}

create_loop_centered_region <- function(gene_idx, anchor_info, loops_gi, buf = buffer) {
  idx       <- anchor_info[[gene_idx]]$all_loops[1]
  lp        <- loops_gi[idx]
  a1        <- anchors(lp, "first")
  a2        <- anchors(lp, "second")
  mid       <- round((start(a1) + end(a2)) / 2)
  list(chrom = as.character(seqnames(a1)),
       start = mid - buf, end = mid + buf, main_loop_idx = idx)
}

anchor_info      <- lapply(seq_along(genes_24h_peak), find_overlapping_anchors,
                           loops_gi = gained_gi)
centered_regions <- lapply(seq_along(genes_24h_peak), create_loop_centered_region,
                           anchor_info = anchor_info, loops_gi = gained_gi)

# Visualization ----
tp_y      <- c(4.4, 5.0, 5.6, 6.2, 6.8, 7.4, 8.0)
tp_labels <- c("0h", "1h", "3h", "6h", "9h", "12h", "24h")
tp_fwd    <- c(paste0(rna_signal_dir, "STRS_HEK293_WT_cont_0h_fwd.bw"),
               paste0(rna_signal_dir, "STRS_HEK293_WT_sorb_1h_fwd.bw"),
               paste0(rna_signal_dir, "STRS_HEK293_WT_sorb_3h_fwd.bw"),
               paste0(rna_signal_dir, "STRS_HEK293_WT_sorb_6h_fwd.bw"),
               paste0(rna_signal_dir, "STRS_HEK293_WT_sorb_9h_fwd.bw"),
               paste0(rna_signal_dir, "STRS_HEK293_WT_sorb_12h_fwd.bw"),
               paste0(rna_signal_dir, "STRS_HEK293_WT_sorb_24h_fwd.bw"))
tp_rev <- sub("_fwd.bw", "_rev.bw", tp_fwd)

if (length(genes_24h_peak) == 0) {
  message("No 24h-peaking genes found; skipping PDF output.")
} else {
  pdf(file = output_pdf, width = page_width, height = page_height)

  for (i in seq_along(genes_24h_peak)) {

    region <- centered_regions[[i]]
    g      <- genes_24h_peak[i]

    p <- pgParams(assembly   = "hg38",
                  resolution = 10e3,
                  chrom      = region$chrom,
                  chromstart = region$start,
                  chromend   = region$end,
                  zrange     = c(0, 100),
                  norm       = "SCALE",
                  x          = 0.25,
                  width      = 5,
                  length     = 5,
                  height     = 2)

    pageCreate(width = page_width, height = page_height, showGuides = FALSE)

    plotText(
      label    = glue("{g$symbol} | Peaks at 24h | log2FC: 1h={round(g$wald_1h_log2FC, 2)}, 6h={round(g$wald_6h_log2FC, 2)}, 12h={round(g$wald_12h_log2FC, 2)}, 24h={round(g$wald_24h_log2FC, 2)} | padj={format(g$wald_24h_padj, digits = 2, scientific = TRUE)}"),
      x        = 2.875, y = 0.1, just = "center", fontface = "bold", fontsize = 7
    )

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
    signalRange <- calcSignalRange(data       = rna_files,
                                   chrom      = region$chrom,
                                   chromstart = region$start,
                                   chromend   = region$end,
                                   assembly   = "hg38",
                                   negData    = FALSE)

    for (j in seq_along(tp_labels)) {
      col <- tp_colors[tp_labels[j]]
      plotSignal(tp_fwd[j], params = p, x = 0.25, y = tp_y[j], height = 0.25,
                 linecolor = col, fill = col, scale = TRUE, range = signalRange)
      plotSignal(tp_rev[j], params = p, x = 0.25, y = tp_y[j], height = 0.25,
                 linecolor = col, fill = col, scale = TRUE, range = signalRange)
      plotText(tp_labels[j], fontcolor = col, rot = 90,
               y = tp_y[j] + 0.1, x = 0.1)
    }

    ## Genes and genome label ----
    plotGenes(param = p, chrom = p$chrom, x = 0.25, y = 8.25, height = 1.0,
              fontsize = 8,
              geneHighlights = data.frame(gene = g$symbol, color = "#DC3220"))
    plotGenomeLabel(params = p, x = 0.25, y = 9.35)

    ## Anchor highlights ----
    main_loop <- gained_gi[region$main_loop_idx]
    anch1     <- anchors(main_loop, "first")
    anch2     <- anchors(main_loop, "second")

    for (plot_obj in list(control, sorb)) {
      annoHighlight(plot = plot_obj, fill = "lightgrey",
                    chrom = as.character(seqnames(anch1)),
                    chromstart = start(anch1), chromend = end(anch1),
                    default.units = "inches",
                    y = 0.25, x = 0.25, height = page_height, just = c("left", "top"))
      annoHighlight(plot = plot_obj, fill = "lightgrey",
                    chrom = as.character(seqnames(anch2)),
                    chromstart = start(anch2), chromend = end(anch2),
                    default.units = "inches",
                    y = 0.25, x = 0.25, height = page_height, just = c("left", "top"))
    }
  }

  dev.off()
  message("Saved: ", output_pdf)
}

sessionInfo()
