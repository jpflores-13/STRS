# ##############################################################################
# filename:    surveyPlot_gainedLoopGenes_rnaTimecourse_10kb.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Survey plots centered on gained loop midpoints for genes at
#              gained anchors with significant temporal expression (LRT
#              padj < 0.1); two Hi-C panels plus seven RNA-seq signal tracks;
#              ordered by LRT padj
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
output_pdf     <- "plots/surveyPlot_gainedLoopGenes_rnaTimecourse_10kb.pdf"
padj_lrt_cutoff <- 0.1
buffer          <- 200e3
page_width      <- 5.75
page_height     <- 9.1

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
dge_gr       <- results(dds, format = "GRanges")
mcols(dge_gr) <- rowData(dds)
dge_gr_prom  <- promoters(dge_gr) |>
  keepStandardChromosomes(pruning.mode = "coarse")
seqlevelsStyle(dge_gr_prom) <- "UCSC"
dge_gr_prom$padj_lrt <- p.adjust(dge_gr_prom$LRTPvalue, method = "BH")

lostLoops   <- diff_loopCounts[diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange < 0]
gainedLoops <- diff_loopCounts[diff_loopCounts$padj < 0.1 & diff_loopCounts$log2FoldChange > 0]

as_gi <- function(x) x |> as.data.frame() |> as_ginteractions()
gained_gi <- as_gi(gainedLoops)

gained_genes <- subsetByOverlaps(dge_gr_prom, gained_gi) |>
  subset(padj_lrt < padj_lrt_cutoff)
gained_genes$loop_type <- "gained"
gained_genes <- gained_genes[order(gained_genes$padj_lrt)]

message("Found ", length(gained_genes), " genes at gained loop anchors with significant temporal expression changes")
message("Top 5 genes: ", paste(head(gained_genes$symbol, 5), collapse = ", "))

find_overlapping_anchors <- function(gene_gr, loops_gi) {
  anchor1    <- anchors(loops_gi, "first")
  anchor2    <- anchors(loops_gi, "second")
  ov1        <- findOverlaps(gene_gr, anchor1)
  ov2        <- findOverlaps(gene_gr, anchor2)
  loop_idx1  <- subjectHits(ov1)
  loop_idx2  <- subjectHits(ov2)
  list(anchor1_loops = loop_idx1,
       anchor2_loops = loop_idx2,
       all_loops     = unique(c(loop_idx1, loop_idx2)))
}

create_loop_centered_region <- function(gene_idx, anchor_info, loops_gi, buf = buffer) {
  loop_indices <- anchor_info[[gene_idx]]$all_loops
  main_loop    <- loops_gi[loop_indices[1]]
  anchor1      <- anchors(main_loop, "first")
  anchor2      <- anchors(main_loop, "second")
  loop_midpoint <- round((start(anchor1) + end(anchor2)) / 2)
  list(chrom         = as.character(seqnames(anchor1)),
       start         = loop_midpoint - buf,
       end           = loop_midpoint + buf,
       main_loop_idx = loop_indices[1])
}

anchor_info      <- lapply(seq_along(gained_genes), find_overlapping_anchors,
                           loops_gi = gained_gi)
centered_regions <- lapply(seq_along(gained_genes), create_loop_centered_region,
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

pdf(file = output_pdf, width = page_width, height = page_height)

for (i in seq_along(gained_genes)) {

  region <- centered_regions[[i]]

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

  ## Title ----
  plotText(
    label    = glue("{gained_genes[i]$symbol} (gained loop, LRT padj = {format(gained_genes[i]$padj_lrt, digits = 3, scientific = TRUE)})"),
    x        = 2.875,
    y        = 0.1,
    just     = "center",
    fontface = "bold",
    fontsize = 9
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

message("Survey plots created for ", length(gained_genes), " genes at gained loop anchors")

sessionInfo()
