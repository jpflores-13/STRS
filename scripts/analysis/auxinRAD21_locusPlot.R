# ##############################################################################
# filename:    auxinRAD21_locusPlot.R
# author:      JP Flores
# date:        2026-05-29
# project:     STRS
# description: Genome browser panel at the SOX8 locus (chr16) analogous to
#              Figure 5D. Shows two Hi-C contact maps (untreated and sorbitol,
#              HEK293 eGFP-YAP) alongside three RNA-seq signal tracks
#              (untreated, sorbitol 6h, sorbitol+auxin 6h) from HCT116
#              mAID2-RAD21 cells. No CUT&Tag tracks are included.
#              Note: Hi-C is from HEK293 eGFP-YAP cells; RNA-seq is from
#              HCT116 mAID2-RAD21 cells (cross-cell-line visualization).
# input:       HEK293 eGFP-YAP Hi-C .hic maps, AuxinRAD21 merged stranded
#              bigWig files, differential loop RDS
# output:      plots/auxinRAD21_locusPlot.pdf
# ##############################################################################


# Libraries -------------------------------------------------------------------

library(GenomicRanges)
library(GenomeInfoDb)
library(InteractionSet)
library(IRanges)
library(mariner)
library(plotgardener)
library(S4Vectors)
library(stringr)


# Utility scripts -------------------------------------------------------------

source("scripts/utils/ggplot2_pgTheme.R")


# Parameters ------------------------------------------------------------------

project_dir <- "/work/users/j/p/jpflores/projects/STRS"
plots_dir   <- file.path(project_dir, "plots")

## Hi-C maps (HEK293 eGFP-YAP)
hic_control  <- file.path(project_dir,
                          "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic")
hic_sorbitol <- file.path(project_dir,
                          "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic")

## Differential loops
diff_loops_rds <- file.path(project_dir,
                            "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds")

## AuxinRAD21 merged stranded bigWig directory
rna_bw_dir <- file.path(project_dir,
                        "data/processed/rna/AuxinRAD21/output/mergeSignal/stranded")

## Locus — SOX8 locus, chr16 (matching Figure 5D)
locus_chrom      <- "chr16"
locus_chromstart <- 870000
locus_chromend   <- 1260000

## Hi-C display
resolution  <- 10e3
zrange_hic  <- c(0, 100)
norm_method <- "SCALE"

## Track dimensions
page_width   <- 3.5
page_height  <- 5.5
x_start      <- 0.25
plot_width   <- 3.0
hic_height   <- 1.20
hic_gap      <- 0.10
rna_height   <- 0.20
rna_gap      <- 0.05
gene_height  <- 0.25

## Colors
color_ctrl     <- "#999999"
color_sorbitol <- "#F8766D"
color_auxin    <- "#009E73"   # Okabe-Ito colorblind-friendly green
gray_label     <- "#666666"

## Differential loop thresholds
padj_loop   <- 0.1

## Output
output_pdf <- file.path(plots_dir, "auxinRAD21_locusPlot.pdf")


# Load data -------------------------------------------------------------------

## Differential loops
diff_loopCounts <- readRDS(diff_loops_rds) |> interactions()
diff_loopCounts <- keepStandardChromosomes(diff_loopCounts,
                                           pruning.mode = "coarse")

gainedLoops <- diff_loopCounts[diff_loopCounts$padj < padj_loop &
                                 diff_loopCounts$log2FoldChange > 0]
lostLoops   <- diff_loopCounts[diff_loopCounts$padj < padj_loop &
                                 diff_loopCounts$log2FoldChange < 0]

## Find the gained loop(s) at the target locus for anchor highlights
locus_gr    <- GRanges(locus_chrom, IRanges(locus_chromstart, locus_chromend))
target_loop <- gainedLoops[
  overlapsAny(anchors(gainedLoops, "first"),  locus_gr) |
    overlapsAny(anchors(gainedLoops, "second"), locus_gr)
]

## RNA-seq bigWig files (merged across replicates, stranded)
bw_ctrl_fwd  <- file.path(rna_bw_dir, "YAPP_HCT116_mAID2-RAD21_untreated_0h_fwd.bw")
bw_ctrl_rev  <- file.path(rna_bw_dir, "YAPP_HCT116_mAID2-RAD21_untreated_0h_rev.bw")
bw_sorb_fwd  <- file.path(rna_bw_dir, "YAPP_HCT116_mAID2-RAD21_sorbitol_6h_fwd.bw")
bw_sorb_rev  <- file.path(rna_bw_dir, "YAPP_HCT116_mAID2-RAD21_sorbitol_6h_rev.bw")
bw_auxin_fwd <- file.path(rna_bw_dir, "YAPP_HCT116_mAID2-RAD21_auxin_6h_fwd.bw")
bw_auxin_rev <- file.path(rna_bw_dir, "YAPP_HCT116_mAID2-RAD21_auxin_6h_rev.bw")

## Compute a shared signal range across all 3 fwd-strand bigWigs
## so all RNA tracks use the same y-scale
all_fwd_bws <- c(bw_ctrl_fwd, bw_sorb_fwd, bw_auxin_fwd)
signalRange  <- calcSignalRange(
  data       = all_fwd_bws,
  chrom      = locus_chrom,
  chromstart = locus_chromstart,
  chromend   = locus_chromend,
  assembly   = "hg38",
  negData    = FALSE
)


# Visualization ---------------------------------------------------------------

## RNA track metadata — order determines top-to-bottom rendering
rna_tracks <- list(
  list(fwd   = bw_ctrl_fwd,  rev = bw_ctrl_rev,
       color = color_ctrl,     label = "untreated"),
  list(fwd   = bw_sorb_fwd,  rev = bw_sorb_rev,
       color = color_sorbitol, label = "+ 6h sorbitol"),
  list(fwd   = bw_auxin_fwd, rev = bw_auxin_rev,
       color = color_auxin,    label = "+ 6h sorbitol + auxin")
)

pdf(file = output_pdf, width = page_width, height = page_height)

pageCreate(width     = page_width,
           height    = page_height,
           xgrid     = 0,
           ygrid     = 0,
           showGuides = FALSE)

## Shared pgParams for genomic coordinates — passed explicitly to each function
## (pgParams fields are not accessible via $; always pass args directly)
p <- pgParams(
  assembly   = "hg38",
  resolution = resolution,
  chrom      = locus_chrom,
  chromstart = locus_chromstart,
  chromend   = locus_chromend,
  zrange     = zrange_hic,
  norm       = norm_method,
  x          = x_start,
  width      = plot_width
)

## 1) Untreated Hi-C ---------------------------------------------------------
y_hic_ctrl <- 0.30

ctrl_hic <- plotHicRectangle(
  data   = hic_control,
  params = p,
  y      = y_hic_ctrl,
  height = hic_height
)

annoPixels(ctrl_hic, data = gainedLoops,
           shift = 0.5, type = "arrow", col = color_sorbitol)
annoPixels(ctrl_hic, data = lostLoops,
           shift = 0.5, type = "arrow", col = color_ctrl)

plotText("untreated",
         x = x_start, y = y_hic_ctrl,
         just = c("left", "top"), fontsize = 7,
         fontcolor = gray_label)

## 2) Sorbitol Hi-C ----------------------------------------------------------
y_hic_sorb <- y_hic_ctrl + hic_height + hic_gap

sorb_hic <- plotHicRectangle(
  data   = hic_sorbitol,
  params = p,
  y      = y_hic_sorb,
  height = hic_height
)

annoPixels(sorb_hic, data = gainedLoops,
           shift = 0.5, type = "arrow", col = color_sorbitol)
annoPixels(sorb_hic, data = lostLoops,
           shift = 0.5, type = "arrow", col = color_ctrl)

plotText("+ 1h sorbitol",
         x = x_start, y = y_hic_sorb,
         just = c("left", "top"), fontsize = 7,
         fontcolor = gray_label)

## Hi-C color scale legend (right of first map, spans both)
annoHeatmapLegend(ctrl_hic,
                  orientation = "v",
                  fontsize    = 6,
                  fontcolor   = "black",
                  digits      = 2,
                  x           = x_start + plot_width + 0.05,
                  y           = y_hic_ctrl,
                  width       = 0.07,
                  height      = 0.70,
                  just        = c("left", "top"),
                  default.units = "inches")

## 3) RNA-seq signal tracks --------------------------------------------------
## Leave a small gap below the Hi-C maps
y_rna_start <- y_hic_sorb + hic_height + 0.18

for (j in seq_along(rna_tracks)) {
  yj    <- y_rna_start + (j - 1) * (rna_height + rna_gap)
  track <- rna_tracks[[j]]
  
  ## Forward strand — filled polygon
  if (file.exists(track$fwd)) {
    plotSignal(
      data      = track$fwd,
      params    = p,
      x         = x_start,
      y         = yj,
      height    = rna_height,
      linecolor = track$color,
      fill      = track$color,
      scale     = FALSE,
      range     = signalRange
    )
  }
  
  ## Reverse strand — line only (no fill), same range
  if (file.exists(track$rev)) {
    plotSignal(
      data      = track$rev,
      params    = p,
      x         = x_start,
      y         = yj,
      height    = rna_height,
      linecolor = track$color,
      scale     = FALSE,
      range     = signalRange
    )
  }
  
  ## Condition label — placed just inside the top of the track (yj + 0.06)
  ## so it sits in the gap below the previous track rather than overlapping it
  plotText(track$label,
           x         = x_start,
           y         = yj + 0.06,
           just      = c("left", "bottom"),
           fontsize  = 6,
           fontcolor = track$color)
}

## Vertical "RNA-seq" label to the left of the tracks
n_rna     <- length(rna_tracks)
rna_end_y <- y_rna_start + n_rna * (rna_height + rna_gap) - rna_gap
rna_mid_y <- (y_rna_start + rna_end_y) / 2

plotText("RNA-seq",
         rot       = 90,
         x         = x_start - 0.12,
         y         = rna_mid_y,
         fontsize  = 7,
         fontcolor = gray_label)

## Signal range label (top-right, above first track)
range_label <- sprintf("%d\u2013%d",
                       round(signalRange[1]),
                       round(signalRange[2]))
plotText(range_label,
         x         = x_start + plot_width,
         y         = y_rna_start - 0.01,
         just      = c("right", "bottom"),
         fontsize  = 6,
         fontcolor = gray_label)

## 4) Gene track -------------------------------------------------------------
y_genes <- rna_end_y + 0.10

plotGenes(
  params  = p,
  chrom   = locus_chrom,
  x       = x_start,
  y       = y_genes,
  height  = gene_height,
  fontsize = 7,
  geneHighlights = data.frame(gene = "SOX8", color = color_sorbitol)
)

## 5) Genome label -----------------------------------------------------------
y_genome <- y_genes + gene_height + 0.05

plotGenomeLabel(
  params   = p,
  x        = x_start,
  y        = y_genome,
  length   = plot_width,
  fontsize = 7
)

## 6) Loop anchor highlights (grey shading at both anchors) ------------------
if (length(target_loop) >= 1L) {
  loop <- target_loop[1L]      # use the top-ranked loop at this locus
  hl_y <- y_hic_ctrl
  hl_h <- y_genome - hl_y     # span from top of first Hi-C map to genome label
  
  ## Left anchor
  annoHighlight(
    plot       = ctrl_hic,
    fill       = "lightgrey",
    alpha      = 0.5,
    chrom      = as.character(seqnames(anchors(loop, "first"))),
    chromstart = start(anchors(loop, "first")),
    chromend   = start(anchors(loop, "first")) + resolution,
    x             = x_start,
    y             = hl_y,
    width         = plot_width,
    height        = hl_h,
    just          = c("left", "top"),
    default.units = "inches"
  )
  
  ## Right anchor
  annoHighlight(
    plot       = ctrl_hic,
    fill       = "lightgrey",
    alpha      = 0.5,
    chrom      = as.character(seqnames(anchors(loop, "second"))),
    chromstart = end(anchors(loop, "second")) - resolution,
    chromend   = end(anchors(loop, "second")),
    x             = x_start,
    y             = hl_y,
    width         = plot_width,
    height        = hl_h,
    just          = c("left", "top"),
    default.units = "inches"
  )
}


# Save outputs ----------------------------------------------------------------

dev.off()
message("Locus plot saved: ", output_pdf)


# Session info ----------------------------------------------------------------

sessionInfo()