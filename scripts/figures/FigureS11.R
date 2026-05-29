# ##############################################################################
# filename:    FigureS11.R
# author:      JP Flores
# date:        2026-05-29
# project:     STRS
# description: Figure S11 — Cohesin depletion RNA-seq panels:
#                Panel A (left):  Genome browser locus plot at SOX8 (chr16)
#                                 showing Hi-C maps and AuxinRAD21 RNA-seq
#                                 signal tracks (auxinRAD21_locusPlot.R)
#                Panel B (right, top): Density plot of log2FC at sorbitol-
#                                 induced loop anchors for sorbitol and
#                                 sorbitol+auxin conditions
#                                 (auxinRAD21_gainedLoops_histogram.R)
#                Panel C (right, bottom): Grouped barplot of mean log2FC at
#                                 gained vs pre-existing loop anchors
#                                 (auxinRAD21_gainedLoops_barplot.R)
#              NOTE: Run auxinRAD21_gainedLoops_histogram.R and
#              auxinRAD21_gainedLoops_barplot.R first to generate the
#              intermediate ggplot .rds files loaded here.
# input:       Hi-C .hic maps, AuxinRAD21 bigWig files, differential loop RDS,
#              intermediate ggplot RDS objects from histogram and barplot scripts
# output:      figures/FigureS11.pdf
# ##############################################################################


# Libraries -------------------------------------------------------------------

library(GenomicRanges)
library(GenomeInfoDb)
library(ggplot2)
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
figures_dir <- file.path(project_dir, "figures")

## Hi-C maps (Panel A)
hic_control  <- file.path(project_dir,
                          "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic")
hic_sorbitol <- file.path(project_dir,
                          "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic")

## Differential loops
diff_loops_rds <- file.path(project_dir,
                            "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds")

## AuxinRAD21 merged stranded bigWig directory (Panel A)
rna_bw_dir <- file.path(project_dir,
                        "data/processed/rna/AuxinRAD21/output/mergeSignal/stranded")

## Pre-computed ggplot objects from individual plot scripts (Panels B and C)
## Run auxinRAD21_gainedLoops_histogram.R and auxinRAD21_gainedLoops_barplot.R first
gg_histogram_rds <- file.path(project_dir,
                              "data/processed/rna/AuxinRAD21/output/deseqObjs",
                              "auxinRAD21_gg_histogram.rds")

gg_barplot_rds <- file.path(project_dir,
                            "data/processed/rna/AuxinRAD21/output/deseqObjs",
                            "auxinRAD21_gg_barplot.rds")

## Locus parameters — SOX8 locus, chr16 (matching Figure 5D)
locus_chrom      <- "chr16"
locus_chromstart <- 870000
locus_chromend   <- 1260000

## Hi-C display
resolution  <- 10e3
zrange_hic  <- c(0, 100)
norm_method <- "SCALE"

## Track dimensions for Panel A (locus plot)
hic_height  <- 1.20
hic_gap     <- 0.10
rna_height  <- 0.20
rna_gap     <- 0.05
gene_height <- 0.25

## Page and layout dimensions
page_width  <- 11.0    # inches, landscape
page_height <-  5.0    # reduced to match content height

## Panel A (locus — left column)
panelA_x     <- 0.35
panelA_y     <- 0.40
panelA_w     <- 3.50
panelA_hic_x <- panelA_x + 0.20   # inset for plotgardener genomic tracks

## Panel B and C height is set so B + gap + C = Panel A content height (~3.78")
## Panel A ends at y ≈ 4.18 (panelA_y + 2 Hi-C + RNA + genes + genome label)
## 2 × 1.80 + 0.18 gap = 3.78 ✓
panelB_x <- 5.00
panelB_y <- 0.40
panelB_w <- 5.50
panelB_h <- 1.80

## Panel C starts immediately after Panel B + gap
panelC_x <- 5.00
panelC_y <- panelB_y + panelB_h + 0.18   # = 2.38
panelC_w <- 5.50
panelC_h <- 1.80

## Colors
color_ctrl     <- "#999999"
color_sorbitol <- "#F8766D"
color_auxin    <- "#009E73"   # Okabe-Ito colorblind-friendly green
gray_label     <- "#666666"

## Differential loop threshold
padj_loop <- 0.1

## Output
output_pdf <- file.path(figures_dir, "FigureS11.pdf")


# Load data -------------------------------------------------------------------

## Differential loops
diff_loopCounts <- readRDS(diff_loops_rds) |> interactions()
diff_loopCounts <- keepStandardChromosomes(diff_loopCounts,
                                           pruning.mode = "coarse")

gainedLoops <- diff_loopCounts[diff_loopCounts$padj < padj_loop &
                                 diff_loopCounts$log2FoldChange > 0]
lostLoops   <- diff_loopCounts[diff_loopCounts$padj < padj_loop &
                                 diff_loopCounts$log2FoldChange < 0]

locus_gr    <- GRanges(locus_chrom, IRanges(locus_chromstart, locus_chromend))
target_loop <- gainedLoops[
  overlapsAny(anchors(gainedLoops, "first"),  locus_gr) |
    overlapsAny(anchors(gainedLoops, "second"), locus_gr)
]

## RNA-seq bigWigs
bw_ctrl_fwd  <- file.path(rna_bw_dir, "YAPP_HCT116_mAID2-RAD21_untreated_0h_fwd.bw")
bw_ctrl_rev  <- file.path(rna_bw_dir, "YAPP_HCT116_mAID2-RAD21_untreated_0h_rev.bw")
bw_sorb_fwd  <- file.path(rna_bw_dir, "YAPP_HCT116_mAID2-RAD21_sorbitol_6h_fwd.bw")
bw_sorb_rev  <- file.path(rna_bw_dir, "YAPP_HCT116_mAID2-RAD21_sorbitol_6h_rev.bw")
bw_auxin_fwd <- file.path(rna_bw_dir, "YAPP_HCT116_mAID2-RAD21_auxin_6h_fwd.bw")
bw_auxin_rev <- file.path(rna_bw_dir, "YAPP_HCT116_mAID2-RAD21_auxin_6h_rev.bw")

all_fwd_bws <- c(bw_ctrl_fwd, bw_sorb_fwd, bw_auxin_fwd)
signalRange  <- calcSignalRange(
  data       = all_fwd_bws,
  chrom      = locus_chrom,
  chromstart = locus_chromstart,
  chromend   = locus_chromend,
  assembly   = "hg38",
  negData    = FALSE
)

## Pre-computed ggplot objects from Panels B and C
stopifnot(file.exists(gg_histogram_rds),
          file.exists(gg_barplot_rds))

p_hist <- readRDS(gg_histogram_rds)
p_bar  <- readRDS(gg_barplot_rds)

## RNA track metadata
rna_tracks <- list(
  list(fwd   = bw_ctrl_fwd,  rev = bw_ctrl_rev,
       color = color_ctrl,     label = "untreated"),
  list(fwd   = bw_sorb_fwd,  rev = bw_sorb_rev,
       color = color_sorbitol, label = "+ 6h sorbitol"),
  list(fwd   = bw_auxin_fwd, rev = bw_auxin_rev,
       color = color_auxin,    label = "+ 6h sorbitol + auxin")
)


# Visualization ---------------------------------------------------------------

pdf(file = output_pdf, width = page_width, height = page_height)

pageCreate(width      = page_width,
           height     = page_height,
           xgrid      = 0,
           ygrid      = 0,
           showGuides = FALSE)

## Panel labels ---------------------------------------------------------------
plotText("A", x = panelA_x - 0.10, y = panelA_y - 0.10,
         just = c("left", "top"), fontsize = 14, fontface = "bold")
plotText("B", x = panelB_x - 0.30, y = panelB_y - 0.10,
         just = c("left", "top"), fontsize = 14, fontface = "bold")
plotText("C", x = panelC_x - 0.30, y = panelC_y - 0.10,
         just = c("left", "top"), fontsize = 14, fontface = "bold")


## PANEL A: Locus plot --------------------------------------------------------

## Shared pgParams — arguments must always be passed directly to each function
## (pgParams fields are not accessible via $)
p <- pgParams(
  assembly   = "hg38",
  resolution = resolution,
  chrom      = locus_chrom,
  chromstart = locus_chromstart,
  chromend   = locus_chromend,
  zrange     = zrange_hic,
  norm       = norm_method,
  x          = panelA_hic_x,
  width      = panelA_w
)

## 1) Untreated Hi-C
y_hic_ctrl <- panelA_y

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
         x = panelA_hic_x, y = y_hic_ctrl,
         just = c("left", "top"), fontsize = 7, fontcolor = gray_label)

## 2) Sorbitol Hi-C
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
         x = panelA_hic_x, y = y_hic_sorb,
         just = c("left", "top"), fontsize = 7, fontcolor = gray_label)

## Hi-C legend (right of first map)
annoHeatmapLegend(ctrl_hic,
                  orientation   = "v",
                  fontsize      = 6,
                  fontcolor     = "black",
                  digits        = 2,
                  x             = panelA_hic_x + panelA_w + 0.05,
                  y             = y_hic_ctrl,
                  width         = 0.07,
                  height        = 0.65,
                  just          = c("left", "top"),
                  default.units = "inches")

## 3) RNA-seq signal tracks
y_rna_start <- y_hic_sorb + hic_height + 0.18

for (j in seq_along(rna_tracks)) {
  yj    <- y_rna_start + (j - 1) * (rna_height + rna_gap)
  track <- rna_tracks[[j]]
  
  if (file.exists(track$fwd)) {
    plotSignal(data = track$fwd, params = p,
               x = panelA_hic_x, y = yj,
               height = rna_height,
               linecolor = track$color, fill = track$color,
               scale = FALSE, range = signalRange)
  }
  
  if (file.exists(track$rev)) {
    plotSignal(data = track$rev, params = p,
               x = panelA_hic_x, y = yj,
               height = rna_height,
               linecolor = track$color,
               scale = FALSE, range = signalRange)
  }
  
  plotText(track$label,
           x = panelA_hic_x, y = yj + 0.06,
           just = c("left", "bottom"),
           fontsize = 6, fontcolor = track$color)
}

## Vertical RNA-seq label
n_rna     <- length(rna_tracks)
rna_end_y <- y_rna_start + n_rna * (rna_height + rna_gap) - rna_gap
rna_mid_y <- (y_rna_start + rna_end_y) / 2

plotText("RNA-seq",
         rot = 90,
         x = panelA_hic_x - 0.12, y = rna_mid_y,
         fontsize = 7, fontcolor = gray_label)

## Signal range label
range_label <- sprintf("%d\u2013%d",
                       round(signalRange[1]),
                       round(signalRange[2]))
plotText(range_label,
         x = panelA_hic_x + panelA_w, y = y_rna_start - 0.01,
         just = c("right", "bottom"),
         fontsize = 6, fontcolor = gray_label)

## 4) Gene track
y_genes <- rna_end_y + 0.10

plotGenes(
  params   = p,
  chrom    = locus_chrom,
  x        = panelA_hic_x,
  y        = y_genes,
  height   = gene_height,
  fontsize = 7,
  geneHighlights = data.frame(gene = "SOX8", color = color_sorbitol)
)

## 5) Genome label
y_genome <- y_genes + gene_height + 0.05

plotGenomeLabel(
  params   = p,
  x        = panelA_hic_x,
  y        = y_genome,
  length   = panelA_w,
  fontsize = 7
)

## 6) Loop anchor highlights
if (length(target_loop) >= 1L) {
  loop <- target_loop[1L]
  hl_y <- y_hic_ctrl
  hl_h <- y_genome - hl_y
  
  annoHighlight(
    plot = ctrl_hic, fill = "lightgrey", alpha = 0.5,
    chrom      = as.character(seqnames(anchors(loop, "first"))),
    chromstart = start(anchors(loop, "first")),
    chromend   = start(anchors(loop, "first")) + resolution,
    x = panelA_hic_x, y = hl_y,
    width = panelA_w, height = hl_h,
    just = c("left", "top"), default.units = "inches"
  )
  
  annoHighlight(
    plot = ctrl_hic, fill = "lightgrey", alpha = 0.5,
    chrom      = as.character(seqnames(anchors(loop, "second"))),
    chromstart = end(anchors(loop, "second")) - resolution,
    chromend   = end(anchors(loop, "second")),
    x = panelA_hic_x, y = hl_y,
    width = panelA_w, height = hl_h,
    just = c("left", "top"), default.units = "inches"
  )
}

## PANEL B: Histogram ---------------------------------------------------------

plotGG(
  plot   = p_hist,
  x      = panelB_x,
  y      = panelB_y,
  width  = panelB_w,
  height = panelB_h,
  just   = c("left", "top")
)

## PANEL C: Barplot -----------------------------------------------------------

plotGG(
  plot   = p_bar,
  x      = panelC_x,
  y      = panelC_y,
  width  = panelC_w,
  height = panelC_h,
  just   = c("left", "top")
)


# Save outputs ----------------------------------------------------------------

dev.off()
message("Multi-panel figure saved: ", output_pdf)


# Session info ----------------------------------------------------------------

sessionInfo()