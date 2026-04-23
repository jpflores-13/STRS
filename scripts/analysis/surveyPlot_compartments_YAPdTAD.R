# surveyPlot_compartments_YAPvsdTAD.R -------------------------------------
# Author:      JP Flores
# Date:        2026-04-19
# Project:     STRS
# Description: Draws survey plots for 100 randomly sampled autosomal 5 Mb
#              windows. Each page shows a 1x2 grid of Hi-C squares: left
#              column = control (upper triangle eGFP-YAP, lower triangle
#              dTAD), right column = sorbitol (upper triangle eGFP-YAP,
#              lower triangle dTAD). Each column has stacked PC1 eigenvector
#              tracks (eGFP-YAP then dTAD), its own gene track, and genome
#              label. zrange is shared across both squares per page.
# Input:       data/processed/hic/maps/  (.hic files)
#              data/processed/hic/eigenvectors/control.cis.vecs.tsv
#              data/processed/hic/eigenvectors/sorbitol.cis.vecs.tsv
#              data/processed/hic/eigenvectors_YAPdTAD/control.cis.vecs.tsv
#              data/processed/hic/eigenvectors_YAPdTAD/sorbitol.cis.vecs.tsv
# Output:      plots/surveyPlot_compartments_YAPvsdTAD.pdf
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

project_dir   <- "/work/users/j/p/jpflores/projects/STRS"
hic_dir       <- file.path(project_dir, "data/processed/hic/maps")
eigvec_dir    <- file.path(project_dir, "data/processed/hic/eigenvectors")
eigvec_dtad   <- file.path(project_dir, "data/processed/hic/eigenvectors_YAPdTAD")
output_dir    <- file.path(project_dir, "plots")

## eGFP-YAP megaMaps
hic_yap_ctrl  <- file.path(hic_dir, "YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic")
hic_yap_sorb  <- file.path(hic_dir, "YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic")

## YAPdTAD single-rep files
hic_dtad_ctrl <- file.path(hic_dir, "YAPP_HEK293_eGFP-YAPdTAD_Cai_control_1_1_inter_30.hic")
hic_dtad_sorb <- file.path(hic_dir, "YAPP_HEK293_eGFP-YAPdTAD_Cai_sorbitol_1_1_inter_30.hic")

output_pdf    <- file.path(output_dir, "surveyPlot_compartments_YAPvsdTAD.pdf")

## Sampling parameters
n_regions     <- 100L
window_size   <- 5e6        # 5 Mb per survey window
resolution    <- 10e3       # 10 kb -- matches FigureSX.R
rng_seed      <- 42L

## Autosomes to sample from
autosomes     <- paste0("chr", 1:22)

## Hi-C display
norm          <- "SCALE"
zrange        <- c(0, 100)   # shared across both squares per page

## Page dimensions — landscape to accommodate two columns
page_width    <- 11.0
page_height   <- 8.0

## Track dimensions (inches)
sq_width      <- 4.0        # Hi-C square width = height
eig_height    <- 0.35
gene_height   <- 0.45
gap           <- 0.05
legend_width  <- 0.1
legend_height <- 1.5
legend_gap    <- 0.1

## Column x positions
x_col1        <- 0.5                              # left square
x_col2        <- x_col1 + sq_width + 0.6         # right square

## Colors -- consistent with existing STRS figures
col_yap       <- "#4C72B0"   # eGFP-YAP  (control color)
col_dtad      <- "#DD8452"   # YAPdTAD   (sorbitol color)


# Libraries ---------------------------------------------------------------

library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(GenomicRanges)
library(IRanges)
library(plotgardener)
library(readr)


# Load data ---------------------------------------------------------------

## eGFP-YAP eigenvectors
eig_yap_ctrl <- read_tsv(
  file.path(eigvec_dir, "control.cis.vecs.tsv"),
  show_col_types = FALSE
)
eig_yap_sorb <- read_tsv(
  file.path(eigvec_dir, "sorbitol.cis.vecs.tsv"),
  show_col_types = FALSE
)

## YAPdTAD eigenvectors
eig_dtad_ctrl <- read_tsv(
  file.path(eigvec_dtad, "control.cis.vecs.tsv"),
  show_col_types = FALSE
)
eig_dtad_sorb <- read_tsv(
  file.path(eigvec_dtad, "sorbitol.cis.vecs.tsv"),
  show_col_types = FALSE
)


# Sample random genomic windows -------------------------------------------

chrom_lens <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[autosomes]

sampling_gr <- GRanges(
  seqnames = names(chrom_lens),
  ranges   = IRanges(
    start = window_size + 1L,
    end   = as.integer(chrom_lens) - window_size
  )
)
sampling_gr <- sampling_gr[width(sampling_gr) > 0]

chrom_weights <- width(sampling_gr) / sum(width(sampling_gr))

set.seed(rng_seed)

sampled_chroms <- sample(
  as.character(seqnames(sampling_gr)),
  size    = n_regions,
  replace = TRUE,
  prob    = chrom_weights
)

sampled_starts <- vapply(sampled_chroms, function(chrom) {
  gr <- sampling_gr[as.character(seqnames(sampling_gr)) == chrom]
  sample(start(gr):end(gr), size = 1L)
}, integer(1L))

survey_windows <- GRanges(
  seqnames = sampled_chroms,
  ranges   = IRanges(
    start = sampled_starts,
    end   = sampled_starts + window_size - 1L
  )
)


# Helper function ---------------------------------------------------------

## draw_eig_track() --------------------------------------------------------
## Draws a plotSignal eigenvector track for one condition/genotype.

draw_eig_track <- function(eig_win, x, y, width, chrom_val, chrom_start,
                           chrom_end, fill_color, label, eig_range) {
  
  
  plotText(
    label     = label,
    x         = x,
    y         = y,
    just      = c("left", "top"),
    fontsize  = 6,
    fontcolor = fill_color,
    default.units = "inches"
  )
  
  plotText(
    label     = "Eigenvalue",
    x         = x - 0.1,
    y         = y + eig_height / 2,
    rot       = 90,
    just      = "center",
    fontsize  = 4,
    default.units = "inches"
  )
  
  if (nrow(eig_win) == 0L || all(is.na(eig_win$E1))) return(invisible(NULL))
  
  plotSignal(
    data      = eig_win |>
      select(chrom, start, end, score = E1) |>
      filter(!is.na(score)) |>
      mutate(end = end - 1L),
    chrom     = chrom_val,
    chromstart = chrom_start,
    chromend  = chrom_end,
    assembly  = "hg38",
    x         = x,
    width     = width,
    y         = y,
    height    = eig_height,
    range     = eig_range,
    fill      = fill_color,
    linecolor = fill_color,
    baseline  = TRUE,
    default.units = "inches"
  )
}


# Visualization -----------------------------------------------------------

## Page layout per column -- top-to-bottom (inches):
##
##  condition label
##  Hi-C square + heatmap legend      (sq_width x sq_width)
##  PC1 (eGFP-YAP)                    (eig_height)
##  PC1 (YAPdTAD)                      (eig_height)
##  gene track                         (gene_height)
##  genome label

y_label   <- 0.25
y_hic     <- 0.45
y_eig_yap <- y_hic + sq_width + gap
y_eig_dtad <- y_eig_yap + eig_height + gap
y_genes   <- y_eig_dtad + eig_height + gap * 2
y_genome  <- y_genes + gene_height + gap

pdf(file = output_pdf, width = page_width, height = page_height)

for (i in seq_along(survey_windows)) {
  
  chrom     <- as.character(seqnames(survey_windows))[i]
  win_start <- start(survey_windows)[i]
  win_end   <- end(survey_windows)[i]
  
  ## Subset eigenvectors to this window
  sub_eig <- function(eig_df) {
    eig_df |>
      filter(.data$chrom == !!chrom,
             .data$end   >= win_start,
             .data$start <= win_end)
  }
  
  eig_yap_ctrl_win  <- sub_eig(eig_yap_ctrl)
  eig_yap_sorb_win  <- sub_eig(eig_yap_sorb)
  eig_dtad_ctrl_win <- sub_eig(eig_dtad_ctrl)
  eig_dtad_sorb_win <- sub_eig(eig_dtad_sorb)
  
  ## Shared eigenvector y-axis range across all four tracks this page
  all_eig_vals <- c(eig_yap_ctrl_win$E1, eig_yap_sorb_win$E1,
                    eig_dtad_ctrl_win$E1, eig_dtad_sorb_win$E1)
  all_eig_vals <- all_eig_vals[is.finite(all_eig_vals)]
  eig_range    <- if (length(all_eig_vals) > 0) {
    c(min(all_eig_vals), max(all_eig_vals))
  } else {
    c(-1, 1)
  }
  
  pageCreate(
    width      = page_width,
    height     = page_height,
    xgrid      = 0,
    ygrid      = 0,
    showGuides = FALSE
  )
  
  ## -- Column labels ------------------------------------------------------
  plotText(label = "control",    x = x_col1, y = y_label,
           just = c("left", "top"), fontsize = 9,
           default.units = "inches")
  plotText(label = "+ sorbitol", x = x_col2, y = y_label,
           just = c("left", "top"), fontsize = 9,
           default.units = "inches")
  
  ## -- Left square: control -----------------------------------------------
  ## Upper triangle = eGFP-YAP, lower triangle = YAPdTAD
  sq_ctrl <- plotHicSquare(
    data       = hic_yap_ctrl,
    chrom      = chrom, chromstart = win_start, chromend = win_end,
    assembly   = "hg38", resolution = resolution,
    norm       = norm,   zrange     = zrange,
    x = x_col1, y = y_hic, width = sq_width, height = sq_width,
    half = "top", just = c("left", "top"), default.units = "inches"
  )
  
  plotHicSquare(
    data       = hic_dtad_ctrl,
    chrom      = chrom, chromstart = win_start, chromend = win_end,
    assembly   = "hg38", resolution = resolution,
    norm       = norm,   zrange     = zrange,
    x = x_col1, y = y_hic, width = sq_width, height = sq_width,
    half = "bottom", just = c("left", "top"), default.units = "inches"
  )
  
  annoHeatmapLegend(
    sq_ctrl,
    x = x_col1 + sq_width + legend_gap, y = y_hic,
    width = legend_width, height = legend_height,
    orientation = "v", fontsize = 7, fontcolor = "black", digits = 2,
    just = c("left", "top"), default.units = "inches"
  )
  
  ## Genotype labels
  plotText(label = "eGFP-YAP", x = x_col1 + 0.05, y = y_hic + 0.05,
           just = c("left", "top"), fontsize = 6, fontcolor = col_yap,
           default.units = "inches")
  plotText(label = "YAPdTAD",  x = x_col1 + sq_width - 0.05,
           y = y_hic + sq_width - 0.05,
           just = c("right", "bottom"), fontsize = 6, fontcolor = col_dtad,
           default.units = "inches")
  
  ## -- Left column eigenvector tracks ------------------------------------
  draw_eig_track(eig_yap_ctrl_win,
                 x = x_col1, y = y_eig_yap, width = sq_width,
                 chrom_val = chrom, chrom_start = win_start, chrom_end = win_end,
                 fill_color = col_yap, label = "eGFP-YAP",
                 eig_range = eig_range)
  
  draw_eig_track(eig_dtad_ctrl_win,
                 x = x_col1, y = y_eig_dtad, width = sq_width,
                 chrom_val = chrom, chrom_start = win_start, chrom_end = win_end,
                 fill_color = col_dtad, label = "YAPdTAD",
                 eig_range = eig_range)
  
  ## -- Left column genes + genome label ----------------------------------
  plotGenes(
    chrom = chrom, chromstart = win_start, chromend = win_end,
    assembly = "hg38", x = x_col1, y = y_genes,
    width = sq_width, height = gene_height,
    just = c("left", "top"), default.units = "inches"
  )
  
  annoGenomeLabel(
    plot = sq_ctrl, scale = "Mb", axis = "x",
    x = x_col1, y = y_genome,
    just = c("left", "top"), default.units = "inches"
  )
  
  ## -- Right square: sorbitol --------------------------------------------
  sq_sorb <- plotHicSquare(
    data       = hic_yap_sorb,
    chrom      = chrom, chromstart = win_start, chromend = win_end,
    assembly   = "hg38", resolution = resolution,
    norm       = norm,   zrange     = zrange,
    x = x_col2, y = y_hic, width = sq_width, height = sq_width,
    half = "top", just = c("left", "top"), default.units = "inches"
  )
  
  plotHicSquare(
    data       = hic_dtad_sorb,
    chrom      = chrom, chromstart = win_start, chromend = win_end,
    assembly   = "hg38", resolution = resolution,
    norm       = norm,   zrange     = zrange,
    x = x_col2, y = y_hic, width = sq_width, height = sq_width,
    half = "bottom", just = c("left", "top"), default.units = "inches"
  )
  
  annoHeatmapLegend(
    sq_sorb,
    x = x_col2 + sq_width + legend_gap, y = y_hic,
    width = legend_width, height = legend_height,
    orientation = "v", fontsize = 7, fontcolor = "black", digits = 2,
    just = c("left", "top"), default.units = "inches"
  )
  
  ## Genotype labels
  plotText(label = "eGFP-YAP", x = x_col2 + 0.05, y = y_hic + 0.05,
           just = c("left", "top"), fontsize = 6, fontcolor = col_yap,
           default.units = "inches")
  plotText(label = "YAPdTAD",  x = x_col2 + sq_width - 0.05,
           y = y_hic + sq_width - 0.05,
           just = c("right", "bottom"), fontsize = 6, fontcolor = col_dtad,
           default.units = "inches")
  
  ## -- Right column eigenvector tracks -----------------------------------
  draw_eig_track(eig_yap_sorb_win,
                 x = x_col2, y = y_eig_yap, width = sq_width,
                 chrom_val = chrom, chrom_start = win_start, chrom_end = win_end,
                 fill_color = col_yap, label = "eGFP-YAP",
                 eig_range = eig_range)
  
  draw_eig_track(eig_dtad_sorb_win,
                 x = x_col2, y = y_eig_dtad, width = sq_width,
                 chrom_val = chrom, chrom_start = win_start, chrom_end = win_end,
                 fill_color = col_dtad, label = "YAPdTAD",
                 eig_range = eig_range)
  
  ## -- Right column genes + genome label ---------------------------------
  plotGenes(
    chrom = chrom, chromstart = win_start, chromend = win_end,
    assembly = "hg38", x = x_col2, y = y_genes,
    width = sq_width, height = gene_height,
    just = c("left", "top"), default.units = "inches"
  )
  
  annoGenomeLabel(
    plot = sq_sorb, scale = "Mb", axis = "x",
    x = x_col2, y = y_genome,
    just = c("left", "top"), default.units = "inches"
  )
}

dev.off()


# Session info ------------------------------------------------------------

sessionInfo()