# surveyPlot_compartments.R -----------------------------------------------
# Author:      JP Flores
# Date:        2026-04-17
# Project:     STRS
# Description: Draws survey plots for 100 randomly sampled autosomal 2 Mb
#              windows. Each page mirrors the FigureSX.R layout per condition:
#              eigenvector track -> Hi-C rectangle + heatmap legend, repeated
#              for control then sorbitol, followed by gene track and genome
#              label at the bottom.
# Input:       data/processed/hic/maps/  (.hic files)
#              data/processed/hic/eigenvectors/control.cis.vecs.tsv
#              data/processed/hic/eigenvectors/sorbitol.cis.vecs.tsv
# Output:      plots/surveyPlot_compartments.pdf
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

project_dir  <- "/work/users/j/p/jpflores/projects/STRS"
hic_dir      <- file.path(project_dir, "data/processed/hic/maps")
eigvec_dir   <- file.path(project_dir, "data/processed/hic/eigenvectors")
output_dir   <- file.path(project_dir, "plots")

hic_control  <- file.path(
  hic_dir, "YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic"
)
hic_sorbitol <- file.path(
  hic_dir, "YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic"
)

output_pdf   <- file.path(output_dir, "surveyPlot_compartments.pdf")

## Sampling parameters
n_regions    <- 100L
window_size  <- 5e6         # 5 Mb per survey window
resolution   <- 10e3        # 10 kb -- matches FigureSX.R
rng_seed     <- 42L

## Autosomes to sample from
autosomes    <- paste0("chr", 1:22)

## Hi-C display
zrange       <- c(0, 50)
norm         <- "SCALE"

## Page dimensions
page_width   <- 5.75
page_height  <- 8.0

## Track dimensions (inches)
hic_height   <- 2.0
hic_width    <- 5.0
eig_height   <- 0.4
gene_height  <- 0.5
gap          <- 0.05
x_left       <- 0.25

## Colors -- consistent with domain_analysis.R and existing figures
col_control  <- "#4C72B0"
col_sorbitol <- "#DD8452"
col_A        <- "#C44E52"   # positive PC1 -> A compartment
col_B        <- "#4C72B0"   # negative PC1 -> B compartment


# Libraries ---------------------------------------------------------------

library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(GenomicRanges)
library(IRanges)
library(plotgardener)
library(readr)


# Load data ---------------------------------------------------------------

## Eigenvector TSVs from cooltools eigs-cis; columns: chrom, start, end, E1
eig_ctrl <- read_tsv(
  file.path(eigvec_dir, "control.cis.vecs.tsv"),
  show_col_types = FALSE
)

eig_sorb <- read_tsv(
  file.path(eigvec_dir, "sorbitol.cis.vecs.tsv"),
  show_col_types = FALSE
)


# Sample random genomic windows -------------------------------------------

## Chromosome lengths from BSgenome
chrom_lens <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[autosomes]

## Valid sampling range: shrink each end by window_size to avoid edge windows
sampling_gr <- GRanges(
  seqnames = names(chrom_lens),
  ranges   = IRanges(
    start = window_size + 1L,
    end   = as.integer(chrom_lens) - window_size
  )
)
sampling_gr <- sampling_gr[width(sampling_gr) > 0]

## Weight chromosomes proportional to available sampling space
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
## Draws a plotSignal-style eigenvector track using plotgardener's
## plotSignal(), matching FigureSX.R exactly. Positive PC1 -> A compartment,
## negative -> B compartment. Uses a shared eig_range across both conditions
## per page for a fair comparison.

draw_eig_track <- function(eig_win, y, params_shared, fill_color,
                           label, eig_range) {
  
  plotText(
    label     = label,
    x         = x_left,
    y         = y,
    just      = c("left", "top"),
    fontsize  = 6,
    fontcolor = fill_color,
    default.units = "inches"
  )
  
  plotText(
    label     = "Eigenvalue",
    x         = x_left - 0.1,
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
    params    = params_shared,
    x         = x_left,
    y         = y,
    width     = hic_width,
    height    = eig_height,
    range     = eig_range,
    fill      = fill_color,
    linecolor = fill_color,
    baseline  = TRUE,
    default.units = "inches"
  )
}


# Visualization -----------------------------------------------------------

## Page layout -- top-to-bottom (inches), mirroring FigureSX.R:
##
##  y_ctrl_eig    control eigenvector track
##  y_ctrl_hic    control Hi-C rectangle    [heatmap legend]
##  y_sorb_eig    sorbitol eigenvector track
##  y_sorb_hic    sorbitol Hi-C rectangle   [heatmap legend]
##  y_genes       gene track
##  y_genome      genome label

y_ctrl_eig  <- 0.5
y_ctrl_hic  <- y_ctrl_eig + eig_height + gap
y_sorb_eig  <- y_ctrl_hic + hic_height + gap * 4
y_sorb_hic  <- y_sorb_eig + eig_height + gap
y_genes     <- y_sorb_hic + hic_height + gap * 4
y_genome    <- y_genes    + gene_height + gap

pdf(file = output_pdf, width = page_width, height = page_height)

for (i in seq_along(survey_windows)) {
  
  chrom     <- as.character(seqnames(survey_windows))[i]
  win_start <- start(survey_windows)[i]
  win_end   <- end(survey_windows)[i]
  
  ## Shared pgParams -- governs coordinate mapping for all tracks this page
  p <- pgParams(
    assembly   = "hg38",
    resolution = resolution,
    chrom      = chrom,
    chromstart = win_start,
    chromend   = win_end,
    zrange     = zrange,
    norm       = norm,
    x          = x_left,
    width      = hic_width,
    length     = hic_width
  )
  
  ## Subset eigenvectors to this window
  eig_ctrl_win <- eig_ctrl |>
    filter(.data$chrom == !!chrom,
           .data$end   >= win_start,
           .data$start <= win_end)
  
  eig_sorb_win <- eig_sorb |>
    filter(.data$chrom == !!chrom,
           .data$end   >= win_start,
           .data$start <= win_end)
  
  ## Shared y-axis range across both conditions for fair comparison
  eig_vals  <- c(eig_ctrl_win$E1, eig_sorb_win$E1)
  eig_vals  <- eig_vals[is.finite(eig_vals)]
  eig_range <- if (length(eig_vals) > 0) {
    c(min(eig_vals), max(eig_vals))
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
  
  ## -- Control eigenvector track ----------------------------------------
  draw_eig_track(
    eig_win    = eig_ctrl_win,
    y          = y_ctrl_eig,
    params_shared = p,
    fill_color = col_control,
    label      = "untreated",
    eig_range  = eig_range
  )
  
  ## -- Control Hi-C rectangle -------------------------------------------
  ctrl_plot <- plotHicRectangle(
    data   = hic_control,
    params = p,
    height = hic_height,
    y      = y_ctrl_hic
  )
  
  annoHeatmapLegend(
    ctrl_plot,
    orientation   = "v",
    fontsize      = 8,
    fontcolor     = "black",
    digits        = 2,
    x             = 5.5,
    y             = y_ctrl_hic,
    width         = 0.1,
    height        = 1.5,
    just          = c("left", "top"),
    default.units = "inches"
  )
  
  ## -- Sorbitol eigenvector track ---------------------------------------
  draw_eig_track(
    eig_win    = eig_sorb_win,
    y          = y_sorb_eig,
    params_shared = p,
    fill_color = col_sorbitol,
    label      = "+ sorbitol",
    eig_range  = eig_range
  )
  
  ## -- Sorbitol Hi-C rectangle ------------------------------------------
  sorb_plot <- plotHicRectangle(
    data   = hic_sorbitol,
    params = p,
    height = hic_height,
    y      = y_sorb_hic
  )
  
  annoHeatmapLegend(
    sorb_plot,
    orientation   = "v",
    fontsize      = 8,
    fontcolor     = "black",
    digits        = 2,
    x             = 5.5,
    y             = y_sorb_hic,
    width         = 0.1,
    height        = 1.5,
    just          = c("left", "top"),
    default.units = "inches"
  )
  
  ## -- Gene track -------------------------------------------------------
  plotGenes(
    params = p,
    x      = x_left,
    y      = y_genes,
    height = gene_height
  )
  
  ## -- Genome label -----------------------------------------------------
  plotGenomeLabel(params = p, x = x_left, y = y_genome)
}

dev.off()


# Session info ------------------------------------------------------------

sessionInfo()