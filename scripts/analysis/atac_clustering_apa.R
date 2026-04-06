# ATAC Peak Clustering APA -------------------------------------------------
# Author:      JP Flores
# Date:        2026-04-02
# Project:     STRS
# Description: APA analysis asking whether open chromatin regions cluster in
#              3D space after sorbitol-induced hyperosmotic stress, motivated
#              by Cai et al. 2019 (Nat Cell Biol) which showed that accessible
#              chromatin regions cluster following sorbitol treatment. Gained
#              ATAC peaks (significantly increased under sorbitol) are used as
#              anchors, ranked by log2FoldChange with a ceiling of top_n peaks.
# Input:       - ATAC peak counts: data/processed/atac/output/peaks/YAPP_HEK_1_peakCounts.tsv
#              - Hi-C maps:        data/processed/hic/maps/
# Output:      - plots/atac_clustering_apa.pdf
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

top_n        <- 1000         ## Max gained ATAC peaks to use (ceiling; ranked by log2FoldChange)
max_dist     <- 1e6          ## Maximum genomic distance between pairs (bp)
min_dist     <- 25e3         ## Minimum genomic distance between pairs (bp)
bin_size     <- 10e3         ## Hi-C bin size (bp)
buffer       <- 10           ## APA matrix buffer (bins on each side)
min_count    <- 15           ## Pre-filtering threshold: min mean normalized count

apa_width    <- 0.75         ## APA panel width (inches)
apa_height   <- 0.75         ## APA panel height (inches)
col_gap      <- 0.1          ## Gap between APA columns (inches)
left_margin  <- 0.25         ## Left page margin (inches)
top_margin   <- 0.75         ## Top page margin within each panel (inches)

project_dir  <- "/work/users/j/p/jpflores/projects/STRS"

atac_counts_path <- file.path(project_dir,
                              "data/processed/atac/output/peaks/YAPP_HEK_1_peakCounts.tsv")

hic_ctrl <- file.path(project_dir,
                      "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic")
hic_sorb <- file.path(project_dir,
                      "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic")

output_path <- file.path(project_dir, "plots/atac_clustering_apa.pdf")


# Libraries ---------------------------------------------------------------

library(DESeq2)
library(GenomicRanges)
library(InteractionSet)
library(mariner)
library(plotgardener)
library(RColorBrewer)
library(tidyverse)


# Utility scripts ---------------------------------------------------------

source(file.path(project_dir, "scripts/utils/calculate_apa_score.R"))


# Load data ---------------------------------------------------------------

## Raw ATAC peak count table (peaks x samples)
atac_raw <- read_tsv(atac_counts_path, show_col_types = FALSE) |>
  dplyr::select(-starts_with("Unnamed"))

cat("Total ATAC peaks loaded:", nrow(atac_raw), "\n")


# Analysis ----------------------------------------------------------------

## Build count matrix for DESeq2 -------------------------------------------

count_matrix <- atac_raw |>
  dplyr::select(starts_with("YAPP")) |>
  as.matrix()

storage.mode(count_matrix) <- "integer"

rownames(count_matrix) <- paste0(
  atac_raw$chr, ":",
  atac_raw$start, "-",
  atac_raw$stop
)

col_data <- DataFrame(
  treatment = factor(
    c("ctrl", "ctrl", "ctrl", "ctrl", "sorb", "sorb", "sorb", "sorb"),
    levels = c("ctrl", "sorb")
  ),
  replicate = factor(c(1, 2, 3, 4, 1, 2, 3, 4)),
  row.names = colnames(count_matrix)
)

## Run DESeq2 --------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = col_data,
  design    = ~ replicate + treatment
)

## Fix size factors to 1: hyperosmotic stress alters global chromatin
## accessibility, so normalization assuming count stability would mask
## true genome-wide changes
sizeFactors(dds) <- rep(1, ncol(dds))

keep <- rowMeans(counts(dds, normalized = TRUE)) >= min_count
dds  <- dds[keep, ]
cat("Peaks retained after pre-filtering:", sum(keep), "\n")

dds <- DESeq(dds)

res <- lfcShrink(
  dds,
  coef = "treatment_sorb_vs_ctrl",
  type = "apeglm"
)

res_df <- as.data.frame(res) |>
  tibble::rownames_to_column("peak_id")

cat("Total tested peaks:", nrow(res_df), "\n")

## Select gained ATAC peaks ------------------------------------------------

## Motivated by Cai et al. 2019: open chromatin regions cluster *after*
## sorbitol treatment, so we select peaks that are significantly gained
## under sorbitol. Peaks are ranked by log2FoldChange (strongest gainers
## first) with a ceiling of top_n in case of a very large gained set.
top_peaks <- res_df |>
  dplyr::filter(padj < 0.05, log2FoldChange > 0) |>
  dplyr::arrange(desc(log2FoldChange)) |>
  dplyr::slice_head(n = top_n)

cat("Gained ATAC peaks selected:", nrow(top_peaks),
    "(padj < 0.05, log2FC > 1, top", top_n, "by log2FC)\n")

## Convert peak IDs to GRanges ---------------------------------------------

peak_id_to_granges <- function(peak_ids) {
  coords <- str_match(peak_ids, "^(chr[^:]+):(\\d+)-(\\d+)$")
  GRanges(
    seqnames = coords[, 2],
    ranges   = IRanges(
      start = as.integer(coords[, 3]),
      end   = as.integer(coords[, 4])
    )
  )
}

cres <- peak_id_to_granges(top_peaks$peak_id)

## Generate ATAC-ATAC pairs ------------------------------------------------

make_atac_pairs <- function(cres) {
  cres  <- keepStandardChromosomes(cres, pruning.mode = "coarse")
  hits  <- findOverlaps(cres, cres, maxgap = max_dist, ignore.strand = TRUE)
  hits  <- hits[queryHits(hits) < subjectHits(hits)]
  
  pairs <- GInteractions(
    anchor1 = cres[queryHits(hits)],
    anchor2 = cres[subjectHits(hits)]
  )
  
  pairs <- pairs[pairdist(pairs) >= min_dist]
  
  ## Ensure anchor1 is always the upstream anchor
  flipped <- start(anchors(pairs, "first")) > start(anchors(pairs, "second"))
  pairs[flipped] <- swapAnchors(pairs[flipped])
  
  pairs <- assignToBins(pairs, binSize = bin_size)
  
  cat("ATAC-ATAC pairs after distance filtering:", length(pairs), "\n")
  pairs
}

pairs <- make_atac_pairs(cres)

## Pull Hi-C matrices and aggregate (APA) ----------------------------------

run_apa <- function(pairs, hic_file) {
  pairs |>
    pixelsToMatrices(buffer = buffer) |>
    removeShortPairs() |>
    pullHicMatrices(
      binSize = bin_size,
      files   = hic_file,
      half    = "upper",
      norm    = "NONE",
      matrix  = "observed"
    ) |>
    aggHicMatrices(FUN = sum)
}

apa_ctrl <- run_apa(pairs, hic_ctrl) / length(pairs)
apa_sorb <- run_apa(pairs, hic_sorb) / length(pairs)


# Visualization -----------------------------------------------------------

pal_hic <- colorRampPalette(brewer.pal(9, "YlGnBu"))

x_ctrl   <- left_margin
x_sorb   <- x_ctrl + apa_width + col_gap
x_center <- x_ctrl + apa_width + col_gap / 2   ## midpoint between the two matrices

## Helper: add "untreated" / "+ sorbitol" column headers
add_column_headers <- function(y_header) {
  for (lbl in list(
    list(label = "untreated",  x = x_ctrl + apa_width / 2),
    list(label = "+ sorbitol", x = x_sorb + apa_width / 2)
  )) {
    plotText(
      label         = lbl$label,
      fontsize      = 6,
      fontcolor     = "black",
      x             = lbl$x,
      y             = y_header,
      just          = c("center", "bottom"),
      default.units = "inches"
    )
  }
}

## Helper: plot one APA row with a centered title above the matrices
plot_apa_row <- function(apa_ctrl_mat, apa_sorb_mat,
                         y_pos, zrange, row_title, n) {
  
  ## Centered row title
  plotText(
    label         = row_title,
    fontsize      = 7,
    fontface      = "italic",
    fontcolor     = "black",
    x             = x_center,
    y             = y_pos - 0.12,
    just          = c("center", "bottom"),
    default.units = "inches"
  )
  
  ## n label below row title
  plotText(
    label         = paste0("n = ", n),
    fontsize      = 5,
    fontcolor     = "grey40",
    x             = x_center,
    y             = y_pos - 0.02,
    just          = c("center", "bottom"),
    default.units = "inches"
  )
  
  apaParams <- pgParams(
    width    = apa_width,
    height   = apa_height,
    fontsize = 5,
    palette  = pal_hic
  )
  
  plotMatrix(data = apa_ctrl_mat, x = x_ctrl, y = y_pos,
             zrange = zrange, params = apaParams)
  
  plotText(
    label     = sprintf("%.2f", calculate_apa_score(apa_ctrl_mat)),
    fontsize  = 5,
    fontcolor = "black",
    x         = x_ctrl + apa_width - 0.05,
    y         = y_pos + 0.03,
    just      = c("right", "top")
  )
  
  plotMatrix(data = apa_sorb_mat, x = x_sorb, y = y_pos,
             zrange = zrange, params = apaParams) |>
    annoHeatmapLegend(
      x             = x_sorb + apa_width + 0.05,
      y             = y_pos,
      width         = 0.05,
      height        = apa_height,
      fontsize      = 5,
      fontcolor     = "black",
      just          = c("left", "top"),
      default.units = "inches"
    )
  
  plotText(
    label     = sprintf("%.2f", calculate_apa_score(apa_sorb_mat)),
    fontsize  = 5,
    fontcolor = "black",
    x         = x_sorb + apa_width - 0.05,
    y         = y_pos + 0.03,
    just      = c("right", "top")
  )
}

## Page dimensions
panel_height <- top_margin + apa_height + 0.5
page_width   <- left_margin + 2 * apa_width + col_gap + 0.35
page_height  <- panel_height + 0.25

zrange <- c(0, max(apa_ctrl, apa_sorb, na.rm = TRUE))


# Save outputs ------------------------------------------------------------

pdf(output_path, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = FALSE)

add_column_headers(y_header = top_margin - 0.05)

plot_apa_row(
  apa_ctrl_mat = apa_ctrl,
  apa_sorb_mat = apa_sorb,
  y_pos        = top_margin,
  zrange       = zrange,
  row_title    = "Gained ATAC peaks",
  n            = length(pairs)
)

dev.off()


# Session info ------------------------------------------------------------

sessionInfo()