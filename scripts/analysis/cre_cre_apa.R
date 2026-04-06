# CRE-CRE APA Analysis ----------------------------------------------------
# Author:      JP Flores
# Date:        2026-04-02
# Project:     STRS
# Description: APA analysis of Hi-C contact frequency between all pairs of
#              H3K27ac CRE peaks, as requested by reviewers (akin to Friman
#              Genome Research 2023). All confidently detected H3K27ac peaks
#              (baseMean >= min_count) are used — no differential filter —
#              with and without an ATAC accessibility filter. Pairs within
#              1 Mb are aggregated in Hi-C maps +/- sorbitol treatment.
# Input:       - H3K27ac DESeq2 results: data/processed/cutntag/deseq2/diff_H3K27ac_counts.rds
#              - ATAC peak counts:       data/processed/atac/output/peaks/YAPP_HEK_1_peakCounts.tsv
#              - Hi-C maps:              data/processed/hic/maps/
# Output:      - plots/cre_cre_apa.pdf
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

min_count    <- 50           ## Min baseMean to retain a peak as confidently detected
top_n        <- 1000         ## Max CREs to use after threshold filter (ranked by baseMean)
max_dist     <- 1e6          ## Maximum genomic distance between CRE pairs (bp)
min_dist     <- 25e3         ## Minimum genomic distance between CRE pairs (bp)
bin_size     <- 10e3         ## Hi-C bin size (bp)
buffer       <- 10           ## APA matrix buffer (bins on each side)

apa_width    <- 0.75         ## APA panel width (inches)
apa_height   <- 0.75         ## APA panel height (inches)
col_gap      <- 0.1          ## Gap between APA columns (inches)
row_gap      <- 0.15         ## Gap between APA rows (inches)
left_margin  <- 0.25         ## Left page margin (inches) — no left-side labels needed
top_margin   <- 0.75         ## Top page margin within each panel (inches) — room for row titles + col headers

project_dir  <- "/work/users/j/p/jpflores/projects/STRS"

diff_h3k27ac_path <- file.path(project_dir,
                               "data/processed/cutntag/deseq2/diff_H3K27ac_counts.rds")
atac_counts_path  <- file.path(project_dir,
                               "data/processed/atac/output/peaks/YAPP_HEK_1_peakCounts.tsv")

hic_ctrl <- file.path(project_dir,
                      "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic")
hic_sorb <- file.path(project_dir,
                      "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic")

output_path <- file.path(project_dir, "plots/cre_cre_apa.pdf")


# Libraries ---------------------------------------------------------------

library(DESeq2)
library(GenomicRanges)
library(InteractionSet)
library(mariner)
library(plotgardener)
library(RColorBrewer)
library(plyranges)
library(tidyverse)


# Utility scripts ---------------------------------------------------------

source(file.path(project_dir, "scripts/utils/calculate_apa_score.R"))


# Load data ---------------------------------------------------------------

## Pre-computed H3K27ac DESeq2 results (DESeqDataSet or DESeqResults GRanges)
diff_h3k27ac <- readRDS(diff_h3k27ac_path)

## ATAC peak count table — used to optionally restrict CRE selection
atac_raw <- read_tsv(atac_counts_path, show_col_types = FALSE) |>
  dplyr::select(-starts_with("Unnamed"))


# Wrangle data ------------------------------------------------------------

## Build a flat data frame of H3K27ac peaks with genomic coordinates
res_df <- as.data.frame(mcols(diff_h3k27ac)) |>
  mutate(peak_id = paste0(
    as.character(seqnames(diff_h3k27ac)), ":",
    start(diff_h3k27ac), "-",
    end(diff_h3k27ac)
  ))

cat("Total H3K27ac peaks loaded:", nrow(res_df), "\n")

## Build ATAC GRanges for overlap-based filtering
atac_gr <- GRanges(
  seqnames = atac_raw$chr,
  ranges   = IRanges(start = atac_raw$start, end = atac_raw$stop)
)

## Restrict the H3K27ac peak table to those overlapping ATAC peaks
k27ac_gr    <- granges(diff_h3k27ac)
atac_hits   <- findOverlaps(k27ac_gr, atac_gr)
res_df_atac <- res_df[unique(queryHits(atac_hits)), ]

cat("H3K27ac peaks overlapping ATAC:", nrow(res_df_atac), "\n")


# Analysis ----------------------------------------------------------------

## Helper: retain all H3K27ac peaks with mean signal >= min_count, then
## cap at the top_n strongest by baseMean. The threshold excludes noise;
## the cap keeps the pair count tractable while still using the most
## confidently detected CREs — matching the reviewer's intent.
select_cres_by_threshold <- function(df) {
  selected <- df |>
    dplyr::filter(baseMean >= min_count) |>
    dplyr::arrange(desc(baseMean)) |>
    dplyr::slice_head(n = top_n)
  cat("CREs retained (baseMean >=", min_count, ", top", top_n, "):",
      nrow(selected), "out of", nrow(df), "total\n")
  selected
}

## Helper: convert peak_id strings ("chrN:start-end") to GRanges
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

## Helper: generate all valid intra-chromosomal CRE-CRE pairs within
## the specified distance range, snapped to the Hi-C bin grid
make_cre_pairs <- function(cres) {
  cres  <- keepStandardChromosomes(cres, pruning.mode = "coarse")
  hits  <- findOverlaps(cres, cres, maxgap = max_dist, ignore.strand = TRUE)
  hits  <- hits[queryHits(hits) < subjectHits(hits)]   # upper triangle only
  
  pairs <- GInteractions(
    anchor1 = cres[queryHits(hits)],
    anchor2 = cres[subjectHits(hits)]
  )
  
  pairs <- pairs[pairdist(pairs) >= min_dist]
  
  ## Ensure anchor1 is always the upstream anchor
  flipped <- start(anchors(pairs, "first")) > start(anchors(pairs, "second"))
  pairs[flipped] <- swapAnchors(pairs[flipped])
  
  pairs <- assignToBins(pairs, binSize = bin_size)
  
  cat("CRE-CRE pairs after distance filtering:", length(pairs), "\n")
  pairs
}

## Helper: pull Hi-C matrices and aggregate into an APA matrix
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

## Select all confidently detected CREs from the full peak set (no ATAC filter)
top_noatac      <- select_cres_by_threshold(res_df)
cres_noatac     <- peak_id_to_granges(top_noatac$peak_id)
pairs_noatac    <- make_cre_pairs(cres_noatac)

## Select all confidently detected CREs from ATAC-intersected peaks
top_atac        <- select_cres_by_threshold(res_df_atac)
cres_atac       <- peak_id_to_granges(top_atac$peak_id)
pairs_atac      <- make_cre_pairs(cres_atac)

## Run APA for both Hi-C conditions, both peak sets
apa_noatac_ctrl <- run_apa(pairs_noatac, hic_ctrl) / length(pairs_noatac)
apa_noatac_sorb <- run_apa(pairs_noatac, hic_sorb) / length(pairs_noatac)

apa_atac_ctrl   <- run_apa(pairs_atac, hic_ctrl) / length(pairs_atac)
apa_atac_sorb   <- run_apa(pairs_atac, hic_sorb) / length(pairs_atac)


# Visualization -----------------------------------------------------------

pal_hic <- colorRampPalette(brewer.pal(9, "YlGnBu"))

## Column x positions
x_ctrl <- left_margin
x_sorb <- x_ctrl + apa_width + col_gap

## Helper: plot one APA row (control | sorbitol) with a centered row title
## and n label above the matrices. Left-side labels are intentionally omitted.
plot_apa_row <- function(apa_ctrl_mat, apa_sorb_mat,
                         y_pos, zrange, row_title, n) {
  
  ## Row title — centered over both matrices, sitting above them
  x_center <- x_ctrl + apa_width + col_gap / 2   # midpoint between the two matrices
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
  
  ## n label — centered below the row title
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

## Page layout
panel_height <- top_margin + apa_height + 0.2
page_width   <- left_margin + 2 * apa_width + col_gap + 0.35
page_height  <- 2 * panel_height + 0.25


# Save outputs ------------------------------------------------------------

pdf(output_path, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = FALSE)

## ---- Panel 1: No ATAC filter ----

zrange_noatac <- c(0, max(apa_noatac_ctrl, apa_noatac_sorb, na.rm = TRUE))

add_column_headers(y_header = top_margin - 0.05)

plot_apa_row(
  apa_ctrl_mat = apa_noatac_ctrl,
  apa_sorb_mat = apa_noatac_sorb,
  y_pos        = top_margin,
  zrange       = zrange_noatac,
  row_title    = "Top K27ac peaks",
  n            = length(pairs_noatac)
)

## ---- Panel 2: With ATAC filter ----

y_panel2_offset <- panel_height + 0.1

zrange_atac <- c(0, max(apa_atac_ctrl, apa_atac_sorb, na.rm = TRUE))

add_column_headers(y_header = y_panel2_offset + top_margin - 0.05)

plot_apa_row(
  apa_ctrl_mat = apa_atac_ctrl,
  apa_sorb_mat = apa_atac_sorb,
  y_pos        = y_panel2_offset + top_margin,
  zrange       = zrange_atac,
  row_title    = "Top K27ac/ATAC peaks",
  n            = length(pairs_atac)
)

dev.off()


# Session info ------------------------------------------------------------

sessionInfo()