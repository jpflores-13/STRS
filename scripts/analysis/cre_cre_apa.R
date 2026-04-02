# ##############################################################################
# filename:    cre_cre_apa.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: APA analysis between all pairs of H3K27ac CRE peaks (with and
#              without ATAC intersection) under control and sorbitol conditions
# ##############################################################################

# Libraries ----
library(DESeq2)
library(GenomicRanges)
library(InteractionSet)
library(mariner)
library(plotgardener)
library(RColorBrewer)
library(plyranges)
library(tidyverse)

source("scripts/utils/calculate_apa_score.R")

# Parameters ----
diff_h3k27ac_path <- "data/processed/cutntag/deseq2/diff_H3K27ac_counts.rds"
atac_counts_path  <- "data/processed/atac/output/peaks/YAPP_HEK_1_peakCounts.tsv"
hic_ctrl          <- "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic"
hic_sorb          <- "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic"
output_pdf        <- "plots/cre_cre_apa.pdf"

# Data import ----
diff_h3k27ac <- readRDS(diff_h3k27ac_path)

res_df <- as.data.frame(mcols(diff_h3k27ac)) |>
  mutate(peak_id = paste0(
    as.character(seqnames(diff_h3k27ac)), ":",
    start(diff_h3k27ac), "-",
    end(diff_h3k27ac)
  ))

cat("Total H3K27ac peaks loaded:", nrow(res_df), "\n")

atac_raw <- read_tsv(atac_counts_path, show_col_types = FALSE) |>
  dplyr::select(-starts_with("Unnamed"))

atac_gr <- GRanges(
  seqnames = atac_raw$chr,
  ranges   = IRanges(start = atac_raw$start, end = atac_raw$stop)
)

# Analysis ----

## Select top gained and lost CREs ----
select_top_cres <- function(res_df) {

  gained <- res_df |>
    dplyr::filter(padj < 0.05, log2FoldChange > 1) |>
    dplyr::arrange(desc(log2FoldChange))

  lost <- res_df |>
    dplyr::filter(padj < 0.05, log2FoldChange < -1) |>
    dplyr::arrange(log2FoldChange)

  cat("Gained CREs selected:", nrow(gained), "\n")
  cat("Lost CREs selected:  ", nrow(lost), "\n")

  list(gained = gained, lost = lost)
}

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

## No ATAC filter
top_noatac <- select_top_cres(res_df)
cres_noatac_gained <- peak_id_to_granges(top_noatac$gained$peak_id)
cres_noatac_lost   <- peak_id_to_granges(top_noatac$lost$peak_id)

## With ATAC filter
k27ac_gr    <- granges(diff_h3k27ac)
atac_hits   <- findOverlaps(k27ac_gr, atac_gr)
res_df_atac <- res_df[unique(queryHits(atac_hits)), ]
cat("H3K27ac peaks overlapping ATAC:", nrow(res_df_atac), "\n")

top_atac <- select_top_cres(res_df_atac)
cres_atac_gained <- peak_id_to_granges(top_atac$gained$peak_id)
cres_atac_lost   <- peak_id_to_granges(top_atac$lost$peak_id)

## Generate CRE-CRE pairs ----
make_cre_pairs <- function(cres, max_dist = 1e6, min_dist = 25e3, bin_size = 10e3) {
  cres  <- keepStandardChromosomes(cres, pruning.mode = "coarse")
  hits  <- findOverlaps(cres, cres, maxgap = max_dist, ignore.strand = TRUE)
  hits  <- hits[queryHits(hits) < subjectHits(hits)]

  pairs <- GInteractions(
    anchor1 = cres[queryHits(hits)],
    anchor2 = cres[subjectHits(hits)]
  )

  pairs <- pairs[pairdist(pairs) >= min_dist]

  flipped <- start(anchors(pairs, "first")) > start(anchors(pairs, "second"))
  pairs[flipped] <- swapAnchors(pairs[flipped])

  pairs <- assignToBins(pairs, binSize = bin_size)
  cat("CRE-CRE pairs after distance filtering:", length(pairs), "\n")
  pairs
}

pairs_noatac_gained <- make_cre_pairs(cres_noatac_gained)
pairs_noatac_lost   <- make_cre_pairs(cres_noatac_lost)
pairs_atac_gained   <- make_cre_pairs(cres_atac_gained)
pairs_atac_lost     <- make_cre_pairs(cres_atac_lost)

## Pull Hi-C matrices and aggregate ----
run_apa <- function(pairs, hic_file, bin_size = 10e3, buffer = 10) {
  pairs |>
    pixelsToMatrices(buffer = buffer) |>
    removeShortPairs() |>
    pullHicMatrices(binSize = bin_size, files = hic_file,
                    half = "upper", norm = "NONE", matrix = "observed") |>
    aggHicMatrices(FUN = sum)
}

apa_noatac_gained_ctrl <- run_apa(pairs_noatac_gained, hic_ctrl)
apa_noatac_gained_sorb <- run_apa(pairs_noatac_gained, hic_sorb)
apa_noatac_lost_ctrl   <- run_apa(pairs_noatac_lost,   hic_ctrl)
apa_noatac_lost_sorb   <- run_apa(pairs_noatac_lost,   hic_sorb)

apa_atac_gained_ctrl   <- run_apa(pairs_atac_gained, hic_ctrl)
apa_atac_gained_sorb   <- run_apa(pairs_atac_gained, hic_sorb)
apa_atac_lost_ctrl     <- run_apa(pairs_atac_lost,   hic_ctrl)
apa_atac_lost_sorb     <- run_apa(pairs_atac_lost,   hic_sorb)

apa_noatac_gained_ctrl <- apa_noatac_gained_ctrl / length(pairs_noatac_gained)
apa_noatac_gained_sorb <- apa_noatac_gained_sorb / length(pairs_noatac_gained)
apa_noatac_lost_ctrl   <- apa_noatac_lost_ctrl   / length(pairs_noatac_lost)
apa_noatac_lost_sorb   <- apa_noatac_lost_sorb   / length(pairs_noatac_lost)

apa_atac_gained_ctrl   <- apa_atac_gained_ctrl / length(pairs_atac_gained)
apa_atac_gained_sorb   <- apa_atac_gained_sorb / length(pairs_atac_gained)
apa_atac_lost_ctrl     <- apa_atac_lost_ctrl   / length(pairs_atac_lost)
apa_atac_lost_sorb     <- apa_atac_lost_sorb   / length(pairs_atac_lost)

# Visualization ----

apa_width   <- 0.75
apa_height  <- 0.75
col_gap     <- 0.1
row_gap     <- 0.5
left_margin <- 0.5
top_margin  <- 0.75

x_ctrl  <- left_margin
x_sorb  <- x_ctrl + apa_width + col_gap
y_gained <- top_margin
y_lost   <- y_gained + apa_height + row_gap

pal_hic <- colorRampPalette(brewer.pal(9, "YlGnBu"))

plot_apa_row <- function(apa_ctrl_mat, apa_sorb_mat, y_pos, zrange, label, n) {
  apaParams <- pgParams(width = apa_width, height = apa_height,
                        fontsize = 5, palette = pal_hic)

  plotMatrix(data = apa_ctrl_mat, x = x_ctrl, y = y_pos,
             zrange = zrange, params = apaParams)

  plotText(label = sprintf("%.2f", calculate_apa_score(apa_ctrl_mat)),
           fontsize = 5, fontcolor = "black",
           x = x_ctrl + apa_width - 0.05, y = y_pos + 0.03,
           just = c("right", "top"))

  plotMatrix(data = apa_sorb_mat, x = x_sorb, y = y_pos,
             zrange = zrange, params = apaParams) |>
    annoHeatmapLegend(x = x_sorb + apa_width + 0.05, y = y_pos,
                      width = 0.05, height = apa_height,
                      fontsize = 5, fontcolor = "black",
                      just = c("left", "top"), default.units = "inches")

  plotText(label = sprintf("%.2f", calculate_apa_score(apa_sorb_mat)),
           fontsize = 5, fontcolor = "black",
           x = x_sorb + apa_width - 0.05, y = y_pos + 0.03,
           just = c("right", "top"))

  plotText(label = label, fontsize = 6, fontcolor = "black",
           x = x_ctrl - 0.05, y = y_pos + apa_height / 2,
           just = c("right", "center"), default.units = "inches")

  plotText(label = paste0("n = ", n), fontsize = 5, fontcolor = "grey40",
           x = x_ctrl - 0.05, y = y_pos + apa_height / 2 + 0.15,
           just = c("right", "center"), default.units = "inches")
}

add_column_headers <- function(y_header) {
  for (lbl in list(
    list(label = "untreated",  x = x_ctrl + apa_width / 2),
    list(label = "+ sorbitol", x = x_sorb + apa_width / 2)
  )) {
    plotText(label = lbl$label, fontsize = 6, fontcolor = "black",
             x = lbl$x, y = y_header, just = c("center", "bottom"),
             default.units = "inches")
  }
}

panel_height  <- top_margin + 2 * apa_height + row_gap + 0.5
page_width    <- left_margin + 2 * apa_width + col_gap + 0.35
page_height   <- 2 * panel_height + 0.5

# Save outputs ----
pdf(output_pdf, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = FALSE)

## Panel 1: No ATAC filter ----
zrange_noatac <- c(0, max(apa_noatac_gained_ctrl, apa_noatac_gained_sorb,
                           apa_noatac_lost_ctrl,   apa_noatac_lost_sorb, na.rm = TRUE))

plotText(label = "H3K27ac CRE-CRE APA (no ATAC filter)", fontsize = 8, fontface = "bold",
         fontcolor = "black", x = left_margin, y = 0.15,
         just = c("left", "top"), default.units = "inches")
add_column_headers(y_header = top_margin - 0.05)
plot_apa_row(apa_noatac_gained_ctrl, apa_noatac_gained_sorb, y_gained, zrange_noatac,
             "gained\nH3K27ac", length(pairs_noatac_gained))
plot_apa_row(apa_noatac_lost_ctrl, apa_noatac_lost_sorb, y_lost, zrange_noatac,
             "lost\nH3K27ac", length(pairs_noatac_lost))

## Panel 2: With ATAC filter ----
y_panel2_offset <- panel_height + 0.25
zrange_atac <- c(0, max(apa_atac_gained_ctrl, apa_atac_gained_sorb,
                         apa_atac_lost_ctrl,   apa_atac_lost_sorb, na.rm = TRUE))

plotText(label = "H3K27ac CRE-CRE APA (intersected with ATAC)", fontsize = 8, fontface = "bold",
         fontcolor = "black", x = left_margin, y = y_panel2_offset + 0.15,
         just = c("left", "top"), default.units = "inches")
add_column_headers(y_header = y_panel2_offset + top_margin - 0.05)
plot_apa_row(apa_atac_gained_ctrl, apa_atac_gained_sorb, y_panel2_offset + y_gained, zrange_atac,
             "gained\nH3K27ac", length(pairs_atac_gained))
plot_apa_row(apa_atac_lost_ctrl, apa_atac_lost_sorb, y_panel2_offset + y_lost, zrange_atac,
             "lost\nH3K27ac", length(pairs_atac_lost))

dev.off()

sessionInfo()
