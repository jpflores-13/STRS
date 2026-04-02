# ##############################################################################
# filename:    atac_clustering_apa.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: APA analysis asking whether open chromatin regions (ATAC peaks)
#              co-localize in 3D space under sorbitol stress; uses top 1000
#              gained and top 1000 lost differential ATAC peaks as anchors
# ##############################################################################

# Libraries ----
library(DESeq2)
library(GenomicRanges)
library(InteractionSet)
library(mariner)
library(plotgardener)
library(RColorBrewer)
library(tidyverse)

source("scripts/utils/calculate_apa_score.R")

# Parameters ----
atac_counts_path <- "data/processed/atac/output/peaks/YAPP_HEK_1_peakCounts.tsv"
hic_ctrl         <- "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic"
hic_sorb         <- "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic"
output_pdf       <- "plots/atac_clustering_apa.pdf"

# Data import ----
atac_raw <- read_tsv(atac_counts_path, show_col_types = FALSE) |>
  dplyr::select(-starts_with("Unnamed"))

cat("Total ATAC peaks loaded:", nrow(atac_raw), "\n")

# Analysis ----

## Run DESeq2 on ATAC peaks ----
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

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = col_data,
  design    = ~ replicate + treatment
)

# Fix size factors to 1: hyperosmotic stress alters global chromatin
# accessibility, so we avoid normalization assuming count stability
sizeFactors(dds) <- rep(1, ncol(dds))

keep <- rowMeans(counts(dds, normalized = TRUE)) >= 15
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

## Select top 1000 gained and top 1000 lost ATAC peaks ----
gained <- res_df |>
  dplyr::filter(padj < 0.05, log2FoldChange > 1) |>
  dplyr::arrange(desc(log2FoldChange)) |>
  dplyr::slice_head(n = 1000)

lost <- res_df |>
  dplyr::filter(padj < 0.05, log2FoldChange < -1) |>
  dplyr::arrange(log2FoldChange) |>
  dplyr::slice_head(n = 1000)

cat("Gained ATAC peaks selected:", nrow(gained), "\n")
cat("Lost ATAC peaks selected:  ", nrow(lost), "\n")

## Convert peak IDs to GRanges ----
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

cres_gained <- peak_id_to_granges(gained$peak_id)
cres_lost   <- peak_id_to_granges(lost$peak_id)

## Generate ATAC-ATAC pairs ----
make_atac_pairs <- function(cres,
                            max_dist = 1e6,
                            min_dist = 25e3,
                            bin_size = 10e3) {

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

  cat("ATAC-ATAC pairs after distance filtering:", length(pairs), "\n")
  pairs
}

pairs_gained <- make_atac_pairs(cres_gained)
pairs_lost   <- make_atac_pairs(cres_lost)

## Pull Hi-C matrices and aggregate (APA) ----
run_apa <- function(pairs, hic_file, bin_size = 10e3, buffer = 10) {
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

apa_gained_ctrl <- run_apa(pairs_gained, hic_ctrl)
apa_gained_sorb <- run_apa(pairs_gained, hic_sorb)
apa_lost_ctrl   <- run_apa(pairs_lost,   hic_ctrl)
apa_lost_sorb   <- run_apa(pairs_lost,   hic_sorb)

apa_gained_ctrl <- apa_gained_ctrl / length(pairs_gained)
apa_gained_sorb <- apa_gained_sorb / length(pairs_gained)
apa_lost_ctrl   <- apa_lost_ctrl   / length(pairs_lost)
apa_lost_sorb   <- apa_lost_sorb   / length(pairs_lost)

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

panel_height <- top_margin + 2 * apa_height + row_gap + 0.5
page_width   <- left_margin + 2 * apa_width + col_gap + 0.35
page_height  <- panel_height + 0.25

zrange <- c(0, max(apa_gained_ctrl, apa_gained_sorb,
                   apa_lost_ctrl,   apa_lost_sorb, na.rm = TRUE))

# Save outputs ----
pdf(output_pdf, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = FALSE)

plotText(label = "ATAC peak clustering APA", fontsize = 8, fontface = "bold",
         fontcolor = "black", x = left_margin, y = 0.15,
         just = c("left", "top"), default.units = "inches")

add_column_headers(y_header = top_margin - 0.05)

plot_apa_row(apa_ctrl_mat = apa_gained_ctrl, apa_sorb_mat = apa_gained_sorb,
             y_pos = y_gained, zrange = zrange, label = "gained\nATAC",
             n = length(pairs_gained))

plot_apa_row(apa_ctrl_mat = apa_lost_ctrl, apa_sorb_mat = apa_lost_sorb,
             y_pos = y_lost, zrange = zrange, label = "lost\nATAC",
             n = length(pairs_lost))

dev.off()

sessionInfo()
