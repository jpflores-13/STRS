# FigureS10 ----------------------------------------------------------------
# Author:      JP Flores
# Date:        2026-04-23
# Project:     STRS
# Description: Composite supplemental figure with two panels. Panel A: stacked
#              barplot showing the fraction of gained, lost, and static loops
#              that are enhancer-promoter (E-P) loops — legend replaced with
#              direct text labels at the top of each bar. Panel B: APA
#              analysis of Hi-C contact frequency between pairs of confidently
#              detected H3K27ac CRE peaks overlapping ATAC peaks; control and
#              sorbitol matrices are stacked vertically (control top, sorbitol
#              bottom).
# Input:       - H3K27ac DESeq2 results: data/processed/cutntag/deseq2/diff_H3K27ac_counts.rds
#              - ATAC peak counts:       data/processed/atac/output/peaks/YAPP_HEK_1_peakCounts.tsv
#              - Hi-C maps:              data/processed/hic/maps/
#              - Differential loops:    data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds
#              - CUT&Tag peaks dir:     data/processed/cutntag/output/peaks/
# Output:      - figures/FigureS10.pdf
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## E-P loop enrichment
diff_loops_rds    <- file.path("data/processed/hic/diffLoops",
                               "diffLoops_eGFP-YAP_noDroso_10kb.rds")
cutntag_peaks_dir <- "data/processed/cutntag/output/peaks/"

## CRE-CRE APA
min_count <- 50           ## Min baseMean to retain a peak as confidently detected
top_n     <- 1000         ## Max CREs to use after threshold filter (ranked by baseMean)
max_dist  <- 1e6          ## Maximum genomic distance between CRE pairs (bp)
min_dist  <- 25e3         ## Minimum genomic distance between CRE pairs (bp)
bin_size  <- 10e3         ## Hi-C bin size (bp)
buffer    <- 10           ## APA matrix buffer (bins on each side)

## Page / panel geometry
bar_width   <- 2.25       ## Width of panel A barplot (inches)
bar_height  <- 2.5        ## Height of panel A barplot (inches)
apa_width   <- 0.9        ## APA matrix width (inches)
apa_height  <- 0.9        ## APA matrix height (inches)
row_gap     <- 0.15       ## Vertical gap between stacked APA matrices (inches)
panel_gap   <- 0.5        ## Horizontal gap between panel A and panel B (inches)
left_margin <- 0.4        ## Left page margin (inches)
top_margin  <- 0.5        ## Top page margin (inches)
apa_top     <- 0.55       ## Top offset within panel B (for header + first matrix)

## x-axis label colors matching Figure 4 convention
col_gained  <- "#F8766D"
col_lost    <- "#619CFF"
col_static  <- "#999999"

## Paths
project_dir       <- "/work/users/j/p/jpflores/projects/STRS"
diff_h3k27ac_path <- file.path(project_dir,
                               "data/processed/cutntag/deseq2/diff_H3K27ac_counts.rds")
atac_counts_path  <- file.path(project_dir,
                               "data/processed/atac/output/peaks/YAPP_HEK_1_peakCounts.tsv")
hic_ctrl          <- file.path(project_dir,
                               "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic")
hic_sorb          <- file.path(project_dir,
                               "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic")
diff_loops_path   <- file.path(project_dir, diff_loops_rds)
peaks_dir         <- file.path(project_dir, cutntag_peaks_dir)
output_path       <- file.path(project_dir, "figures/FigureS10.pdf")


# Libraries ---------------------------------------------------------------

library(BiocGenerics)
library(DESeq2)
library(GenomeInfoDb)
library(GenomicRanges)
library(ggplot2)
library(ggtext)
library(InteractionSet)
library(mariner)
library(plotgardener)
library(plyranges)
library(RColorBrewer)
library(scales)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


# Utility scripts ---------------------------------------------------------

source(file.path(project_dir, "scripts/utils/calculate_apa_score.R"))


# Load data ---------------------------------------------------------------

## Pre-computed H3K27ac DESeq2 results
diff_h3k27ac <- readRDS(diff_h3k27ac_path)

## ATAC peak count table
atac_raw <- read_tsv(atac_counts_path, show_col_types = FALSE) |>
  dplyr::select(-starts_with("Unnamed"))

## Differential loop calls
noDroso_loops <- readRDS(diff_loops_path) |>
  interactions() |>
  as.data.frame() |>
  as_ginteractions()

## TxDb for promoter definitions
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene |> keepStandardChromosomes()


# Wrangle data ------------------------------------------------------------

## ---- Panel A: E-P loop enrichment ----

mcols(noDroso_loops)$loop_type <- case_when(
  mcols(noDroso_loops)$padj <= 0.1 & mcols(noDroso_loops)$log2FoldChange > 0 ~ "gained",
  mcols(noDroso_loops)$padj <= 0.1 & mcols(noDroso_loops)$log2FoldChange < 0 ~ "lost",
  mcols(noDroso_loops)$padj >  0.1 ~ "static",
  TRUE ~ "other"
)

seqlevels(noDroso_loops) <- seqlevels(txdb)
seqinfo(noDroso_loops)   <- seqinfo(txdb)

gainedLoops <- noDroso_loops |> subset(loop_type == "gained")
lostLoops   <- noDroso_loops |> subset(loop_type == "lost")
staticLoops <- noDroso_loops |> subset(loop_type == "static")

genes            <- genes(txdb)
promoter_regions <- promoters(genes)

k27ac_peaks <- list.files(peaks_dir, full.names = TRUE, pattern = ".narrowPeak") |>
  stringr::str_subset("H3K27ac") |>
  lapply(plyranges::read_narrowpeaks) |>
  lapply(keepStandardChromosomes, pruning.mode = "coarse")

k27ac_consensus   <- GenomicRanges::reduce(unlist(GRangesList(k27ac_peaks)))
overlapping_peaks <- subsetByOverlaps(k27ac_consensus, promoter_regions)
enhancers         <- GenomicRanges::setdiff(k27ac_consensus, overlapping_peaks)

## ---- Panel B: CRE-CRE APA ----

res_df <- as.data.frame(mcols(diff_h3k27ac)) |>
  mutate(peak_id = paste0(
    as.character(seqnames(diff_h3k27ac)), ":",
    start(diff_h3k27ac), "-",
    end(diff_h3k27ac)
  ))

cat("Total H3K27ac peaks loaded:", nrow(res_df), "\n")

atac_gr     <- GRanges(
  seqnames = atac_raw$chr,
  ranges   = IRanges(start = atac_raw$start, end = atac_raw$stop)
)
k27ac_gr    <- granges(diff_h3k27ac)
atac_hits   <- findOverlaps(k27ac_gr, atac_gr)
res_df_atac <- res_df[unique(queryHits(atac_hits)), ]

cat("H3K27ac peaks overlapping ATAC:", nrow(res_df_atac), "\n")


# Analysis ----------------------------------------------------------------

## ---- Panel A helpers ----

analyze_ep_loops <- function(loops, promoters, enhancers) {
  ep_loops <- loops |>
    linkOverlaps(promoters, enhancers) |>
    as.data.frame() |>
    dplyr::select(query) |>
    dplyr::distinct()
  
  total_loops <- length(loops)
  ep_count    <- nrow(ep_loops)
  
  list(total       = total_loops,
       ep_count    = ep_count,
       ep_fraction = ep_count / total_loops)
}

ep_results <- list(
  gained = analyze_ep_loops(gainedLoops, promoter_regions, enhancers),
  lost   = analyze_ep_loops(lostLoops,   promoter_regions, enhancers),
  static = analyze_ep_loops(staticLoops, promoter_regions, enhancers)
)

ep_viz_data <- data.frame(
  loop_type       = c("Gained", "Lost", "Static"),
  ep_fraction     = c(ep_results$gained$ep_fraction,
                      ep_results$lost$ep_fraction,
                      ep_results$static$ep_fraction),
  non_ep_fraction = c(1 - ep_results$gained$ep_fraction,
                      1 - ep_results$lost$ep_fraction,
                      1 - ep_results$static$ep_fraction)
) |>
  pivot_longer(cols      = c(ep_fraction, non_ep_fraction),
               names_to  = "category",
               values_to = "fraction") |>
  group_by(loop_type) |>
  mutate(pos = cumsum(fraction) - 0.5 * fraction) |>
  ungroup()

ep_viz_data$category <- factor(ep_viz_data$category,
                               levels = c("non_ep_fraction", "ep_fraction"))

## Only label the E-P fraction segment
ep_label_df <- ep_viz_data |>
  dplyr::filter(category == "ep_fraction")

## ---- Panel B helpers ----

select_cres_by_threshold <- function(df) {
  selected <- df |>
    dplyr::filter(baseMean >= min_count) |>
    dplyr::arrange(desc(baseMean)) |>
    dplyr::slice_head(n = top_n)
  cat("CREs retained (baseMean >=", min_count, ", top", top_n, "):",
      nrow(selected), "out of", nrow(df), "total\n")
  selected
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

make_cre_pairs <- function(cres) {
  cres  <- keepStandardChromosomes(cres, pruning.mode = "coarse")
  hits  <- findOverlaps(cres, cres, maxgap = max_dist, ignore.strand = TRUE)
  hits  <- hits[queryHits(hits) < subjectHits(hits)]
  
  pairs <- GInteractions(
    anchor1 = cres[queryHits(hits)],
    anchor2 = cres[subjectHits(hits)]
  )
  
  pairs <- pairs[pairdist(pairs) >= min_dist]
  
  flipped        <- start(anchors(pairs, "first")) > start(anchors(pairs, "second"))
  pairs[flipped] <- swapAnchors(pairs[flipped])
  pairs          <- assignToBins(pairs, binSize = bin_size)
  
  cat("CRE-CRE pairs after distance filtering:", length(pairs), "\n")
  pairs
}

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

top_atac      <- select_cres_by_threshold(res_df_atac)
cres_atac     <- peak_id_to_granges(top_atac$peak_id)
pairs_atac    <- make_cre_pairs(cres_atac)

apa_atac_ctrl <- run_apa(pairs_atac, hic_ctrl) / length(pairs_atac)
apa_atac_sorb <- run_apa(pairs_atac, hic_sorb) / length(pairs_atac)


# Visualization -----------------------------------------------------------

## ---- Panel A: barplot ----

p_bar <- ggplot(ep_viz_data, aes(x = loop_type, y = fraction, fill = category)) +
  geom_bar(position = "stack", stat = "identity", width = 0.85) +
  geom_text(
    data     = ep_label_df,
    aes(y = pos, label = scales::percent(fraction, accuracy = 0.1) |>
          stringr::str_remove("%")),
    color    = "black",
    size     = 2.8,
    fontface = "plain"
  ) +
  scale_fill_manual(values = c("non_ep_fraction" = "lightgrey",
                               "ep_fraction"     = "#B8A9D9")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels = c("Gained" = "Gained",
                              "Lost"   = "Lost",
                              "Static" = "Static")) +
  labs(x = "", y = "% E-P Loops") +
  theme_minimal() +
  theme(
    axis.text.x        = element_text(size = 7),
    axis.text.y        = element_text(size = 8.5),
    axis.title.y       = element_text(size = 9.5),
    axis.line.x        = element_blank(),
    axis.line.y        = element_blank(),
    axis.ticks         = element_blank(),
    legend.position    = "none",
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin        = margin(2, 5, 2, 2, "pt")
  )

## ---- Panel B: plotgardener APA (stacked: ctrl top, sorb bottom) ----

pal_hic    <- colorRampPalette(brewer.pal(9, "YlGnBu"))
zrange_b   <- c(0, max(apa_atac_ctrl, apa_atac_sorb, na.rm = TRUE))

x_b_origin <- left_margin + bar_width + panel_gap
y_ctrl     <- top_margin + apa_top
y_sorb     <- y_ctrl + apa_height + row_gap
x_legend   <- x_b_origin + apa_width + 0.05

## Expand page width to give legend and labels plenty of room
page_width  <- x_b_origin + apa_width + 0.8
page_height <- top_margin + apa_top + 2 * apa_height + row_gap + 0.5

apaParams <- pgParams(
  width    = apa_width,
  height   = apa_height,
  fontsize = 5,
  palette  = pal_hic
)

add_apa_score <- function(mat, x_pos, y_pos) {
  plotText(
    label     = sprintf("%.2f", calculate_apa_score(mat)),
    fontsize  = 5,
    fontcolor = "black",
    x         = x_pos + apa_width - 0.05,
    y         = y_pos + 0.03,
    just      = c("right", "top")
  )
}


# Save outputs ------------------------------------------------------------

pdf(output_path, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = FALSE)

## ---- Panel A label ----
plotText(
  label         = "A",
  fontsize      = 10,
  fontface      = "bold",
  x             = left_margin,
  y             = top_margin,
  just          = c("left", "top"),
  default.units = "inches"
)

## ---- Panel A: barplot ----
plotGG(
  plot          = p_bar,
  x             = left_margin,
  y             = top_margin + 0.15,
  width         = bar_width,
  height        = bar_height,
  just          = c("left", "top"),
  default.units = "inches"
)

## ---- Panel B label ----
plotText(
  label         = "B",
  fontsize      = 10,
  fontface      = "bold",
  x             = x_b_origin,
  y             = top_margin,
  just          = c("left", "top"),
  default.units = "inches"
)

## ---- Panel B title ----
plotText(
  label         = "Top Overlapping H3K27ac & ATAC Peaks",
  fontsize      = 7,
  fontface      = "bold",
  fontcolor     = "black",
  x             = x_b_origin + apa_width / 2,
  y             = top_margin + 0.15,
  just          = c("center", "top"),
  default.units = "inches"
)

## ---- n label ----
plotText(
  label         = paste0("n = ", length(pairs_atac), " pairs"),
  fontsize      = 5,
  fontcolor     = "grey40",
  x             = x_b_origin + apa_width / 2,
  y             = top_margin + 0.35,
  just          = c("center", "top"),
  default.units = "inches"
)

## ---- Control matrix (top) ----
plotMatrix(
  data          = apa_atac_ctrl,
  x             = x_b_origin,
  y             = y_ctrl,
  zrange        = zrange_b,
  params        = apaParams
)
add_apa_score(apa_atac_ctrl, x_b_origin, y_ctrl)

plotText(
  label         = "untreated",
  fontsize      = 5,
  fontcolor     = "black",
  x             = x_b_origin + 0.03,
  y             = y_ctrl + 0.03,
  just          = c("left", "top"),
  default.units = "inches"
)

## ---- Sorbitol matrix (bottom) with shared legend ----
plotMatrix(
  data          = apa_atac_sorb,
  x             = x_b_origin,
  y             = y_sorb,
  zrange        = zrange_b,
  params        = apaParams
) |>
  annoHeatmapLegend(
    x             = x_legend,
    y             = y_ctrl,
    width         = 0.05,
    height        = apa_height,    ## single-matrix height, not full span
    fontsize      = 5,
    fontcolor     = "black",
    just          = c("left", "top"),
    default.units = "inches"
  )
add_apa_score(apa_atac_sorb, x_b_origin, y_sorb)

plotText(
  label         = "+ sorbitol",
  fontsize      = 5,
  fontcolor     = "black",
  x             = x_b_origin + 0.03,
  y             = y_sorb + 0.03,
  just          = c("left", "top"),
  default.units = "inches"
)

dev.off()


# Session info ------------------------------------------------------------

sessionInfo()