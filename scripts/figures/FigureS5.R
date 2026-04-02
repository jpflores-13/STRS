# ##############################################################################
# filename:    FigureS5.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Supplementary Figure S5 — absence of Downstream-of-Gene (DoG)
#              transcription; RNA-seq signal tracks for SAMD4A and KDM6A with
#              gene body and downstream region highlights
# ##############################################################################

# Libraries ----
library(plotgardener)
library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Parameters ----
bed_file               <- "data/processed/rna/Amat2019/Supplemental_Table_S2.bed"
rna_1h                 <- "data/processed/rna/timecourse/output/mergeSignal/unstranded/STRS_HEK293_WT_sorb_1h.bw"
rna_0h                 <- "data/processed/rna/timecourse/output/mergeSignal/unstranded/STRS_HEK293_WT_cont_0h.bw"
output_pdf             <- "figures/FigureS5.pdf"
extension_downstream_kb <- 20
extension_upstream_kb   <- 5
target_genes           <- c("SAMD4A", "KDM6A")
rna_color_1h           <- "#41B6C4"
rna_color_0h           <- "#7FCDBB"
page_width             <- 8.5
page_height            <- 11
margin_left            <- 0.8
margin_right           <- 0.5

# Data import ----
bed_data <- read.table(bed_file, header = FALSE, sep = "\t",
                       col.names = c("chr", "start", "end", "gene_id",
                                     "score", "strand", "symbol", "biotype")) |>
  dplyr::mutate(gene_length = end - start)

selected_genes <- bed_data |>
  dplyr::filter(symbol %in% target_genes)

stopifnot(nrow(selected_genes) > 0)

# Analysis ----
title_h          <- 0.15
rna_track_h      <- 0.55
gap_between_0h_1h <- 0.12
gap_before_gene  <- 0.12
gene_track_h     <- 0.55
gap_after_gene   <- 0.08
genome_label_h   <- 0.25
gap_after_genome <- 0.50

one_gene_height <- title_h + (rna_track_h * 2) + gap_between_0h_1h +
  gap_before_gene + gene_track_h + gap_after_gene +
  genome_label_h + gap_after_genome

rna_files_all   <- c(rna_0h, rna_1h)
rna_files_exist <- rna_files_all[file.exists(rna_files_all)]

if (length(rna_files_exist) > 0) {
  gene_max_signals       <- numeric(nrow(selected_genes))
  downstream_max_signals <- numeric(nrow(selected_genes))

  for (i in seq_len(nrow(selected_genes))) {
    gene_row <- selected_genes[i, ]

    gene_body_range <- calcSignalRange(
      data        = rna_files_exist,
      chrom       = gene_row$chr,
      chromstart  = gene_row$start,
      chromend    = gene_row$end,
      assembly    = "hg38",
      negData     = FALSE
    )
    gene_max_signals[i] <- gene_body_range[2]

    if (gene_row$strand == "+") {
      downstream_start <- gene_row$end
      downstream_end   <- gene_row$end + (extension_downstream_kb * 1000)
    } else {
      downstream_start <- max(0, gene_row$start - (extension_downstream_kb * 1000))
      downstream_end   <- gene_row$start
    }

    downstream_range <- calcSignalRange(
      data        = rna_files_exist,
      chrom       = gene_row$chr,
      chromstart  = downstream_start,
      chromend    = downstream_end,
      assembly    = "hg38",
      negData     = FALSE
    )
    downstream_max_signals[i] <- downstream_range[2]
  }

  downstream_max <- max(downstream_max_signals, na.rm = TRUE)
  signalRange    <- c(0, downstream_max * 3)
} else {
  signalRange <- NULL
  message("Warning: no RNA-seq files found; signal range will be calculated per region")
}

create_dog_panel <- function(gene_row, panel_x, panel_y, gene_number) {
  chrom          <- gene_row$chr
  gene_start     <- gene_row$start
  gene_end       <- gene_row$end
  strand         <- gene_row$strand
  symbol         <- gene_row$symbol

  if (strand == "+") {
    region_start    <- max(0, gene_start - (extension_upstream_kb * 1000))
    region_end      <- gene_end + (extension_downstream_kb * 1000)
    downstream_start <- gene_end
    downstream_end   <- region_end
  } else {
    region_start    <- max(0, gene_start - (extension_downstream_kb * 1000))
    region_end      <- gene_end + (extension_upstream_kb * 1000)
    downstream_start <- region_start
    downstream_end   <- gene_start
  }

  plotText(label = symbol, x = panel_x, y = panel_y,
           fontsize = 10, fontface = "bold", just = c("left", "top"))

  y           <- panel_y + title_h
  panel_width <- page_width - margin_left - margin_right
  params      <- pgParams(chrom = chrom, chromstart = region_start, chromend = region_end,
                          assembly = "hg38", x = panel_x,
                          width = panel_width, length = panel_width)

  if (file.exists(rna_0h)) {
    plotSignal(data = rna_0h, params = params, y = y, height = rna_track_h,
               fill = rna_color_0h, linecolor = rna_color_0h,
               baseline = TRUE, baseline.color = "grey50", baseline.lwd = 0.5,
               scale = FALSE, range = signalRange)
    plotText("0h", x = panel_x - 0.15, y = y + rna_track_h / 2,
             just = c("right", "center"), fontsize = 7)
    y <- y + rna_track_h + gap_between_0h_1h
  }

  if (file.exists(rna_1h)) {
    plotSignal(data = rna_1h, params = params, y = y, height = rna_track_h,
               fill = rna_color_1h, linecolor = rna_color_1h,
               baseline = TRUE, baseline.color = "grey50", baseline.lwd = 0.5,
               scale = FALSE, range = signalRange)
    plotText("1h", x = panel_x - 0.15, y = y + rna_track_h / 2,
             just = c("right", "center"), fontsize = 7)
    y <- y + rna_track_h + gap_before_gene
  }

  genes_plot <- plotGenes(params = params, y = y, height = gene_track_h,
                          fontsize = 6, strandLabels = TRUE)
  y <- y + gene_track_h + gap_after_gene

  plotGenomeLabel(params = params, y = y, scale = "bp", fontsize = 7)
  y <- y + genome_label_h

  highlight_top    <- panel_y + title_h
  highlight_height <- y - highlight_top

  annoHighlight(plot = genes_plot, chrom = chrom,
                chromstart = gene_start, chromend = gene_end,
                fill = "lightblue", alpha = 0.2, linecolor = "blue",
                lwd = 0.8, lty = 2,
                x = panel_x, y = highlight_top,
                width = panel_width, height = highlight_height,
                just = c("left", "top"))

  annoHighlight(plot = genes_plot, chrom = chrom,
                chromstart = downstream_start, chromend = downstream_end,
                fill = "pink", alpha = 0.15, linecolor = "red",
                lwd = 0.8, lty = 3,
                x = panel_x, y = highlight_top,
                width = panel_width, height = highlight_height,
                just = c("left", "top"))

  return(y)
}

# Visualization ----
n_panels <- nrow(selected_genes)

pdf(output_pdf, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = FALSE)

legend_y        <- 0.3
legend_center_x <- page_width / 2

plotText(label = "Legend:", x = legend_center_x - 1.8, y = legend_y,
         just = c("left", "top"), fontsize = 8, fontface = "bold")
plotRect(x = legend_center_x - 1.2, y = legend_y + 0.02,
         width = 0.25, height = 0.10,
         fill = "lightblue", alpha = 0.2, linecolor = "blue", lty = 2,
         just = c("left", "top"))
plotText(label = "Gene body", x = legend_center_x - 0.9, y = legend_y + 0.07,
         just = c("left", "center"), fontsize = 7)
plotRect(x = legend_center_x + 0.2, y = legend_y + 0.02,
         width = 0.25, height = 0.10,
         fill = "pink", alpha = 0.15, linecolor = "red", lty = 3,
         just = c("left", "top"))
plotText(label = sprintf("Downstream (%d kb)", extension_downstream_kb),
         x = legend_center_x + 0.5, y = legend_y + 0.07,
         just = c("left", "center"), fontsize = 7)

y_position <- 0.6

for (i in seq_len(n_panels)) {
  gene_row <- selected_genes[i, ]
  tryCatch({
    final_y    <- create_dog_panel(gene_row, margin_left, y_position, i)
    y_position <- final_y + gap_after_genome
  }, error = function(e) {
    y_position <<- y_position + 3.0
  })
}

# Save outputs ----
dev.off()

sessionInfo()
