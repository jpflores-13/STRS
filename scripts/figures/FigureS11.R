# FigureS11 ----------------------------------------------------------------
# Author:      JP Flores
# Date:        2026-04-23
# Project:     STRS
# Description: Three-row supplementary figure. Row 1 (3 panels): HOMER motif
#              QQ scatter plots for CTCF, RAD21, and ATAC loop anchors.
#              Row 2 (3 panels): SP/KLF zoom-in of the same three HOMER
#              comparisons (CTCF/BORIS removed so axis rescales).
#              Row 3 (3 panels): ENCODE HEK293/HEK293T ChIP-seq overlap
#              scatter plots for the same three comparisons.
#              Zoom-in dashed rectangles on each Row 1 panel connect via
#              dotted lines to the corresponding Row 2 panel directly below,
#              highlighting the low-enrichment / SP/KLF region. Matching
#              dashed rectangles on each Row 2 panel connect via dotted lines
#              to the corresponding Row 3 panel below.
# Input:       - data/processed/cutntag/homer_motifs/<condition>/knownResults.txt
#              - data/processed/atac/homer_input/<condition>/knownResults.txt
#              - data/external/encode_chipseq/metadata_filtered.tsv
#              - data/external/encode_chipseq/beds/*.bed.gz
#              - data/processed/cutntag/output/peaks/*.narrowPeak
#              - data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds
# Output:      - figures/FigureS11.pdf
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

project_dir    <- "/work/users/j/p/jpflores/projects/STRS"
homer_base     <- file.path(project_dir, "data/processed/cutntag/homer_motifs")
homer_atac     <- file.path(project_dir, "data/processed/atac/homer_input")
encode_dir     <- file.path(project_dir, "data/external/encode_chipseq")
beds_dir       <- file.path(encode_dir, "beds")
metadata_file  <- file.path(encode_dir, "metadata_filtered.tsv")
peaks_dir      <- file.path(project_dir, "data/processed/cutntag/output/peaks")
diff_loops_rds <- file.path(project_dir, "data/processed/hic/diffLoops",
                            "diffLoops_eGFP-YAP_noDroso_10kb.rds")
output_pdf     <- file.path(project_dir, "figures/FigureS11.pdf")

## Peak files
ctcf_narrowpeak  <- file.path(peaks_dir, "STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak")
rad21_narrowpeak <- file.path(peaks_dir, "STRS_HEK293_eGFP-YAP_RAD21_cont_0h_peaks.narrowPeak")
ctcf_sorb_peak   <- file.path(peaks_dir, "STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h_peaks.narrowPeak")
rad21_sorb_peak  <- file.path(peaks_dir, "STRS_HEK293_eGFP-YAP_RAD21_sorbitol_1h_peaks.narrowPeak")

## Analysis thresholds
padj_threshold <- 0.1
peak_buffer    <- 1000
min_pct        <- 1
top_n_label    <- 15
min_log10p     <- 2

## Page layout (inches) — three rows of three equal panels
page_width  <- 10
panel_w     <- 2.8    ## width of each panel (same for all rows)
panel_h     <- 2.8    ## height of each panel (same for all rows)
row_gap     <- 0.6    ## vertical gap between rows (space for connectors)
col_gap     <- 0.35   ## horizontal gap between panels within a row

top_margin  <- 0.5
left_margin <- 0.5

## Derived x origins for each column
x_col <- left_margin + (0:2) * (panel_w + col_gap)

## Derived y origins for each row
y_row1 <- top_margin
y_row2 <- y_row1 + panel_h + row_gap
y_row3 <- y_row2 + panel_h + row_gap

page_height <- y_row3 + panel_h + 0.3

## Zoom box geometry — fraction of the data area covered by the zoom box.
## The box marks the dense / low-enrichment (origin) corner of each QQ plot.
zoom_frac <- 0.35


# Libraries ---------------------------------------------------------------

library(dplyr)
library(GenomeInfoDb)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
library(InteractionSet)
library(mariner)
library(plotgardener)
library(plyranges)
library(rtracklayer)
library(tidyr)


# Utility functions -------------------------------------------------------

read_narrowpeaks <- function(file) {
  peaks <- read.table(file, sep = "\t", header = FALSE,
                      col.names = c("seqnames", "start", "end", "name",
                                    "score", "strand", "signalValue",
                                    "pValue", "qValue", "peak"))
  GRanges(seqnames = peaks$seqnames,
          ranges   = IRanges(start = peaks$start + 1, end = peaks$end),
          strand   = "*",
          score = peaks$score, signalValue = peaks$signalValue,
          pValue = peaks$pValue, qValue = peaks$qValue, peak = peaks$peak)
}

parse_homer <- function(homer_dir) {
  results_file <- file.path(homer_dir, "knownResults.txt")
  if (!file.exists(results_file)) stop("knownResults.txt not found in: ", homer_dir)
  read.delim(results_file, header = TRUE, stringsAsFactors = FALSE) |>
    dplyr::mutate(
      neg_log10p  = -Log.P.value / log(10),
      pct_target  = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif)),
      pct_bg      = as.numeric(gsub("%", "", X..of.Background.Sequences.with.Motif)),
      motif_short = gsub("/.*", "", Motif.Name) |> gsub("\\(.*", "", x = _)
    ) |>
    dplyr::select(motif_short, neg_log10p, pct_target, pct_bg)
}

assign_tf_family <- function(motif_name) {
  dplyr::case_when(
    grepl("CTCF|BORIS|CTFL|CTCFL", motif_name, ignore.case = TRUE)    ~ "CTCF/BORIS",
    grepl("^KLF|^Klf|^SP[0-9]|^Sp[0-9]", motif_name)                 ~ "KLF/SP",
    grepl("GATA|TAL|SCL|E-box|MyoD|NeuroD|Atoh|Ngn",
          motif_name, ignore.case = TRUE)                              ~ "bHLH",
    grepl("AP-1|FOS|JUN|ATF|CREB|bZIP|Fra|Fosl|MafB|MafA|NRL",
          motif_name, ignore.case = TRUE)                              ~ "bZIP/AP-1",
    grepl("NR[0-9]|RAR|RXR|ER[^G]|AR[^E]|GR[^C]|COUP|THR|PPAR",
          motif_name, ignore.case = TRUE)                              ~ "Nuclear receptor",
    grepl("ETS|ELF|ERG|FLI|ETV|PU.1|GABP|ELK",
          motif_name, ignore.case = TRUE)                              ~ "ETS",
    grepl("RUNX|CBF", motif_name, ignore.case = TRUE)                 ~ "RUNX",
    TRUE                                                               ~ "Other"
  )
}

build_qq_df <- function(dir_x, dir_y) {
  df_x <- parse_homer(dir_x) |>
    dplyr::rename(neg_log10p_x = neg_log10p, pct_target_x = pct_target, pct_bg_x = pct_bg)
  df_y <- parse_homer(dir_y) |>
    dplyr::rename(neg_log10p_y = neg_log10p, pct_target_y = pct_target, pct_bg_y = pct_bg)
  dplyr::inner_join(df_x, df_y, by = "motif_short") |>
    dplyr::filter(pmax(neg_log10p_x, neg_log10p_y) >= min_log10p) |>
    dplyr::mutate(
      tf_family = assign_tf_family(motif_short),
      diag_dist = neg_log10p_y - neg_log10p_x,
      abs_dist  = abs(diag_dist)
    )
}

assign_factor_class <- function(target, assay) {
  dplyr::case_when(
    assay == "Histone ChIP-seq" &
      grepl("H3K27ac|H3K4me1|H3K4me3|H3K36me3|H3K9ac", target,
            ignore.case = TRUE)                                    ~ "Active mark",
    assay == "Histone ChIP-seq" &
      grepl("H3K27me3|H3K9me3|H3K9me2", target,
            ignore.case = TRUE)                                    ~ "Repressive mark",
    grepl("^CTCF$|^RAD21$|^SMC|^STAG", target, ignore.case = TRUE) ~ "CTCF/Cohesin",
    grepl("^SP[0-9]|^KLF", target, ignore.case = TRUE)            ~ "SP/KLF",
    TRUE                                                           ~ "TF"
  )
}

identify_peak_changes <- function(control_peaks, sorbitol_peaks,
                                  buffer = peak_buffer) {
  overlaps <- countOverlaps(control_peaks, sorbitol_peaks + buffer) > 0
  list(retained = control_peaks[overlaps], lost = control_peaks[!overlaps])
}

extract_anchors <- function(gi) {
  unique(c(anchors(gi, "first"), anchors(gi, "second")))
}

## Draw a dashed zoom box + two dotted connector lines from its bottom
## corners down to the top-left and top-right corners of the panel below.
draw_zoom_connectors <- function(px, y_top_row, y_bot_row) {
  
  data_x0 <- px          + panel_w * 0.14   ## left edge of data area
  data_y0 <- y_top_row   + panel_h * 0.08   ## top edge of data area
  data_w  <- panel_w * 0.78
  data_h  <- panel_h * 0.78
  
  box_w <- data_w * zoom_frac
  box_h <- data_h * zoom_frac
  box_x <- data_x0                          ## left-aligned to data area
  box_y <- data_y0 + data_h - box_h         ## bottom-aligned (origin corner)
  
  ## Dashed zoom rectangle
  plotRect(
    x     = box_x + box_w / 2,
    y     = box_y + box_h / 2,
    width = box_w, height = box_h,
    just  = "center",
    lty   = 2, lwd = 0.8,
    linecolor = "grey40",
    fill  = "transparent",
    default.units = "inches"
  )
  
  ## Left connector: bottom-left of box → top-left of panel below
  plotSegments(
    x0 = box_x,        y0 = box_y + box_h,
    x1 = px,           y1 = y_bot_row,
    lty = 3, lwd = 0.6, linecolor = "grey50",
    default.units = "inches"
  )
  
  ## Right connector: bottom-right of box → top-right of panel below
  plotSegments(
    x0 = box_x + box_w, y0 = box_y + box_h,
    x1 = px + panel_w,  y1 = y_bot_row,
    lty = 3, lwd = 0.6, linecolor = "grey50",
    default.units = "inches"
  )
}


# Load data ---------------------------------------------------------------

loops <- readRDS(diff_loops_rds) |>
  interactions() |>
  as.data.frame() |>
  as_ginteractions()

ctcf_ctrl_peaks  <- read_narrowpeaks(ctcf_narrowpeak)
ctcf_sorb_peaks  <- read_narrowpeaks(ctcf_sorb_peak)
rad21_ctrl_peaks <- read_narrowpeaks(rad21_narrowpeak)
rad21_sorb_peaks <- read_narrowpeaks(rad21_sorb_peak)

metadata <- read.delim(metadata_file, header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE)


# Wrangle data ------------------------------------------------------------

## ---- Row 1: HOMER QQ data ----
ctcf_qq_df  <- build_qq_df(file.path(homer_base, "ctcf_lost_peaks"),
                           file.path(homer_base, "ctcf_retained_peaks"))
rad21_qq_df <- build_qq_df(file.path(homer_base, "rad21_lost_peaks"),
                           file.path(homer_base, "rad21_retained_peaks"))
atac_qq_df  <- build_qq_df(file.path(homer_atac,  "atac_lost_anchors"),
                           file.path(homer_atac,  "atac_gained_anchors"))

## ---- Row 3: ENCODE overlap data ----
ctcf_changes  <- identify_peak_changes(ctcf_ctrl_peaks,  ctcf_sorb_peaks)
rad21_changes <- identify_peak_changes(rad21_ctrl_peaks, rad21_sorb_peaks)

gained_loops <- loops[loops$padj < padj_threshold & loops$log2FoldChange > 0]
lost_loops   <- loops[loops$padj < padj_threshold & loops$log2FoldChange < 0]

query_sets <- list(
  ctcf_retained  = ctcf_changes$retained,
  ctcf_lost      = ctcf_changes$lost,
  rad21_retained = rad21_changes$retained,
  rad21_lost     = rad21_changes$lost,
  atac_gained    = extract_anchors(gained_loops),
  atac_lost      = extract_anchors(lost_loops)
)

bed_files    <- list.files(beds_dir, pattern = "\\.bed\\.gz$", full.names = TRUE)
overlap_list <- lapply(bed_files, function(bf) {
  enc <- tryCatch(rtracklayer::import(bf, format = "narrowPeak"),
                  error = function(e) NULL)
  if (is.null(enc) || length(enc) == 0) return(NULL)
  accession <- gsub("\\.bed\\.gz$", "", basename(bf))
  pcts <- vapply(query_sets, function(qs) {
    sum(countOverlaps(qs, enc) > 0) / length(qs) * 100
  }, numeric(1))
  data.frame(file_accession = accession, as.list(pcts), check.names = FALSE)
})

meta_slim <- metadata |>
  dplyr::select(file_accession = `File accession`,
                target_raw     = `Experiment target`,
                assay          = `Assay`) |>
  dplyr::mutate(target = gsub("-human$", "", target_raw))

overlap_collapsed <- overlap_list[!sapply(overlap_list, is.null)] |>
  dplyr::bind_rows() |>
  dplyr::left_join(meta_slim, by = "file_accession") |>
  dplyr::group_by(target, assay) |>
  dplyr::summarise(dplyr::across(c(ctcf_retained, ctcf_lost,
                                   rad21_retained, rad21_lost,
                                   atac_gained, atac_lost),
                                 mean, na.rm = TRUE),
                   n_experiments = dplyr::n(), .groups = "drop") |>
  dplyr::mutate(
    factor_class = assign_factor_class(target, assay),
    factor_class = factor(factor_class,
                          levels = c("Active mark", "Repressive mark",
                                     "CTCF/Cohesin", "SP/KLF", "TF"))
  )


# Visualization -----------------------------------------------------------

okabe_ito <- c("CTCF/BORIS" = "#E69F00", "KLF/SP" = "#0072B2",
               "bHLH" = "#CCCCCC", "bZIP/AP-1" = "#CCCCCC",
               "Nuclear receptor" = "#CCCCCC", "ETS" = "#CCCCCC",
               "RUNX" = "#CCCCCC", "Other TF" = "#CCCCCC", "Other" = "#CCCCCC")

encode_colors <- c("Active mark"     = "#E69F00",
                   "Repressive mark" = "#56B4E9",
                   "CTCF/Cohesin"    = "#009E73",
                   "SP/KLF"          = "#0072B2",
                   "TF"              = "#CCCCCC")
encode_sizes  <- c("Active mark" = 2.0, "Repressive mark" = 2.0,
                   "CTCF/Cohesin" = 2.0, "SP/KLF" = 2.0, "TF" = 1.2)
encode_alpha  <- c("Active mark" = 0.9, "Repressive mark" = 0.9,
                   "CTCF/Cohesin" = 0.9, "SP/KLF" = 0.9, "TF" = 0.3)

theme_s8 <- theme_classic(base_size = 8) +
  theme(plot.title   = element_text(size = 8, face = "bold", hjust = 0.5,
                                    margin = margin(b = 3)),
        axis.title   = element_text(size = 7),
        axis.text    = element_text(size = 6),
        legend.text  = element_text(size = 6),
        legend.title = element_text(size = 6, face = "bold"))

## Row 1: full HOMER QQ (all motifs, CTCF/BORIS + KLF/SP highlighted)
make_qq <- function(df, x_label, y_label, title) {
  ax_max     <- max(c(df$neg_log10p_x, df$neg_log10p_y), na.rm = TRUE) * 1.05
  top_motifs <- df |> dplyr::slice_max(abs_dist, n = top_n_label, with_ties = FALSE)
  df_grey    <- df |> dplyr::filter(!tf_family %in% c("CTCF/BORIS", "KLF/SP"))
  df_hi      <- df |> dplyr::filter( tf_family %in% c("CTCF/BORIS", "KLF/SP"))
  ggplot(df, aes(x = neg_log10p_x, y = neg_log10p_y)) +
    geom_abline(slope = 1, intercept = 0, color = "grey60",
                linetype = "dashed", linewidth = 0.4) +
    geom_point(data = df_grey, aes(color = tf_family), alpha = 0.3,  size = 1.4) +
    geom_point(data = df_hi,   aes(color = tf_family), alpha = 0.85, size = 1.8) +
    ggrepel::geom_text_repel(
      data = top_motifs, aes(label = motif_short, color = tf_family),
      size = 2.2, fontface = "bold", max.overlaps = 25,
      segment.size = 0.3, segment.color = "grey50", show.legend = FALSE
    ) +
    scale_color_manual(values = okabe_ito, name = "TF family",
                       breaks = c("CTCF/BORIS", "KLF/SP")) +
    coord_fixed(xlim = c(0, ax_max), ylim = c(0, ax_max)) +
    labs(title = title,
         x = paste0("-log10(p)  [", x_label, "]"),
         y = paste0("-log10(p)  [", y_label, "]")) +
    theme_s8 + theme(legend.position = "none")
}

## Row 2: SP/KLF zoom (CTCF/BORIS excluded so axis rescales)
make_zoom <- function(df, x_label, y_label, title) {
  df_zoom    <- df |> dplyr::filter(tf_family != "CTCF/BORIS")
  ax_max     <- max(c(df_zoom$neg_log10p_x, df_zoom$neg_log10p_y),
                    na.rm = TRUE) * 1.05
  top_motifs <- df_zoom |>
    dplyr::filter(tf_family == "KLF/SP") |>
    dplyr::slice_max(abs_dist, n = top_n_label, with_ties = FALSE)
  df_grey    <- df_zoom |> dplyr::filter(tf_family != "KLF/SP")
  df_hi      <- df_zoom |> dplyr::filter(tf_family == "KLF/SP")
  ggplot(df_zoom, aes(x = neg_log10p_x, y = neg_log10p_y)) +
    geom_abline(slope = 1, intercept = 0, color = "grey60",
                linetype = "dashed", linewidth = 0.4) +
    geom_point(data = df_grey, aes(color = tf_family), alpha = 0.3,  size = 1.4) +
    geom_point(data = df_hi,   aes(color = tf_family), alpha = 0.85, size = 1.8) +
    ggrepel::geom_text_repel(
      data = top_motifs, aes(label = motif_short, color = tf_family),
      size = 2.2, fontface = "bold", max.overlaps = 25,
      segment.size = 0.3, segment.color = "grey50", show.legend = FALSE
    ) +
    scale_color_manual(values = okabe_ito, name = "TF family",
                       breaks = c("CTCF/BORIS", "KLF/SP")) +
    coord_fixed(xlim = c(0, ax_max), ylim = c(0, ax_max)) +
    labs(title = paste0(title, "  \u2014  SP/KLF zoom"),
         x = paste0("-log10(p)  [", x_label, "]"),
         y = paste0("-log10(p)  [", y_label, "]")) +
    theme_s8 + theme(legend.position = "none")
}

## Row 3: ENCODE ChIP-seq overlap scatter
make_encode <- function(df, x_col, y_col, x_label, y_label, title) {
  df <- df |> dplyr::mutate(x = .data[[x_col]], y = .data[[y_col]],
                            abs_dist = abs(y - x))
  ax_max     <- max(c(df$x, df$y), na.rm = TRUE) * 1.05
  top_points <- df |> dplyr::filter(pmax(x, y) >= min_pct) |>
    dplyr::slice_max(abs_dist, n = top_n_label, with_ties = FALSE)
  df_bg <- df |> dplyr::filter(factor_class == "TF")
  df_fg <- df |> dplyr::filter(factor_class != "TF")
  ggplot(df, aes(x = x, y = y)) +
    geom_abline(slope = 1, intercept = 0, color = "grey60",
                linetype = "dashed", linewidth = 0.4) +
    geom_point(data = df_bg,
               aes(color = factor_class, size = factor_class, alpha = factor_class)) +
    geom_point(data = df_fg,
               aes(color = factor_class, size = factor_class, alpha = factor_class)) +
    ggrepel::geom_text_repel(
      data = top_points, aes(label = target, color = factor_class),
      size = 2.2, fontface = "bold", max.overlaps = 25,
      segment.size = 0.3, segment.color = "grey50", show.legend = FALSE
    ) +
    scale_color_manual(values = encode_colors, name = "Factor class") +
    scale_size_manual(values  = encode_sizes,  name = "Factor class") +
    scale_alpha_manual(values = encode_alpha,  name = "Factor class") +
    coord_fixed(xlim = c(0, ax_max), ylim = c(0, ax_max)) +
    labs(title = title,
         x = paste0("% overlap  [", x_label, "]"),
         y = paste0("% overlap  [", y_label, "]")) +
    theme_s8 + theme(legend.position = "none")
}

## ---- Row 1: HOMER QQ panels ----
p_qq_ctcf  <- make_qq(ctcf_qq_df,
                      x_label = "CTCF Lost",         y_label = "CTCF Retained",
                      title   = "CTCF: Retained vs. Lost")
p_qq_rad21 <- make_qq(rad21_qq_df,
                      x_label = "RAD21 Lost",        y_label = "RAD21 Retained",
                      title   = "RAD21: Retained vs. Lost")
p_qq_atac  <- make_qq(atac_qq_df,
                      x_label = "ATAC Lost anchors", y_label = "ATAC Gained anchors",
                      title   = "ATAC: Gained vs. Lost anchors")

## ---- Row 2: SP/KLF zoom panels ----
p_zoom_ctcf  <- make_zoom(ctcf_qq_df,
                          x_label = "CTCF Lost",         y_label = "CTCF Retained",
                          title   = "CTCF")
p_zoom_rad21 <- make_zoom(rad21_qq_df,
                          x_label = "RAD21 Lost",        y_label = "RAD21 Retained",
                          title   = "RAD21")
p_zoom_atac  <- make_zoom(atac_qq_df,
                          x_label = "ATAC Lost anchors", y_label = "ATAC Gained anchors",
                          title   = "ATAC")

## ---- Row 3: ENCODE scatter panels ----
p_enc_ctcf  <- make_encode(overlap_collapsed,
                           x_col = "ctcf_lost",  y_col = "ctcf_retained",
                           x_label = "CTCF lost",  y_label = "CTCF retained",
                           title   = "ENCODE HEK293 ChIP-seq:\nCTCF Retained vs. Lost")
p_enc_rad21 <- make_encode(overlap_collapsed,
                           x_col = "rad21_lost", y_col = "rad21_retained",
                           x_label = "RAD21 lost", y_label = "RAD21 retained",
                           title   = "ENCODE HEK293 ChIP-seq:\nRAD21 Retained vs. Lost")
p_enc_atac  <- make_encode(overlap_collapsed,
                           x_col = "atac_lost",  y_col = "atac_gained",
                           x_label = "ATAC lost", y_label = "ATAC gained",
                           title   = "ENCODE HEK293 ChIP-seq:\nGained vs. Lost Loop Anchors")


# Save outputs ------------------------------------------------------------

pdf(output_pdf, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = FALSE)

## ---- Panel labels --------------------------------------------------------
label_params <- list(fontsize = 12, fontface = "bold", default.units = "inches")

plotText("A", x = x_col[1] - 0.2, y = y_row1 - 0.2, params = label_params)
plotText("B", x = x_col[2] - 0.2, y = y_row1 - 0.2, params = label_params)
plotText("C", x = x_col[3] - 0.2, y = y_row1 - 0.2, params = label_params)
plotText("D", x = x_col[1] - 0.2, y = y_row2 - 0.2, params = label_params)
plotText("E", x = x_col[2] - 0.2, y = y_row2 - 0.2, params = label_params)
plotText("F", x = x_col[3] - 0.2, y = y_row2 - 0.2, params = label_params)
plotText("G", x = x_col[1] - 0.2, y = y_row3 - 0.2, params = label_params)
plotText("H", x = x_col[2] - 0.2, y = y_row3 - 0.2, params = label_params)
plotText("I", x = x_col[3] - 0.2, y = y_row3 - 0.2, params = label_params)

## ---- Row 1: HOMER QQ (rendered first so connectors overlap their edges) ----
for (pnl in list(list(p_qq_ctcf, x_col[1]),
                 list(p_qq_rad21, x_col[2]),
                 list(p_qq_atac,  x_col[3]))) {
  plotGG(plot = pnl[[1]], x = pnl[[2]], y = y_row1,
         width = panel_w, height = panel_h,
         just = c("left", "top"), default.units = "inches")
}

## ---- Zoom connectors: Row 1 → Row 2 -------------------------------------
for (px in x_col) draw_zoom_connectors(px, y_top_row = y_row1, y_bot_row = y_row2)

## ---- Row 2: SP/KLF zoom (rendered after connectors so panels sit on top) ----
for (pnl in list(list(p_zoom_ctcf, x_col[1]),
                 list(p_zoom_rad21, x_col[2]),
                 list(p_zoom_atac,  x_col[3]))) {
  plotGG(plot = pnl[[1]], x = pnl[[2]], y = y_row2,
         width = panel_w, height = panel_h,
         just = c("left", "top"), default.units = "inches")
}

## ---- Row 3: ENCODE scatter (rendered last) --------------------------------
for (pnl in list(list(p_enc_ctcf, x_col[1]),
                 list(p_enc_rad21, x_col[2]),
                 list(p_enc_atac,  x_col[3]))) {
  plotGG(plot = pnl[[1]], x = pnl[[2]], y = y_row3,
         width = panel_w, height = panel_h,
         just = c("left", "top"), default.units = "inches")
}

dev.off()

message("Saved: ", output_pdf)


# Session info ------------------------------------------------------------

sessionInfo()