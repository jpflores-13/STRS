# HOMER Motif QQ-Plot (Retained vs. Lost / Gained vs. Lost) ------------------
# Author:      JP Flores
# Date:        2026-04-06
# Project:     STRS
# Description: Three-panel QQ-style scatter plot comparing HOMER known-motif
#              -log10(p-values) between paired peak sets: CTCF retained vs.
#              lost, RAD21 retained vs. lost, and ATAC gained vs. lost loop
#              anchors. Points above the diagonal are more enriched in the
#              "surviving" peak set (retained/gained); points below are more
#              enriched in the lost set. Top motifs by distance from diagonal
#              are labeled per panel. Colorblind-friendly (Okabe-Ito palette)
#              coloring by TF family.
#              Row 1 (A–C): full QQ overview with CTCF/BORIS + KLF/SP highlighted.
#              Row 2 (D–F): SP/KLF zoom — CTCF/BORIS removed so axis rescales
#                to reveal SP/KLF family separation. Dashed zoom rectangles on
#                each Row 1 panel connect via dotted lines to the corresponding
#                Row 2 panel directly below.
# Input:       data/processed/cutntag/homer_motifs/<condition>/knownResults.txt
#              data/processed/atac/homer_input/<condition>/knownResults.txt
# Output:      plots/homerQQ.pdf
# -----------------------------------------------------------------------------


# Parameters ------------------------------------------------------------------

project_dir  <- "/work/users/j/p/jpflores/projects/STRS"
homer_base   <- file.path(project_dir, "data/processed/cutntag/homer_motifs")
homer_atac   <- file.path(project_dir, "data/processed/atac/homer_input")
output_pdf   <- file.path(project_dir, "plots/homerQQ.pdf")

## Minimum -log10(p) in at least one condition to keep a motif
min_log10p   <- 2

## Number of top motifs (by distance from diagonal) to label per panel
top_n_label  <- 10

## Page and panel dimensions (inches)
page_width   <- 14
panel_width  <- 4.0
panel_height <- 4.5
col_gap      <- 0.5
row_gap      <- 0.8    ## vertical gap between rows (space for connectors)
left_margin  <- 0.5
top_margin   <- 0.5
label_offset_x <- 0.1
label_offset_y <- 0.15

## Zoom box geometry — fraction of the data area covered by the zoom box
zoom_frac <- 0.35

## Derived positions
x_positions <- c(
  left_margin,
  left_margin + panel_width + col_gap,
  left_margin + (panel_width + col_gap) * 2
)
row1_y      <- top_margin
row2_y      <- row1_y + panel_height + row_gap
page_height <- row2_y + panel_height + 0.3


# Libraries -------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggrepel)
library(plotgardener)


# Utility functions -----------------------------------------------------------

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
    grepl("CTCF|BORIS|CTFL|CTCFL", motif_name, ignore.case = TRUE)          ~ "CTCF/BORIS",
    grepl("^KLF|^Klf|^SP[0-9]|^Sp[0-9]", motif_name)                       ~ "KLF/SP",
    grepl("GATA|TAL|SCL|E-box|MyoD|NeuroD|Atoh|Ngn", motif_name,
          ignore.case = TRUE)                                                ~ "bHLH",
    grepl("AP-1|FOS|JUN|ATF|CREB|bZIP|Fra|Fosl|MafB|MafA|NRL",
          motif_name, ignore.case = TRUE)                                    ~ "bZIP/AP-1",
    grepl("NR[0-9]|RAR|RXR|ER[^G]|AR[^E]|GR[^C]|COUP|THR|PPAR|NF-[IE]",
          motif_name, ignore.case = TRUE)                                    ~ "Nuclear receptor",
    grepl("ETS|ELF|ERG|FLI|ETV|PU.1|GABP|ELK", motif_name,
          ignore.case = TRUE)                                                ~ "ETS",
    grepl("RUNX|CBF", motif_name, ignore.case = TRUE)                       ~ "RUNX",
    grepl("SOX|OCT|POU|NANOG|Hox|HOX|PAX|FOX|IRF|STAT|NF-?[Kk][Bb]",
          motif_name, ignore.case = TRUE)                                    ~ "Other TF",
    TRUE                                                                     ~ "Other"
  )
}

build_qq_df <- function(dir_x, dir_y, label_x, label_y) {
  df_x <- parse_homer(dir_x) |>
    dplyr::rename(neg_log10p_x = neg_log10p, pct_target_x = pct_target, pct_bg_x = pct_bg)
  df_y <- parse_homer(dir_y) |>
    dplyr::rename(neg_log10p_y = neg_log10p, pct_target_y = pct_target, pct_bg_y = pct_bg)
  dplyr::inner_join(df_x, df_y, by = "motif_short") |>
    dplyr::filter(pmax(neg_log10p_x, neg_log10p_y) >= min_log10p) |>
    dplyr::mutate(
      tf_family    = assign_tf_family(motif_short),
      diag_dist    = neg_log10p_y - neg_log10p_x,
      abs_dist     = abs(diag_dist),
      label_x_name = label_x,
      label_y_name = label_y
    )
}

homer_colors <- c(
  "CTCF/BORIS"       = "#E69F00",
  "KLF/SP"           = "#0072B2",
  "bHLH"             = "#CCCCCC",
  "bZIP/AP-1"        = "#CCCCCC",
  "Nuclear receptor"  = "#CCCCCC",
  "ETS"              = "#CCCCCC",
  "RUNX"             = "#CCCCCC",
  "Other TF"         = "#CCCCCC",
  "Other"            = "#CCCCCC"
)

make_qq_panel <- function(df, x_label, y_label, title) {
  ax_max     <- max(c(df$neg_log10p_x, df$neg_log10p_y), na.rm = TRUE) * 1.05
  top_motifs <- df |> dplyr::slice_max(abs_dist, n = top_n_label, with_ties = FALSE)
  df_grey    <- df |> dplyr::filter(!tf_family %in% c("CTCF/BORIS", "KLF/SP"))
  df_hi      <- df |> dplyr::filter(tf_family  %in% c("CTCF/BORIS", "KLF/SP"))
  ggplot(df, aes(x = neg_log10p_x, y = neg_log10p_y)) +
    geom_abline(slope = 1, intercept = 0,
                color = "grey60", linetype = "dashed", linewidth = 0.5) +
    geom_point(data = df_grey, aes(color = tf_family), alpha = 0.35, size = 1.6) +
    geom_point(data = df_hi,   aes(color = tf_family), alpha = 0.85, size = 2.0) +
    ggrepel::geom_text_repel(
      data = top_motifs, aes(label = motif_short, color = tf_family),
      size = 2.6, fontface = "bold", max.overlaps = 25,
      segment.size = 0.3, segment.color = "grey50", show.legend = FALSE
    ) +
    scale_color_manual(values = homer_colors, name = "TF family",
                       breaks = c("CTCF/BORIS", "KLF/SP")) +
    coord_fixed(xlim = c(0, ax_max), ylim = c(0, ax_max)) +
    labs(
      title = title,
      x     = paste0("\u2212log\u2081\u2080(p)  [", x_label, "]"),
      y     = paste0("\u2212log\u2081\u2080(p)  [", y_label, "]")
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title      = element_text(size = 10, face = "bold",
                                     hjust = 0.5, margin = margin(b = 4)),
      axis.title      = element_text(size = 9),
      axis.text       = element_text(size = 8),
      legend.position = "none"
    )
}

make_zoom_panel <- function(df, x_label, y_label, title) {
  df_zoom    <- df |> dplyr::filter(tf_family != "CTCF/BORIS")
  ax_max     <- max(c(df_zoom$neg_log10p_x, df_zoom$neg_log10p_y), na.rm = TRUE) * 1.05
  top_motifs <- df_zoom |>
    dplyr::filter(tf_family == "KLF/SP") |>
    dplyr::slice_max(abs_dist, n = top_n_label, with_ties = FALSE)
  df_grey    <- df_zoom |> dplyr::filter(tf_family != "KLF/SP")
  df_hi      <- df_zoom |> dplyr::filter(tf_family == "KLF/SP")
  ggplot(df_zoom, aes(x = neg_log10p_x, y = neg_log10p_y)) +
    geom_abline(slope = 1, intercept = 0,
                color = "grey60", linetype = "dashed", linewidth = 0.5) +
    geom_point(data = df_grey, aes(color = tf_family), alpha = 0.35, size = 1.6) +
    geom_point(data = df_hi,   aes(color = tf_family), alpha = 0.85, size = 2.0) +
    ggrepel::geom_text_repel(
      data = top_motifs, aes(label = motif_short, color = tf_family),
      size = 2.6, fontface = "bold", max.overlaps = 25,
      segment.size = 0.3, segment.color = "grey50", show.legend = FALSE
    ) +
    scale_color_manual(values = homer_colors, name = "TF family",
                       breaks = c("CTCF/BORIS", "KLF/SP")) +
    coord_fixed(xlim = c(0, ax_max), ylim = c(0, ax_max)) +
    labs(
      title = paste0(title, "  \u2014  SP/KLF zoom"),
      x     = paste0("\u2212log\u2081\u2080(p)  [", x_label, "]"),
      y     = paste0("\u2212log\u2081\u2080(p)  [", y_label, "]")
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title      = element_text(size = 10, face = "bold",
                                     hjust = 0.5, margin = margin(b = 4)),
      axis.title      = element_text(size = 9),
      axis.text       = element_text(size = 8),
      legend.position = "none"
    )
}

## Draw a dashed zoom box on the Row 1 panel + two dotted connector lines
## from its bottom corners down to the top corners of the Row 2 panel below.
draw_zoom_connectors <- function(px, y_top_row, y_bot_row) {
  data_x0 <- px        + panel_width  * 0.14   ## left edge of data area
  data_y0 <- y_top_row + panel_height * 0.08   ## top edge of data area
  data_w  <- panel_width  * 0.78
  data_h  <- panel_height * 0.78
  
  box_w <- data_w * zoom_frac
  box_h <- data_h * zoom_frac
  box_x <- data_x0                        ## left-aligned to data area
  box_y <- data_y0 + data_h - box_h       ## bottom-aligned (origin corner)
  
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
    x1 = px + panel_width, y1 = y_bot_row,
    lty = 3, lwd = 0.6, linecolor = "grey50",
    default.units = "inches"
  )
}

place_panel <- function(plot, label, x, y) {
  plotText(
    label = label, x = x - label_offset_x, y = y - label_offset_y,
    fontsize = 14, fontface = "bold", default.units = "inches"
  )
  plotGG(
    plot = plot, x = x, y = y,
    width = panel_width, height = panel_height,
    just = c("left", "top"), default.units = "inches"
  )
}


# Load data -------------------------------------------------------------------

ctcf_df <- build_qq_df(
  dir_x   = file.path(homer_base, "ctcf_lost_peaks"),
  dir_y   = file.path(homer_base, "ctcf_retained_peaks"),
  label_x = "Lost",
  label_y = "Retained"
)

rad21_df <- build_qq_df(
  dir_x   = file.path(homer_base, "rad21_lost_peaks"),
  dir_y   = file.path(homer_base, "rad21_retained_peaks"),
  label_x = "Lost",
  label_y = "Retained"
)

atac_df <- build_qq_df(
  dir_x   = file.path(homer_atac, "atac_lost_anchors"),
  dir_y   = file.path(homer_atac, "atac_gained_anchors"),
  label_x = "Lost anchors",
  label_y = "Gained anchors"
)


# Visualization ---------------------------------------------------------------

p_ctcf_qq <- make_qq_panel(
  ctcf_df,
  x_label = "CTCF lost peaks",
  y_label = "CTCF retained peaks",
  title   = "CTCF motif enrichment:\nRetained vs. Lost peaks"
)

p_rad21_qq <- make_qq_panel(
  rad21_df,
  x_label = "RAD21 lost peaks",
  y_label = "RAD21 retained peaks",
  title   = "RAD21 motif enrichment:\nRetained vs. Lost peaks"
)

p_atac_qq <- make_qq_panel(
  atac_df,
  x_label = "Lost loop anchors",
  y_label = "Gained loop anchors",
  title   = "Loop anchor motif enrichment:\nGained vs. Lost anchors"
)

p_ctcf_zoom <- make_zoom_panel(
  ctcf_df,
  x_label = "CTCF lost peaks",
  y_label = "CTCF retained peaks",
  title   = "CTCF motif enrichment"
)

p_rad21_zoom <- make_zoom_panel(
  rad21_df,
  x_label = "RAD21 lost peaks",
  y_label = "RAD21 retained peaks",
  title   = "RAD21 motif enrichment"
)

p_atac_zoom <- make_zoom_panel(
  atac_df,
  x_label = "Lost loop anchors",
  y_label = "Gained loop anchors",
  title   = "Loop anchor motif enrichment"
)


# Save outputs ----------------------------------------------------------------

pdf(output_pdf, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = FALSE)

## ---- Row 1: QQ overview (rendered first so connectors overlap their edges) --
place_panel(p_ctcf_qq,  "A", x_positions[1], row1_y)
place_panel(p_rad21_qq, "B", x_positions[2], row1_y)
place_panel(p_atac_qq,  "C", x_positions[3], row1_y)

## ---- Zoom connectors: Row 1 → Row 2 -----------------------------------------
for (px in x_positions) draw_zoom_connectors(px, y_top_row = row1_y, y_bot_row = row2_y)

## ---- Row 2: SP/KLF zoom (rendered after connectors so panels sit on top) ----
place_panel(p_ctcf_zoom,  "D", x_positions[1], row2_y)
place_panel(p_rad21_zoom, "E", x_positions[2], row2_y)
place_panel(p_atac_zoom,  "F", x_positions[3], row2_y)

dev.off()

message("Figure saved to: ", output_pdf)


# Session info ----------------------------------------------------------------

sessionInfo()