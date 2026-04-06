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
# Input:       data/processed/cutntag/homer_motifs/<condition>/knownResults.txt
#              data/processed/atac/homer_motifs/<condition>/knownResults.txt
# Output:      plots/homerQQ.pdf
# -----------------------------------------------------------------------------


# Parameters ------------------------------------------------------------------

project_dir  <- "/work/users/j/p/jpflores/projects/STRS"
homer_base   <- file.path(project_dir, "data/processed/cutntag/homer_motifs")
homer_atac   <- file.path(project_dir, "data/processed/atac/homer_input")
output_pdf   <- file.path(project_dir, "plots/homerQQ.pdf")

## Minimum -log10(p) in at least one condition to keep a motif
min_log10p   <- 2      # corresponds to p <= 0.01

## Number of top motifs (by distance from diagonal) to label per panel
top_n_label  <- 10

## Page and panel dimensions (inches)
page_width   <- 14
page_height  <- 11.5   # two rows of panels + margins
panel_width  <- 4.0
panel_height <- 4.5
row_gap      <- 0.8    # vertical gap between overview and zoom rows


# Libraries -------------------------------------------------------------------

library(plotgardener)
library(ggplot2)
library(ggrepel)
library(dplyr)


# Utility: parse knownResults.txt ---------------------------------------------

## Reads a HOMER knownResults.txt and returns a tidy data frame with
## motif short name, -log10(p-value), and % target / % background columns.
## HOMER stores natural-log p-values in column "Log.P.value"; we convert.
parse_homer <- function(homer_dir) {
  results_file <- file.path(homer_dir, "knownResults.txt")
  if (!file.exists(results_file)) {
    stop("knownResults.txt not found in: ", homer_dir)
  }
  read.delim(results_file, header = TRUE, stringsAsFactors = FALSE) |>
    dplyr::mutate(
      neg_log10p  = -Log.P.value / log(10),   # convert ln(p) -> -log10(p)
      pct_target  = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif)),
      pct_bg      = as.numeric(gsub("%", "", X..of.Background.Sequences.with.Motif)),
      motif_short = gsub("/.*", "", Motif.Name) |> gsub("\\(.*", "", x = _)
    ) |>
    dplyr::select(motif_short, neg_log10p, pct_target, pct_bg)
}


# Utility: assign TF family ---------------------------------------------------

## Assigns a broad TF family label to a motif short name for coloring.
## Uses simple regex matching; "Other" catches everything unclassified.
assign_tf_family <- function(motif_name) {
  dplyr::case_when(
    grepl("CTCF|BORIS|CTFL|CTCFL", motif_name, ignore.case = TRUE)             ~ "CTCF/BORIS",
    grepl("^KLF|^Klf|^SP[0-9]|^Sp[0-9]", motif_name)                          ~ "KLF/SP",
    grepl("GATA|TAL|SCL|E-box|MyoD|NeuroD|Atoh|Ngn", motif_name,
          ignore.case = TRUE)                                                   ~ "bHLH",
    grepl("AP-1|FOS|JUN|ATF|CREB|bZIP|Fra|Fosl|MafB|MafA|NRL",
          motif_name, ignore.case = TRUE)                                       ~ "bZIP/AP-1",
    grepl("NR[0-9]|RAR|RXR|ER[^G]|AR[^E]|GR[^C]|COUP|THR|PPAR|NF-[IE]",
          motif_name, ignore.case = TRUE)                                       ~ "Nuclear receptor",
    grepl("ETS|ELF|ERG|FLI|ETV|PU.1|GABP|ELK", motif_name,
          ignore.case = TRUE)                                                   ~ "ETS",
    grepl("RUNX|CBF", motif_name, ignore.case = TRUE)                          ~ "RUNX",
    grepl("SOX|OCT|POU|NANOG|Hox|HOX|PAX|FOX|IRF|STAT|NF-?[Kk][Bb]",
          motif_name, ignore.case = TRUE)                                       ~ "Other TF",
    TRUE                                                                        ~ "Other"
  )
}


# Utility: build QQ data frame ------------------------------------------------

## Joins two parsed HOMER results by motif short name, computes diagonal
## distance, and applies minimum enrichment filter.
build_qq_df <- function(dir_x, dir_y, label_x, label_y) {
  df_x <- parse_homer(dir_x) |>
    dplyr::rename(neg_log10p_x = neg_log10p,
                  pct_target_x = pct_target,
                  pct_bg_x     = pct_bg)
  
  df_y <- parse_homer(dir_y) |>
    dplyr::rename(neg_log10p_y = neg_log10p,
                  pct_target_y = pct_target,
                  pct_bg_y     = pct_bg)
  
  dplyr::inner_join(df_x, df_y, by = "motif_short") |>
    dplyr::filter(
      pmax(neg_log10p_x, neg_log10p_y) >= min_log10p
    ) |>
    dplyr::mutate(
      tf_family    = assign_tf_family(motif_short),
      diag_dist    = neg_log10p_y - neg_log10p_x,  # + = above diagonal
      abs_dist     = abs(diag_dist),
      label_x_name = label_x,
      label_y_name = label_y
    )
}


# Utility: build one QQ panel -------------------------------------------------

## Returns a ggplot object for one comparison. Axes are linear and share a
## common limit set by the true maximum of both conditions, so CTCF/BORIS
## outliers remain fully visible and the diagonal reference line is exact.
## Top motifs by absolute distance from the diagonal are labeled.
make_qq_panel <- function(df, x_label, y_label, title) {
  
  ## Only CTCF/BORIS and KLF/SP receive distinct Okabe-Ito colors.
  ## All other TF families are collapsed to the same light grey so they
  ## recede visually without being removed from the plot.
  ## Source for highlighted colors: Okabe & Ito 2008; Color Universal Design
  okabe_ito <- c(
    "CTCF/BORIS"       = "#E69F00",   # orange  (highlighted)
    "KLF/SP"           = "#0072B2",   # blue    (highlighted)
    "bHLH"             = "#CCCCCC",   # grey
    "bZIP/AP-1"        = "#CCCCCC",   # grey
    "Nuclear receptor"  = "#CCCCCC",   # grey
    "ETS"              = "#CCCCCC",   # grey
    "RUNX"             = "#CCCCCC",   # grey
    "Other TF"         = "#CCCCCC",   # grey
    "Other"            = "#CCCCCC"    # grey
  )
  
  ## Shared axis limit: true max of both axes, with a little padding
  ax_max <- max(c(df$neg_log10p_x, df$neg_log10p_y), na.rm = TRUE) * 1.05
  
  ## Top motifs to label (largest absolute distance from diagonal)
  top_motifs <- df |>
    dplyr::slice_max(order_by = abs_dist, n = top_n_label, with_ties = FALSE)
  
  ## Draw grey (background) points first, then colored (foreground) points on
  ## top, so CTCF/BORIS and KLF/SP are never occluded by the grey cloud.
  df_grey      <- df |> dplyr::filter(!tf_family %in% c("CTCF/BORIS", "KLF/SP"))
  df_highlight <- df |> dplyr::filter(tf_family  %in% c("CTCF/BORIS", "KLF/SP"))
  
  ggplot(df, aes(x = neg_log10p_x, y = neg_log10p_y)) +
    geom_abline(slope = 1, intercept = 0,
                color = "grey60", linetype = "dashed", linewidth = 0.5) +
    geom_point(data  = df_grey,
               aes(color = tf_family), alpha = 0.35, size = 1.6) +
    geom_point(data  = df_highlight,
               aes(color = tf_family), alpha = 0.85, size = 2.0) +
    ggrepel::geom_text_repel(
      data          = top_motifs,
      aes(label = motif_short, color = tf_family),
      size          = 2.6,
      fontface      = "bold",
      max.overlaps  = 25,
      segment.size  = 0.3,
      segment.color = "grey50",
      show.legend   = FALSE
    ) +
    scale_color_manual(values = okabe_ito, name = "TF family",
                       breaks = c("CTCF/BORIS", "KLF/SP")) +
    coord_fixed(xlim = c(0, ax_max), ylim = c(0, ax_max)) +
    labs(
      title = title,
      x     = paste0("-log10(pval)  [", x_label, "]"),
      y     = paste0("-log10(pval)  [", y_label, "]")
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title       = element_text(size = 10, face = "bold",
                                      hjust = 0.5, margin = margin(b = 4)),
      axis.title       = element_text(size = 9),
      axis.text        = element_text(size = 8),
      legend.title     = element_text(size = 8, face = "bold"),
      legend.text      = element_text(size = 7),
      legend.key.size  = unit(0.45, "cm"),
      legend.position  = "right"
    )
}


# Utility: build one KLF/SP zoom panel ----------------------------------------

## Identical to make_qq_panel except:
##   - CTCF/BORIS points are excluded so they don't set the axis range
##   - ax_max is computed from KLF/SP + grey cloud only
##   - Only KLF/SP members are labeled
##   - Legend is suppressed (overview row above carries it)
##   - Title gets a "(KLF/SP zoom)" suffix so the row is self-explanatory
make_zoom_panel <- function(df, x_label, y_label, title) {
  
  okabe_ito <- c(
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
  
  ## Drop CTCF/BORIS so they don't dominate the axis limit
  df_zoom <- df |> dplyr::filter(tf_family != "CTCF/BORIS")
  
  ax_max <- max(c(df_zoom$neg_log10p_x, df_zoom$neg_log10p_y),
                na.rm = TRUE) * 1.05
  
  ## Label only KLF/SP members, ranked by diagonal distance
  top_motifs <- df_zoom |>
    dplyr::filter(tf_family == "KLF/SP") |>
    dplyr::slice_max(order_by = abs_dist, n = top_n_label, with_ties = FALSE)
  
  df_grey      <- df_zoom |> dplyr::filter(tf_family != "KLF/SP")
  df_highlight <- df_zoom |> dplyr::filter(tf_family == "KLF/SP")
  
  ggplot(df_zoom, aes(x = neg_log10p_x, y = neg_log10p_y)) +
    geom_abline(slope = 1, intercept = 0,
                color = "grey60", linetype = "dashed", linewidth = 0.5) +
    geom_point(data  = df_grey,
               aes(color = tf_family), alpha = 0.35, size = 1.6) +
    geom_point(data  = df_highlight,
               aes(color = tf_family), alpha = 0.85, size = 2.0) +
    ggrepel::geom_text_repel(
      data          = top_motifs,
      aes(label = motif_short, color = tf_family),
      size          = 2.6,
      fontface      = "bold",
      max.overlaps  = 25,
      segment.size  = 0.3,
      segment.color = "grey50",
      show.legend   = FALSE
    ) +
    scale_color_manual(values = okabe_ito, name = "TF family",
                       breaks = c("CTCF/BORIS", "KLF/SP")) +
    coord_fixed(xlim = c(0, ax_max), ylim = c(0, ax_max)) +
    labs(
      title = paste0(title, "  \u2014 Zoom"),
      x     = paste0("-log10(pval)  [", x_label, "]"),
      y     = paste0("-log10(pval)  [", y_label, "]")
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title       = element_text(size = 10, face = "bold",
                                      hjust = 0.5, margin = margin(b = 4)),
      axis.title       = element_text(size = 9),
      axis.text        = element_text(size = 8),
      legend.position  = "none"   # suppressed — overview row carries the legend
    )
}


# Load data -------------------------------------------------------------------

## CTCF: retained peaks (y) vs. lost peaks (x)
ctcf_df <- build_qq_df(
  dir_x   = file.path(homer_base, "ctcf_lost_peaks"),
  dir_y   = file.path(homer_base, "ctcf_retained_peaks"),
  label_x = "Lost",
  label_y = "Retained"
)

## RAD21: retained peaks (y) vs. lost peaks (x)
rad21_df <- build_qq_df(
  dir_x   = file.path(homer_base, "rad21_lost_peaks"),
  dir_y   = file.path(homer_base, "rad21_retained_peaks"),
  label_x = "Lost",
  label_y = "Retained"
)

## ATAC: gained anchors (y) vs. lost anchors (x)
atac_df <- build_qq_df(
  dir_x   = file.path(homer_atac, "atac_lost_anchors"),
  dir_y   = file.path(homer_atac, "atac_gained_anchors"),
  label_x = "Lost anchors",
  label_y = "Gained anchors"
)


# Visualization ---------------------------------------------------------------

p_ctcf  <- make_qq_panel(ctcf_df,
                         x_label = "CTCF Lost",
                         y_label = "CTCF Retained",
                         title   = "CTCF: Retained vs. Lost peaks")

p_rad21 <- make_qq_panel(rad21_df,
                         x_label = "RAD21 Lost",
                         y_label = "RAD21 Retained",
                         title   = "RAD21: Retained vs. Lost peaks")

p_atac  <- make_qq_panel(atac_df,
                         x_label = "ATAC Lost anchors",
                         y_label = "ATAC Gained anchors",
                         title   = "ATAC: Gained vs. Lost peaks")

## Zoom panels: same data, CTCF/BORIS excluded, axis rescaled to KLF/SP range
p_ctcf_zoom  <- make_zoom_panel(ctcf_df,
                                x_label = "CTCF Lost",
                                y_label = "CTCF Retained",
                                title   = "CTCF")

p_rad21_zoom <- make_zoom_panel(rad21_df,
                                x_label = "RAD21 Lost",
                                y_label = "RAD21 Retained",
                                title   = "RAD21")

p_atac_zoom  <- make_zoom_panel(atac_df,
                                x_label = "ATAC Lost anchors",
                                y_label = "ATAC Gained anchors",
                                title   = "ATAC")


# Save outputs ----------------------------------------------------------------

pdf(output_pdf, width = page_width, height = page_height)

pageCreate(width  = page_width,
           height = page_height,
           showGuides = FALSE)

## Shared x positions for both rows
x_positions <- c(0.5, 0.5 + panel_width + 0.5, 0.5 + (panel_width + 0.5) * 2)

## Row 1 (overview): y origin at 0.5 in
## Row 2 (zoom):     y origin below row 1 + gap
row1_y <- 0.5
row2_y <- row1_y + panel_height + row_gap

panel_labels <- c("A", "B", "C", "D", "E", "F")
panels       <- list(p_ctcf, p_rad21, p_atac,
                     p_ctcf_zoom, p_rad21_zoom, p_atac_zoom)
y_positions  <- c(rep(row1_y, 3), rep(row2_y, 3))

for (i in seq_along(panels)) {
  plotText(
    label    = panel_labels[i],
    x        = x_positions[((i - 1) %% 3) + 1] - 0.1,
    y        = y_positions[i] - 0.15,
    fontsize = 14,
    fontface = "bold"
  )
  plotGG(
    plot   = panels[[i]],
    x      = x_positions[((i - 1) %% 3) + 1],
    y      = y_positions[i],
    width  = panel_width,
    height = panel_height,
    just   = c("left", "top")
  )
}

dev.off()


# Session info ----------------------------------------------------------------

sessionInfo()