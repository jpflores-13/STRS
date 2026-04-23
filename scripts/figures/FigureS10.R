# FigureS10_auxinDegron.R -----------------------------------------------------
# Author:      JP Flores
# Date:        2026-04-23
# Project:     STRS
# Description: Supplementary figure S10 showing CTCF/RAD21 auxin-inducible
#              degron validation. Panel A: representative microscopy images
#              (DAPI, GFP, MERGE) for each condition across CTCF and RAD21.
#              Panel B: violin plot of mean nuclear GFP intensity per condition
#              with pairwise t-test significance brackets.
# Input:       data/raw/imaging/auxin_images_flat/*.tif
#              data/processed/imaging/cellprofiler_output/CTCF_RAD21_auxin_Nuclei.csv
# Output:      figures/FigureS10.pdf
# -----------------------------------------------------------------------------


# Parameters ------------------------------------------------------------------

project_dir  <- "/work/users/j/p/jpflores/projects/STRS"
image_dir    <- file.path(project_dir, "data/raw/imaging/auxin_images_flat")
cp_output    <- file.path(project_dir,
                          "data/processed/imaging/cellprofiler_output",
                          "CTCF_RAD21_auxin_Nuclei.csv")
output_file  <- file.path(project_dir, "figures/FigureS10.pdf")

## Representative field: fld_3 from row C for each condition
## Col 1-3 = CTCF (control, sorbitol, sorbitol+auxin)
## Col 4-6 = RAD21 (control, sorbitol, sorbitol+auxin)
rep_images <- list(
  ctcf_control   = list(row = "C", col = "1", fld = "3"),
  ctcf_sorbitol  = list(row = "C", col = "2", fld = "3"),
  ctcf_auxin     = list(row = "C", col = "3", fld = "3"),
  rad21_control  = list(row = "C", col = "4", fld = "3"),
  rad21_sorbitol = list(row = "C", col = "5", fld = "3"),
  rad21_auxin    = list(row = "C", col = "6", fld = "3")
)

## Column labels and row labels for Panel A
col_labels <- c("untreated", "sorbitol", "sorbitol + auxin",
                "untreated", "sorbitol", "sorbitol + auxin")
row_labels <- c("DAPI", "GFP", "MERGE")

## Plate layout for Panel B
plate_layout <- data.frame(
  Col       = as.character(1:6),
  Protein   = c("CTCF",  "CTCF",  "CTCF",
                "RAD21", "RAD21", "RAD21"),
  Treatment = c("Control", "Sorbitol", "Sorbitol + Auxin",
                "Control", "Sorbitol", "Sorbitol + Auxin")
)

## QC filters
min_area_px <- 500
max_area_px <- 15000

## Pairwise comparisons for stats
stat_comparisons_ctcf <- list(
  c("CTCF\nControl", "CTCF\nSorbitol"),
  c("CTCF\nControl", "CTCF\nSorbitol + Auxin"),
  c("CTCF\nSorbitol", "CTCF\nSorbitol + Auxin")
)
stat_comparisons_rad21 <- list(
  c("RAD21\nControl", "RAD21\nSorbitol"),
  c("RAD21\nControl", "RAD21\nSorbitol + Auxin"),
  c("RAD21\nSorbitol", "RAD21\nSorbitol + Auxin")
)

## Page layout (inches) — landscape, panels side by side
page_width   <- 14
page_height  <- 4.5

## Panel A geometry — left half of page
img_size     <- 0.9   # width and height of each image tile
img_margin   <- 0.04  # gap between tiles
label_left   <- 0.15  # width reserved for row labels on left
label_top    <- 0.50  # height reserved for col labels + protein headers
panel_a_x    <- 0.1   # left edge of Panel A
panel_a_y    <- 0.3   # top of Panel A

## Panel A total height — computed so Panel B can align to it
panel_a_h    <- label_top + 3 * img_size + 2 * img_margin

## Panel B geometry — right half, aligned top and bottom to Panel A
panel_b_x    <- 7.2
panel_b_y    <- panel_a_y
panel_b_w    <- page_width - panel_b_x - 0.3
panel_b_h    <- panel_a_h   # exact same height = square-ish given the width

## Global display multipliers — EBImage already returns normalized floats,
## so no 65535 division needed. FITC is boosted globally (same factor for
## all conditions) to compensate for its inherently dimmer signal.
## This preserves relative intensities across conditions.
dapi_boost   <- 1.5   # DAPI is bright enough; slight boost for display
fitc_boost   <- 10.0  # FITC max ~0.038; boost to ~0.38 for visibility


# Libraries -------------------------------------------------------------------

library(EBImage)
library(plotgardener)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(forcats)
library(readr)


# Load data -------------------------------------------------------------------

## Helper: build filename from well metadata
make_path <- function(row, col, fld, channel) {
  fname <- sprintf("%s_%s_fld_%s_wv_%s_%s.tif", row, col, fld, channel, channel)
  file.path(image_dir, fname)
}

## Helper: load TIFF — EBImage already returns values in [0,1] for 16-bit
## files, so no additional scaling needed
load_tif <- function(path) {
  readImage(path)
}

## Helper: apply false color to a grayscale EBImage object
## color_vec: named numeric vector c(r=, g=, b=) each in [0,1]
colorize <- function(img_gray, r = 0, g = 0, b = 0) {
  ## img_gray is a 2D matrix (single channel); build an RGB Image
  ch_r <- img_gray * r
  ch_g <- img_gray * g
  ch_b <- img_gray * b
  rgbImage(red = ch_r, green = ch_g, blue = ch_b)
}

## Load all 6 conditions × 2 channels (DAPI + FITC)
conditions <- names(rep_images)

imgs <- lapply(conditions, function(cond) {
  meta  <- rep_images[[cond]]
  dapi  <- load_tif(make_path(meta$row, meta$col, meta$fld, "DAPI"))
  fitc  <- load_tif(make_path(meta$row, meta$col, meta$fld, "FITC"))
  list(dapi = dapi, fitc = fitc)
})
names(imgs) <- conditions

## CellProfiler data for Panel B
nuclei_raw <- read_csv(cp_output, show_col_types = FALSE)


# Wrangle data ----------------------------------------------------------------

imgs_display <- lapply(imgs, function(pair) {
  ## Apply global boost multipliers, then clip to [0,1].
  ## Same multiplier applied to all conditions preserves relative intensities.
  dapi_boosted <- pair$dapi * dapi_boost
  fitc_boosted <- pair$fitc * fitc_boost
  dapi_boosted[dapi_boosted > 1] <- 1
  fitc_boosted[fitc_boosted > 1] <- 1
  
  dapi_rgb  <- colorize(dapi_boosted, r = 0, g = 0, b = 1)
  fitc_rgb  <- colorize(fitc_boosted, r = 0, g = 1, b = 0)
  
  ## Merge: additive blend, clip to [0,1]
  merge_rgb <- dapi_rgb + fitc_rgb
  merge_rgb[merge_rgb > 1] <- 1
  
  list(dapi = dapi_rgb, fitc = fitc_rgb, merge = merge_rgb)
})

## Panel B data wrangling
nuclei <- nuclei_raw |>
  dplyr::rename(Row = Metadata_Row, Col = Metadata_Col) |>
  mutate(Col = as.character(Col)) |>
  left_join(plate_layout, by = "Col") |>
  filter(AreaShape_Area >= min_area_px,
         AreaShape_Area <= max_area_px) |>
  mutate(
    Condition = paste(Protein, Treatment, sep = "\n"),
    Condition = fct_relevel(
      Condition,
      "CTCF\nControl", "CTCF\nSorbitol", "CTCF\nSorbitol + Auxin",
      "RAD21\nControl", "RAD21\nSorbitol", "RAD21\nSorbitol + Auxin"
    )
  )

replicate_means <- nuclei |>
  group_by(Row, Col, Condition, Protein, Treatment) |>
  summarize(mean_gfp = mean(Intensity_MeanIntensity_GFP),
            n_nuclei = dplyr::n(),
            .groups  = "drop")


# Visualization ---------------------------------------------------------------

## Panel B ggplot object
protein_colors <- c("CTCF" = "#E69F00", "RAD21" = "#0072B2")

panel_b <- ggplot(
  nuclei,
  aes(x = Condition, y = Intensity_MeanIntensity_GFP, fill = Protein)
) +
  geom_violin(color = NA, alpha = 0.4, scale = "width", trim = TRUE) +
  geom_point(
    data     = replicate_means,
    aes(y    = mean_gfp),
    shape    = 21, size = 2.5, stroke = 0.6, color = "black",
    position = position_jitter(width = 0.05, seed = 42)
  ) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.4, linewidth = 0.5, color = "black") +
  stat_compare_means(
    comparisons   = stat_comparisons_ctcf,
    data          = nuclei |> filter(Protein == "CTCF"),
    method        = "t.test", label = "p.signif",
    tip.length    = 0.01, bracket.size = 0.4, size = 4,
    y.position    = c(0.038, 0.042, 0.046)
  ) +
  stat_compare_means(
    comparisons   = stat_comparisons_rad21,
    data          = nuclei |> filter(Protein == "RAD21"),
    method        = "t.test", label = "p.signif",
    tip.length    = 0.01, bracket.size = 0.4, size = 4,
    y.position    = c(0.038, 0.042, 0.046)
  ) +
  facet_wrap(~ Protein, scales = "free_x") +
  scale_fill_manual(values = protein_colors) +
  scale_y_continuous(limits = c(NA, 0.07)) +
  labs(x = NULL, y = "Mean Nuclear GFP Intensity (a.u.)") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold", size = 11),
    axis.text.x     = element_text(size = 8),
    panel.spacing   = unit(1.5, "lines")
  )


# Save outputs ----------------------------------------------------------------

pdf(output_file, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height,
           showGuides = FALSE, default.units = "inches")

## ---- Panel A label ----------------------------------------------------------
plotText("A", x = panel_a_x, y = panel_a_y,
         fontsize = 14, fontface = "bold",
         just = c("left", "top"), default.units = "inches")

## ---- Protein group headers (CTCF / RAD21) -----------------------------------
a_img_left   <- panel_a_x + label_left
ctcf_center  <- a_img_left + (3 * img_size + 2 * img_margin) / 2
rad21_center <- a_img_left + 3 * (img_size + img_margin) +
  (3 * img_size + 2 * img_margin) / 2 + img_margin

plotText("CTCF",  x = ctcf_center,  y = panel_a_y + 0.05,
         fontsize = 10, fontface = "bold",
         just = c("center", "top"), default.units = "inches")
plotText("RAD21", x = rad21_center, y = panel_a_y + 0.05,
         fontsize = 10, fontface = "bold",
         just = c("center", "top"), default.units = "inches")

## ---- Column labels ----------------------------------------------------------
for (j in seq_along(col_labels)) {
  x_center <- a_img_left + (j - 1) * (img_size + img_margin) +
    img_size / 2 +
    if (j > 3) img_margin else 0
  plotText(col_labels[j],
           x = x_center,
           y = panel_a_y + 0.25,
           fontsize = 7,
           just = c("center", "top"),
           default.units = "inches")
}

## ---- Row labels + images ----------------------------------------------------
channel_keys <- c("dapi", "fitc", "merge")

for (i in seq_along(row_labels)) {
  y_top <- panel_a_y + label_top + (i - 1) * (img_size + img_margin)
  
  plotText(row_labels[i],
           x = panel_a_x + label_left - 0.05,
           y = y_top + img_size / 2,
           fontsize = 8, fontface = "bold",
           rot = 90,
           just = c("center", "bottom"),
           default.units = "inches")
  
  for (j in seq_along(conditions)) {
    x_left <- a_img_left + (j - 1) * (img_size + img_margin) +
      if (j > 3) img_margin else 0
    
    img_to_plot <- imgs_display[[conditions[j]]][[channel_keys[i]]]
    
    plotRaster(img_to_plot,
               x      = x_left,
               y      = y_top,
               width  = img_size,
               height = img_size,
               just   = c("left", "top"),
               default.units = "inches")
  }
}

## ---- Panel B label + plot ---------------------------------------------------
plotText("B", x = panel_b_x, y = panel_b_y,
         fontsize = 14, fontface = "bold",
         just = c("left", "top"), default.units = "inches")

plotGG(panel_b,
       x      = panel_b_x + 0.2,
       y      = panel_b_y + 0.1,
       width  = panel_b_w,
       height = panel_b_h,
       just   = c("left", "top"),
       default.units = "inches")

dev.off()


# Session info ----------------------------------------------------------------

sessionInfo()