# ##############################################################################
# filename:    hic_timecourse.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: plotgardener APA matrix grid for gained and lost loops across
#              the hyperosmotic stress timecourse (0h–24h); each matrix uses
#              its own background-scaled zrange
# ##############################################################################

# Libraries ----
library(plotgardener)
library(tidyverse)
library(RColorBrewer)

# Parameters ----
apa_dir    <- "data/processed/hic/normalizedAPA"
output_pdf <- "plots/hic_timecourse.pdf"
timepoints <- c("0h", "10m", "30m", "1h", "3h", "6h", "12h", "24h")
arb_val    <- 2.75
page_width  <- 8
page_height <- 2

# Data import ----
time_order <- setNames(seq_along(timepoints), timepoints)

sort_by_timepoint <- function(files) {
  files[order(sapply(files, function(x) {
    tp <- str_extract(x, "0h|10m|30m|1h|3h|6h|12h|24h")
    time_order[tp]
  }))]
}

gained_files   <- list.files(apa_dir, full.names = TRUE,
                              pattern = "tc_gainedLoops.*normalized.rds") |>
  sort_by_timepoint()

lost_files     <- list.files(apa_dir, full.names = TRUE,
                              pattern = "tc_lostLoops.*normalized.rds") |>
  sort_by_timepoint()

gained_matrices <- lapply(gained_files, readRDS)
lost_matrices   <- lapply(lost_files,   readRDS)

# Analysis ----
calculate_background_median <- function(mat) {
  top_left     <- mat[1:5, 1:5]
  bottom_right <- mat[17:21, 17:21]
  bg_vals      <- c(as.vector(top_left), as.vector(bottom_right))
  bg_vals      <- bg_vals[!is.na(bg_vals)]
  if (length(bg_vals) == 0 || all(bg_vals == 0)) 0.01 else median(bg_vals)
}

# Visualization ----
apaParams <- pgParams(
  assembly = "hg38",
  height   = 0.75,
  width    = 0.75,
  fontsize = 5,
  palette  = colorRampPalette(brewer.pal(n = 9, "YlGnBu"))
)

# Save outputs ----
pdf(output_pdf, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height,
           default.units = "inches", showGuides = FALSE)

for (i in seq_along(timepoints)) {
  x_pos    <- 0.2 + (i - 1) * 0.9
  gained_bg <- calculate_background_median(gained_matrices[[i]])
  lost_bg   <- calculate_background_median(lost_matrices[[i]])

  plotMatrix(gained_matrices[[i]],
             x = x_pos, y = 0.2,
             zrange = c(0, arb_val) * gained_bg,
             params = apaParams) |>
    annoHeatmapLegend(x = x_pos + 0.8, y = 0.2,
                      width = 0.05, height = 0.75,
                      fontcolor = "black", fontsize = 5,
                      zrange = c(0, arb_val) * gained_bg)

  plotMatrix(lost_matrices[[i]],
             x = x_pos, y = 1.0,
             zrange = c(0, arb_val) * lost_bg,
             params = apaParams) |>
    annoHeatmapLegend(x = x_pos + 0.8, y = 1.0,
                      width = 0.05, height = 0.75,
                      fontcolor = "black", fontsize = 5,
                      zrange = c(0, arb_val) * lost_bg)

  plotText(label = timepoints[i],
           x = x_pos + 0.375, y = 0.1,
           params = apaParams, fontsize = 8, just = "center")
}

plotText(label = "gained", x = 0.1, y = 0.6,
         fontcolor = "#F8766D", fontsize = 10, rot = 90)
plotText(label = "lost",   x = 0.1, y = 1.4,
         fontcolor = "#619CFF", fontsize = 10, rot = 90)

dev.off()

sessionInfo()
