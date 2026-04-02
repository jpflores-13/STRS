# ##############################################################################
# filename:    hic_timecourse_loop_enrichment_line_plot.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Dual-axis line plot of APA loop enrichment scores over the
#              hyperosmotic stress timecourse; gained (left axis) and lost
#              (right axis) loops use independent y-axis scales
# ##############################################################################

# Libraries ----
library(tidyverse)
library(scales)
library(grid)
library(gridExtra)

# Parameters ----
enrichment_rds <- "data/processed/hic/apa_loop_enrichment_summary_original_bg.rds"
output_pdf     <- "plots/hic_timecourse_loop_enrichment_line_plot.pdf"
timepoints     <- c("0h", "10m", "30m", "1h", "3h", "6h", "12h", "24h")
color_gained   <- "#F8766D"
color_lost     <- "#619CFF"
page_width     <- 8
page_height    <- 3

# Data import ----
enrichment_summary <- readRDS(enrichment_rds)

gained_data <- enrichment_summary |> filter(loop_type == "Gained")
lost_data   <- enrichment_summary |> filter(loop_type == "Lost")

# Analysis ----
y_range_gained <- c(1, max(gained_data$mean_score, gained_data$upper_ci, na.rm = TRUE))
y_range_lost   <- c(1, max(lost_data$mean_score,   lost_data$upper_ci,   na.rm = TRUE))

lost_to_gained_transform <- function(x) {
  gained_span <- y_range_gained[2] - y_range_gained[1]
  lost_span   <- y_range_lost[2]   - y_range_lost[1]
  (x - y_range_lost[1]) * (gained_span / lost_span) + y_range_gained[1]
}

gained_to_lost_transform <- function(x) {
  gained_span <- y_range_gained[2] - y_range_gained[1]
  lost_span   <- y_range_lost[2]   - y_range_lost[1]
  (x - y_range_gained[1]) * (lost_span / gained_span) + y_range_lost[1]
}

matrix_positions <- 0.2 + (seq_along(timepoints) - 1) * 0.9 + 0.375

gained_data <- gained_data |>
  mutate(x_pos     = match(timepoint, timepoints),
         x_aligned = matrix_positions[x_pos])

lost_data_transformed <- lost_data |>
  mutate(
    mean_score_transformed  = lost_to_gained_transform(mean_score),
    lower_ci_transformed    = lost_to_gained_transform(lower_ci),
    upper_ci_transformed    = lost_to_gained_transform(upper_ci),
    x_pos                   = match(timepoint, timepoints),
    x_aligned               = matrix_positions[x_pos]
  )

# Visualization ----
p_dual <- ggplot() +
  geom_line(data = gained_data,
            aes(x = x_aligned, y = mean_score, group = 1),
            linewidth = 1.2, color = color_gained) +
  geom_point(data = gained_data,
             aes(x = x_aligned, y = mean_score),
             size = 3, color = color_gained) +
  geom_line(data = lost_data_transformed,
            aes(x = x_aligned, y = mean_score_transformed, group = 1),
            linewidth = 1.2, color = color_lost) +
  geom_point(data = lost_data_transformed,
             aes(x = x_aligned, y = mean_score_transformed),
             size = 3, color = color_lost) +
  geom_text(data = gained_data |> filter(timepoint == timepoints[length(timepoints)]),
            aes(x = x_aligned, y = mean_score),
            label = "Gained", nudge_x = 0.3, nudge_y = -0.15,
            size = 3.5, fontface = "bold", color = color_gained) +
  geom_text(data = lost_data_transformed |> filter(timepoint == timepoints[length(timepoints)]),
            aes(x = x_aligned, y = mean_score_transformed),
            label = "Lost", nudge_x = 0.3, nudge_y = 0.01,
            size = 3.5, fontface = "bold", color = color_lost) +
  scale_x_continuous(
    breaks = matrix_positions,
    labels = timepoints,
    limits = c(min(matrix_positions) - 0.2, max(matrix_positions) + 0.5)
  ) +
  scale_y_continuous(
    name    = "APA Score",
    limits  = y_range_gained,
    breaks  = pretty(y_range_gained, n = 5),
    sec.axis = sec_axis(
      trans  = ~gained_to_lost_transform(.),
      name   = NULL,
      breaks = pretty(y_range_lost, n = 6)
    )
  ) +
  labs(x = "Time after hyperosmotic stress") +
  theme_minimal() +
  theme(
    legend.position          = "none",
    panel.grid               = element_blank(),
    panel.border             = element_blank(),
    axis.line.x              = element_line(color = "black",       linewidth = 0.5),
    axis.line.y.left         = element_line(color = color_gained,  linewidth = 0.5),
    axis.line.y.right        = element_line(color = color_lost,    linewidth = 0.5),
    axis.ticks.x             = element_line(color = "black",       linewidth = 0.5),
    axis.ticks.y.left        = element_line(color = color_gained,  linewidth = 0.5),
    axis.ticks.y.right       = element_line(color = color_lost,    linewidth = 0.5),
    axis.ticks.length.y.left  = unit(-0.15, "cm"),
    axis.ticks.length.y.right = unit(-0.15, "cm"),
    axis.text.y.left         = element_text(color = color_gained, size = 10, margin = margin(r = 5)),
    axis.text.y.right        = element_text(color = color_lost,   size = 10, margin = margin(l = 5)),
    axis.title.y.left        = element_text(color = "black", angle = 90, size = 11),
    axis.text.x              = element_text(size = 10),
    axis.title.x             = element_text(size = 11, margin = margin(t = 15)),
    plot.margin              = margin(10, 10, 10, 10)
  )

# Save outputs ----
pdf(output_pdf, width = page_width, height = page_height)
print(p_dual)
dev.off()

sessionInfo()
