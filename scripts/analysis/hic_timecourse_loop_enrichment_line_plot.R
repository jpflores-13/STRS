## Hi-C Timecourse Loop Enrichment Line Plot with Dual Y-Axes
## This script creates a line plot visualizing loop enrichment calculated with the Original Background
## (Top-Left + Bottom-Right corner selection), aligned with matrix centers
## Modified to have independent y-axis scales for Gained and Lost loops

## Load required libraries
library(tidyverse)
library(scales)
library(grid)
library(gridExtra)

# Define timepoints
timepoints <- c("0h", "10m", "30m", "1h", "3h", "6h", "12h", "24h")

# Convert timepoints to numeric hours for plotting
timepoint_hours <- c(0, 10/60, 30/60, 1, 3, 6, 12, 24)
names(timepoint_hours) <- timepoints

# Load the enrichment data for Original Background selection
enrichment_summary <- readRDS("data/processed/hic/apa_loop_enrichment_summary_original_bg.rds")

# Separate datasets for gained and lost loops
gained_data <- enrichment_summary |> filter(loop_type == "Gained")
lost_data <- enrichment_summary |> filter(loop_type == "Lost")

# Calculate y-axis ranges based on actual data maxima
# Set minimum to 1, and maximum to the max value in each dataset
y_range_gained <- c(1, max(gained_data$mean_score, gained_data$upper_ci, na.rm = TRUE))
y_range_lost <- c(1, max(lost_data$mean_score, lost_data$upper_ci, na.rm = TRUE))

# Define colors for consistency
color_gained <- "#F8766D"
color_lost <- "#619CFF"

# Transformation functions to map between the two y-axis scales
# Maps values from the Lost scale to the Gained scale for plotting
lost_to_gained_transform <- function(x) {
  # Linear transformation: map [y_range_lost[1], y_range_lost[2]] to [y_range_gained[1], y_range_gained[2]]
  gained_span <- y_range_gained[2] - y_range_gained[1]
  lost_span <- y_range_lost[2] - y_range_lost[1]
  
  (x - y_range_lost[1]) * (gained_span / lost_span) + y_range_gained[1]
}

# Inverse transformation for the right y-axis labels
# Maps values from the Gained scale back to the Lost scale
gained_to_lost_transform <- function(x) {
  gained_span <- y_range_gained[2] - y_range_gained[1]
  lost_span <- y_range_lost[2] - y_range_lost[1]
  
  (x - y_range_gained[1]) * (lost_span / gained_span) + y_range_lost[1]
}

# Transform lost data for plotting on the gained scale
lost_data_transformed <- lost_data |>
  mutate(
    mean_score_transformed = lost_to_gained_transform(mean_score),
    lower_ci_transformed = lost_to_gained_transform(lower_ci),
    upper_ci_transformed = lost_to_gained_transform(upper_ci)
  )

# Calculate matrix positions
# The first matrix is positioned at x=0.2, and each subsequent matrix is 0.9 units to the right
# Each matrix has a width of 0.75, so its center is at x+0.375
matrix_positions <- numeric(length(timepoints))
for (i in seq_along(timepoints)) {
  matrix_positions[i] <- 0.2 + (i-1) * 0.9 + 0.375
}

# Adjust x-coordinate for each data point to align with matrix positions
gained_data <- gained_data |>
  mutate(x_pos = match(timepoint, timepoints),
         x_aligned = matrix_positions[x_pos])

lost_data_transformed <- lost_data_transformed |>
  mutate(x_pos = match(timepoint, timepoints),
         x_aligned = matrix_positions[x_pos])

# Create the plot with dual y-axes
p_dual <- ggplot() +
  # First layer: Gained loops (left axis)
  geom_line(data = gained_data, 
            aes(x = x_aligned, y = mean_score, group = 1), 
            linewidth = 1.2, color = color_gained) +
  geom_point(data = gained_data, 
             aes(x = x_aligned, y = mean_score), 
             size = 3, color = color_gained) +
  
  # Second layer: Lost loops (right axis) - using transformed values
  geom_line(data = lost_data_transformed, 
            aes(x = x_aligned, y = mean_score_transformed, group = 1), 
            linewidth = 1.2, color = color_lost) +
  geom_point(data = lost_data_transformed, 
             aes(x = x_aligned, y = mean_score_transformed), 
             size = 3, color = color_lost) +
  
  # Add direct labels at end of lines
  geom_text(
    data = gained_data |> filter(timepoint == timepoints[length(timepoints)]),
    aes(x = x_aligned, y = mean_score),
    label = "Gained",
    nudge_x = 0.3,
    nudge_y = -0.15,
    size = 3.5,
    fontface = "bold",
    color = color_gained,
    check_overlap = FALSE
  ) +
  geom_text(
    data = lost_data_transformed |> filter(timepoint == timepoints[length(timepoints)]),
    aes(x = x_aligned, y = mean_score_transformed),
    label = "Lost",
    nudge_x = 0.3,
    nudge_y = 0.01,
    size = 3.5,
    fontface = "bold",
    color = color_lost,
    check_overlap = FALSE
  ) +
  
  # Set x-axis scale with fixed breaks at matrix positions
  scale_x_continuous(
    breaks = matrix_positions,
    labels = timepoints,
    limits = c(min(matrix_positions) - 0.2, max(matrix_positions) + 0.5)
  ) +
  
  # Set up left y-axis (primary) for Gained data
  scale_y_continuous(
    name = "APA Score",
    limits = y_range_gained,
    breaks = pretty(y_range_gained, n = 5),  # Automatic pretty breaks
    # Add secondary axis for Lost data with proper transformation
    sec.axis = sec_axis(
      trans = ~ gained_to_lost_transform(.),  # Transform from Gained scale to Lost scale
      name = NULL,  # No label for right axis
      breaks = pretty(y_range_lost, n = 6)  # Pretty breaks for the Lost scale (increased for more labels)
    )
  ) +
  
  # Add only necessary labels
  labs(
    x = "Time after hyperosmotic stress"
  ) +
  
  # Theme with axis lines and ticks, no gridlines
  theme_minimal() +
  theme(
    # Remove legend since we're using direct labels
    legend.position = "none",
    
    # Clean up the plot by removing grid lines
    panel.grid = element_blank(),
    panel.border = element_blank(),
    
    # Style the axis lines - color-coded for each y-axis
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y.left = element_line(color = color_gained, linewidth = 0.5),
    axis.line.y.right = element_line(color = color_lost, linewidth = 0.5),
    
    # Style the axis ticks - inward pointing and color-coded
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y.left = element_line(color = color_gained, linewidth = 0.5),
    axis.ticks.y.right = element_line(color = color_lost, linewidth = 0.5),
    axis.ticks.length.y.left = unit(-0.15, "cm"),   # Negative for inward ticks
    axis.ticks.length.y.right = unit(-0.15, "cm"),  # Negative for inward ticks
    
    # Style the left y-axis text in Gained color
    axis.text.y.left = element_text(color = color_gained, size = 10, margin = margin(r = 5)),
    
    # Style the right y-axis text in Lost color
    axis.text.y.right = element_text(color = color_lost, size = 10, margin = margin(l = 5)),
    
    # Style the axis titles - left axis label in black
    axis.title.y.left = element_text(color = "black", angle = 90, size = 11),
    
    # X-axis styling
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11, margin = margin(t = 15)),
    
    # Adjust plot margins
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

# Set up the output plots
pdf("plots/hic_timecourse_loop_enrichment_line_plot.pdf", width = 8, height = 3)
print(p_dual)
dev.off()

# Alternative save using ggsave with matching dimensions
ggsave("plots/hic_timecourse_loop_enrichment_line_plot_aligned.pdf", 
       p_dual, 
       width = 8, 
       height = 3, 
       units = "in")

# Display the plot in the R environment
print(p_dual)