# ##############################################################################
# filename:    Figure2.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Figure 2 — APA heatmaps for gained and lost loops over the
#              sorbitol timecourse; aggregated TAD plots; APA enrichment
#              line plot with dual y-axes; matchRanges null set comparison
# ##############################################################################

# Libraries ----
library(InteractionSet)
library(strawr)
library(tidyverse)
library(glue)
library(nullranges)
library(data.table)
library(mariner)
library(raster)
library(reshape2)
library(RColorBrewer)
library(plotgardener)
library(scales)
library(grid)
library(gridExtra)

# Define consistent gray color for styling
gray_color <- "#666666"  # Light but readable gray

# Define pixel to inches conversion helper
pix_to_inch <- function(pixels) {
  return(pixels * 0.25)
}

# Helper function for timecourse background calculations
calculate_background_median <- function(matrix) {
  # Extract background regions (top-left and bottom-right 5x5 corners)
  top_left <- matrix[1:5, 1:5]
  bottom_right <- matrix[17:21, 17:21]
  
  # Combine values
  background_values <- c(as.vector(top_left), as.vector(bottom_right))
  
  # Remove NA values
  background_values <- background_values[!is.na(background_values)]
  
  # Calculate median or return a small positive value if all values are 0 or no values remain
  if(length(background_values) == 0 || all(background_values == 0)) {
    return(0.01)  # Small default value to avoid zrange issues
  } else {
    return(median(background_values))
  }
}

# Source helper functions
source("scripts/utils/aggregateLoops.R")
source("scripts/utils/aggregateTAD.R")
source("scripts/utils/plotAggTAD.R")

# Parameters ----
hic_rep_dir        <- "/proj/phanstiel_lab/Data/processed/YAPP/hic/hg38/220715_dietJuicerCore/output/"
hic_merge_dir      <- "/proj/phanstiel_lab/Data/processed/YAPP/hic/hg38/220716_dietJuicerMerge_condition/output/"
diff_loops_rds     <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
normalized_apa_dir <- "data/processed/hic/normalizedAPA"
enrichment_rds     <- "data/processed/hic/apa_loop_enrichment_summary_original_bg.rds"
output_pdf         <- "figures/Figure2.pdf"
page_width         <- 7.5
page_height        <- 7.0

# Data import ----

# Define timepoints for time course data
timepoints <- c("0h", "10m", "30m", "1h", "3h", "6h", "12h", "24h")

hicFiles <- list.files(hic_rep_dir,
                       full.names = TRUE,
                       recursive  = TRUE,
                       pattern    = "inter_30.hic")

merged_hicFiles <- list.files(hic_merge_dir,
                              full.names = TRUE,
                              recursive  = TRUE,
                              pattern    = "inter_30.hic")

# Get and sort time course files
gained_files <- list.files(normalized_apa_dir,
                           full.names = TRUE,
                           pattern    = "tc_gainedLoops.*normalized.rds")

lost_files <- list.files(normalized_apa_dir,
                         full.names = TRUE,
                         pattern    = "tc_lostLoops.*normalized.rds")

# Sort files by time point order
sort_files <- function(files) {
  time_order <- setNames(seq_along(timepoints), timepoints)
  files[order(sapply(files, function(x) {
    time <- str_extract(x, "0h|10m|30m|1h|3h|6h|12h|24h")
    time_order[time]
  }))]
}

gained_files <- sort_files(gained_files)
lost_files <- sort_files(lost_files)

# Load differential loops and process data
noDroso_loops <- readRDS(diff_loops_rds) |>
  interactions() |> 
  pullHicPixels(binSize = 10e3,
                files = hicFiles,
                half = "both",
                norm = "VC_SQRT",
                matrix = "observed")

# Load matrices for time course
gained_matrices <- lapply(gained_files, readRDS)
lost_matrices <- lapply(lost_files, readRDS)

# Load the enrichment data for timecourse line plot
enrichment_summary <- readRDS(enrichment_rds)

# DATA TRANSFORMATIONS AND CALCULATIONS -----------------------------------

# Create an mcol for loop size
mcols(noDroso_loops)$loop_size <- pairdist(noDroso_loops)

# Create an mcol for loop type
mcols(noDroso_loops)$loop_type <- case_when(
  mcols(noDroso_loops)$padj < 0.05 & mcols(noDroso_loops)$log2FoldChange > 1 &
    mcols(noDroso_loops)$loop_size >= 150000 ~ "truegained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange > 0 &
    mcols(noDroso_loops)$loop_size >= 150000 ~ "gained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange < 0 ~ "lost",
  mcols(noDroso_loops)$padj > 0.1 ~ "static",
  TRUE ~ "other")

# Create aggregated sorb counts mcol
mcols(noDroso_loops)$sorb_contacts <- counts(noDroso_loops)[,"YAPP_HEK_sorbitol_4_2_inter_30.hic"] +
  counts(noDroso_loops)[,"YAPP_HEK_sorbitol_5_2_inter_30.hic"] + 
  counts(noDroso_loops)[,"YAPP_HEK_sorbitol_6_2_inter_30.hic"]

# Create aggregated cont counts mcol
mcols(noDroso_loops)$cont_contacts <- counts(noDroso_loops)[,"YAPP_HEK_control_1_2_inter_30.hic"] +
  counts(noDroso_loops)[,"YAPP_HEK_control_2_2_inter_30.hic"] + 
  counts(noDroso_loops)[,"YAPP_HEK_control_3_2_inter_30.hic"]

# Create aggregated sorb+cont counts mcol
mcols(noDroso_loops)$agg_contacts <- mcols(noDroso_loops)$sorb_contacts + mcols(noDroso_loops)$cont_contacts

# Log transform loop_size to get a normal distribution for matchRanges
mcols(noDroso_loops)$loop_size <- log(mcols(noDroso_loops)$loop_size)

# Convert all `0` agg_contact values to NAs and log transform for matchRanges
mcols(noDroso_loops)$agg_contacts[mcols(noDroso_loops)$agg_contacts == 0] <- NA
noDroso_loops <- interactions(noDroso_loops) |> 
  as.data.frame() |>
  na.omit() |>
  as_ginteractions()

mcols(noDroso_loops)$agg_contacts <- log((mcols(noDroso_loops)$agg_contacts + 1))

# Create matched set for gained YAPP loops
focal <- noDroso_loops[!noDroso_loops$loop_type %in% c("static", "lost","other")] 
pool <- noDroso_loops[noDroso_loops$loop_type %in% c("static", "lost","other")] 

nullSet <- matchRanges(focal = focal,
                       pool = pool,
                       covar = ~ agg_contacts + loop_size, 
                       method = 'stratified',
                       replace = FALSE)

# Calculate sample sizes for labels
n_gained <- length(focal)
n_matched <- length(nullSet)
n_lost <- sum(noDroso_loops$loop_type == "lost")

gainedLoops <- focal
gained_bed <- gainedLoops |> 
  as.data.frame() |> 
  dplyr::select(c(1,2,3,6,7,8))

nullSet_df <- nullSet |> 
  as.data.frame() |> 
  dplyr::select(c(1,2,3,6,7,8))

# Create aggregated TADs
aggtad_gain <- aggregateTAD(loops = gained_bed,
                            hic = merged_hicFiles[2],
                            res = 10e3,
                            buffer = 0.5,
                            norm = "VC_SQRT")

aggtad_match <- aggregateTAD(loops = nullSet_df,
                             hic = merged_hicFiles[1],
                             res = 10e3,
                             buffer = 0.5,
                             norm = "VC_SQRT")

# Create data list for enrichment plots with correct keys
data_list <- list("matchedSet" = aggtad_match, "Gained" = aggtad_gain)
plot_colors <- c("matchedSet" = "#4a5568", "Gained" = "#F8766D")

# TIMECOURSE DATA PREPARATION (UPDATED FOR DUAL Y-AXES) -------------------

# Separate datasets for gained and lost loops from enrichment summary
gained_data <- enrichment_summary |> filter(loop_type == "Gained")
lost_data <- enrichment_summary |> filter(loop_type == "Lost")

# Calculate y-axis ranges based on actual data maxima
y_max_gained <- max(gained_data$mean_score, gained_data$upper_ci, na.rm = TRUE)
y_max_lost <- max(lost_data$mean_score, lost_data$upper_ci, na.rm = TRUE)

# Print actual data maxima for debugging
message("Gained data max: ", round(y_max_gained, 3))
message("Lost data max: ", round(y_max_lost, 3))

# Set specific maximum values: 2.4 for Gained, 3.0 for Lost
y_range_gained <- c(1, 2.4)
y_range_lost <- c(1, 3.0)

# Print ranges for debugging
message("Gained range: ", y_range_gained[1], " to ", y_range_gained[2])
message("Lost range: ", y_range_lost[1], " to ", y_range_lost[2])
message("Transformation needed: ", ifelse(y_range_gained[2] != y_range_lost[2], "YES", "NO"))

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

# MODIFIED PLOTTING FUNCTIONS WITH GRAY STYLING ---------------------------

# Function to plot differential (log2FC) heatmaps with divergent color palette
plotDifferentialHeatmap <- function(gained_matrix, 
                                    preexisting_matrix, 
                                    pseudocount = 1,
                                    zrange = c(-1, 1),
                                    cols = NULL,
                                    title = "",
                                    show_legend = FALSE) {
  
  # Set default Red-White-Blue divergent palette if not provided
  if (is.null(cols)) {
    # Create Red-White-Blue palette
    # Red for positive (higher in Gained), Blue for negative (higher in Pre-Existing)
    cols <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
                               "white", 
                               "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(100)
  }
  
  # Calculate log2 fold change with pseudocount
  # log2((gained + pseudocount) / (preexisting + pseudocount))
  log2fc_matrix <- log2((gained_matrix + pseudocount) / (preexisting_matrix + pseudocount))
  
  # Convert to long format for ggplot
  log2fc_long <- setNames(reshape2::melt(log2fc_matrix), c('x', 'y', 'log2fc'))
  
  # Create the plot
  p <- ggplot(data = log2fc_long, mapping = aes(x = x, y = y, fill = log2fc)) + 
    geom_tile() + 
    theme_void() + 
    theme(
      # Core aspect ratio setting
      aspect.ratio = 1,
      
      # Remove ALL outer margins
      plot.margin = margin(0, 0, 0, 0, unit = "pt"),
      
      # Handle title elements
      plot.title = if(title == "") element_blank() else element_text(margin = margin(0, 0, 0, 0)),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      
      # Remove axis elements
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.ticks.length = unit(0, "pt"),
      axis.line = element_blank(),
      
      # Remove panel elements
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.spacing = unit(0, "pt"),
      panel.background = element_rect(fill = "transparent", color = NA),
      
      # Handle plot background
      plot.background = element_rect(fill = "transparent", color = NA),
      
      # Legend settings - ultra-compact
      legend.position = if(!show_legend) "none" else "right",
      legend.title = element_text(size = 3, margin = margin(0, 0, 5, 0)),  # Increased bottom margin more
      legend.text = element_text(size = 2),
      legend.key.width = unit(0.04, "inches"),  # Ultra thin
      legend.key.height = unit(0.08, "inches"),  # Ultra short
      legend.margin = margin(0, 0, 0, 5, "pt"),  # Space to the left (pushes right)
      legend.box.margin = margin(0, 0, 0, 0, "pt"),
      legend.spacing = unit(0, "pt"),
      legend.box.spacing = unit(0, "pt"),
      legend.justification = "center"  # Center the legend vertically
    )
  
  # Add title if provided
  if(title != "") {
    p <- p + ggtitle(title)
  }
  
  # Add divergent color scale centered at 0
  p <- p + scale_fill_gradientn(
    colours = cols,
    limits = zrange,
    oob = scales::squish,  # Squish values outside range
    na.value = "gray80",
    name = "log2fc(gained/pre-existing)",  # Updated legend title
    breaks = seq(zrange[1], zrange[2], by = 0.5),  # Fewer breaks
    labels = function(x) format(x, nsmall = 1)
  )
  
  # Add coordinate system with zero expansion
  p <- p + coord_fixed(expand = FALSE)
  
  return(p)
}

# Modified create_enrichment_plot function with gray styling
create_enrichment_plot <- function(data_list, 
                                   plot_type = c("tad", "loop"),
                                   colors = NULL,
                                   title = NULL,
                                   legend_position = c(0.75, 0.9),
                                   line_size = 0.5,
                                   gray_color = "#666666") {
  
  plot_type <- match.arg(plot_type)
  
  if (is.null(colors)) {
    colors <- c("#4a5568", "#F8766D")
    names(colors) <- names(data_list)
  }
  
  plot_data <- data.frame()
  
  if (plot_type == "tad") {
    vlines <- c(26, 51)
    x_labels <- c("left", "right")
    y_limits <- c(0.9, 2)
    y_breaks <- seq(1, 2.5, 0.5)
    
    for (cond_name in names(data_list)) {
      mat <- data_list[[cond_name]]
      data_values <- mat[cbind(1:76, 25:100)]
      bg <- median(data_values[c(1:20, 56:76)])
      
      cond_df <- data.frame(
        position = 1:76,
        enrichment = data_values/bg,
        condition = cond_name
      )
      plot_data <- rbind(plot_data, cond_df)
    }
    
  } else {
    vlines <- 26
    x_labels <- "loop pixel"
    y_limits <- c(0.9, 4.25)
    y_breaks <- seq(1, 5.0, 1.0)
    
    mat <- data_list[["matchedSet"]]
    match_data <- mat[cbind(1:51, 50:100)]
    match_bg <- median(match_data[c(1:10, 40:50)])
    match_df <- data.frame(
      position = 1:51,
      enrichment = match_data/match_bg,
      condition = "matchedSet"
    )
    
    mat <- data_list[["Gained"]]
    gain_data <- mat[cbind(1:50, 50:99)]
    gain_bg <- median(gain_data[c(1:10, 40:50)])
    gain_df <- data.frame(
      position = 1:50,
      enrichment = gain_data/gain_bg,
      condition = "Gained"
    )
    
    plot_data <- rbind(match_df, gain_df)
  }
  
  p <- ggplot(plot_data, aes(x = position, y = enrichment, color = condition)) +
    geom_vline(xintercept = vlines, linetype = "dashed", color = "grey90") +
    geom_line(size = line_size) +
    scale_color_manual(values = colors) +
    scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = c(0, 0)) +
    scale_x_continuous(breaks = vlines, labels = x_labels, expand = c(0, 0),
                       limits = range(plot_data$position)) +
    labs(x = NULL, y = NULL, color = NULL, title = title) +
    coord_fixed(ratio = 1/diff(y_limits) * diff(range(plot_data$position))) +
    theme_bw() +
    theme(
      plot.margin = margin(0, 0, 0, 0, unit = "pt"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = gray_color, fill = NA, linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.3, color = gray_color),
      axis.ticks.length = unit(3, "pt"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      plot.title = if(is.null(title)) element_blank() else 
        element_text(hjust = 0.5, size = 12, margin = margin(0, 0, 0, 0)),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      legend.position = legend_position,
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.8, "lines"),
      legend.margin = margin(0, 0, 0, 0, "pt"),
      legend.box.margin = margin(0, 0, 0, 0, "pt"),
      legend.spacing = unit(0, "pt"),
      legend.background = element_rect(fill = "white", color = NA),
      aspect.ratio = 1
    )
  
  return(p)
}

# Modified add_enrichment_labels function with gray styling
add_enrichment_labels <- function(plot_type = c("tad", "loop"),
                                  x_pos,
                                  y_pos,
                                  plot_width,
                                  plot_height,
                                  gray_color = "#666666",
                                  debug = FALSE) {
  
  plot_type <- match.arg(plot_type)
  
  if (plot_type == "tad") {
    vlines <- c(26, 51)
    x_labels <- c("left", "right")
    y_limits <- c(0.9, 2)
    y_breaks <- c(1.0, 1.5, 2)
    data_range <- 76
  } else {
    vlines <- 26
    x_labels <- "loop pixel"
    y_limits <- c(0.9, 4.25)
    y_breaks <- c(1.0, 2.0, 3.0, 4.0)
    data_range <- 51
  }
  
  plot_bottom <- y_pos + plot_height
  plot_right <- x_pos + plot_width
  
  if (debug) {
    cat("Plot Type:", plot_type, "\n")
    cat("Plot Position (x, y):", x_pos, y_pos, "\n")
    cat("Plot Size (w, h):", plot_width, plot_height, "\n")
    cat("Plot Boundaries (bottom, right):", plot_bottom, plot_right, "\n")
    cat("vlines:", vlines, "\n")
    cat("data_range:", data_range, "\n")
  }
  
  usable_width <- plot_width * 0.90
  x_offset <- (plot_width - usable_width) / 2
  
  if (length(vlines) == 1) {
    norm_x <- (vlines - min(1:data_range)) / (max(1:data_range) - min(1:data_range))
    label_x <- x_pos + x_offset + (usable_width * norm_x)
    
    plotText(
      label = x_labels,
      x = label_x,
      y = plot_bottom + 0.07,
      fontsize = 4,
      fontcolor = gray_color,
      just = c("center", "top")
    )
    
  } else {
    for (i in seq_along(vlines)) {
      norm_x <- (vlines[i] - min(1:data_range)) / (max(1:data_range) - min(1:data_range))
      label_x <- x_pos + x_offset + (usable_width * norm_x)
      
      plotText(
        label = x_labels[i],
        x = label_x,
        y = plot_bottom + 0.07,
        fontsize = 4,
        fontcolor = gray_color,
        just = c("center", "top")
      )
    }
  }
  
  usable_height <- plot_height * 0.90
  y_offset <- (plot_height - usable_height) / 2
  
  for (y_val in y_breaks) {
    norm_y <- (y_val - y_limits[1]) / (y_limits[2] - y_limits[1])
    label_y <- plot_bottom - y_offset - (usable_height * norm_y)
    
    plotText(
      label = format(y_val, nsmall = 1),
      x = x_pos - 0.07,
      y = label_y,
      fontsize = 4,
      fontcolor = gray_color,
      just = c("right", "center")
    )
  }
  
  invisible(NULL)
}

# VISUALIZATION PREPARATION -----------------------------------------------

page_width <- 7.5
right_margin <- pix_to_inch(1.05)
effective_width <- page_width - right_margin

# Row 1: Gained (left) and Pre-Existing (right) domain matrices
x_pos <- pix_to_inch(1.5)
y_pos <- pix_to_inch(0.25)  # Back to original position
plot_width <- pix_to_inch(3)
plot_height <- pix_to_inch(3)
buffer_height <- pix_to_inch(1)
buffer_width <- pix_to_inch(2.5)

# Row 1: Gained (left) and Pre-Existing (right) - removing loop matrices
x_pos_preexisting <- x_pos + plot_width + buffer_width

# Row 2: Enrichment plots (domain boundaries left, peaks right)
y_pos_2 <- y_pos + plot_height + buffer_height

# Row 3: Centered differential matrix
y_pos_diff <- y_pos_2 + plot_height + buffer_height
total_space_between <- x_pos_preexisting - (x_pos + plot_width)
x_pos_centered_diff <- x_pos + plot_width + (total_space_between - plot_width) / 2

# Timecourse section starts after Pre-Existing matrix
x_pos_timecourse <- x_pos_preexisting + plot_width + buffer_width
y_pos_timecourse <- y_pos

right_section_max_width <- effective_width - x_pos_timecourse
y_axis_label_space <- 0.3
right_section_plot_width <- right_section_max_width - y_axis_label_space

num_timepoints <- length(timepoints)
min_buffer <- pix_to_inch(0.05)
total_buffer_space <- (num_timepoints - 1) * min_buffer
max_matrix_width <- (right_section_plot_width - total_buffer_space) / num_timepoints

matrix_width <- min(pix_to_inch(1.5), max_matrix_width * 0.9)
matrix_buffer <- (right_section_plot_width - (num_timepoints * matrix_width)) / (num_timepoints - 1)

matrix_positions <- numeric(num_timepoints)
for (i in 1:num_timepoints) {
  matrix_positions[i] <- x_pos_timecourse + (i-1) * (matrix_width + matrix_buffer) + (matrix_width / 2)
}

gained_data <- gained_data |>
  mutate(x_pos = match(timepoint, timepoints),
         x_aligned = matrix_positions[x_pos])

lost_data_transformed <- lost_data |>
  mutate(
    mean_score_transformed = lost_to_gained_transform(mean_score),
    lower_ci_transformed = lost_to_gained_transform(lower_ci),
    upper_ci_transformed = lost_to_gained_transform(upper_ci),
    x_pos = match(timepoint, timepoints),
    x_aligned = matrix_positions[x_pos]
  )


# VISUALIZATION AND PLOTTING ----------------------------------------------
# Increased height to accommodate differential row
pdf(output_pdf, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = F)

apaParams <- pgParams(
  assembly = "hg38",
  height = matrix_width,
  width = matrix_width,
  fontsize = 5,
  palette = colorRampPalette(brewer.pal(n = 9, "YlGnBu"))
)

# PLOT LEFT SECTION: TAD Plots --------------------------------------------

plotText(
  label = "A",
  x = x_pos - 0.3,
  y = y_pos - 0.3,
  fontsize = 12,
  fontface = "bold",
  just = c("left", "top")
)

aggTAD_gainedPlot <- plotAggTAD(
  aggtad_gain, 
  maxval = aggtad_gain[26,75]*1.2,
  title = ""
) 

plotGG(
  aggTAD_gainedPlot,
  x = x_pos,
  y = y_pos,
  width = plot_width,
  height = plot_height,
  just = c("left", "top")
)

# Add first red diagonal line (from original Panel A, Gained)
plotSegments(
  x0 = x_pos + (plot_width * 0.21),
  y0 = y_pos + (plot_height * 0.51),
  x1 = x_pos + (plot_width * 0.51),
  y1 = y_pos + (plot_height * 0.21),
  linewidth = 1,
  linecolor = "#F8766D"
)

# Label for first segment (domains)
plotText(
  label = "domains",
  x = x_pos + (plot_width * 0.21) - 0.03,
  y = y_pos + (plot_height * 0.51) + 0.03,
  fontsize = 3,
  fontcolor = "#F8766D",
  just = c("right", "top")
)

# Add second red diagonal line (from original Panel B, Gained)
plotSegments(
  x0 = x_pos + (plot_width * 0.16),
  y0 = y_pos + (plot_height * 0.43),
  x1 = x_pos + (plot_width * 0.46),
  y1 = y_pos + (plot_height * 0.13),
  linewidth = 1,
  linecolor = "#F8766D"
)

# Label for second segment (loops)
plotText(
  label = "loops",
  x = x_pos + (plot_width * 0.16) - 0.03,
  y = y_pos + (plot_height * 0.43) + 0.03,
  fontsize = 3,
  fontcolor = "#F8766D",
  just = c("right", "top")
)

# PLOT PRE-EXISTING MATRIX (Row 1, Right) -------------------------------------

# PLOT PRE-EXISTING MATRIX (Row 1, Right) -------------------------------------

aggTAD_matchedPlot <- plotAggTAD(
  aggtad_match, 
  maxval = aggtad_match[26,75]*1.2,
  title = ""
)

plotGG(
  aggTAD_matchedPlot,
  x = x_pos_preexisting,
  y = y_pos,
  width = plot_width,
  height = plot_height,
  just = c("left", "top")
)

# Add first black diagonal line (from original Panel A, Pre-Existing)
plotSegments(
  x0 = x_pos_preexisting + (plot_width * 0.21),
  y0 = y_pos + (plot_height * 0.51),
  x1 = x_pos_preexisting + (plot_width * 0.51),
  y1 = y_pos + (plot_height * 0.21),
  linewidth = 1,
  linecolor = "black"
)

# Label for first segment (domains)
plotText(
  label = "domains",
  x = x_pos_preexisting + (plot_width * 0.21) - 0.03,
  y = y_pos + (plot_height * 0.51) + 0.03,
  fontsize = 3,
  fontcolor = "black",
  just = c("right", "top")
)

# Add second black diagonal line (from original Panel B, Pre-Existing)
plotSegments(
  x0 = x_pos_preexisting + (plot_width * 0.16),
  y0 = y_pos + (plot_height * 0.43),
  x1 = x_pos_preexisting + (plot_width * 0.46),
  y1 = y_pos + (plot_height * 0.13),
  linewidth = 1,
  linecolor = "black"
)

# Label for second segment (loops)
plotText(
  label = "loops",
  x = x_pos_preexisting + (plot_width * 0.16) - 0.03,
  y = y_pos + (plot_height * 0.43) + 0.03,
  fontsize = 3,
  fontcolor = "black",
  just = c("right", "top")
)

# PLOT ENRICHMENT LINE PLOTS (Row 2) ------------------------------------------

# Position enrichment plots in Row 2 with same spacing as Row 1

tad_plot <- create_enrichment_plot(
  data_list = data_list,
  plot_type = "tad",
  colors = plot_colors,
  line_size = 0.3,
  gray_color = gray_color) +
  theme(legend.position = "none")

plotGG(
  tad_plot,
  x = x_pos,
  y = y_pos_2,
  width = plot_width,
  height = plot_height,
  just = c("left", "top")
)

add_enrichment_labels(
  plot_type = "tad",
  x_pos = x_pos,
  y_pos = y_pos_2,
  plot_width = plot_width,
  plot_height = plot_height,
  gray_color = gray_color
)

# Legend labels for domain boundaries enrichment plot (top-left inside plot)
plotText(
  label = "Gained",
  x = x_pos + 0.05,
  y = y_pos_2 + 0.02,
  fontsize = 3.5,
  fontcolor = "#F8766D",
  just = c("left", "top")
)

plotText(
  label = "Pre-Existing",
  x = x_pos + 0.05,
  y = y_pos_2 + 0.08,
  fontsize = 3.5,
  fontcolor = "#4a5568",
  just = c("left", "top")
)

# PLOT ENRICHMENT LINE PLOT FOR PEAKS (Row 2, Right) --------------------------

loop_plot <- create_enrichment_plot(
  data_list = data_list,
  plot_type = "loop",
  colors = plot_colors,
  line_size = 0.3,
  gray_color = gray_color) +
  theme(legend.position = "none")

plotGG(
  loop_plot,
  x = x_pos_preexisting,
  y = y_pos_2,
  width = plot_width,
  height = plot_height,
  just = c("left", "top")
)

add_enrichment_labels(
  plot_type = "loop", 
  x_pos = x_pos_preexisting, 
  y_pos = y_pos_2,
  plot_width = plot_width, 
  plot_height = plot_height,
  gray_color = gray_color
)

# Legend labels for peaks enrichment plot (top-left inside plot)
plotText(
  label = "Gained",
  x = x_pos_preexisting + 0.05,
  y = y_pos_2 + 0.02,
  fontsize = 3.5,
  fontcolor = "#F8766D",
  just = c("left", "top")
)

plotText(
  label = "Pre-Existing",
  x = x_pos_preexisting + 0.05,
  y = y_pos_2 + 0.08,
  fontsize = 3.5,
  fontcolor = "#4a5568",
  just = c("left", "top")
)

# PLOT CENTERED DIFFERENTIAL MATRIX (Row 3) ------------------------------------

# Create and plot centered differential matrix WITH LEGEND (no title)
differential_plot <- plotDifferentialHeatmap(
  gained_matrix = aggtad_gain,
  preexisting_matrix = aggtad_match,
  pseudocount = 1,
  zrange = c(-1, 1),
  title = "",  # No title
  show_legend = TRUE
)

# Increased width to allow matrix portion to expand to match aggTAD matrices
plotGG(
  differential_plot,
  x = x_pos_centered_diff,
  y = y_pos_diff,
  width = plot_width + 0.6,  # More total space for matrix to expand
  height = plot_height,
  just = c("left", "top")
)

# PLOT RIGHT SECTION: Hi-C Timecourse Matrices ----------------------------

plotText(
  label = "C",
  x = x_pos_timecourse - 0.3,
  y = y_pos - 0.3,
  fontsize = 12,
  fontface = "bold",
  just = c("left", "top")
)

matrix_edges_left <- numeric(num_timepoints)
matrix_edges_right <- numeric(num_timepoints)

matrix_plot <- NULL

for(i in seq_along(timepoints)) {
  x_pos_matrix <- x_pos_timecourse + (i-1) * (matrix_width + matrix_buffer)
  
  matrix_edges_left[i] <- x_pos_matrix
  matrix_edges_right[i] <- x_pos_matrix + matrix_width
  
  gainedBg <- calculate_background_median(gained_matrices[[i]])
  lostBg <- calculate_background_median(lost_matrices[[i]])
  
  arbVal <- 2.75
  
  # Get APA scores from enrichment_summary for this timepoint
  gained_apa_score <- enrichment_summary |> 
    filter(loop_type == "Gained", timepoint == timepoints[i]) |> 
    pull(mean_score)
  
  lost_apa_score <- enrichment_summary |> 
    filter(loop_type == "Lost", timepoint == timepoints[i]) |> 
    pull(mean_score)
  
  if(i == 1) {
    matrix_plot <- plotMatrix(
      gained_matrices[[i]],
      x = x_pos_matrix,
      y = y_pos_timecourse,
      zrange = (c(0, arbVal) * gainedBg),
      params = apaParams
    )
  } else {
    plotMatrix(
      gained_matrices[[i]],
      x = x_pos_matrix,
      y = y_pos_timecourse,
      zrange = (c(0, arbVal) * gainedBg),
      params = apaParams
    )
  }
  
  # Add APA score for gained loop matrix
  if(length(gained_apa_score) > 0) {
    plotText(
      label = sprintf("%.2f", gained_apa_score),
      fontcolor = "black",
      fontsize = 5,
      x = x_pos_matrix + matrix_width - 0.05,
      y = y_pos_timecourse + 0.03,
      just = c("right", "top")
    )
  }
  
  plotMatrix(
    lost_matrices[[i]],
    x = x_pos_matrix,
    y = (y_pos*2) + pix_to_inch(1.9),
    zrange = (c(0, arbVal) * lostBg),
    params = apaParams
  )
  
  # Add APA score for lost loop matrix
  if(length(lost_apa_score) > 0) {
    plotText(
      label = sprintf("%.2f", lost_apa_score),
      fontcolor = "black",
      fontsize = 5,
      x = x_pos_matrix + matrix_width - 0.05,
      y = (y_pos*2) + pix_to_inch(1.9) + 0.03,
      just = c("right", "top")
    )
  }
}

leftmost_edge <- min(matrix_edges_left)
rightmost_edge <- max(matrix_edges_right)
total_timecourse_width <- rightmost_edge - leftmost_edge

legend_buffer <- 0.2
legend_x <- rightmost_edge + legend_buffer
legend_width <- 0.05
legend_height <- matrix_width

legend_y_center <- y_pos_timecourse + ((y_pos*2) + pix_to_inch(1.9) + matrix_width - y_pos_timecourse) / 2

right_offset <- 0.05

annoHeatmapLegend(
  plot = matrix_plot,
  x = (legend_x - 0.1) + right_offset,
  y = legend_y_center - (legend_height / 2),
  width = legend_width,
  height = legend_height,
  just = c("left", "top"),
  fontsize = 0,
  border = FALSE,
  baseline = FALSE,
  orientation = "v",
  at = c(0, 2.75),
  col = colorRampPalette(brewer.pal(n = 9, "YlGnBu"))(100)
)

legend_left <- (legend_x - 0.1) + right_offset
legend_top <- legend_y_center - (legend_height / 2)
legend_bottom <- legend_top + legend_height
legend_middle_x <- legend_left + (legend_width / 2)

plotText(
  label = "2.75x",
  x = legend_middle_x,
  y = legend_top - 0.1,
  fontsize = 4,
  fontcolor = gray_color,
  just = c("center", "bottom")
)

plotText(
  label = "background",
  x = legend_middle_x,
  y = legend_top - 0.02,
  fontsize = 4,
  fontcolor = gray_color,
  just = c("center", "bottom")
)

plotText(
  label = "0",
  x = legend_middle_x,
  y = legend_bottom + 0.05,
  fontsize = 4,
  fontcolor = gray_color,
  just = c("center", "top")
)

label_x_pos <- leftmost_edge - 0.2

plotText(
  label = "Gained",
  x = label_x_pos + 0.1,
  y = y_pos_timecourse + (matrix_width / 2),
  fontsize = 6,
  fontcolor = "#F8766D",
  fontface = "bold",
  rot = 90,
  just = c("center", "center")
)

plotText(
  label = "Lost",
  x = label_x_pos + 0.1,
  y = (y_pos*2) + pix_to_inch(1.9) + (matrix_width / 2),
  fontsize = 6,
  fontcolor = "#619CFF",
  fontface = "bold",
  rot = 90,
  just = c("center", "center")
)

# PLOT PANEL D: TIMECOURSE LINE PLOT WITH DUAL Y-AXES ---------------------
# Line plot stays in original position (aligned with timecourse)
y_pos_lineplot <- y_pos_timecourse + pix_to_inch(5.05)

plotText(
  label = "D",
  x = leftmost_edge - 0.3,
  y = y_pos_lineplot - 0.3,
  fontsize = 12,
  fontface = "bold",
  just = c("left", "top")
)

timeplot <- ggplot() +
  geom_line(data = gained_data, 
            aes(x = x_aligned, y = mean_score, group = 1), 
            linewidth = 0.5, color = "#F8766D") +
  geom_point(data = gained_data, 
             aes(x = x_aligned, y = mean_score), 
             size = 0.5, color = "#F8766D") +
  geom_line(data = lost_data_transformed, 
            aes(x = x_aligned, y = mean_score_transformed, group = 1), 
            linewidth = 0.5, color = "#619CFF") +
  geom_point(data = lost_data_transformed, 
             aes(x = x_aligned, y = mean_score_transformed), 
             size = 0.5, color = "#619CFF") +
  geom_text(
    data = gained_data |> filter(timepoint == timepoints[length(timepoints)]),
    aes(x = x_aligned, y = mean_score),
    label = "Gained",
    nudge_x = -0.1,
    nudge_y = -0.15,
    size = 2,
    fontface = "bold",
    color = "#F8766D",
    check_overlap = FALSE
  ) +
  geom_text(
    data = lost_data_transformed |> filter(timepoint == timepoints[length(timepoints)]),
    aes(x = x_aligned, y = mean_score_transformed),
    label = "Lost",
    nudge_x = -0.1,
    nudge_y = -0.15,
    size = 2,
    fontface = "bold",
    color = "#619CFF",
    check_overlap = FALSE
  ) +
  scale_x_continuous(
    limits = c(leftmost_edge, rightmost_edge + 0.3),
    expand = c(0, 0),
    breaks = NULL,
    labels = NULL
  ) +
  scale_y_continuous(
    limits = y_range_gained,
    expand = c(0, 0),
    breaks = NULL,
    labels = NULL,
    sec.axis = sec_axis(
      trans = ~ .,
      breaks = NULL,
      labels = NULL
    )
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

plotGG(
  timeplot,
  x = leftmost_edge,
  y = y_pos_lineplot,
  width = total_timecourse_width + 0.3,
  height = pix_to_inch(5.05),
  just = c("left", "top")
)

# Calculate actual plot boundaries for axis placement
plot_left_edge <- leftmost_edge
plot_right_edge <- leftmost_edge + total_timecourse_width + 0.3

# PLOT LINE PLOT AXIS ELEMENTS (WITH DUAL Y-AXES) -------------------------

tick_length <- 0.05

for (i in seq_along(matrix_positions)) {
  plotSegments(
    x0 = matrix_positions[i],
    y0 = y_pos_lineplot + pix_to_inch(5.05),
    x1 = matrix_positions[i],
    y1 = y_pos_lineplot + pix_to_inch(5.05) + tick_length,
    linewidth = 0.5,
    linecolor = "black"
  )
  
  plotText(
    label = timepoints[i],
    x = matrix_positions[i],
    y = y_pos_lineplot + pix_to_inch(5.05) + tick_length + 0.08,
    fontsize = 8,
    just = c("center", "center")
  )
}

plotSegments(
  x0 = plot_left_edge,
  y0 = y_pos_lineplot + pix_to_inch(5.05),
  x1 = plot_right_edge,
  y1 = y_pos_lineplot + pix_to_inch(5.05),
  linewidth = 0.5,
  linecolor = "black"
)

plotText(
  label = "Time after hyperosmotic stress",
  x = plot_left_edge + ((plot_right_edge - plot_left_edge)/2),
  y = y_pos_lineplot + pix_to_inch(5.05) + tick_length + 0.25,
  fontsize = 8,
  just = c("center", "center")
)

plotSegments(
  x0 = plot_left_edge,
  y0 = y_pos_lineplot,
  x1 = plot_left_edge,
  y1 = y_pos_lineplot + pix_to_inch(5.05),
  linewidth = 0.5,
  linecolor = "#F8766D"
)

plotText(
  label = "APA Score",
  x = plot_left_edge - 0.3,
  y = y_pos_lineplot + pix_to_inch(2.5),
  fontsize = 6,
  fontcolor = "black",
  rot = 90,
  just = c("center", "center")
)

# Generate tick marks with even 0.2 intervals for left axis
y_ticks <- seq(from = 1.0, to = 2.4, by = 0.2)

plot_y_bottom <- y_pos_lineplot + pix_to_inch(5.05)
plot_y_top <- y_pos_lineplot

value_to_y_pos <- function(value, data_range, plot_y_bottom, plot_y_top) {
  normalized <- (value - data_range[1]) / (data_range[2] - data_range[1])
  y_pos <- plot_y_bottom - (normalized * (plot_y_bottom - plot_y_top))
  return(y_pos)
}

for (tick_value in y_ticks) {
  tick_y <- value_to_y_pos(tick_value, y_range_gained, plot_y_bottom, plot_y_top)
  
  plotSegments(
    x0 = plot_left_edge,
    y0 = tick_y,
    x1 = plot_left_edge + tick_length,
    y1 = tick_y,
    linewidth = 0.5,
    linecolor = "#F8766D"
  )
  
  plotText(
    label = format(tick_value, nsmall = 1),
    x = plot_left_edge - 0.05,
    y = tick_y,
    fontsize = 6,
    fontcolor = "#F8766D",
    just = c("right", "center")
  )
}

# Add right y-axis tick marks & labels for lost APA Score data ------------

plotSegments(
  x0 = plot_right_edge,
  y0 = y_pos_lineplot,
  x1 = plot_right_edge,
  y1 = y_pos_lineplot + pix_to_inch(5.05),
  linewidth = 0.5,
  linecolor = "#619CFF"
)

# Generate tick values for the Lost scale (independent from Gained)
y_ticks_lost <- seq(from = 1, to = y_range_lost[2], by = 0.4)

for (tick_value in y_ticks_lost) {
  # Transform the Lost tick value to the Gained scale to get the correct y position
  tick_value_gained <- lost_to_gained_transform(tick_value)
  tick_y <- value_to_y_pos(tick_value_gained, y_range_gained, plot_y_bottom, plot_y_top)
  
  plotSegments(
    x0 = plot_right_edge,
    y0 = tick_y,
    x1 = plot_right_edge - tick_length,
    y1 = tick_y,
    linewidth = 0.5,
    linecolor = "#619CFF"
  )
  
  plotText(
    label = format(tick_value, nsmall = 1),
    x = plot_right_edge + 0.05,
    y = tick_y,
    fontsize = 6,
    fontcolor = "#619CFF",
    just = c("left", "center")
  )
}

# Add labels for Row 1 matrices ---------------------------------------------

# Label for Gained matrix (Row 1, left)
plotText(
  label = "Gained",
  x = x_pos + (plot_width / 2),
  y = y_pos - 0.05,
  fontsize = 6,
  fontcolor = "#F8766D",
  just = c("center", "top")
)

# Label for Pre-Existing matrix (Row 1, right)
plotText(
  label = "Pre-Existing",
  x = x_pos_preexisting + (plot_width / 2),
  y = y_pos - 0.05,
  fontsize = 6,
  fontcolor = "#4a5568",
  just = c("center", "top")
)

# Add y-axis labels for Row 2 enrichment plots --------------------------------

# Y-axis label for domain boundaries plot (left)
plotText(
  label = "Relative",
  x = x_pos - 0.35,
  y = y_pos_2 + (plot_height / 2),
  fontsize = 5,
  fontcolor = "black",
  rot = 90,
  just = c("center", "center")
)

plotText(
  label = "Contact Frequency",
  x = x_pos - 0.25,
  y = y_pos_2 + (plot_height / 2),
  fontsize = 5,
  fontcolor = "black",
  rot = 90,
  just = c("center", "center")
)

# Y-axis label for peaks plot (right)
plotText(
  label = "Relative",
  x = x_pos_preexisting - 0.35,
  y = y_pos_2 + (plot_height / 2),
  fontsize = 5,
  fontcolor = "black",
  rot = 90,
  just = c("center", "center")
)

plotText(
  label = "Contact Frequency",
  x = x_pos_preexisting - 0.25,
  y = y_pos_2 + (plot_height / 2),
  fontsize = 5,
  fontcolor = "black",
  rot = 90,
  just = c("center", "center")
)

# Bottom labels for enrichment plots (x-axis labels)
plotText(
  label = "Aggregated Domains",
  x = x_pos + (plot_width / 2),
  y = y_pos_2 + plot_height + 0.15,
  fontsize = 5,
  fontcolor = "black",
  just = c("center", "top")
)

plotText(
  label = "Aggregated Loops",
  x = x_pos_preexisting + (plot_width / 2),
  y = y_pos_2 + plot_height + 0.15,
  fontsize = 5,
  fontcolor = "black",
  just = c("center", "top")
)

# Save outputs ----
dev.off()

sessionInfo()