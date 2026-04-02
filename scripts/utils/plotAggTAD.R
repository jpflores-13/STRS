# plotAggTAD.R ----------------------------------------------------------- ## Plot aggregated TAD heatmaps and enrichment line plots

# Libraries ----
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(scales)
library(plotgardener)

#' Plot aggregated TAD heatmap (original version)
#'
#' @param AggTAD matrix. Aggregated TAD data to visualize
#' @param maxval numeric. Maximum value for color scaling (default: 10000)
#' @param cols character. Color palette vector (default: YlGnBu 6-color)
#' @param title character. Plot title (default: "")
#' @param show_legend logical. Whether to display the color legend (default: FALSE)
#' @return ggplot. Heatmap of the aggregated TAD matrix
plotAggTAD_og <- function(AggTAD, maxval = 10000,
                       cols = RColorBrewer::brewer.pal(6, "YlGnBu"),
                       title = "", show_legend = FALSE) {
  # Convert to long format
  AggTAD_long <- setNames(reshape2::melt(AggTAD), c('x', 'y', 'counts'))

  # Create the base plot
  p <- ggplot(data = AggTAD_long, mapping = aes(x = x, y = y, fill = counts)) +
    geom_tile() +
    theme_void() +
    theme(
      aspect.ratio = 1,
      # Remove all margins
      plot.margin = margin(0, 0, 0, 0),
      # Remove legend if not needed
      sizelegend.position = if(!show_legend) "none" else "right",
      # Remove other spacing elements as needed
      panel.spacing = unit(0, "pt"),
      panel.border = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)) +
    ggtitle(title) +
    scale_fill_gradientn(colours = cols,
                         na.value = cols[maxval],
                         limits = c(0, maxval),
                         oob = scales::squish)

  return(p)
}

#' Create enrichment plot for TAD or loop analysis (original version)
#'
#' @param data_list list. Named list of matrices, each representing a condition
#' @param plot_type character. Either "tad" or "loop"
#' @param colors character. Named vector of colors for each condition (optional)
#' @param title character. Optional plot title
#' @param legend_position numeric. Legend position as c(x, y) (default: c(0.75, 0.9))
#' @param line_size numeric. Thickness of data lines (default: 0.5)
#' @return ggplot. Enrichment line plot
create_enrichment_plot_og <- function(data_list,
                                   plot_type = c("tad", "loop"),
                                   colors = NULL,
                                   title = NULL,
                                   legend_position = c(0.75, 0.9),
                                   line_size = 0.5) {

  # Input validation
  plot_type <- match.arg(plot_type)

  # Set default colors if not provided
  if (is.null(colors)) {
    colors <- c("#619CFF", "#F8766D")
    names(colors) <- names(data_list)
  }

  # Create combined dataframe for plotting
  plot_data <- data.frame()

  # Process each condition based on plot type
  if (plot_type == "tad") {
    # TAD enrichment plot parameters
    vlines <- c(26, 51)
    x_labels <- c("left", "right")
    y_limits <- c(0.9, 2.0)
    y_breaks <- seq(1, 2, 0.5)

    # Process each condition for TAD plot
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

  } else {  # loop enrichment plot
    # Loop enrichment plot parameters
    vlines <- 26
    x_labels <- "loop pixel"
    y_limits <- c(0.9, 4.25)
    y_breaks <- seq(1, 5.0, 1.0)

    # Extract matchedSet data
    mat <- data_list[["Pre-Existing"]]
    match_data <- mat[cbind(1:51, 50:100)]
    match_bg <- median(match_data[c(1:10, 40:50)])
    match_df <- data.frame(
      position = 1:51,
      enrichment = match_data/match_bg,
      condition = "Pre-Existing"
    )

    # Extract gained data
    mat <- data_list[["Gained"]]
    gain_data <- mat[cbind(1:50, 50:99)]
    gain_bg <- median(gain_data[c(1:10, 40:50)])
    gain_df <- data.frame(
      position = 1:50,
      enrichment = gain_data/gain_bg,
      condition = "Gained"
    )

    # Combine data frames
    plot_data <- rbind(match_df, gain_df)
  }

  # Create ggplot with a square aspect ratio and optimized legend
  p <- ggplot(plot_data, aes(x = position, y = enrichment, color = condition)) +
    geom_vline(xintercept = vlines, linetype = "dashed", color = "grey90") +
    geom_line(size = line_size) +
    scale_color_manual(values = colors) +
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +
    scale_x_continuous(breaks = vlines, labels = x_labels) +
    labs(x = NULL, y = NULL, color = NULL, title = title) +
    # Force square aspect ratio
    coord_fixed(ratio = 1/diff(y_limits) * diff(range(plot_data$position))) +
    theme_bw() +
    theme(
      # Remove grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      # Add border
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),

      # Regular text (not bold)
      axis.text.x = element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10),
      plot.title = element_text(hjust = 0.5, size = 12),

      # Optimize legend - smaller text, inside plot area
      sizelegend.position = legend_position,
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.8, "lines"),
      legend.margin = margin(1, 1, 1, 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.box.background = element_blank(),

      # Make plot square
      aspect.ratio = 1
    )

  return(p)
}


#' Plot aggregated TAD heatmap with zero margins
#'
#' Creates a heatmap of an aggregated TAD matrix with no extra margins,
#' suitable for precise positioning in plotgardener layouts.
#'
#' @param AggTAD matrix. Aggregated TAD data to visualize
#' @param maxval numeric. Maximum value for color scaling (default: 10000)
#' @param cols character. Color palette vector (default: YlGnBu 6-color)
#' @param title character. Plot title; omit or pass "" to suppress (default: "")
#' @param show_legend logical. Whether to display the color legend (default: FALSE)
#' @return ggplot. Zero-margin heatmap ready for plotgardener placement
plotAggTAD <- function(AggTAD, maxval = 10000,
                       cols = RColorBrewer::brewer.pal(6, "YlGnBu"),
                       title = "", show_legend = FALSE) {
  # Convert to long format
  AggTAD_long <- setNames(reshape2::melt(AggTAD), c('x', 'y', 'counts'))

  # Create the base plot
  p <- ggplot(data = AggTAD_long, mapping = aes(x = x, y = y, fill = counts)) +
    geom_tile() +
    theme_void() +
    theme(
      # Core aspect ratio setting
      aspect.ratio = 1,

      # Remove ALL outer margins with explicit 0 units
      plot.margin = margin(0, 0, 0, 0, unit = "pt"),

      # Handle title elements - remove space even when empty
      plot.title = if(title == "") element_blank() else element_text(margin = margin(0, 0, 0, 0)),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),

      # Ensure all axis elements are removed
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.ticks.length = unit(0, "pt"),
      axis.line = element_blank(),

      # Remove any remaining panel elements
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.spacing = unit(0, "pt"),
      panel.background = element_rect(fill = "transparent", color = NA),

      # Handle the plot background
      plot.background = element_rect(fill = "transparent", color = NA),

      # Completely remove legend space if not shown
      legend.position = if(!show_legend) "none" else "right",
      legend.margin = margin(0, 0, 0, 0, "pt"),
      legend.box.margin = margin(0, 0, 0, 0, "pt"),
      legend.spacing = unit(0, "pt"),
      legend.box.spacing = unit(0, "pt")
    )

  # Only add title if it's not empty
  if(title != "") {
    p <- p + ggtitle(title)
  }

  # Add the color scale with zero expansion
  p <- p + scale_fill_gradientn(
    colours = cols,
    na.value = cols[1],
    limits = c(0, maxval),
    oob = scales::squish,
    expand = c(0, 0)
  )

  return(p)
}

#' Create zero-margin enrichment plot for TAD or loop analysis
#'
#' Generates a publication-ready enrichment line plot with no extra margins
#' but visible tick marks, suitable for plotgardener placement.
#'
#' @param data_list list. Named list of matrices, each representing a condition
#' @param plot_type character. Either "tad" or "loop"
#' @param colors character. Named vector of colors for each condition (optional)
#' @param title character. Optional plot title; NULL suppresses title space
#' @param legend_position numeric. Legend position as c(x, y) (default: c(0.75, 0.9))
#' @param line_size numeric. Thickness of data lines (default: 0.5)
#' @return ggplot. Zero-margin enrichment plot ready for plotgardener placement
create_enrichment_plot <- function(data_list,
                                   plot_type = c("tad", "loop"),
                                   colors = NULL,
                                   title = NULL,
                                   legend_position = c(0.75, 0.9),
                                   line_size = 0.5) {

  # Input validation
  plot_type <- match.arg(plot_type)

  # Set default colors if not provided
  if (is.null(colors)) {
    colors <- c("#619CFF", "#F8766D")
    names(colors) <- names(data_list)
  }

  # Create combined dataframe for plotting
  plot_data <- data.frame()

  # Process each condition based on plot type
  if (plot_type == "tad") {
    # TAD enrichment plot parameters
    vlines <- c(26, 51)
    x_labels <- c("left", "right")
    y_limits <- c(0.9, 2)
    y_breaks <- seq(1, 2.5, 0.5)

    # Process each condition for TAD plot
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

  } else {  # loop enrichment plot
    # Loop enrichment plot parameters
    vlines <- 26
    x_labels <- "loop pixel"
    y_limits <- c(0.9, 4.25)
    y_breaks <- seq(1, 5.0, 1.0)

    # Extract matchedSet data
    mat <- data_list[["Pre-Existing"]]
    match_data <- mat[cbind(1:51, 50:100)]
    match_bg <- median(match_data[c(1:10, 40:50)])
    match_df <- data.frame(
      position = 1:51,
      enrichment = match_data/match_bg,
      condition = "Pre-Existing"
    )

    # Extract gained data
    mat <- data_list[["Gained"]]
    gain_data <- mat[cbind(1:50, 50:99)]
    gain_bg <- median(gain_data[c(1:10, 40:50)])
    gain_df <- data.frame(
      position = 1:50,
      enrichment = gain_data/gain_bg,
      condition = "Gained"
    )

    # Combine data frames
    plot_data <- rbind(match_df, gain_df)
  }

  # Create ggplot with proper tickmarks
  p <- ggplot(plot_data, aes(x = position, y = enrichment, color = condition)) +
    geom_vline(xintercept = vlines, linetype = "dashed", color = "grey90") +
    geom_line(size = line_size) +
    scale_color_manual(values = colors) +
    # Perfect bounds with ticks at specified locations
    scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = c(0, 0)) +
    scale_x_continuous(breaks = vlines, labels = x_labels, expand = c(0, 0),
                       limits = range(plot_data$position)) +
    labs(x = NULL, y = NULL, color = NULL, title = title) +
    # Square aspect ratio
    coord_fixed(ratio = 1/diff(y_limits) * diff(range(plot_data$position))) +
    theme_bw() +
    theme(
      # Remove ALL outer margins
      plot.margin = margin(0, 0, 0, 0, unit = "pt"),

      # Remove grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      # Keep border
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),

      # Keep ticks but remove tick labels
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),

      # Handle title with no extra space
      plot.title = if(is.null(title)) element_blank() else
        element_text(hjust = 0.5, size = 12, margin = margin(0, 0, 0, 0)),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),

      # Keep ticks visible
      axis.title = element_blank(),
      axis.ticks = element_line(linewidth = 0.3),
      axis.ticks.length = unit(3, "pt"),

      # Compact legend
      legend.position = legend_position,
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.8, "lines"),
      legend.margin = margin(0, 0, 0, 0, "pt"),
      legend.box.margin = margin(0, 0, 0, 0, "pt"),
      legend.spacing = unit(0, "pt"),
      legend.background = element_rect(fill = "white", color = NA),

      # Maintain square shape
      aspect.ratio = 1
    )

  return(p)
}

#' Add tick-mark labels to enrichment plots with precise alignment
#'
#' Adds x- and y-axis tick labels that align with the existing tick marks
#' in a plotgardener-placed enrichment plot created by \code{create_enrichment_plot()}.
#'
#' @param plot_type character. Either "tad" or "loop"
#' @param x_pos numeric. x position of the plot in inches
#' @param y_pos numeric. y position of the plot in inches
#' @param plot_width numeric. Width of the plot in inches
#' @param plot_height numeric. Height of the plot in inches
#' @param debug logical. Whether to print alignment debug information (default: FALSE)
#' @return NULL (invisibly)
add_enrichment_labels <- function(plot_type = c("tad", "loop"),
                                  x_pos,
                                  y_pos,
                                  plot_width,
                                  plot_height,
                                  debug = FALSE) {

  # Input validation
  plot_type <- match.arg(plot_type)

  # Define label parameters based on plot type
  if (plot_type == "tad") {
    # TAD boundary enrichment plot
    vlines <- c(26, 51)
    x_labels <- c("left", "right")
    y_limits <- c(0.9, 2)
    y_breaks <- c(1.0, 1.5, 2)
    data_range <- 76
  } else {
    # Loop pixel enrichment plot
    vlines <- 26
    x_labels <- "loop pixel"
    y_limits <- c(0.9, 4.25)
    y_breaks <- c(1.0, 2.0, 3.0, 4.0)
    data_range <- 51
  }

  # Calculate plot boundaries
  plot_bottom <- y_pos + plot_height
  plot_right <- x_pos + plot_width

  # Debug output if requested
  if (debug) {
    cat("Plot Type:", plot_type, "\n")
    cat("Plot Position (x, y):", x_pos, y_pos, "\n")
    cat("Plot Size (w, h):", plot_width, plot_height, "\n")
    cat("Plot Boundaries (bottom, right):", plot_bottom, plot_right, "\n")
    cat("vlines:", vlines, "\n")
    cat("data_range:", data_range, "\n")
  }

  # Account for ggplot padding in x-axis
  usable_width <- plot_width * 0.90
  x_offset <- (plot_width - usable_width) / 2

  # Add x-axis labels with precise alignment
  if (length(vlines) == 1) {
    # Single label for loop plot
    norm_x <- (vlines - min(1:data_range)) / (max(1:data_range) - min(1:data_range))
    label_x <- x_pos + x_offset + (usable_width * norm_x)

    plotText(
      label = x_labels,
      x = label_x,
      y = plot_bottom + 0.07,
      fontsize = 4,
      just = c("center", "top")
    )

  } else {
    # Multiple labels for TAD plot
    for (i in seq_along(vlines)) {
      norm_x <- (vlines[i] - min(1:data_range)) / (max(1:data_range) - min(1:data_range))
      label_x <- x_pos + x_offset + (usable_width * norm_x)

      plotText(
        label = x_labels[i],
        x = label_x,
        y = plot_bottom + 0.07,
        fontsize = 4,
        just = c("center", "top")
      )
    }
  }

  # Account for ggplot padding in y-axis
  usable_height <- plot_height * 0.90
  y_offset <- (plot_height - usable_height) / 2

  # Add y-axis labels with precise alignment
  for (y_val in y_breaks) {
    norm_y <- (y_val - y_limits[1]) / (y_limits[2] - y_limits[1])
    label_y <- plot_bottom - y_offset - (usable_height * norm_y)

    plotText(
      label = format(y_val, nsmall = 1),
      x = x_pos - 0.07,
      y = label_y,
      fontsize = 4,
      just = c("right", "center")
    )
  }

  invisible(NULL)
}

# Footer ----
message("plotAggTAD.R loaded. Available functions:")
message("  plotAggTAD(AggTAD, maxval, cols, title, show_legend)")
message("  plotAggTAD_og(AggTAD, maxval, cols, title, show_legend)  -- original version")
message("  create_enrichment_plot(data_list, plot_type, colors, title, legend_position, line_size)")
message("  create_enrichment_plot_og(...)  -- original version")
message("  add_enrichment_labels(plot_type, x_pos, y_pos, plot_width, plot_height, debug)")
message("Example: p <- plotAggTAD(myTAD, maxval = 5000)")
message("         plotGG(plot = p, x = 0.5, y = 0.5, width = 1.5, height = 1.5)")
