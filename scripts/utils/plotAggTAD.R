#' @param AggTAD Aggregated TAD object to plot
#' @param maxval integer representing the maximum value to plot (sets the top of the zrange)
#' @param cols color palette for plotting
#' @param title title of the plot
#' @param show_legend logical indicating whether to display the legend (default: FALSE)
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

#' Create a square enrichment plot with optimized legend and thinner lines
#'
#' @param data_list List of matrices to plot (named list)
#' @param plot_type Character, either "tad" or "loop" to determine plot parameters
#' @param colors Named vector of colors for each condition
#' @param title Optional plot title
#' @param legend_position Position for the legend
#' @param line_size Thickness of the data lines
#' @return A ggplot object representing the enrichment plot
create_enrichment_plot_og <- function(data_list, 
                                   plot_type = c("tad", "loop"),
                                   colors = NULL,
                                   title = NULL,
                                   legend_position = c(0.75, 0.9),
                                   line_size = 0.5) {  # Reduced line thickness
  
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
    geom_line(size = line_size) + # Thinner lines as requested
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


#' Plot Aggregated TAD Heatmap with Zero Margins
#'
#' Creates a heatmap visualization of an aggregated Topologically Associating Domain (TAD)
#' with absolutely no extra margins, allowing for precise positioning in grid-based layouts
#' like plotgardener. The function handles color scaling, removes all unnecessary space,
#' and optimizes the visualization for publication-quality figures.
#'
#' @param AggTAD Matrix containing the aggregated TAD data to visualize
#' @param maxval Numeric value representing the maximum value for color scaling (default: 10000)
#' @param cols Color palette to use for the heatmap, defaults to YlGnBu color scheme
#' @param title Optional title for the plot (default: ""). If empty, no space is allocated
#' @param show_legend Logical indicating whether to display the color legend (default: FALSE)
#'
#' @return A ggplot object with zero margins, ready for exact positioning
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
      
      # Completely remove legend space if not shown - fixed typo
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
    na.value = cols[1],  # Using first color to avoid potential index issues
    limits = c(0, maxval),
    oob = scales::squish,
    expand = c(0, 0)  # Added zero expansion for exact bounds
  )
  
  return(p)
}

#' Create Zero-Margin Enrichment Plot for TAD or Loop Analysis
#'
#' Generates a publication-ready enrichment plot for either TAD boundary or loop pixel 
#' analysis with absolutely no extra margins, but maintains visible tickmarks.
#'
#' @param data_list Named list of matrices to plot, each representing a condition
#' @param plot_type Character, either "tad" or "loop" to determine plot parameters
#' @param colors Named vector of colors for each condition (optional)
#' @param title Optional plot title, no space allocated if NULL
#' @param legend_position Position for the legend as a numeric vector c(x, y)
#' @param line_size Thickness of the data lines
#'
#' @return A ggplot object with zero margins, ready for exact positioning
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
    theme_bw() + # Use theme_bw to get borders and ticks
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
      axis.ticks = element_line(linewidth = 0.3), # Make ticks visible
      axis.ticks.length = unit(3, "pt"), # Keep ticks a useful length
      
      
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

#' Add Tickmark Labels to Enrichment Plots with Precise Alignment
#'
#' Adds labels that perfectly align with the existing tickmarks in the plot.
#'
#' @param plot_type Character, either "tad" or "loop" to determine layout
#' @param x_pos Numeric, x position of the plot in inches
#' @param y_pos Numeric, y position of the plot in inches
#' @param plot_width Numeric, width of the plot in inches
#' @param plot_height Numeric, height of the plot in inches
#' @param debug Logical, whether to print debugging information
#'
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
    data_range <- 76 # Total data points for x-axis
  } else {
    # Loop pixel enrichment plot
    vlines <- 26
    x_labels <- "loop pixel"
    y_limits <- c(0.9, 4.25)
    y_breaks <- c(1.0, 2.0, 3.0, 4.0)
    data_range <- 51 # Total data points for x-axis
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
  
  # ADJUSTED: Account for ggplot padding in x-axis
  # Calculate usable width within the plot (adjust if needed)
  usable_width <- plot_width * 0.90  # 90% of width to account for ggplot padding
  x_offset <- (plot_width - usable_width) / 2  # Offset from edge
  
  # Add x-axis labels with precise alignment
  if (length(vlines) == 1) {
    # Single label for loop plot - calculate exact position
    norm_x <- (vlines - min(1:data_range)) / (max(1:data_range) - min(1:data_range))
    label_x <- x_pos + x_offset + (usable_width * norm_x)
    
    # Add label
    plotText(
      label = x_labels,
      x = label_x,
      y = plot_bottom + 0.07,  # Increased spacing from bottom
      fontsize = 4,
      just = c("center", "top")
    )
    
  } else {
    # Multiple labels for TAD plot
    for (i in seq_along(vlines)) {
      # Calculate exact normalized position
      norm_x <- (vlines[i] - min(1:data_range)) / (max(1:data_range) - min(1:data_range))
      label_x <- x_pos + x_offset + (usable_width * norm_x)
      
      # Add label
      plotText(
        label = x_labels[i],
        x = label_x,
        y = plot_bottom + 0.07,  # Increased spacing from bottom
        fontsize = 4,
        just = c("center", "top")
      )
    }
  }
  
  # ADJUSTED: Account for ggplot padding in y-axis
  # Calculate usable height within the plot
  usable_height <- plot_height * 0.90  # 90% of height to account for ggplot padding
  y_offset <- (plot_height - usable_height) / 2  # Offset from edge
  
  # Add y-axis labels with precise alignment
  for (y_val in y_breaks) {
    # Calculate exact normalized position
    norm_y <- (y_val - y_limits[1]) / (y_limits[2] - y_limits[1])
    label_y <- plot_bottom - y_offset - (usable_height * norm_y)
    
    # Add label
    plotText(
      label = format(y_val, nsmall = 1),  # Format with one decimal place
      x = x_pos - 0.07,  # Increased spacing from left
      y = label_y,
      fontsize = 4,
      just = c("right", "center")
    )
  }
  
  # Return invisibly
  invisible(NULL)
}