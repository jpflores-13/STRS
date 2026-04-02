# ggplot2_pgTheme.R ------------------------------------------------------ ## ggplot2 theme and axis helpers compatible with plotgardener

# Libraries ----
library(ggplot2)
library(plotgardener)

# theme_pg: ggplot2 theme with axes at the edge of the plotting area ----
# Note: add axis.ticks.length = unit(-.1, "cm") to the theme of your plot
theme_pg <- theme_classic() + theme(
  axis.line = element_line(color = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  plot.background = element_blank(),
  axis.text = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_line(color = "black"),
  plot.margin = margin(0, 0, 0, 0, "cm")
)

#' Add x-axis labels to a plotgardener-compatible ggplot
#'
#' Strips geom layers from the plot and returns a version showing only the
#' x-axis tick labels, suitable for overlaying via plotGG().
#'
#' @param ggplotobject ggplot. A ggplot object built with theme_pg
#' @param fontsize numeric. Font size for axis text (default: 12)
#' @return ggplot. The input plot with all geoms removed and x-axis text visible
addXaxis <- function(ggplotobject, fontsize = 12) {
  plot.x <- ggplotobject
  # remove the geoms
  plot.x$layers <- NULL
  # add the axis info
  plot.x <- plot.x +
    theme(
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(color = "black", size = fontsize),
    )
  return(plot.x)
}

#' Add y-axis labels to a plotgardener-compatible ggplot
#'
#' Strips geom layers from the plot and returns a version showing only the
#' y-axis tick labels, suitable for overlaying via plotGG().
#'
#' @param ggplotobject ggplot. A ggplot object built with theme_pg
#' @param fontsize numeric. Font size for axis text (default: 12)
#' @return ggplot. The input plot with all geoms removed and y-axis text visible
addYaxis <- function(ggplotobject, fontsize = 12) {
  plot.y <- ggplotobject
  # remove the geoms
  plot.y$layers <- NULL
  # add the axis info
  plot.y <- plot.y +
    theme(
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(color = "black", size = fontsize)
    )
  return(plot.y)
}

# Footer ----
message("ggplot2_pgTheme.R loaded. Available objects and functions:")
message("  theme_pg          -- ggplot2 theme with axes at plot edge")
message("  addXaxis(ggplotobject, fontsize)  -- overlay x-axis labels via plotGG()")
message("  addYaxis(ggplotobject, fontsize)  -- overlay y-axis labels via plotGG()")
message("Example:")
message("  p <- ggplot(...) + theme_pg + theme(axis.ticks.length = unit(-.1, 'cm'))")
message("  plotGG(plot = p, x = 0.25, y = 0.25, width = 1, height = 1)")
message("  plotGG(plot = addXaxis(p, fontsize = 5), x = 0.25, y = 0.35, width = 1, height = 1)")
