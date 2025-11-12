## Create ggplot2 + plotgardener compatible theme
library(ggplot2)
library(plotgardener)

# Define a theme that will put axes at edge of plotting area
# note-remember to add "axis.ticks.length = unit(-.1, "cm")" to the theme of your plot
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

# Define a function that add x axes to a plot
addXaxis <- function(ggplotobject,fontsize=12)
{
  plot.x = ggplotobject
  # remove the geoms
  plot.x$layers <- NULL
  # add the axis info
  plot.x= plot.x +
    theme(
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(color = "black",size=fontsize),
    )
  return(plot.x)
}

# Define a function that add x axes to a plot
addYaxis <- function(ggplotobject,fontsize=12)
{
  plot.y = ggplotobject
  # remove the geoms
  plot.y$layers <- NULL
  # add the axis info
  plot.y= plot.y +
    theme(
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(color = "black",size=fontsize)
    )
  return(plot.y)
}


# Example -----------------------------------------------------------------
# 
# # Create the scatter plot
# plot <- ggplot(atac, aes(x = cpm, y = log2FoldChange,color = class)) +
#   geom_hline(yintercept = 0, linetype = "dashed", size = 0.1) +
#   geom_point(size = .25) +
#   scale_color_manual(values = c("static" = "grey", "lost" = "#619CFF", "gained" = "#F8766D")) +
#   scale_x_log10(labels=scales::comma) +
#   ylim(-8, 8) +
#   theme_pg +
#   theme(
#     legend.position = "none",
#     axis.ticks = element_line(color = "black", linewidth = 0.25),
#     axis.line = element_line(size = 0.25),
#     axis.ticks.length = unit(-.1, "cm"),
#     axis.line.y = element_blank()
#   )
# # add the plot to the plotgardener page
# plotGG(
#   plot = plot,
#   x = 0.25, y = 0.25,
#   width = plotwidth, height = plotwidth, just = c("left", "top")
# )
# # Add the x axis
# plotGG(
#   plot = addXaxis(plot,fontsize = 5),
#   x = 0.25, y = 0.35,
#   width = plotwidth, height = plotwidth, just = c("left", "top")
# )
# # Add the y axis
# plotGG(
#   plot = addYaxis(plot,fontsize = 5),
#   x = 0.15, y = 0.25,
#   width = plotwidth, height = plotwidth, just = c("left", "top")
# )
# # add axes labels
# plotText(label="cpm",x=1.0,y=1.95,just = "center",fontsize = 7)
# plotText(label="ATAC log2 (LPS+IFNg/control)",x=0.05,y=1.0,just = "center",rot = 90,fontsize = 7)