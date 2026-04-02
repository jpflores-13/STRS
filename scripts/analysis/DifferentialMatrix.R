## Standalone Differential Matrix Plot
## Just the log2(Gained/Pre-Existing) matrix with legend and title

# SETUP AND LIBRARIES -----------------------------------------------------

library(InteractionSet)
library(strawr)
library(tidyverse)
library(nullranges)
library(data.table)
library(mariner)
library(raster)
library(reshape2)
library(plotgardener)
library(scales)

# Source helper functions
source("scripts/utils/aggregateTAD.R")

# FUNCTION DEFINITION -----------------------------------------------------

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
      
      # Legend settings - very short height
      legend.position = if(!show_legend) "none" else "right",
      legend.title = element_text(size = 4, margin = margin(0, 0, 1, 0)),
      legend.text = element_text(size = 3),
      legend.key.width = unit(0.08, "inches"),  # Thin width
      legend.key.height = unit(0.3, "inches"),  # Very short - ~15% of 2" matrix
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
    name = "log2fc",  # Legend title
    breaks = seq(zrange[1], zrange[2], by = 0.5),  # Fewer breaks
    labels = function(x) format(x, nsmall = 1)
  )
  
  # Add coordinate system with zero expansion
  p <- p + coord_fixed(expand = FALSE)
  
  return(p)
}

# DATA LOADING AND PROCESSING ---------------------------------------------

# Load necessary Hi-C files
hicFiles <- list.files("/proj/phanstiel_lab/Data/processed/YAPP/hic/hg38/220715_dietJuicerCore/output/",
                       full.names = TRUE,
                       recursive = TRUE,
                       pattern = "inter_30.hic")

merged_hicFiles <- list.files("/proj/phanstiel_lab/Data/processed/YAPP/hic/hg38/220716_dietJuicerMerge_condition/output/",
                              full.names = TRUE,
                              recursive = TRUE,
                              pattern = "inter_30.hic")

# Load differential loops
noDroso_loops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  as_ginteractions() |>
  pullHicPixels(binSize = 10e3,
                files = hicFiles,
                half = "both",
                norm = "VC_SQRT",
                matrix = "observed")

# Create loop size and type annotations
mcols(noDroso_loops)$loop_size <- pairdist(noDroso_loops)

mcols(noDroso_loops)$loop_type <- case_when(
  mcols(noDroso_loops)$padj < 0.05 & mcols(noDroso_loops)$log2FoldChange > 1 &
    mcols(noDroso_loops)$loop_size >= 150000 ~ "truegained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange > 0 &
    mcols(noDroso_loops)$loop_size >= 150000 ~ "gained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange < 0 ~ "lost",
  mcols(noDroso_loops)$padj > 0.1 ~ "static",
  is.character("NA") ~ "other")

# Calculate aggregate contacts
mcols(noDroso_loops)$sorb_contacts <- counts(noDroso_loops)[,"YAPP_HEK_sorbitol_4_2_inter_30.hic"] +
  counts(noDroso_loops)[,"YAPP_HEK_sorbitol_5_2_inter_30.hic"] + 
  counts(noDroso_loops)[,"YAPP_HEK_sorbitol_6_2_inter_30.hic"]

mcols(noDroso_loops)$cont_contacts <- counts(noDroso_loops)[,"YAPP_HEK_control_1_2_inter_30.hic"] +
  counts(noDroso_loops)[,"YAPP_HEK_control_2_2_inter_30.hic"] + 
  counts(noDroso_loops)[,"YAPP_HEK_control_3_2_inter_30.hic"]

mcols(noDroso_loops)$agg_contacts <- mcols(noDroso_loops)$sorb_contacts + mcols(noDroso_loops)$cont_contacts

# Log transform for matching
mcols(noDroso_loops)$loop_size <- log(mcols(noDroso_loops)$loop_size)
mcols(noDroso_loops)$agg_contacts[mcols(noDroso_loops)$agg_contacts == 0] <- NA

noDroso_loops <- interactions(noDroso_loops) |> 
  as.data.frame() |>
  na.omit() |>
  as_ginteractions()

mcols(noDroso_loops)$agg_contacts <- log((mcols(noDroso_loops)$agg_contacts + 1))

# Create matched set
focal <- noDroso_loops[!noDroso_loops$loop_type %in% c("static", "lost","other")] 
pool <- noDroso_loops[noDroso_loops$loop_type %in% c("static", "lost","other")] 

nullSet <- matchRanges(focal = focal,
                       pool = pool,
                       covar = ~ agg_contacts + loop_size, 
                       method = 'stratified',
                       replace = FALSE)

# Prepare loop sets for aggregation
gained_bed <- focal |> 
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

# PLOTTING ----------------------------------------------------------------

# Define plot dimensions for 3x4 inch page
page_width <- 3
page_height <- 3
plot_width <- 2  # 2 inches for matrix
plot_height <- 2
legend_space <- 0.5  # Space for legend

# Create PDF
pdf("plots/DifferentialMatrix.pdf",
    width = page_width,
    height = page_height)

# Create page
pageCreate(width = page_width, 
           height = page_height, 
           showGuides = FALSE)

# Create differential plot
differential_plot <- plotDifferentialHeatmap(
  gained_matrix = aggtad_gain,
  preexisting_matrix = aggtad_match,
  pseudocount = 1,
  zrange = c(-1, 1),
  title = "",
  show_legend = TRUE
)

# Center the plot on the page
x_pos <- (page_width - plot_width - legend_space) / 2
y_pos <- (page_height - plot_height) / 2 + 0.3  # Slight offset for title

plotGG(
  differential_plot,
  x = x_pos,
  y = y_pos,
  width = plot_width + legend_space,
  height = plot_height,
  just = c("left", "top")
)

# Add title above matrix - centered with matrix, not the whole plot
plotText(
  label = "log2(Gained/Pre-Existing)",
  x = x_pos + (plot_width / 2),  # Center of matrix only, not including legend space
  y = y_pos - 0.2,
  fontsize = 10,
  fontcolor = "black",
  just = c("center", "bottom")
)

dev.off()

message("Differential matrix plot saved to: plots/DifferentialMatrix.pdf")