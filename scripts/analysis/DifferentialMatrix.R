# ##############################################################################
# filename:    DifferentialMatrix.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Aggregate TAD differential heatmap — log2(Gained/Pre-Existing)
#              contact matrix using matchRanges-matched null set; composed
#              with plotgardener
# ##############################################################################

# Libraries ----
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

# Parameters ----
hic_replicates_dir <- "/proj/phanstiel_lab/Data/processed/YAPP/hic/hg38/220715_dietJuicerCore/output/"
hic_merged_dir     <- "/proj/phanstiel_lab/Data/processed/YAPP/hic/hg38/220716_dietJuicerMerge_condition/output/"
diff_loops_rds     <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
output_pdf         <- "plots/DifferentialMatrix.pdf"
page_width         <- 3
page_height        <- 3

# Source helper functions ----
source("scripts/utils/aggregateTAD.R")

# Helper functions ----

plotDifferentialHeatmap <- function(gained_matrix,
                                    preexisting_matrix,
                                    pseudocount = 1,
                                    zrange = c(-1, 1),
                                    cols = NULL,
                                    title = "",
                                    show_legend = FALSE) {

  if (is.null(cols)) {
    cols <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                               "white",
                               "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(100)
  }

  log2fc_matrix <- log2((gained_matrix + pseudocount) / (preexisting_matrix + pseudocount))
  log2fc_long   <- setNames(reshape2::melt(log2fc_matrix), c("x", "y", "log2fc"))

  p <- ggplot(data = log2fc_long, mapping = aes(x = x, y = y, fill = log2fc)) +
    geom_tile() +
    theme_void() +
    theme(
      aspect.ratio   = 1,
      plot.margin    = margin(0, 0, 0, 0, unit = "pt"),
      plot.title     = if (title == "") element_blank() else element_text(margin = margin(0, 0, 0, 0)),
      plot.subtitle  = element_blank(),
      plot.caption   = element_blank(),
      axis.text      = element_blank(),
      axis.title     = element_blank(),
      axis.ticks     = element_blank(),
      axis.ticks.length = unit(0, "pt"),
      axis.line      = element_blank(),
      panel.grid     = element_blank(),
      panel.border   = element_blank(),
      panel.spacing  = unit(0, "pt"),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background  = element_rect(fill = "transparent", color = NA),
      legend.position  = if (!show_legend) "none" else "right",
      legend.title     = element_text(size = 4, margin = margin(0, 0, 1, 0)),
      legend.text      = element_text(size = 3),
      legend.key.width  = unit(0.08, "inches"),
      legend.key.height = unit(0.3,  "inches"),
      legend.margin     = margin(0, 0, 0, 5, "pt"),
      legend.box.margin = margin(0, 0, 0, 0, "pt"),
      legend.spacing    = unit(0, "pt"),
      legend.box.spacing = unit(0, "pt"),
      legend.justification = "center"
    )

  if (title != "") p <- p + ggtitle(title)

  p <- p + scale_fill_gradientn(
    colours  = cols,
    limits   = zrange,
    oob      = scales::squish,
    na.value = "gray80",
    name     = "log2fc",
    breaks   = seq(zrange[1], zrange[2], by = 0.5),
    labels   = function(x) format(x, nsmall = 1)
  )

  p + coord_fixed(expand = FALSE)
}

# Data import ----
hicFiles <- list.files(hic_replicates_dir, full.names = TRUE,
                       recursive = TRUE, pattern = "inter_30.hic")

merged_hicFiles <- list.files(hic_merged_dir, full.names = TRUE,
                               recursive = TRUE, pattern = "inter_30.hic")

noDroso_loops <- readRDS(diff_loops_rds) |>
  as_ginteractions() |>
  pullHicPixels(binSize = 10e3, files = hicFiles,
                half = "both", norm = "VC_SQRT", matrix = "observed")

# Analysis ----

## Annotate loop size and type ----
mcols(noDroso_loops)$loop_size <- pairdist(noDroso_loops)

mcols(noDroso_loops)$loop_type <- case_when(
  mcols(noDroso_loops)$padj < 0.05 & mcols(noDroso_loops)$log2FoldChange > 1 &
    mcols(noDroso_loops)$loop_size >= 150000 ~ "truegained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange > 0 &
    mcols(noDroso_loops)$loop_size >= 150000 ~ "gained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange < 0 ~ "lost",
  mcols(noDroso_loops)$padj > 0.1 ~ "static",
  is.character("NA") ~ "other"
)

mcols(noDroso_loops)$sorb_contacts <-
  counts(noDroso_loops)[, "YAPP_HEK_sorbitol_4_2_inter_30.hic"] +
  counts(noDroso_loops)[, "YAPP_HEK_sorbitol_5_2_inter_30.hic"] +
  counts(noDroso_loops)[, "YAPP_HEK_sorbitol_6_2_inter_30.hic"]

mcols(noDroso_loops)$cont_contacts <-
  counts(noDroso_loops)[, "YAPP_HEK_control_1_2_inter_30.hic"] +
  counts(noDroso_loops)[, "YAPP_HEK_control_2_2_inter_30.hic"] +
  counts(noDroso_loops)[, "YAPP_HEK_control_3_2_inter_30.hic"]

mcols(noDroso_loops)$agg_contacts <-
  mcols(noDroso_loops)$sorb_contacts + mcols(noDroso_loops)$cont_contacts

mcols(noDroso_loops)$loop_size <- log(mcols(noDroso_loops)$loop_size)
mcols(noDroso_loops)$agg_contacts[mcols(noDroso_loops)$agg_contacts == 0] <- NA

noDroso_loops <- interactions(noDroso_loops) |>
  as.data.frame() |>
  na.omit() |>
  as_ginteractions()

mcols(noDroso_loops)$agg_contacts <- log((mcols(noDroso_loops)$agg_contacts + 1))

## Create matched null set ----
focal   <- noDroso_loops[!noDroso_loops$loop_type %in% c("static", "lost", "other")]
pool    <- noDroso_loops[ noDroso_loops$loop_type %in% c("static", "lost", "other")]

nullSet <- matchRanges(focal = focal, pool = pool,
                       covar   = ~ agg_contacts + loop_size,
                       method  = "stratified",
                       replace = FALSE)

gained_bed  <- focal   |> as.data.frame() |> dplyr::select(c(1, 2, 3, 6, 7, 8))
nullSet_df  <- nullSet |> as.data.frame() |> dplyr::select(c(1, 2, 3, 6, 7, 8))

aggtad_gain  <- aggregateTAD(loops = gained_bed,  hic = merged_hicFiles[2],
                              res = 10e3, buffer = 0.5, norm = "VC_SQRT")
aggtad_match <- aggregateTAD(loops = nullSet_df,  hic = merged_hicFiles[1],
                              res = 10e3, buffer = 0.5, norm = "VC_SQRT")

# Save outputs ----
plot_width   <- 2
plot_height  <- 2
legend_space <- 0.5

differential_plot <- plotDifferentialHeatmap(
  gained_matrix     = aggtad_gain,
  preexisting_matrix = aggtad_match,
  pseudocount = 1,
  zrange      = c(-1, 1),
  title       = "",
  show_legend = TRUE
)

x_pos <- (page_width  - plot_width  - legend_space) / 2
y_pos <- (page_height - plot_height) / 2 + 0.3

pdf(output_pdf, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = FALSE)

plotGG(differential_plot,
       x = x_pos, y = y_pos,
       width  = plot_width + legend_space,
       height = plot_height,
       just   = c("left", "top"))

plotText(label     = "log2(Gained/Pre-Existing)",
         x         = x_pos + (plot_width / 2),
         y         = y_pos - 0.2,
         fontsize  = 10,
         fontcolor = "black",
         just      = c("center", "bottom"))

dev.off()

message("Differential matrix plot saved to: ", output_pdf)

sessionInfo()
