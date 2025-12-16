## Figure 4

pdf("figures/Figure4.pdf",
    width = 12,
    height = 8.5)

library(DESeq2)
library(InteractionSet)
library(glue)
library(apeglm)
library(RColorBrewer)
library(plotgardener)
library(tidyverse)
library(dbscan)
library(data.table)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(mariner)
library(plyranges)
library(ggplot2)
library(cowplot)
library(rsvg)
library(png)
library(grid)
source("scripts/utils/ggplot2_pgTheme.R")

# Load data ---------------------------------------------------------------

## Load APA data for Panel A
apa_list <- list.files("data/processed/hic/normalizedAPA", full.names = TRUE) |>
  sapply(readRDS, USE.NAMES = TRUE)
names(apa_list) <- list.files("data/processed/hic/normalizedAPA/", recursive = TRUE)

## Load differential loops for Panel B
diff_loopCounts <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions()
diff_loopCounts <- keepStandardChromosomes(diff_loopCounts, pruning.mode = "coarse")
mcols(diff_loopCounts)$loop_size <- pairdist(diff_loopCounts)

noDroso_loops <- diff_loopCounts

# Panel B helper function (promoter enrichment) ---------------------------

create_promoter_enrichment_plot <- function(loops_object) {
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene |> keepStandardChromosomes()
  seqlevels(loops_object) <- seqlevels(txdb)
  seqinfo(loops_object) <- seqinfo(txdb)
  
  mcols(loops_object)$loop_type <- dplyr::case_when(
    mcols(loops_object)$padj < 0.05 & mcols(loops_object)$log2FoldChange >  1 ~ "Gained",
    mcols(loops_object)$padj < 0.05 & mcols(loops_object)$log2FoldChange < -1 ~ "Lost",
    mcols(loops_object)$padj > 0.05 ~ "Static",
    TRUE ~ "Other"
  )
  
  genes_gr <- genes(txdb)
  promoter_regions <- promoters(genes_gr, upstream=2000, downstream=500)
  
  analyze_promoter_overlap <- function(loop_subset) {
    if (length(loop_subset) == 0) return(0)
    anchors_first  <- anchors(loop_subset, "first")
    anchors_second <- anchors(loop_subset, "second")
    all_anchors <- c(anchors_first, anchors_second)
    mean(countOverlaps(all_anchors, promoter_regions) > 0)
  }
  
  loop_categories <- c("Gained","Lost","Static")
  promoter_fractions <- purrr::map_dbl(loop_categories, function(cat) {
    loop_subset <- loops_object[mcols(loops_object)$loop_type == cat]
    analyze_promoter_overlap(loop_subset)
  })
  
  promoter_data <- tibble(
    loop_type = factor(loop_categories, levels = loop_categories),
    promoter_fraction = promoter_fractions,
    non_promoter_fraction = 1 - promoter_fractions
  ) |>
    pivot_longer(cols=c(promoter_fraction, non_promoter_fraction),
                 names_to="category", values_to="fraction") |>
    mutate(category=factor(category,
                           levels=c("non_promoter_fraction","promoter_fraction")))
  
  color_upregulated <- "#F8766D"
  color_downregulated <- "#619CFF"
  color_static <- "#999999"
  
  ggplot(promoter_data,
         aes(x = loop_type, y = fraction, fill = interaction(category, loop_type))) +
    geom_col(position = "stack", width = 0.7) +
    geom_text(aes(label = round(fraction * 100, 1)),
              position = position_stack(vjust=.5), size = 2.5) +
    scale_fill_manual(values = c(
      "non_promoter_fraction.Gained" = "lightgrey",
      "non_promoter_fraction.Lost"   = "lightgrey",
      "non_promoter_fraction.Static" = "lightgrey",
      "promoter_fraction.Gained" = color_upregulated,
      "promoter_fraction.Lost"   = color_downregulated,
      "promoter_fraction.Static" = color_static
    )) +
    labs(x="Loop Anchor Type", y="Promoter Overlap (%)") +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(
      axis.text = element_text(size=8.5),
      axis.title = element_text(size=9.5),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(0,0,0,0,"pt")
    )
}

# Panel C helper functions (motif analysis) -------------------------------

parse_homer_known_results <- function(homer_dir) {
  results_file <- file.path(homer_dir, "knownResults.txt")
  
  if (!file.exists(results_file)) {
    stop(paste("Cannot find", results_file))
  }
  
  results <- read.delim(results_file, 
                        header = TRUE, 
                        stringsAsFactors = FALSE)
  
  results$neg_log_pval <- -log10(results$P.value)
  results$motif_short <- gsub("/.*", "", results$Motif.Name)
  results$motif_short <- gsub("\\(.*", "", results$motif_short)
  
  return(results)
}

get_known_logo_files <- function(homer_dir, motif_names) {
  knownResults_dir <- file.path(homer_dir, "knownResults")
  
  logo_files <- character(length(motif_names))
  
  for (i in seq_along(motif_names)) {
    pattern <- paste0("^known", i, "\\.logo\\.svg$")
    logo_file <- list.files(knownResults_dir, pattern = pattern, full.names = TRUE)
    
    if (length(logo_file) > 0) {
      logo_files[i] <- logo_file[1]
    } else {
      logo_files[i] <- NA
    }
  }
  
  return(logo_files)
}

svg_to_raster <- function(svg_file) {
  if (is.na(svg_file) || !file.exists(svg_file)) {
    return(NULL)
  }
  
  temp_png <- tempfile(fileext = ".png")
  rsvg::rsvg_png(svg_file, temp_png, width = 1200, height = 300)
  img <- png::readPNG(temp_png)
  unlink(temp_png)
  
  return(img)
}

create_logo_column <- function(logo_files) {
  n_logos <- length(logo_files)
  
  p <- ggplot() + 
    xlim(0, 1) + 
    ylim(0.5, n_logos + 0.5) +
    theme_void()
  
  for (i in seq_along(logo_files)) {
    if (!is.na(logo_files[i]) && file.exists(logo_files[i])) {
      tryCatch({
        logo_raster <- svg_to_raster(logo_files[i])
        
        if (!is.null(logo_raster)) {
          p <- p + annotation_raster(logo_raster,
                                     xmin = 0, xmax = 1,
                                     ymin = i - 0.4, ymax = i + 0.4,
                                     interpolate = TRUE)
        }
      }, error = function(e) {
        message(paste("Could not read logo:", logo_files[i]))
      })
    }
  }
  
  return(p)
}

create_known_motif_plot <- function(homer_dir, top_n = 15) {
  results <- parse_homer_known_results(homer_dir)
  
  if (nrow(results) == 0) {
    stop("No motifs found in knownResults.txt")
  }
  
  n_available <- min(top_n, nrow(results))
  top_motifs <- results[1:n_available, ]
  
  logo_files <- get_known_logo_files(homer_dir, 1:n_available)
  logo_files_reversed <- rev(logo_files)
  
  top_motifs$motif_factor <- factor(top_motifs$motif_short, 
                                    levels = rev(top_motifs$motif_short))
  top_motifs$y_pos <- seq(n_available, 1, by = -1)
  
  # Create barplot matching original script
  bar_plot <- ggplot(top_motifs, aes(x = y_pos, y = neg_log_pval)) +
    geom_col(fill = "#3182bd") +
    geom_text(aes(y = 0, label = motif_short), 
              hjust = 0, 
              nudge_y = 0.3,
              size = 2.8,
              color = "white",
              fontface = "bold") +
    geom_text(aes(label = sprintf("%.1f", neg_log_pval)), 
              hjust = -0.2, 
              size = 2.5,
              color = "white") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), 
                       breaks = c(0, 10, 20)) +
    scale_x_continuous(expand = expansion(add = c(0.5, 0.5))) +
    labs(x = NULL, 
         y = "-log10(P-value)") +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 10, 5, 5)
    )
  
  logo_column <- create_logo_column(logo_files_reversed)
  
  # Use original rel_widths from standalone script
  combined <- plot_grid(logo_column, bar_plot, 
                        ncol = 2, 
                        rel_widths = c(1.2, 2.5),
                        align = "h",
                        axis = "tb")
  
  return(combined)
}

# Create page -------------------------------------------------------------

pageCreate(width = 12,
           height = 8.5,
           showGuides = FALSE)

# Panel A: APA Matrices ---------------------------------------------------

## Panel A positioning - enlarged matrices, moved down to align with other panels
apa_height <- 0.85
apa_x_col1 <- 0.5
apa_x_col2 <- apa_x_col1 + apa_height + 0.25
apa_x_col3 <- apa_x_col2 + apa_height + 0.25
apa_y_row1 <- 1.0  # Start lower to align with B and C
apa_y_row2 <- apa_y_row1 + apa_height + 0.15

## Calculate total height of Panel A for matching B and C
total_apa_height <- (apa_y_row2 - apa_y_row1) + apa_height

plotText(label="A", 
         x=apa_x_col1 - 0.15, 
         y=apa_y_row1 - 0.3, 
         fontsize=12, 
         fontface="bold")

## Define APA parameters
apaParams <- pgParams(
  assembly="hg38", 
  height=apa_height, 
  width=apa_height,
  fontsize=5, 
  palette=colorRampPalette(brewer.pal(9,"YlGnBu"))
)

## Calculate shared z-ranges for each column
## Column 1: Use separate z-ranges for top and bottom due to different signal intensities
col1_top_max <- max(apa_list$oe_gainedLoops_sorbAPA_normalized.rds, na.rm = TRUE)
col1_bottom_max <- max(apa_list$dtad_gainedLoops_sorbAPA_normalized.rds, na.rm = TRUE)

## Columns 2 and 3: Use combined max of both matrices for shared z-range
ctcf_apa_combined <- c(
  apa_list$ctcf_gainedLoops_sorbAPA_normalized.rds,
  apa_list$ctcf_gainedLoops_sorb_auxAPA_normalized.rds
)

rad21_apa_combined <- c(
  apa_list$rad21_gainedLoops_sorbAPA_normalized.rds,
  apa_list$rad21_gainedLoops_sorb_auxAPA_normalized.rds
)

col2_max <- max(ctcf_apa_combined, na.rm = TRUE)
col3_max <- max(rad21_apa_combined, na.rm = TRUE)

## Adjust z-ranges for columns 2 and 3 to increase contrast
col2_max_adjusted <- col2_max * 1.35
col3_max_adjusted <- col3_max * 1.35

legend_width <- 0.03
legend_spacing <- 0.05
legend_right_shift <- 0.10

## Column 1 TOP: HEK293T eGFP-YAP + sorbitol (own z-range)
plotMatrix(
  apa_list$oe_gainedLoops_sorbAPA_normalized.rds, 
  x=apa_x_col1, 
  y=apa_y_row1,
  zrange=c(0, col1_top_max), 
  params=apaParams
) |>
  annoHeatmapLegend(
    x=apa_x_col1 + apa_height + legend_spacing,
    y=apa_y_row1, 
    width=legend_width, 
    height=apa_height, 
    fontsize=5
  )

## Column 1 BOTTOM: HEK293T eGFP-YAPdTAD + sorbitol (own z-range)
plotMatrix(
  apa_list$dtad_gainedLoops_sorbAPA_normalized.rds, 
  x=apa_x_col1, 
  y=apa_y_row2,
  zrange=c(0, col1_bottom_max), 
  params=apaParams
) |>
  annoHeatmapLegend(
    x=apa_x_col1 + apa_height + legend_spacing,
    y=apa_y_row2, 
    width=legend_width, 
    height=apa_height, 
    fontsize=5
  )

## Column 2: mAID2-CTCF HCT116
plotMatrix(
  apa_list$ctcf_gainedLoops_sorbAPA_normalized.rds, 
  x=apa_x_col2, 
  y=apa_y_row1,
  zrange=c(0, col2_max_adjusted), 
  params=apaParams
) |>
  annoHeatmapLegend(
    x=apa_x_col2 + apa_height + legend_spacing,
    y=apa_y_row1, 
    width=legend_width, 
    height=apa_height, 
    fontsize=5
  )

plotMatrix(
  apa_list$ctcf_gainedLoops_sorb_auxAPA_normalized.rds, 
  x=apa_x_col2, 
  y=apa_y_row2,
  zrange=c(0, col2_max_adjusted), 
  params=apaParams
)

## Column 3: mAID2-RAD21 HCT116
plotMatrix(
  apa_list$rad21_gainedLoops_sorbAPA_normalized.rds, 
  x=apa_x_col3, 
  y=apa_y_row1,
  zrange=c(0, col3_max_adjusted), 
  params=apaParams
) |>
  annoHeatmapLegend(
    x=apa_x_col3 + apa_height + legend_spacing,
    y=apa_y_row1, 
    width=legend_width, 
    height=apa_height, 
    fontsize=5
  )

plotMatrix(
  apa_list$rad21_gainedLoops_sorb_auxAPA_normalized.rds, 
  x=apa_x_col3, 
  y=apa_y_row2,
  zrange=c(0, col3_max_adjusted), 
  params=apaParams
)

## Add column labels below matrices
label_y <- apa_y_row2 + apa_height + 0.1

plotText(
  "HEK293T", 
  x=apa_x_col1 + (apa_height/2), 
  y=label_y, 
  fontsize=6,
  just=c("center","top")
)

plotText(
  "HCT116", 
  x=apa_x_col2 + (apa_height/2), 
  y=label_y, 
  fontsize=6,
  just=c("center","top")
)

plotText(
  "HCT116", 
  x=apa_x_col3 + (apa_height/2), 
  y=label_y, 
  fontsize=6,
  just=c("center","top")
)

## Add condition labels in upper left corners of each matrix
## Column 1 labels - simplified
plotText(
  "eGFP-YAP", 
  x=apa_x_col1 + 0.03, 
  y=apa_y_row1 + 0.03, 
  fontsize=5,
  just=c("left","top"),
  fontcolor="black"
)

plotText(
  "eGFP-YAPdTAD", 
  x=apa_x_col1 + 0.03, 
  y=apa_y_row2 + 0.03, 
  fontsize=5,
  just=c("left","top"),
  fontcolor="black"
)

## Column 2 labels
plotText(
  "+ CTCF", 
  x=apa_x_col2 + 0.03, 
  y=apa_y_row1 + 0.03, 
  fontsize=5,
  just=c("left","top"),
  fontcolor="black"
)

plotText(
  "- CTCF", 
  x=apa_x_col2 + 0.03, 
  y=apa_y_row2 + 0.03, 
  fontsize=5,
  just=c("left","top"),
  fontcolor="black"
)

## Column 3 labels
plotText(
  "+ RAD21", 
  x=apa_x_col3 + 0.03, 
  y=apa_y_row1 + 0.03, 
  fontsize=5,
  just=c("left","top"),
  fontcolor="black"
)

plotText(
  "- RAD21", 
  x=apa_x_col3 + 0.03, 
  y=apa_y_row2 + 0.03, 
  fontsize=5,
  just=c("left","top"),
  fontcolor="black"
)

## Add "gained loops" line and label (like in Figure 1) - text above line
line_y <- apa_y_row1 - 0.05
line_x_start <- apa_x_col1
line_x_end <- apa_x_col3 + apa_height

plotText(
  label = "gained loops",
  x = (line_x_start + line_x_end) / 2,
  y = line_y - 0.12,
  fontcolor = "#F8766D",
  fontsize = 10,
  just = c("center", "bottom")
)

plotSegments(
  x0 = line_x_start,
  y0 = line_y,
  x1 = line_x_end,
  y1 = line_y,
  default.units = "inches",
  linecolor = "#F8766D",
  lwd = 1,
  lty = 1
)

# Panel B: Promoter Enrichment --------------------------------------------

## Position to the right of Panel A, aligned at top and bottom with A's labels
panel_b_x <- apa_x_col3 + apa_height + 0.5
panel_b_y <- apa_y_row1  # Align with top of Panel A
label_y <- apa_y_row2 + apa_height + 0.1  # This is where cell type labels are
panel_b_width <- 2.2  # Wider to avoid squishing
panel_b_height <- label_y - panel_b_y + 0.2  # Extend to align with bottom labels

plotText(label="B", 
         x=panel_b_x - 0.15, 
         y=panel_b_y - 0.3, 
         fontsize=12, 
         fontface="bold")

promoter_plot <- create_promoter_enrichment_plot(noDroso_loops)

plotGG(promoter_plot,
       x = panel_b_x,
       y = panel_b_y,
       width = panel_b_width,
       height = panel_b_height,
       just = c("left","top"))

# Panel C: Motif Analysis -------------------------------------------------

## Position to the right of Panel B
panel_c_x <- panel_b_x + panel_b_width + 0.3

## Calculate title position to align with "gained loops" text
line_y_ref <- apa_y_row1 - 0.05  # Same as the gained loops line
title_y <- line_y_ref - 0.07  # Same offset as "gained loops" text

## Start plot right below where matrices start in Panel A
panel_c_y <- apa_y_row1
panel_c_width <- 5.0
## Extend down to match Panel A's label position
panel_c_height <- (label_y + 0.3) - panel_c_y

## Add Panel C label aligned with A and B
plotText(label="C", 
         x=panel_c_x - 0.15, 
         y=apa_y_row1 - 0.3,  # Same as A and B
         fontsize=12, 
         fontface="bold")

## Add title using plotText (aligned with "gained loops")
plotText(
  label = "Motifs at Promoters of Gained Loop Anchors",
  x = panel_c_x + panel_c_width/2,
  y = title_y,
  fontsize = 10,
  fontface = "bold",
  just = c("center", "bottom")
)

homer_dir <- "data/processed/cutntag/homer_motifs/gained_anchor_promoters"

if (file.exists(file.path(homer_dir, "knownResults.txt"))) {
  motif_plot <- create_known_motif_plot(
    homer_dir = homer_dir,
    top_n = 12
  )
  
  plotGG(motif_plot,
         x = panel_c_x,
         y = panel_c_y,
         width = panel_c_width,
         height = panel_c_height,
         just = c("left","top"))
} else {
  plotText(
    label = "HOMER motif results not found",
    x = panel_c_x + panel_c_width/2,
    y = panel_c_y + panel_c_height/2,
    fontsize = 10,
    just = c("center","center"),
    fontcolor = "red"
  )
}

dev.off()