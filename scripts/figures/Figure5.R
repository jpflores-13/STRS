## Figure 5 


# libraries ---------------------------------------------------------------

library(plotgardener)
library(DESeq2)
library(ComplexHeatmap)
library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(scales)
library(InteractionSet)
library(mariner)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(stringr)
library(grid)
library(rstatix)
library(cowplot)
library(forcats)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomeInfoDb)
library(S4Vectors)
library(IRanges)
library(rtracklayer)


# utils -------------------------------------------------------------------

source("scripts/utils/make_norm_matrix.R")
source("scripts/utils/ggplot2_pgTheme.R")


# color palette set-up ----------------------------------------------------

color_upregulated <- "#F8766D"
color_downregulated <- "#619CFF"
color_static <- "#999999"


# page & grid -------------------------------------------------------------

figure_width  <- 10.5
figure_height <- 11.0  

panel_buffer  <- 0.25
row_buffer    <- 0.15  
margin_left   <- 0.5
margin_top    <- 0.4
margin_right  <- 0.3  

# Three columns for survey plots
hic_panel_width <- 2.85  
legend_width <- 0.10     
survey_spacing <- 0.20   

# Top row: Panel A (promoter - square) and Panel B (bar+heatmap spans remaining width)
available_top_width <- figure_width - margin_left - margin_right - hic_panel_width - (2 * panel_buffer)
square_width <- available_top_width / 2

# Top row height (for panels A, B, C)
top_row_height <- 3.2  

# Bottom row height (for survey plots D, E, F) 
survey_plot_height <- figure_height - margin_top - top_row_height - row_buffer

# Panel B dimensions (stacked bar + heatmap) 
panelB_bar_height  <- top_row_height * 0.60   # Bar plot portion
panelB_heat_height <- top_row_height * 0.28   # Heatmap portion
panelB_total_height <- panelB_bar_height + panelB_heat_height  # Total used space
panelB_vertical_offset <- (top_row_height - panelB_total_height) * 0.4  

# Make Panel B width match its height for square aspect ratio
panelB_width <- panelB_total_height  

# Squares for panels A and C
square_side <- min(square_width, top_row_height) * 0.95  

# Column anchors (for top row)
col1_x <- margin_left
col2_x <- margin_left + square_width + panel_buffer
col3_x <- margin_left + square_width + panelB_width + (2 * panel_buffer) + 0.25  

# Survey plot column anchors (for bottom row) - align with top panels
survey_col1_x <- col1_x + (square_side - hic_panel_width) / 2  
survey_col2_x <- col2_x + (panelB_width - hic_panel_width) / 2  
survey_col3_x <- col3_x + (square_side - hic_panel_width) / 2  

# Row anchors
row1_y <- margin_top
row2_y <- margin_top + top_row_height + row_buffer


# data prep ---------------------------------------------------------------

dds <- readRDS("data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")
clustering_data <- readRDS("data/processed/rna/timecourse/output/clustering_results.rds")

diff_loopCounts <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> interactions()
diff_loopCounts <- keepStandardChromosomes(diff_loopCounts, pruning.mode = "coarse")
mcols(diff_loopCounts)$loop_size <- pairdist(diff_loopCounts)

gainedLoops <- diff_loopCounts[mcols(diff_loopCounts)$padj < 0.1 & mcols(diff_loopCounts)$log2FoldChange > 0]
lostLoops   <- diff_loopCounts[mcols(diff_loopCounts)$padj < 0.1 & mcols(diff_loopCounts)$log2FoldChange < 0]
noDroso_loops <- diff_loopCounts


# Panel B helpers (bar plot + heatmap) ------------------------------------

create_differential_gene_plot <- function(dds_object) {
  coef <- resultsNames(dds_object) |> str_subset("Time")
  numDiff <- lapply(coef, function(cf) {
    results(dds_object, name = cf) |>
      as.data.frame() |>
      filter(padj < 0.05 & (log2FoldChange > 2 | log2FoldChange < -2)) |>
      mutate(timepoint = cf)
  })
  numDiff_combined <- bind_rows(numDiff) |>
    mutate(
      timepoint = recode(timepoint,
                         "Time_1h_vs_0h"="1h","Time_3h_vs_0h"="3h","Time_6h_vs_0h"="6h",
                         "Time_9h_vs_0h"="9h","Time_12h_vs_0h"="12h","Time_24h_vs_0h"="24h"),
      regulation = case_when(
        padj < 0.05 & log2FoldChange >  2 ~ "upregulated",
        padj < 0.05 & log2FoldChange < -2 ~ "downregulated"
      )
    )
  
  ggplot(numDiff_combined |> count(timepoint, regulation),
         aes(x = timepoint, y = n, fill = fct_rev(regulation))) +
    geom_col(width = 0.75) +
    geom_text(aes(label = n), vjust = -0.7, size = 2.2) +
    labs(x="", y="# of Differential Genes", title="", fill="") +
    scale_x_discrete(limits = c("1h","3h","6h","9h","12h","24h")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme(
      axis.text.x = element_text(size = 7, color = "black"),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.title.y = element_text(size = 8),
      legend.position = "top",
      legend.text = element_text(size = 7),
      legend.title = element_blank(),
      legend.key.size = unit(0.35,"cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color="gray90", linewidth=0.3),
      plot.margin = margin(5,5,0,5,"pt"),
      panel.background = element_rect(fill="white", color=NA)
    ) +
    scale_fill_manual(values = c(color_upregulated, color_downregulated))
}

create_clustering_heatmap <- function(clustering_results) {
  ordered_matrix <- clustering_results$ordered_matrix
  row_order <- clustering_results$row_order
  km_object <- clustering_results$kmeans_object
  Heatmap(
    ordered_matrix,
    col = circlize::colorRamp2(c(-2,0,2), c("deepskyblue2","black","gold")),
    column_order = colnames(ordered_matrix),
    heatmap_legend_param = list(
      title = "", 
      at = c(-2,0,2),
      title_gp = gpar(fontsize=8), 
      labels_gp = gpar(fontsize=7),
      legend_width = unit(0.8, "cm"),
      legend_height = unit(2, "cm")
    ),
    show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = FALSE,
    column_names_gp = gpar(fontsize = 8), column_names_rot = 0, column_names_centered = TRUE,
    column_labels = c("0h","1h","3h","6h","9h","12h","24h"),
    row_split = factor(km_object$cluster[row_order], levels = c("1","2", "3")),
    row_gap = unit(2,"mm"), row_title_rot = 90, row_title_gp = gpar(fontsize=9, fontface="bold"),
    width  = unit(panelB_width * 2.54, "cm"),
    height = unit(panelB_heat_height * 2.54, "cm")
  )
}


# Panel A (promoter enrichment) -------------------------------------------

create_promoter_enrichment_plot <- function(loops_object) {
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene |> keepStandardChromosomes()
  seqlevels(loops_object) <- seqlevels(txdb); seqinfo(loops_object) <- seqinfo(txdb)
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
  
  ggplot(promoter_data,
         aes(x = loop_type, y = fraction, fill = interaction(category, loop_type))) +
    geom_col(position = "stack", width = 0.7) +
    geom_text(aes(label = round(fraction * 100, 1)),
              position = position_stack(vjust=.5), size = 2.8) +
    scale_fill_manual(values = c(
      "non_promoter_fraction.Gained" = "lightgrey",
      "non_promoter_fraction.Lost"   = "lightgrey",
      "non_promoter_fraction.Static" = "lightgrey",
      "promoter_fraction.Gained" = color_upregulated,
      "promoter_fraction.Lost"   = color_downregulated,
      "promoter_fraction.Static" = color_static
    )) +
    labs(x="Loop Type", y="Promoter Overlap (%)") +
    scale_y_continuous(labels = percent_format()) +
    theme_minimal() +
    theme(
      axis.text = element_text(size=8.5),
      axis.title = element_text(size=9.5),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(4,4,4,4,"pt")
    )
}


# Panel C (transcription/loop intersection) -------------------------------

create_intersection_plot <- function(dds_object, loops_object) {
  dge_gr <- results(dds_object, format = "GRanges")
  mcols(dge_gr) <- rowData(dds_object)
  dge_gr_prom <- promoters(dge_gr) |> keepStandardChromosomes(pruning.mode="coarse")
  seqlevelsStyle(dge_gr_prom) <- "UCSC"
  
  as_gi <- function(x) x |> as.data.frame() |> as_ginteractions()
  gained <- as_gi(loops_object) |> subset(padj <= .1 & log2FoldChange > 0)
  lost   <- as_gi(loops_object) |> subset(padj <= .1 & log2FoldChange < 0)
  static <- as_gi(loops_object) |> subset(padj > .1)
  
  gained_df <- as.data.frame(subsetByOverlaps(dge_gr_prom, gained)) %>% mutate(type="gained")
  lost_df   <- as.data.frame(subsetByOverlaps(dge_gr_prom, lost))   %>% mutate(type="lost")
  static_df <- as.data.frame(subsetByOverlaps(dge_gr_prom, static)) %>% mutate(type="static")
  
  combined <- bind_rows(gained_df, lost_df, static_df) |>
    dplyr::select(starts_with("Time_"), type, symbol, gene_id) |>
    pivot_longer(cols = starts_with("Time"), values_to="log2FoldChange", names_to="timepoint") |>
    mutate(timepoint = recode(timepoint,
                              "Time_1h_vs_0h"="1","Time_3h_vs_0h"="3","Time_6h_vs_0h"="6",
                              "Time_9h_vs_0h"="9","Time_12h_vs_0h"="12","Time_24h_vs_0h"="24"),
           timepoint = factor(timepoint, levels = c("1","3","6","9","12","24")),
           type = factor(type, levels = c("static","gained","lost"))
    )
  
  loop_colors <- c(static=color_static, gained=color_upregulated, lost=color_downregulated)
  
  make_panel <- function(loop_type, show_x=FALSE, show_xtitle=FALSE) {
    ggplot(filter(combined, type==loop_type), aes(x=timepoint, y=log2FoldChange)) +
      geom_boxplot(fill="white", color=loop_colors[loop_type], outlier.shape=NA, width=.7, linewidth=.5) +
      geom_hline(yintercept=0, color="gray30", linetype="dashed", linewidth=.4) +
      stat_summary(fun=median, geom="line", aes(group=1), color=loop_colors[loop_type], linewidth=.8) +
      stat_summary(fun=median, geom="point", color=loop_colors[loop_type], size=2, shape=18) +
      annotate("text", x=-Inf, y=Inf, label=str_to_title(loop_type), hjust=-.1, vjust=1.5,
               size=3, fontface="bold", color=loop_colors[loop_type]) +
      scale_y_continuous(breaks=c(-1,0,1), limits=c(-1.5,1.5), expand=c(0,0),
                         labels = function(x) sprintf("%.0f", x)) +
      coord_cartesian(ylim=c(-1.5,1.5), clip="off") +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill="white", color=NA),
        panel.grid.major.y = element_line(color="gray90", linewidth=.3),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = if (show_x) element_text(size=7, color="black") else element_blank(),
        axis.text.y = element_text(size=7, color="black"),
        axis.title.x = if (show_xtitle) element_text(size=8, face="bold") else element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(5,5,5,5,"pt")
      ) +
      labs(x = if (show_xtitle) "Hours after hyperosmotic stress" else "", y = "")
  }
  
  p_static <- make_panel("static", show_x=FALSE, show_xtitle=FALSE)
  p_gained <- make_panel("gained", show_x=FALSE, show_xtitle=FALSE)
  p_lost   <- make_panel("lost",   show_x=TRUE,  show_xtitle=TRUE)
  
  plot_grid(p_static, p_gained, p_lost, ncol=1, nrow=3, align="v", axis="lr",
            rel_heights = c(1,1,1.3))
}


# Survey plot helper (for loops) ------------------------------------------

create_hic_loop_visualization <- function(target_loop_index = 313) {
  target_loop <- gainedLoops[target_loop_index]
  loop_region <- GRanges(
    seqnames = seqnames(anchors(target_loop, "first")),
    ranges = IRanges(start = start(anchors(target_loop, "first")),
                     end   = end  (anchors(target_loop, "second")))
  )
  loop_region_buffed <- resize(loop_region, width = width(loop_region) + 2*100e3, fix="center")
  
  list(
    chrom      = as.character(seqnames(loop_region_buffed)),
    chromstart = start(loop_region_buffed),
    chromend   = end(loop_region_buffed),
    assembly   = "hg38",
    resolution = 10e3,
    zrange     = c(0,100),
    norm       = "SCALE",
    target_loop = target_loop,
    gained_loops = gainedLoops,
    lost_loops   = lostLoops,
    rna_colors = c("#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"),
    rna_timepoints = c("0h","1h","3h","6h","9h","12h","24h"),
    rna_files = c("STRS_HEK293_WT_cont_0h","STRS_HEK293_WT_sorb_1h","STRS_HEK293_WT_sorb_3h",
                  "STRS_HEK293_WT_sorb_6h","STRS_HEK293_WT_sorb_9h","STRS_HEK293_WT_sorb_12h",
                  "STRS_HEK293_WT_sorb_24h")
  )
}


# Survey plot helper (for ERRFI1 gene) ------------------------------------

create_hic_errfi1_visualization <- function() {
  # Get ERRFI1 coordinates
  gene_symbol <- "ERRFI1"
  gene_ids <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_symbol,
    keytype = "SYMBOL",
    columns = c("ENSEMBL", "ENTREZID")
  ) |>
    distinct(SYMBOL, .keep_all = TRUE)
  
  # Get gene coordinates from TxDb
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene |>
    keepStandardChromosomes()
  
  gene_range <- genes(txdb, filter = list(gene_id = gene_ids$ENTREZID))
  
  # Expand region by 200kb buffer (matching the loop visualization approach)
  buffer <- 200e3
  loop_region <- GRanges(
    seqnames = seqnames(gene_range),
    ranges = IRanges(start = start(gene_range) - buffer,
                     end = end(gene_range) + buffer)
  )
  
  # Find a gained loop near ERRFI1 for highlighting
  # Check which loops overlap with the ERRFI1 region
  anchor1_overlaps <- overlapsAny(anchors(gainedLoops, "first"), loop_region)
  anchor2_overlaps <- overlapsAny(anchors(gainedLoops, "second"), loop_region)
  loops_in_region <- anchor1_overlaps & anchor2_overlaps
  
  if (any(loops_in_region)) {
    # Use the first gained loop in the region
    target_loop <- gainedLoops[which(loops_in_region)[1]]
  } else {
    # If no gained loops, try any differential loop in the region
    anchor1_overlaps_all <- overlapsAny(anchors(diff_loopCounts, "first"), loop_region)
    anchor2_overlaps_all <- overlapsAny(anchors(diff_loopCounts, "second"), loop_region)
    loops_in_region_all <- anchor1_overlaps_all & anchor2_overlaps_all
    
    if (any(loops_in_region_all)) {
      target_loop <- diff_loopCounts[which(loops_in_region_all)[1]]
    } else {
      # Create dummy loop at gene body for highlighting
      target_loop <- GInteractions(
        anchor1 = GRanges(seqnames = seqnames(gene_range),
                          ranges = IRanges(start = start(gene_range), width = 10e3)),
        anchor2 = GRanges(seqnames = seqnames(gene_range),
                          ranges = IRanges(start = end(gene_range) - 10e3, width = 10e3))
      )
    }
  }
  
  # Return parameters in same format as create_hic_loop_visualization
  list(
    chrom      = as.character(seqnames(loop_region)),
    chromstart = start(loop_region),
    chromend   = end(loop_region),
    assembly   = "hg38",
    resolution = 10e3,
    zrange     = c(0, 100),
    norm       = "SCALE",
    target_loop = target_loop,
    gained_loops = gainedLoops,
    lost_loops   = lostLoops,
    rna_colors = c("#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58"),
    rna_timepoints = c("0h", "1h", "3h", "6h", "9h", "12h", "24h"),
    rna_files = c("STRS_HEK293_WT_cont_0h", "STRS_HEK293_WT_sorb_1h", "STRS_HEK293_WT_sorb_3h",
                  "STRS_HEK293_WT_sorb_6h", "STRS_HEK293_WT_sorb_9h", "STRS_HEK293_WT_sorb_12h",
                  "STRS_HEK293_WT_sorb_24h")
  )
}


# Create plots ------------------------------------------------------------
diff_gene_plot <- create_differential_gene_plot(dds)
clustering_heatmap <- create_clustering_heatmap(clustering_data)
promoter_plot <- create_promoter_enrichment_plot(noDroso_loops)
intersection_plot <- create_intersection_plot(dds, noDroso_loops)

# Survey plot parameters for each loop
hic_params_275 <- create_hic_loop_visualization(275)    # Panel D - NEW loop
hic_params_313 <- create_hic_loop_visualization(313)  # Panel E - was Panel D
hic_params_153 <- create_hic_loop_visualization(153)  # Panel F - was Panel E


# plotgardener visualization ----------------------------------------------

pdf("figures/Figure5.pdf", width = figure_width, height = figure_height)
pageCreate(width = figure_width, height = figure_height, showGuides = FALSE)

## ROW 1 (TOP)

# Panel A - Promoter enrichment
plotText("A", x = col1_x - 0.16, y = row1_y - 0.10, fontsize=12, fontface="bold",
         just=c("left","top"))
plotGG(promoter_plot,
       x = col1_x, y = row1_y,
       width = square_side, height = square_side,
       just = c("left","top"))

# Panel B - Barplot + Heatmap
plotText("B", x = col2_x - 0.12, y = row1_y - 0.10, fontsize=12, fontface="bold",
         just=c("left","top"))
plotGG(diff_gene_plot,
       x = col2_x, y = row1_y,
       width = panelB_width, height = panelB_bar_height,
       just = c("left","top"))
plotGG(plot = grid.grabExpr(draw(clustering_heatmap)),
       x = col2_x, y = row1_y + panelB_bar_height,
       width = panelB_width, height = panelB_heat_height,
       just = c("left","top"))
plotText("Time after hyperosmotic stress",
         x = col2_x + panelB_width/2, y = row1_y + panelB_total_height + 0.08,
         fontsize = 9, just = c("center","top"))

# Panel C - Transcription/Loop intersection
plotText("C", x = col3_x - 0.22, y = row1_y - 0.10, fontsize=12, fontface="bold",
         just=c("left","top"))
plotGG(intersection_plot,
       x = col3_x, y = row1_y,
       width = square_side, height = square_side,
       just = c("left","top"))
plotText("log2FoldChange(sorbitol/control)",
         x = col3_x - 0.12, y = row1_y + square_side/2,
         rot = 90, fontsize = 9, just = c("center","center"))

## ROW 2 (BOTTOM) - Survey plots 

# Helper function for creating survey plot at specified position
create_survey_plot <- function(hic_params, x_pos, y_pos, panel_label) {
  # Use the full survey_plot_height for each survey plot
  survey_total_h <- survey_plot_height - 0.05
  
  # Layout fractions
  frac_top_pad    <- 0.02
  frac_hic_each   <- 0.19
  frac_hic_gap    <- 0.02
  frac_rna_total  <- 0.34
  frac_gene       <- 0.07
  frac_genome     <- 0.04
  frac_bottom_pad <- 0.03
  
  # Convert to inches
  to_in <- function(frac) survey_total_h * frac
  hic_h        <- to_in(frac_hic_each)
  hic_gap_h    <- to_in(frac_hic_gap)
  rna_block_h  <- to_in(frac_rna_total)
  gene_h       <- to_in(frac_gene)
  genome_h     <- to_in(frac_genome)
  top_pad_h    <- to_in(frac_top_pad)
  
  # Plot label
  plotText(panel_label, x = x_pos - 0.16, y = y_pos - 0.10, 
           fontsize=12, fontface="bold", just=c("left","top"))
  
  # pgParams for this survey plot
  p <- pgParams(
    assembly   = hic_params$assembly,
    resolution = hic_params$resolution,
    chrom      = hic_params$chrom,
    chromstart = hic_params$chromstart,
    chromend   = hic_params$chromend,
    zrange     = hic_params$zrange,
    norm       = hic_params$norm,
    x = x_pos, 
    width = hic_panel_width, 
    length = hic_panel_width,
    height = hic_h,
    fontsize = 5
  )
  
  # 1) untreated Hi-C
  y0 <- y_pos + top_pad_h
  control_hic <- plotHicRectangle(
    data = "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic",
    params = p, y = y0
  )
  annoPixels(control_hic, data = hic_params$lost_loops,   shift=.5, type="arrow", col=color_downregulated)
  annoPixels(control_hic, data = hic_params$gained_loops, shift=.5, type="arrow", col=color_upregulated)
  plotText("untreated", x = x_pos, y = y0, just=c("left","top"), fontsize=8)
  
  # Add heatmap legend to the right of each survey plot
  annoHeatmapLegend(control_hic, orientation="v", fontcolor="black", digits=1,
                    x = x_pos + hic_panel_width + 0.05,
                    y = y0, width = 0.05, height = hic_h * 0.85,
                    default.units="inches", fontsize=6)
  
  # 2) +sorbitol Hi-C
  y1 <- y0 + hic_h + hic_gap_h
  sorb_hic <- plotHicRectangle(
    data = "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic",
    params = p, y = y1
  )
  annoPixels(sorb_hic, data = hic_params$lost_loops,   shift=.5, type="arrow", col=color_downregulated)
  annoPixels(sorb_hic, data = hic_params$gained_loops, shift=.5, type="arrow", col=color_upregulated)
  plotText("+ sorbitol", x = x_pos, y = y1, just=c("left","top"), fontsize=8)
  
  # 3) RNA block
  y2 <- y1 + hic_h
  rna_n <- 7
  rna_gap <- rna_block_h * 0.02
  rna_track_h <- (rna_block_h - (rna_n - 1) * rna_gap) / rna_n
  
  rna_files_all <- list.files("data/processed/rna/timecourse/output/mergeSignal/stranded/",
                              full.names=TRUE, pattern="\\.bw$")
  signalRange <- NULL
  if (length(rna_files_all) > 0) {
    signalRange <- calcSignalRange(
      data = rna_files_all,
      chrom = hic_params$chrom, chromstart = hic_params$chromstart, chromend = hic_params$chromend,
      assembly = "hg38", negData = FALSE
    )
  }
  
  for (j in seq_len(rna_n)) {
    yj <- y2 + (j-1) * (rna_track_h + rna_gap)
    fwd <- sprintf("data/processed/rna/timecourse/output/mergeSignal/stranded/%s_fwd.bw",
                   hic_params$rna_files[j])
    if (file.exists(fwd)) {
      plotSignal(data=fwd, params=p, x=x_pos, y=yj,
                 height=rna_track_h, linecolor=hic_params$rna_colors[j],
                 fill=hic_params$rna_colors[j], scale=TRUE, range=signalRange)
    }
    rev <- sprintf("data/processed/rna/timecourse/output/mergeSignal/stranded/%s_rev.bw",
                   hic_params$rna_files[j])
    if (file.exists(rev)) {
      plotSignal(data=rev, params=p, x=x_pos, y=yj,
                 height=rna_track_h, linecolor=hic_params$rna_colors[j],
                 scale=TRUE, range=signalRange)
    }
    # timepoint label
    plotText(hic_params$rna_timepoints[j], fontcolor=hic_params$rna_colors[j],
             rot=90, y=yj + rna_track_h/2, x=x_pos - 0.12, fontsize=6)
  }
  
  # 4) Genes
  y_genes <- y2 + rna_block_h + 0.01
  plotGenes(
    params = p, 
    chrom = p$chrom, 
    x = x_pos, 
    y = y_genes, 
    height = survey_total_h * 0.08
  )
  
  # 5) Genome label
  genome_label_y <- y_genes + survey_total_h * 0.08 + 0.005
  plotGenomeLabel(
    params = p, 
    x = x_pos, 
    y = genome_label_y
  )
  
  # 6) Vertical highlights across full stack
  loop <- hic_params$target_loop
  hl_y <- y0
  hl_h <- genome_label_y - hl_y + 0.2
  
  annoHighlight(
    plot = control_hic,
    fill = "lightgrey",
    alpha = 0.5,
    chrom = as.character(seqnames(anchors(loop,"first"))),
    chromstart = start(anchors(loop,"first")),
    chromend   = start(anchors(loop,"first")) + 10e3,
    x = x_pos, 
    y = hl_y, 
    width = hic_panel_width,
    height = hl_h, 
    just = c("left","top"),
    default.units = "inches"
  )
  
  annoHighlight(
    plot = control_hic,
    fill = "lightgrey",
    alpha = 0.5,
    chrom = as.character(seqnames(anchors(loop,"second"))),
    chromstart = end(anchors(loop,"second")) - 10e3,
    chromend   = end(anchors(loop,"second")),
    x = x_pos, 
    y = hl_y, 
    width = hic_panel_width,
    height = hl_h, 
    just = c("left","top"),
    default.units = "inches"
  )
}

# Panel D - Survey plot for loop #275
create_survey_plot(hic_params_275, survey_col1_x, row2_y, "D")

# Panel E - Survey plot for loop #313
create_survey_plot(hic_params_313, survey_col2_x, row2_y, "E")

# Panel F - Survey plot for loop #153
create_survey_plot(hic_params_153, survey_col3_x, row2_y, "F")

# Finish
dev.off()
