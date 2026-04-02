# ##############################################################################
# filename:    Figure5.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Figure 5 — RNA-seq timecourse analysis; differential gene
#              expression, clustering heatmap, GO enrichment, EISA line plot,
#              and Hi-C survey at gained loop anchor with RNA-seq tracks
# ##############################################################################

# Libraries ----
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
gray_color <- "#666666"  # Matching Figure 3's gray color
h3k27ac_color <- "#807DBA"  # H3K27ac track color (purple)


# Parameters ----
dds_rds           <- "data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds"
clustering_rds    <- "data/processed/rna/timecourse/output/clustering_results.rds"
diff_loops_rds    <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
output_pdf        <- "figures/Figure5.pdf"
figure_width      <- 10.5
figure_height     <- 11.6  

panel_buffer  <- 0.25
row_buffer    <- 0.20  
margin_left   <- 0.5
margin_top    <- 0.4
margin_right  <- 0.3  

# Layout dimensions
# Top row: Panel A (barplot+heatmap) and Panel B (GO enrichment)
top_row_height <- 3.0

# Panel A dimensions (barplot + heatmap stacked)
panelA_width <- 2.5
panelA_heatmap_width <- 2.2  # Narrower heatmap width to avoid Panel B
panelA_bar_height <- top_row_height * 0.65  # Increased from 0.60
panelA_heat_height <- top_row_height * 0.30  # Increased from 0.26 to elongate downward
panelA_total_height <- panelA_bar_height + panelA_heat_height

# Panel B dimensions (GO enrichment) - ADJUSTED for single panel
panelB_width <- 3.5
panelB_height <- top_row_height  # Use full height to align with Panel A

# Panel C dimensions (combined line plot - spans width of A+B)
panelC_width <- panelA_width + panel_buffer + panelB_width
panelC_height <- 2.6  # Increased from 2.2 to make plot taller

# Panel D dimensions (survey plot - aligns with Panels A-C)
survey_width <- 2.85
survey_height <- top_row_height + row_buffer + panelC_height

# Calculate positions
col1_x <- margin_left + 0.1  # Panel A (shifted right)
col2_x <- margin_left + 0.1 + panelA_width + panel_buffer + 0.15  # Panel B (shifted right)
col3_x <- margin_left + panelC_width + panel_buffer + 0.3  # Panel D (survey)

row1_y <- margin_top + 0.1  # Panels A & B - moved down to align labels
row2_y <- margin_top + top_row_height + row_buffer - 0.05  # Panel C - moved up slightly
row3_y <- margin_top + 0.1  # Panel D - moved down to align with A and B

# Update Panel C position to match shifted A and B
panelC_x <- col1_x  # Panel C starts at same x as Panel A


# data prep ---------------------------------------------------------------

dds             <- readRDS(dds_rds)
clustering_data <- readRDS(clustering_rds)

diff_loopCounts <- readRDS(diff_loops_rds) |>
  interactions()
diff_loopCounts <- keepStandardChromosomes(diff_loopCounts, pruning.mode = "coarse")
mcols(diff_loopCounts)$loop_size <- pairdist(diff_loopCounts)

gainedLoops <- diff_loopCounts[mcols(diff_loopCounts)$padj < 0.1 & 
                                 mcols(diff_loopCounts)$log2FoldChange > 0]
lostLoops   <- diff_loopCounts[mcols(diff_loopCounts)$padj < 0.1 & 
                                 mcols(diff_loopCounts)$log2FoldChange < 0]
noDroso_loops <- diff_loopCounts

# Load EISA results for Panel C
eisa_results <- readRDS("data/processed/rna/timecourse/output/EISA/eisa_results_with_loops.rds")


# Panel A: Barplot + Heatmap (was Panel B) -------------------------------

create_differential_gene_plot <- function(dds_object) {
  coef <- resultsNames(dds_object) |> 
    str_subset("Time")
  
  numDiff <- lapply(coef, function(cf) {
    results(dds_object, name = cf) |>
      as.data.frame() |>
      filter(padj < 0.05 & (log2FoldChange > 2 | log2FoldChange < -2)) |>
      mutate(timepoint = cf)
  })
  
  numDiff_combined <- bind_rows(numDiff) |>
    mutate(
      timepoint = recode(timepoint,
                         "Time_1h_vs_0h"="1h", "Time_3h_vs_0h"="3h", "Time_6h_vs_0h"="6h",
                         "Time_9h_vs_0h"="9h", "Time_12h_vs_0h"="12h", "Time_24h_vs_0h"="24h"),
      regulation = case_when(
        padj < 0.05 & log2FoldChange >  2 ~ "upregulated",
        padj < 0.05 & log2FoldChange < -2 ~ "downregulated"
      )
    )
  
  ggplot(numDiff_combined |> count(timepoint, regulation),
         aes(x = timepoint, y = n, fill = fct_rev(regulation))) +
    geom_col(width = 0.75) +
    geom_text(aes(label = n), vjust = -0.7, size = 2.2) +
    labs(x = "", y = "# of Differential Genes", title = "", fill = "") +
    scale_x_discrete(limits = c("1h", "3h", "6h", "9h", "12h", "24h")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme(
      axis.text.x = element_text(size = 7, color = "black"),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.title.y = element_text(size = 8),
      legend.position = "top",
      legend.text = element_text(size = 7),
      legend.title = element_blank(),
      legend.key.size = unit(0.35, "cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      plot.margin = margin(5, 5, 0, 5, "pt"),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    scale_fill_manual(values = c(color_upregulated, color_downregulated))
}

create_clustering_heatmap <- function(clustering_results) {
  ordered_matrix <- clustering_results$ordered_matrix
  row_order <- clustering_results$row_order
  km_object <- clustering_results$kmeans_object
  
  Heatmap(
    ordered_matrix,
    col = circlize::colorRamp2(c(-2, 0, 2), c("deepskyblue2", "black", "gold")),
    column_order = colnames(ordered_matrix),
    heatmap_legend_param = list(
      title = "", 
      at = c(-2, 0, 2),
      labels_gp = gpar(fontsize = 7),
      grid_width = unit(1.5, "mm"),  # Even thinner - reduced from 2mm to 1.5mm
      grid_height = unit(3, "mm"), # Height of each grid cell
      legend_height = unit(2, "cm"),
      legend_direction = "vertical"
    ),
    show_column_dend = FALSE, 
    show_row_dend = FALSE, 
    show_row_names = FALSE,
    column_names_gp = gpar(fontsize = 8), 
    column_names_rot = 0, 
    column_names_centered = TRUE,
    column_labels = c("0h", "1h", "3h", "6h", "9h", "12h", "24h"),
    width  = unit(panelA_width * 2.54, "cm"),  # Use full barplot width
    height = unit(panelA_heat_height * 2.54, "cm")
  )
}


# Panel B: GO Enrichment (Upregulated Only) ------------------------------

create_go_enrichment_plot <- function(dds_object, top_n_terms = 10) {
  # Parameters
  PADJ_THRESHOLD <- 0.05
  LFC_THRESHOLD <- 2
  
  # Get all time coefficients
  time_coefficients <- resultsNames(dds_object) |> 
    str_subset("Time")
  
  # Extract upregulated genes
  extract_upreg_genes <- function(coefficient) {
    res <- results(dds_object, name = coefficient) |>
      as.data.frame() |>
      filter(!is.na(padj))
    
    upreg_genes <- res |>
      filter(padj < PADJ_THRESHOLD & log2FoldChange > LFC_THRESHOLD) |>
      rownames()
    
    return(upreg_genes)
  }
  
  # Extract genes for all timepoints
  upreg_list <- map(time_coefficients, extract_upreg_genes)
  
  # Combine all genes
  all_upreg_genes <- upreg_list |>
    unlist() |>
    unique()
  
  # Get gene info
  gene_info <- rowData(dds_object) |>
    as.data.frame()
  
  if (!"gene_id" %in% colnames(gene_info)) {
    gene_info <- gene_info |>
      rownames_to_column("gene_id")
  }
  
  gene_info <- gene_info |>
    dplyr::select(gene_id, symbol)
  
  # Convert to Entrez IDs - Upregulated
  upreg_symbols <- gene_info |>
    filter(gene_id %in% all_upreg_genes) |>
    pull(symbol)
  
  upreg_entrez <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = upreg_symbols,
    keytype = "SYMBOL",
    columns = "ENTREZID"
  ) |>
    pull(ENTREZID) |>
    na.omit() |>
    unique()
  
  # Background
  background_symbols <- gene_info |>
    pull(symbol) |>
    na.omit()
  
  background_entrez <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = background_symbols,
    keytype = "SYMBOL",
    columns = "ENTREZID"
  ) |>
    pull(ENTREZID) |>
    na.omit() |>
    unique()
  
  # GO enrichment - Upregulated
  go_upreg <- enrichGO(
    gene = upreg_entrez,
    universe = background_entrez,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  # Manually select top 10 most relevant terms for osmotic stress
  selected_upreg_ids <- c(
    "GO:0006814",  # sodium ion transport (THE key osmotic pathway!)
    "GO:0009612",  # response to mechanical stimulus (osmotic = mechanical stress)
    "GO:0007159",  # leukocyte cell-cell adhesion (cell reorganization)
    "GO:0050727",  # regulation of inflammatory response (stress response)
    "GO:0001525",  # angiogenesis (tissue adaptation)
    "GO:0030198",  # extracellular matrix organization
    "GO:0006936",  # muscle contraction
    "GO:0042127",  # regulation of cell proliferation
    "GO:0001666",  # response to hypoxia
    "GO:0071356"   # cellular response to tumor necrosis factor
  )
  
  # Prepare plot data - Upregulated only
  go_upreg_data <- go_upreg@result |>
    as.data.frame() |>
    filter(ID %in% selected_upreg_ids) |>
    arrange(p.adjust) |>
    head(top_n_terms) |>
    mutate(
      log_padj = -log10(p.adjust),
      Description = str_trunc(Description, width = 70),
      Description = fct_reorder(Description, log_padj)
    )
  
  # Create plot with improved aesthetics
  ggplot(go_upreg_data, aes(x = log_padj, y = Description)) +
    geom_col(width = 0.75, fill = color_upregulated, color = NA) +
    geom_text(aes(x = 0.5, label = Description), 
              hjust = 0, size = 3, color = "black", fontface = "plain") +
    labs(
      x = expression(-log[10]("adjusted p-value")),
      y = NULL,
      title = "Upregulated GO Terms"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.3),
      axis.line.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10, "pt"),
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      legend.position = "none"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
}


# Panel C: EISA Line Plot (Introns solid + Exons dotted) -----------------

create_eisa_lineplot <- function(eisa_df) {
  
  # Loop category colors
  loop_colors <- c(
    "static" = "#999999",
    "gained" = "#F8766D",
    "lost" = "#619CFF"
  )
  
  # Prepare intron data (Δintron - transcriptional regulation)
  data_intron <- eisa_df |>
    filter(loop_category != "none") |>
    dplyr::select(gene_id, symbol, loop_category, starts_with("delta_intron_")) |>
    pivot_longer(
      cols = starts_with("delta_intron_"),
      names_to = "timepoint",
      values_to = "value"
    ) |>
    mutate(
      timepoint = str_extract(timepoint, "[0-9]+h") |> str_remove("h"),
      timepoint = factor(timepoint, levels = c("1", "3", "6", "9", "12", "24")),
      loop_category = factor(loop_category, levels = c("static", "gained", "lost")),
      metric = "intron"
    )
  
  # Prepare exon data (Δexon - total expression changes)
  data_exon <- eisa_df |>
    filter(loop_category != "none") |>
    dplyr::select(gene_id, symbol, loop_category, starts_with("delta_exon_")) |>
    pivot_longer(
      cols = starts_with("delta_exon_"),
      names_to = "timepoint",
      values_to = "value"
    ) |>
    mutate(
      timepoint = str_extract(timepoint, "[0-9]+h") |> str_remove("h"),
      timepoint = factor(timepoint, levels = c("1", "3", "6", "9", "12", "24")),
      loop_category = factor(loop_category, levels = c("static", "gained", "lost")),
      metric = "exon"
    )
  
  # Combine data
  combined_data <- bind_rows(data_intron, data_exon)
  
  # Calculate summary statistics
  data_summary <- combined_data |>
    group_by(timepoint, loop_category, metric) |>
    summarise(
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      se = sd(value, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    ) |>
    mutate(timepoint_num = as.numeric(as.character(timepoint)))
  
  # Create plot
  p <- ggplot(data_summary, aes(x = timepoint_num, y = median, color = loop_category)) +
    # Zero reference line
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.4) +
    
    # Intron lines FIRST (solid) - will be underneath
    geom_line(
      data = data_summary |> filter(metric == "intron"),
      aes(linetype = "intron"),
      linewidth = 0.9,  # Thinner
      alpha = 0.9
    ) +
    
    # Intron points with jitter
    geom_point(
      data = data_summary |> filter(metric == "intron"),
      size = 2.5,  # Smaller
      alpha = 0.9,
      position = position_jitter(width = 0.15, height = 0, seed = 123)
    ) +
    
    # Exon lines SECOND (dotted) - will be on top
    geom_line(
      data = data_summary |> filter(metric == "exon"),
      aes(linetype = "exon"),
      linewidth = 0.85,  # Thinner
      alpha = 0.85
    ) +
    
    # Exon points on top with jitter
    geom_point(
      data = data_summary |> filter(metric == "exon"),
      size = 2.0,  # Smaller
      alpha = 0.85,
      shape = 1,  # Open circles
      stroke = 1.0,  # Adjusted stroke for smaller size
      position = position_jitter(width = 0.15, height = 0, seed = 123)
    ) +
    
    # Direct text labels for loop categories in upper left (aligned at x=1)
    annotate("text", x = 1, y = 0.45, label = "Gained", 
             color = loop_colors["gained"], fontface = "bold", size = 3.5, hjust = 0) +
    annotate("text", x = 1, y = 0.35, label = "Static", 
             color = loop_colors["static"], fontface = "bold", size = 3.5, hjust = 0) +
    annotate("text", x = 1, y = 0.25, label = "Lost", 
             color = loop_colors["lost"], fontface = "bold", size = 3.5, hjust = 0) +
    
    # Line type labels with line segments (moved to the right)
    annotate("segment", x = 5, xend = 7, y = 0.45, yend = 0.45, 
             linetype = "solid", color = "black", linewidth = 0.8) +
    annotate("text", x = 7.5, y = 0.45, label = "nascent RNA", 
             color = "black", size = 3.5, hjust = 0) +
    
    annotate("segment", x = 5, xend = 7, y = 0.35, yend = 0.35, 
             linetype = "dotted", color = "black", linewidth = 0.8) +
    annotate("text", x = 7.5, y = 0.35, label = "mature mRNA", 
             color = "black", size = 3.5, hjust = 0) +
    
    # Manual line types
    scale_linetype_manual(
      name = "",
      values = c("intron" = "solid", "exon" = "dotted"),
      labels = c("intron" = "nascent RNA", "exon" = "mature mRNA")
    ) +
    
    # Colors
    scale_color_manual(values = loop_colors, name = "Loop category") +
    scale_fill_manual(values = loop_colors, guide = "none") +
    
    # X-axis
    scale_x_continuous(
      breaks = c(1, 3, 6, 9, 12, 24),
      labels = c("1", "3", "6", "9", "12", "24")
    ) +
    
    # Y-axis - labels at -0.5, 0, 0.5 but tick marks at all quarters
    scale_y_continuous(
      breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
      labels = c("-0.5", "", "0", "", "0.5"),
      minor_breaks = NULL
    ) +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    
    # Labels
    labs(
      x = "Hours after hyperosmotic stress",
      y = "log2FoldChange (treated/untreated)",
      title = ""
    ) +
    
    # Theme
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.major.x = element_line(color = "gray95", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 9),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      legend.position = "none",  # Remove legend
      plot.margin = margin(10, 5, 5, 5, "pt")
    ) +
    guides(
      color = "none",
      linetype = "none"
    )
  
  return(p)
}

# Panel D: Survey Plot ----------------------------------------------------

create_hic_loop_visualization <- function(target_loop_index = 313) {
  target_loop <- gainedLoops[target_loop_index]
  loop_region <- GRanges(
    seqnames = seqnames(anchors(target_loop, "first")),
    ranges = IRanges(start = start(anchors(target_loop, "first")),
                     end   = end(anchors(target_loop, "second")))
  )
  loop_region_buffed <- resize(loop_region, 
                               width = width(loop_region) + 2*100e3, 
                               fix = "center")
  
  list(
    chrom      = as.character(seqnames(loop_region_buffed)),
    chromstart = start(loop_region_buffed),
    chromend   = end(loop_region_buffed),
    assembly   = "hg38",
    resolution = 10e3,
    zrange     = c(0, 100),
    norm       = "SCALE",
    target_loop = target_loop,
    gained_loops = gainedLoops,
    lost_loops   = lostLoops,
    rna_colors = c("#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", 
                   "#225EA8", "#253494", "#081D58"),
    rna_timepoints = c("0h", "1h", "3h", "6h", "9h", "12h", "24h"),
    rna_files = c("STRS_HEK293_WT_cont_0h", "STRS_HEK293_WT_sorb_1h", 
                  "STRS_HEK293_WT_sorb_3h", "STRS_HEK293_WT_sorb_6h", 
                  "STRS_HEK293_WT_sorb_9h", "STRS_HEK293_WT_sorb_12h",
                  "STRS_HEK293_WT_sorb_24h"),
    h3k27ac_control_bw = "data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_H3K27ac_cont_0h.bw",
    h3k27ac_sorb_bw = "data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_H3K27ac_sorbitol_1h.bw"
  )
}


# Create all plots --------------------------------------------------------

diff_gene_plot <- create_differential_gene_plot(dds)
clustering_heatmap <- create_clustering_heatmap(clustering_data)
go_plot <- create_go_enrichment_plot(dds, top_n_terms = 10)
eisa_lineplot <- create_eisa_lineplot(eisa_results)
hic_params <- create_hic_loop_visualization(313)


# plotgardener visualization ----------------------------------------------

pdf(output_pdf, width = figure_width, height = figure_height)
pageCreate(width = figure_width, height = figure_height, showGuides = FALSE)

## PANEL A: Barplot + Heatmap (top-left)
plotText("A", x = col1_x - 0.16, y = row1_y - 0.15, 
         fontsize = 12, fontface = "bold", just = c("left", "top"))

# Move Panel A CONTENT UP to align x-axis labels with Panel B
# The "Time after hyperosmotic stress" should align with "-log10(adjusted p-value)"
# Also aligns the legend with "Upregulated Genes" title
# Negative offset moves content up
panelA_y_offset <- -0.15  # Reduced from -0.35 to move content down a bit
panelA_x_offset <- 0.15   # Shift barplot content to the right

# Barplot
plotGG(diff_gene_plot,
       x = col1_x + panelA_x_offset, y = row1_y + panelA_y_offset,
       width = panelA_width, height = panelA_bar_height,
       just = c("left", "top"))

# Heatmap - positioned to the right so "0h" starts where barplot x-axis begins
# Make it slightly narrower than full barplot width to account for plot margins
heatmap_x_offset <- 0.35  # Positive value moves right, towards Panel B
heatmap_plot_width <- 2.3  # Slightly less than panelA_width to match barplot plotting area
plotGG(plot = grid.grabExpr(draw(clustering_heatmap)),
       x = col1_x + heatmap_x_offset, y = row1_y + panelA_bar_height + panelA_y_offset,
       width = heatmap_plot_width, height = panelA_heat_height,
       just = c("left", "top"))

# X-axis label for Panel A - centered under heatmap
plotText("Time after hyperosmotic stress",
         x = col1_x + heatmap_x_offset + heatmap_plot_width/2, 
         y = row1_y + panelA_total_height + 0.08 + panelA_y_offset,
         fontsize = 9, just = c("center", "top"))


## PANEL B: GO Enrichment (top-right) - NOW UPREGULATED ONLY
plotText("B", x = col2_x - 0.12, y = row1_y - 0.15, 
         fontsize = 12, fontface = "bold", just = c("left", "top"))

# Panel B now has no facets, so we can adjust positioning for better alignment
panelB_y_offset <- 0.05

plotGG(go_plot,
       x = col2_x, y = row1_y + panelB_y_offset,
       width = panelB_width, height = top_row_height - panelB_y_offset,
       just = c("left", "top"))


## PANEL C: Combined Line Plot (middle, full width)
plotText("C", x = panelC_x - 0.16, y = row2_y - 0.15, 
         fontsize = 12, fontface = "bold", just = c("left", "top"))

# Add vertical offset for Panel C if needed
panelC_y_offset <- 0  # Adjust this value: positive moves down, negative moves up

plotGG(eisa_lineplot,
       x = panelC_x, y = row2_y + panelC_y_offset,
       width = panelC_width, height = panelC_height,
       just = c("left", "top"))


## PANEL D: Survey Plot (right side, full height)
plotText("D", x = col3_x - 0.16, y = row3_y - 0.15, 
         fontsize = 12, fontface = "bold", just = c("left", "top"))

# Survey plot dimensions
survey_total_h <- survey_height

# Layout fractions
frac_top_pad    <- 0.01
frac_hic_each   <- 0.20  # Reduced from 0.25 to make room
frac_hic_gap    <- 0.015
frac_rna_total  <- 0.35  # Keep RNA spacing
frac_gene       <- 0.12  # Increased for more space
frac_genome     <- 0.10  # Increased to push down
frac_bottom_pad <- 0.01

# Convert to inches
to_in <- function(frac) survey_total_h * frac
hic_h       <- to_in(frac_hic_each)
hic_gap_h   <- to_in(frac_hic_gap)
rna_block_h <- to_in(frac_rna_total)
gene_h      <- to_in(frac_gene)
genome_h    <- to_in(frac_genome)
top_pad_h   <- to_in(frac_top_pad)

# pgParams for survey plot
p <- pgParams(
  assembly   = hic_params$assembly,
  resolution = hic_params$resolution,
  chrom      = hic_params$chrom,
  chromstart = hic_params$chromstart,
  chromend   = hic_params$chromend,
  zrange     = hic_params$zrange,
  norm       = hic_params$norm,
  x = col3_x, 
  width = survey_width, 
  length = survey_width,
  height = hic_h,
  fontsize = 5
)

# 1) Untreated Hi-C
y0 <- row3_y + top_pad_h
control_hic <- plotHicRectangle(
  data = "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic",
  params = p, y = y0
)

## Add hg38 label to top-right corner
plotText(label = "hg38",
         x = col3_x + survey_width - 0.05,
         y = y0 + 0.05,
         fontsize = 6,
         fontcolor = "black",
         just = c("right", "top"),
         default.units = "inches")

annoPixels(control_hic, data = hic_params$lost_loops, 
           shift = 0.5, type = "arrow", col = color_downregulated)
annoPixels(control_hic, data = hic_params$gained_loops, 
           shift = 0.5, type = "arrow", col = color_upregulated)
plotText("untreated", x = col3_x, y = y0, just = c("left", "top"), fontsize = 8)

# Heatmap legend
annoHeatmapLegend(control_hic, orientation = "v", fontcolor = "black", digits = 1,
                  x = col3_x + survey_width + 0.05,
                  y = y0, width = 0.05, height = hic_h * 0.85,
                  default.units = "inches", fontsize = 6)

# 2) +Sorbitol Hi-C
y1 <- y0 + hic_h + hic_gap_h
sorb_hic <- plotHicRectangle(
  data = "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic",
  params = p, y = y1
)
annoPixels(sorb_hic, data = hic_params$lost_loops, 
           shift = 0.5, type = "arrow", col = color_downregulated)
annoPixels(sorb_hic, data = hic_params$gained_loops, 
           shift = 0.5, type = "arrow", col = color_upregulated)
plotText("+ sorbitol", x = col3_x, y = y1, just = c("left", "top"), fontsize = 8)

# H3K27ac tracks (control and sorbitol)
y_h3k27ac_start <- y1 + hic_h + 0.10  # Spacing after Hi-C maps

# Set H3K27ac signal range to fixed 0-125
h3k27ac_signalRange <- c(0, 125)

h3k27ac_track_h <- 0.22  # Increased to space out more
h3k27ac_gap <- 0.10       # Increased to space out more

# H3K27ac control track
plotSignal(
  data = hic_params$h3k27ac_control_bw,
  params = p,
  x = col3_x,
  y = y_h3k27ac_start,
  height = h3k27ac_track_h,
  linecolor = h3k27ac_color,
  fill = h3k27ac_color,
  scale = FALSE,
  range = h3k27ac_signalRange
)

# Label for control track - horizontal, above track
plotText("untreated",
         fontcolor = gray_color,
         x = col3_x,
         y = y_h3k27ac_start - 0.01,
         just = c("left", "bottom"),
         fontsize = 6)

# H3K27ac sorbitol track
y_h3k27ac_sorb <- y_h3k27ac_start + h3k27ac_track_h + h3k27ac_gap

plotSignal(
  data = hic_params$h3k27ac_sorb_bw,
  params = p,
  x = col3_x,
  y = y_h3k27ac_sorb,
  height = h3k27ac_track_h,
  linecolor = h3k27ac_color,
  fill = h3k27ac_color,
  scale = FALSE,
  range = h3k27ac_signalRange
)

# Label for sorbitol track - horizontal, above track
plotText("+ sorbitol",
         fontcolor = gray_color,
         x = col3_x,
         y = y_h3k27ac_sorb - 0.01,
         just = c("left", "bottom"),
         fontsize = 6)

# Add H3K27ac signal range label
if (!is.null(h3k27ac_signalRange)) {
  h3k27ac_min <- h3k27ac_signalRange[1]
  h3k27ac_max <- h3k27ac_signalRange[2]
  h3k27ac_range_text <- sprintf("%d-%d", round(h3k27ac_min), round(h3k27ac_max))
  
  plotText(
    label = h3k27ac_range_text,
    x = col3_x + survey_width,
    y = y_h3k27ac_start - 0.02,
    just = c("right", "bottom"),
    fontsize = 6,
    fontcolor = gray_color
  )
}

# Add H3K27ac label (vertical, aligned with RNA timepoint labels)
h3k27ac_mid_y <- y_h3k27ac_start + (h3k27ac_track_h + h3k27ac_gap + h3k27ac_track_h) / 2
plotText("H3K27ac",
         fontcolor = h3k27ac_color,
         rot = 90,
         y = h3k27ac_mid_y,
         x = col3_x - 0.12,
         fontsize = 8)

# 3) RNA tracks with custom signal range labels
y2 <- y_h3k27ac_sorb + h3k27ac_track_h + 0.02  # REDUCED from 0.10 to move RNA tracks (and genes/genome label) up to avoid overlap with timepoint labels
rna_n <- 7
rna_track_h <- 0.22  # Keep track height the same
rna_gap <- 0.02      # FURTHER REDUCED from 0.05 to squish tracks more

rna_files_all <- list.files(
  "data/processed/rna/timecourse/output/mergeSignal/stranded/",
  full.names = TRUE, pattern = "\\.bw$"
)

signalRange <- NULL
if (length(rna_files_all) > 0) {
  signalRange <- calcSignalRange(
    data = rna_files_all,
    chrom = hic_params$chrom, 
    chromstart = hic_params$chromstart, 
    chromend = hic_params$chromend,
    assembly = "hg38", 
    negData = FALSE
  )
}

# Plot RNA signal tracks
for (j in seq_len(rna_n)) {
  yj <- y2 + (j - 1) * (rna_track_h + rna_gap)
  
  fwd <- sprintf("data/processed/rna/timecourse/output/mergeSignal/stranded/%s_fwd.bw",
                 hic_params$rna_files[j])
  if (file.exists(fwd)) {
    plotSignal(data = fwd, params = p, x = col3_x, y = yj,
               height = rna_track_h, linecolor = hic_params$rna_colors[j],
               fill = hic_params$rna_colors[j], scale = FALSE, range = signalRange)
  }
  
  rev <- sprintf("data/processed/rna/timecourse/output/mergeSignal/stranded/%s_rev.bw",
                 hic_params$rna_files[j])
  if (file.exists(rev)) {
    plotSignal(data = rev, params = p, x = col3_x, y = yj,
               height = rna_track_h, linecolor = hic_params$rna_colors[j],
               scale = FALSE, range = signalRange)
  }
  
  # Timepoint label - horizontal, ABOVE track, matching track color
  plotText(hic_params$rna_timepoints[j], 
           fontcolor = hic_params$rna_colors[j],
           rot = 0, 
           y = yj + 0.08,  # Moved down even more - changed from yj + 0.02 to yj + 0.08
           x = col3_x, 
           just = c("left", "bottom"),
           fontsize = 6)
}

# Add vertical "RNA-seq" label (aligned with H3K27ac label)
rna_start_y <- y2
rna_end_y <- y2 + (rna_n * (rna_track_h + rna_gap)) - rna_gap
rna_mid_y <- rna_start_y + (rna_end_y - rna_start_y) / 2

plotText("RNA-seq",
         fontcolor = "#F4A6D7",  # Lighter pastel pink/magenta
         rot = 90,
         y = rna_mid_y,
         x = col3_x - 0.12,
         fontsize = 8)

# Add custom signal range label in top right (above first RNA track)
# This replaces the bracketed [0-1583] style labels
if (!is.null(signalRange)) {
  signal_min <- signalRange[1]
  signal_max <- signalRange[2]
  
  # Format the range text
  range_text <- sprintf("%d-%d", round(signal_min), round(signal_max))
  
  # Calculate position: align with "0h" label (first RNA track at y2 + 0.08)
  range_x <- col3_x + survey_width
  range_y <- y2 + 0.08  # Align with "0h" label position
  
  plotText(
    label = range_text,
    x = range_x,
    y = range_y,
    just = c("right", "bottom"),
    fontsize = 6,  # Reduced from 7 to match other labels
    fontcolor = gray_color  # Using gray_color to match Figure 3
  )
}

# 4) Genes - positioned right below last RNA track
last_rna_track_y <- y2 + (rna_n - 1) * (rna_track_h + rna_gap)
last_rna_track_bottom <- last_rna_track_y + rna_track_h
gene_track_height <- 0.25
rna_to_gene_gap <- 0.05  # Small gap between RNA and genes
y_genes <- last_rna_track_bottom + rna_to_gene_gap

plotGenes(
  params = p, 
  chrom = p$chrom, 
  x = col3_x, 
  y = y_genes, 
  height = gene_track_height,
  fontsize = 7,
  geneHighlights = data.frame(
    gene = "SOX8",
    color = color_upregulated
  )
)

# 5) Genome label - positioned right below genes
gene_to_genome_gap <- 0.05
genome_label_y <- y_genes + gene_track_height + gene_to_genome_gap

plotGenomeLabel(
  params = p, 
  x = col3_x, 
  y = genome_label_y,
  fontsize = 8
)

# 6) Vertical highlights
loop <- hic_params$target_loop
hl_y <- y0
hl_h <- genome_label_y - hl_y  # Removed the + 0.2 to stop at genome label

annoHighlight(
  plot = control_hic,
  fill = "lightgrey",
  alpha = 0.5,
  chrom = as.character(seqnames(anchors(loop, "first"))),
  chromstart = start(anchors(loop, "first")),
  chromend   = start(anchors(loop, "first")) + 10e3,
  x = col3_x, 
  y = hl_y, 
  width = survey_width,
  height = hl_h, 
  just = c("left", "top"),
  default.units = "inches"
)

annoHighlight(
  plot = control_hic,
  fill = "lightgrey",
  alpha = 0.5,
  chrom = as.character(seqnames(anchors(loop, "second"))),
  chromstart = end(anchors(loop, "second")) - 10e3,
  chromend   = end(anchors(loop, "second")),
  x = col3_x, 
  y = hl_y, 
  width = survey_width,
  height = hl_h, 
  just = c("left", "top"),
  default.units = "inches"
)

# Save outputs ----
dev.off()

sessionInfo()