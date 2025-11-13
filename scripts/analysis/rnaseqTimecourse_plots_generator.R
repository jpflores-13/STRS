## Generate MA plots and GO enrichment plots for RNA-seq timecourse
## This script creates differential expression visualizations and functional enrichment
## analyses for each timepoint in the sorbitol stress timecourse

## Load required libraries ----
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)

## Set parameters ----
# Thresholds for differential expression
LOG2FC_THRESHOLD <- 2
PADJ_THRESHOLD <- 0.05

# Number of genes to label on MA plots
TOP_GENES_TO_LABEL <- 30

# GO enrichment parameters
GO_MIN_GENESET <- 10
GO_MAX_GENESET <- 1000
GO_PADJ_CUTOFF <- 0.05

# Background gene set parameters
# This combines "expressed genes" with "sufficiently expressed genes"
# Only genes meeting this baseMean threshold will be included in GO background
MIN_BASEMEAN_FOR_BACKGROUND <- 10  # Adjust this as needed (common values: 5, 10, 20)

## Define timepoints ----
TIMEPOINTS <- c("1h", "3h", "6h", "9h", "12h", "24h")

## Load DESeq2 object ----
message("Loading DESeq2 object...")
dds <- readRDS("data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")

## Create output directories ----
dir.create("plots/rnaseqTimecourse", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed/rna/timecourse/output/GOobjs", showWarnings = FALSE, recursive = TRUE)

## Function to create MA plot ----
create_ma_plot <- function(res_df, timepoint, log2fc_thresh = LOG2FC_THRESHOLD, 
                           padj_thresh = PADJ_THRESHOLD, n_label = TOP_GENES_TO_LABEL) {
  
  # Classify genes
  res_df <- res_df |>
    mutate(
      regulation = case_when(
        padj < padj_thresh & log2FoldChange > log2fc_thresh ~ "upregulated",
        padj < padj_thresh & log2FoldChange < -log2fc_thresh ~ "downregulated",
        TRUE ~ "not significant"
      ),
      regulation = factor(regulation, 
                          levels = c("upregulated", "downregulated", "not significant"))
    ) |>
    arrange(regulation)  # Plot non-significant first so colored points are on top
  
  # Select top genes to label (most extreme log2FC among significant genes)
  genes_to_label <- res_df |>
    filter(regulation != "not significant") |>
    arrange(desc(abs(log2FoldChange))) |>
    head(n_label)
  
  # Count differential genes
  n_up <- sum(res_df$regulation == "upregulated", na.rm = TRUE)
  n_down <- sum(res_df$regulation == "downregulated", na.rm = TRUE)
  
  # Create color mapping
  colors <- c("upregulated" = "#F8766D",      
              "downregulated" = "#619CFF",    
              "not significant" = "grey80")   
  
  # Create MA plot combining clean styling with gene labels
  p <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = regulation)) +
    # Plot points with full alpha (no transparency)
    geom_point(alpha = 1, size = 1.0, stroke = 0) +
    scale_color_manual(values = colors) +
    # Dashed horizontal line at y=0
    geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
    # Log scale on x-axis
    scale_x_log10(
      breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000)
    ) +
    # Y-axis limits
    ylim(c(-4, 4)) +
    # Add gene labels with ggrepel
    geom_text_repel(
      data = genes_to_label,
      aes(label = SYMBOL),
      size = 3.2,
      color = "black",
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.3,
      segment.size = 0.25,
      segment.color = "black",
      min.segment.length = 0,
      seed = 42
    ) +
    # Axis labels and title
    labs(
      x = "mean of normalized counts",
      y = "log2FoldChange",
      title = paste0("+", timepoint, "sorbitol")
    ) +
    # Clean theme
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0, size = 11, face = "plain")
    ) +
    # Add text annotations for upregulated genes (bottom right of upper quadrant)
    annotate(
      geom = "text",
      label = "Upregulated",
      x = 10000,
      y = 0.5,
      color = "#F8766D",
      size = 3.5,
      hjust = 1
    ) +
    annotate(
      geom = "text",
      label = paste0("n = ", n_up),
      x = 10000,
      y = 0.1,
      color = "#F8766D",
      size = 3.5,
      hjust = 1
    ) +
    # Add text annotations for downregulated genes (bottom right of lower quadrant)
    annotate(
      geom = "text",
      label = "Downregulated",
      x = 10000,
      y = -0.5,
      color = "#619CFF",
      size = 3.5,
      hjust = 1
    ) +
    annotate(
      geom = "text",
      label = paste0("n = ", n_down),
      x = 10000,
      y = -0.9,
      color = "#619CFF",
      size = 3.5,
      hjust = 1
    )
  
  return(p)
}

## Function to perform GO enrichment ----
perform_go_enrichment <- function(gene_list, universe, timepoint, direction) {
  
  message(sprintf("  Running GO enrichment for %s genes at %s...", direction, timepoint))
  
  # Convert ENSEMBL to ENTREZ
  entrez_ids <- mapIds(
    org.Hs.eg.db,
    keys = gene_list,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  message(sprintf("    Mapped %d out of %d genes to ENTREZ IDs", 
                  length(entrez_ids), length(gene_list)))
  
  # Skip if too few genes
  if (length(entrez_ids) < 10) {
    message(sprintf("    Skipping: fewer than 10 genes mapped"))
    return(NULL)
  }
  
  # Perform GO enrichment
  ego <- enrichGO(
    gene = entrez_ids,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = GO_MIN_GENESET,
    maxGSSize = GO_MAX_GENESET,
    pAdjustMethod = "BH",
    pvalueCutoff = GO_PADJ_CUTOFF
  )
  
  # Check if results exist
  if (is.null(ego) || nrow(ego@result) == 0) {
    message(sprintf("    No enriched GO terms found"))
    return(NULL)
  }
  
  # Simplify terms
  ego_simplified <- clusterProfiler::simplify(
    ego,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  message(sprintf("    Found %d enriched GO terms", nrow(ego_simplified@result)))
  
  return(ego_simplified)
}

## Function to create GO barplot ----
create_go_barplot <- function(ego, timepoint, direction, n_terms = 8) {
  
  if (is.null(ego) || nrow(ego@result) == 0) {
    return(NULL)
  }
  
  # Get top terms
  top_terms <- ego@result |>
    slice_min(order_by = p.adjust, n = n_terms) |>
    as.data.frame()
  
  if (nrow(top_terms) == 0) {
    return(NULL)
  }
  
  # Calculate qscore (enrichment score)
  # qscore is defined as -log10(p.adjust) * (Count/GeneRatio denominator)
  top_terms <- top_terms |>
    mutate(
      qscore = -log10(p.adjust) * (Count / as.numeric(str_extract(GeneRatio, "(?<=/)\\d+"))),
      Description = str_wrap(Description, width = 40)
    )
  
  # Determine color based on direction
  bar_color <- if (direction == "upregulated") "#F8766D" else "#619CFF"
  
  # Create plot
  p <- ggplot(top_terms, aes(x = qscore, y = reorder(Description, qscore))) +
    geom_col(fill = bar_color, alpha = 0.8) +
    scale_fill_gradient(low = bar_color, high = bar_color) +
    labs(
      x = "qscore",
      y = "",
      title = paste0("+", timepoint, "sorbitol")
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 9, color = if (direction == "upregulated") "#F8766D" else "#619CFF"),
      plot.title = element_text(hjust = 0, size = 11),
      axis.title.x = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

## Main analysis loop ----
message("\n=== Starting timecourse analysis ===\n")

# Get universe for GO enrichment using COMBINED APPROACH:
# 1. Genes that passed DESeq2 filtering (expressed in HEK293s)
# 2. Genes with baseMean above threshold (sufficient expression for DE detection)
#
# This is the most biologically relevant and statistically sound approach because:
# - Only includes genes actually expressed in your system (HEK293s)
# - Only includes genes with sufficient counts to detect differential expression
# - Removes very lowly expressed genes that lack statistical power
# - Provides the most appropriate comparison set for enrichment testing

message("Preparing background gene set for GO enrichment...")
message(sprintf("  Strategy: Expressed genes with baseMean >= %.1f", 
                MIN_BASEMEAN_FOR_BACKGROUND))

# Get all genes that passed filtering
expressed_ensembl <- rownames(dds)

# Calculate baseMean for each gene
basemeans <- rowMeans(counts(dds, normalized = TRUE))

# Filter to genes above expression threshold
well_expressed_ensembl <- names(basemeans[basemeans >= MIN_BASEMEAN_FOR_BACKGROUND])

message(sprintf("  Total genes in DESeq2 object: %d", length(expressed_ensembl)))
message(sprintf("  Genes with baseMean >= %.1f: %d", 
                MIN_BASEMEAN_FOR_BACKGROUND, length(well_expressed_ensembl)))

# Convert to ENTREZ IDs for GO analysis
expressed_entrez <- mapIds(
  org.Hs.eg.db,
  keys = well_expressed_ensembl,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

# Remove NAs and duplicates
expressed_entrez <- expressed_entrez[!is.na(expressed_entrez)]
expressed_entrez <- unique(expressed_entrez)

message(sprintf("  Final background set: %d genes (ENTREZ IDs)", 
                length(expressed_entrez)))
message(sprintf("  (%.1f%% of total genes in object)\n", 
                100 * length(well_expressed_ensembl) / length(expressed_ensembl)))

# Use well-expressed genes as universe
all_genes <- expressed_entrez

# Get coefficient names
coef_names <- resultsNames(dds) |>
  str_subset("Time")

message(sprintf("Found %d timepoint comparisons\n", length(coef_names)))

## Process each timepoint ----
for (i in seq_along(TIMEPOINTS)) {
  
  timepoint <- TIMEPOINTS[i]
  coef <- coef_names[i]
  
  message(sprintf("Processing timepoint: %s (coefficient: %s)", timepoint, coef))
  
  # Extract results
  res <- results(dds, name = coef)
  
  # Convert to dataframe with gene IDs
  res_df <- res |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "ENSEMBL") |>
    filter(!is.na(padj))  # Remove genes with NA padj
  
  # Map ENSEMBL to SYMBOL
  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = res_df$ENSEMBL,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  res_df$SYMBOL <- gene_symbols[res_df$ENSEMBL]
  
  # Identify upregulated and downregulated genes
  upreg_genes <- res_df |>
    filter(padj < PADJ_THRESHOLD, log2FoldChange > LOG2FC_THRESHOLD) |>
    pull(ENSEMBL)
  
  downreg_genes <- res_df |>
    filter(padj < PADJ_THRESHOLD, log2FoldChange < -LOG2FC_THRESHOLD) |>
    pull(ENSEMBL)
  
  message(sprintf("  Found %d upregulated and %d downregulated genes", 
                  length(upreg_genes), length(downreg_genes)))
  
  ## Create and save MA plot ----
  message("  Creating MA plot...")
  ma_plot <- create_ma_plot(res_df, timepoint)
  
  pdf(paste0("plots/rnaseqTimecourse/", timepoint, "_annotatedMA.pdf"), 
      width = 8, height = 6)
  print(ma_plot)
  dev.off()
  
  ## GO enrichment for upregulated genes ----
  if (length(upreg_genes) >= 10) {
    ego_up <- perform_go_enrichment(upreg_genes, all_genes, timepoint, "upregulated")
    
    if (!is.null(ego_up)) {
      # Save GO object
      saveRDS(ego_up, 
              paste0("data/processed/rna/timecourse/output/GOobjs/", 
                     timepoint, "_upregulated_GO_results.rds"))
      
      # Create and save barplot
      go_plot_up <- create_go_barplot(ego_up, timepoint, "upregulated")
      
      if (!is.null(go_plot_up)) {
        pdf(paste0("plots/rnaseqTimecourse/", timepoint, "_upregGO.pdf"), 
            width = 8, height = 6)
        print(go_plot_up)
        dev.off()
      }
    }
  } else {
    message(sprintf("  Skipping GO enrichment for upregulated genes: only %d genes", 
                    length(upreg_genes)))
  }
  
  ## GO enrichment for downregulated genes ----
  if (length(downreg_genes) >= 10) {
    ego_down <- perform_go_enrichment(downreg_genes, all_genes, timepoint, "downregulated")
    
    if (!is.null(ego_down)) {
      # Save GO object
      saveRDS(ego_down, 
              paste0("data/processed/rna/timecourse/output/GOobjs/", 
                     timepoint, "_downregulated_GO_results.rds"))
      
      # Create and save barplot
      go_plot_down <- create_go_barplot(ego_down, timepoint, "downregulated")
      
      if (!is.null(go_plot_down)) {
        pdf(paste0("plots/rnaseqTimecourse/", timepoint, "_downregGO.pdf"), 
            width = 8, height = 6)
        print(go_plot_down)
        dev.off()
      }
    }
  } else {
    message(sprintf("  Skipping GO enrichment for downregulated genes: only %d genes", 
                    length(downreg_genes)))
  }
  
  message(sprintf("Completed timepoint: %s\n", timepoint))
}

message("\n=== Analysis complete! ===")
message("\nGenerated files:")
message("  MA plots: plots/rnaseqTimecourse/*_annotatedMA.pdf")
message("  GO plots: plots/rnaseqTimecourse/*_upregGO.pdf and *_downregGO.pdf")
message("  GO objects: data/processed/rna/timecourse/output/GOobjs/*_GO_results.rds")