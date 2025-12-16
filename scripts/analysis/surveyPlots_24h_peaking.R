# Survey Plots for Genes with 24h-Peaking Expression at Gained Loop Anchors

# Load packages -----------------------------------------------------------

library(plotgardener)
library(InteractionSet)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DESeq2)
library(RColorBrewer)
library(GenomicRanges)
library(glue)

# Load data ---------------------------------------------------------------

# Load DESeq2 object
dds <- readRDS("data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")

# Load differential loops
diff_loopCounts <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions()
diff_loopCounts <- keepStandardChromosomes(diff_loopCounts, pruning.mode = "coarse")

# Define gained loops -----------------------------------------------------

gainedLoops <- diff_loopCounts[which(diff_loopCounts$padj < 0.1 & 
                                       diff_loopCounts$log2FoldChange > 0)]

lostLoops <- diff_loopCounts[which(diff_loopCounts$padj < 0.1 & 
                                     diff_loopCounts$log2FoldChange < 0)]

# Convert gained loops to GInteractions
as_gi <- function(x) x |> as.data.frame() |> as_ginteractions()
gained_gi <- as_gi(gainedLoops)

# Prepare gene data -------------------------------------------------------

# Get DGE results as GRanges with gene coordinates
dge_gr <- results(dds, format = "GRanges")
mcols(dge_gr) <- rowData(dds)

# Get promoter regions
dge_gr_prom <- promoters(dge_gr) |> 
  keepStandardChromosomes(pruning.mode = "coarse")
seqlevelsStyle(dge_gr_prom) <- "UCSC"

# Find genes overlapping GAINED loops -------------------------------------

genes_at_gained_loops <- subsetByOverlaps(dge_gr_prom, gained_gi)

message("Found ", length(genes_at_gained_loops), " genes at gained loop anchors")

# Extract Wald test results for all timepoints ---------------------------

# Get Wald test results for each timepoint
wald_1h <- results(dds, name = "Time_1h_vs_0h")
wald_3h <- results(dds, name = "Time_3h_vs_0h")
wald_6h <- results(dds, name = "Time_6h_vs_0h")
wald_9h <- results(dds, name = "Time_9h_vs_0h")
wald_12h <- results(dds, name = "Time_12h_vs_0h")
wald_24h <- results(dds, name = "Time_24h_vs_0h")

# Add Wald results to genes_at_gained_loops
mcols(genes_at_gained_loops)$wald_1h_log2FC <- wald_1h[names(genes_at_gained_loops), "log2FoldChange"]
mcols(genes_at_gained_loops)$wald_3h_log2FC <- wald_3h[names(genes_at_gained_loops), "log2FoldChange"]
mcols(genes_at_gained_loops)$wald_6h_log2FC <- wald_6h[names(genes_at_gained_loops), "log2FoldChange"]
mcols(genes_at_gained_loops)$wald_9h_log2FC <- wald_9h[names(genes_at_gained_loops), "log2FoldChange"]
mcols(genes_at_gained_loops)$wald_12h_log2FC <- wald_12h[names(genes_at_gained_loops), "log2FoldChange"]
mcols(genes_at_gained_loops)$wald_24h_log2FC <- wald_24h[names(genes_at_gained_loops), "log2FoldChange"]
mcols(genes_at_gained_loops)$wald_24h_padj <- wald_24h[names(genes_at_gained_loops), "padj"]

# Selection Strategy: Genes peaking at 24h --------------------------------

# Filter for genes that:
# 1. Are significantly upregulated at 24h (padj < 0.1, log2FC > 0)
# 2. Peak at 24h (24h log2FC is highest across all timepoints)

genes_24h_peak <- genes_at_gained_loops[
  # Must be significant at 24h
  !is.na(genes_at_gained_loops$wald_24h_padj) &
    !is.na(genes_at_gained_loops$wald_24h_log2FC) &
    genes_at_gained_loops$wald_24h_padj < 0.1 &
    genes_at_gained_loops$wald_24h_log2FC > 0 &
    
    # 24h must be higher than all earlier timepoints (continuously rising)
    !is.na(genes_at_gained_loops$wald_1h_log2FC) &
    !is.na(genes_at_gained_loops$wald_3h_log2FC) &
    !is.na(genes_at_gained_loops$wald_6h_log2FC) &
    !is.na(genes_at_gained_loops$wald_9h_log2FC) &
    !is.na(genes_at_gained_loops$wald_12h_log2FC) &
    genes_at_gained_loops$wald_24h_log2FC > genes_at_gained_loops$wald_1h_log2FC &
    genes_at_gained_loops$wald_24h_log2FC > genes_at_gained_loops$wald_3h_log2FC &
    genes_at_gained_loops$wald_24h_log2FC > genes_at_gained_loops$wald_6h_log2FC &
    genes_at_gained_loops$wald_24h_log2FC > genes_at_gained_loops$wald_9h_log2FC &
    genes_at_gained_loops$wald_24h_log2FC > genes_at_gained_loops$wald_12h_log2FC
]

# Rank by 24h log2FC (highest first)
genes_24h_peak <- genes_24h_peak[order(-genes_24h_peak$wald_24h_log2FC)]

# Limit to 150 genes maximum
genes_24h_peak <- head(genes_24h_peak, 150)

message("\nGenes peaking at 24h: ", length(genes_24h_peak), " genes")

if(length(genes_24h_peak) > 0) {
  message("\nTop 10 genes (by 24h log2FC):")
  message("  ", paste(head(genes_24h_peak$symbol, 10), collapse = ", "))
  message("\nlog2FC at 24h for top gene: ", round(genes_24h_peak$wald_24h_log2FC[1], 2))
}

# Helper functions --------------------------------------------------------

# For each gene, find which loop anchors it overlaps
find_overlapping_anchors <- function(gene_gr, loops_gi) {
  # Extract anchors as GRanges
  anchor1 <- anchors(loops_gi, "first")
  anchor2 <- anchors(loops_gi, "second")
  
  # Find overlaps with each anchor
  ov1 <- findOverlaps(gene_gr, anchor1)
  ov2 <- findOverlaps(gene_gr, anchor2)
  
  # Get the loop indices that overlap
  loop_idx1 <- subjectHits(ov1)
  loop_idx2 <- subjectHits(ov2)
  
  # Store which anchors overlap for this gene
  list(
    anchor1_loops = loop_idx1,
    anchor2_loops = loop_idx2,
    all_loops = unique(c(loop_idx1, loop_idx2))
  )
}

# Function to calculate loop center and create buffered region
create_loop_centered_region <- function(gene_idx, anchor_info, loops_gi, buffer = 400e3) {
  # Get the loop(s) this gene overlaps
  loop_indices <- anchor_info[[gene_idx]]$all_loops
  
  # If multiple loops, use the first one
  main_loop <- loops_gi[loop_indices[1]]
  
  # Get anchor positions
  anchor1 <- anchors(main_loop, "first")
  anchor2 <- anchors(main_loop, "second")
  
  # Calculate loop midpoint
  loop_start <- start(anchor1)
  loop_end <- end(anchor2)
  loop_midpoint <- round((loop_start + loop_end) / 2)
  
  # Create centered region
  region_start <- loop_midpoint - buffer
  region_end <- loop_midpoint + buffer
  
  # Return as a list with region coordinates and the main loop
  list(
    chrom = as.character(seqnames(anchor1)),
    start = region_start,
    end = region_end,
    main_loop_idx = loop_indices[1]
  )
}

# Function to create survey plots -----------------------------------------

create_survey_plots <- function(genes, output_file) {
  
  message("\nCreating survey plots for 24h-peaking genes...")
  
  if(length(genes) == 0) {
    message("  No genes found. Skipping.")
    return(invisible(NULL))
  }
  
  # Store anchor information for each gene
  anchor_info <- lapply(seq_along(genes), function(i) {
    find_overlapping_anchors(genes[i], gained_gi)
  })
  
  # Create centered regions for each gene
  centered_regions <- lapply(seq_along(genes), function(i) {
    create_loop_centered_region(i, anchor_info, gained_gi, buffer = 400e3)
  })
  
  # RNA-seq timecourse colors (YlGnBu palette)
  timepoint_colors <- c(
    "0h" = "#C7E9B4",
    "1h" = "#7FCDBB",
    "3h" = "#41B6C4",
    "6h" = "#1D91C0",
    "9h" = "#225EA8",
    "12h" = "#253494",
    "24h" = "#081D58"
  )
  
  # Define bigWig files explicitly (only .bw files)
  rna_files <- list.files(
    "data/processed/rna/timecourse/output/mergeSignal/stranded/",
    pattern = "\\.bw$",
    full.names = TRUE
  )
  
  # Create PDF with increased height
  pdf(file = output_file, width = 5.75, height = 10.5)
  
  # Loop through each gene
  for(i in seq_along(genes)) {
    
    # Get the centered region for this gene
    region <- centered_regions[[i]]
    
    # Define parameters for this region (centered on loop)
    p <- pgParams(
      assembly = "hg38",
      resolution = 10e3,
      chrom = region$chrom,
      chromstart = region$start,
      chromend = region$end,
      zrange = c(0, 100),
      norm = "SCALE",
      x = 0.25,
      width = 5,
      length = 5,
      height = 2
    )
    
    # Create page with increased height
    pageCreate(width = 5.75, height = 10.5, showGuides = FALSE)
    
    # Add info as title at the very top
    gene_symbol <- genes[i]$symbol
    gene_padj <- genes[i]$wald_24h_padj
    gene_lfc_24h <- genes[i]$wald_24h_log2FC
    
    # Show the trajectory pattern in the title
    gene_lfc_1h <- genes[i]$wald_1h_log2FC
    gene_lfc_6h <- genes[i]$wald_6h_log2FC
    gene_lfc_12h <- genes[i]$wald_12h_log2FC
    
    plotText(
      label = glue("{gene_symbol} | Peaks at 24h | log2FC: 1h={round(gene_lfc_1h, 2)}, 6h={round(gene_lfc_6h, 2)}, 12h={round(gene_lfc_12h, 2)}, 24h={round(gene_lfc_24h, 2)} | padj={format(gene_padj, digits=2, scientific=TRUE)}"),
      x = 2.875,
      y = 0.1,
      just = "center",
      fontface = "bold",
      fontsize = 7
    )
    
    # Plot control Hi-C rectangle
    control <- plotHicRectangle(
      data = "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic",
      params = p,
      y = 0.25
    )
    
    plotText(
      label = "untreated",
      x = 0.25,
      y = 0.25,
      just = c("top", "left")
    )
    
    annoHeatmapLegend(
      control,
      orientation = "v",
      fontsize = 8,
      fontcolor = "black",
      digits = 2,
      x = 5.5,
      y = 0.25,
      width = 0.1,
      height = 1.5,
      just = c("left", "top"),
      default.units = "inches"
    )
    
    # Annotate loops on control
    annoPixels(control, data = lostLoops, shift = 0.5, 
               type = "arrow", col = "#005AB5")
    annoPixels(control, data = gainedLoops, shift = 0.5, 
               type = "arrow", col = "#DC3220")
    
    # Plot sorbitol Hi-C rectangle
    sorb <- plotHicRectangle(
      data = "data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic",
      params = p,
      y = 2.35
    )
    
    plotText(
      label = "+ sorbitol",
      x = 0.25,
      y = 2.35,
      just = c("top", "left")
    )
    
    annoHeatmapLegend(
      sorb,
      orientation = "v",
      fontsize = 8,
      fontcolor = "black",
      digits = 2,
      x = 5.5,
      y = 2.35,
      width = 0.1,
      height = 1.5,
      just = c("left", "top"),
      default.units = "inches"
    )
    
    # Annotate loops on sorbitol
    annoPixels(sorb, data = lostLoops, shift = 0.5, 
               type = "arrow", col = "#005AB5")
    annoPixels(sorb, data = gainedLoops, shift = 0.5, 
               type = "arrow", col = "#DC3220")
    
    # Calculate signal range for RNA-seq tracks
    signalRange <- calcSignalRange(
      data = rna_files,
      chrom = region$chrom,
      chromstart = region$start,
      chromend = region$end,
      assembly = "hg38",
      negData = FALSE
    )
    
    # Plot RNA-seq timecourse tracks with fill
    y_positions <- c(4.4, 5.0, 5.6, 6.2, 6.8, 7.4, 8.0)
    timepoints <- c("0h", "1h", "3h", "6h", "9h", "12h", "24h")
    
    for(j in seq_along(timepoints)) {
      tp <- timepoints[j]
      y_pos <- y_positions[j]
      color <- timepoint_colors[tp]
      
      # Construct file paths
      if(tp == "0h") {
        fwd_file <- "data/processed/rna/timecourse/output/mergeSignal/stranded/STRS_HEK293_WT_cont_0h_fwd.bw"
        rev_file <- "data/processed/rna/timecourse/output/mergeSignal/stranded/STRS_HEK293_WT_cont_0h_rev.bw"
      } else {
        fwd_file <- glue("data/processed/rna/timecourse/output/mergeSignal/stranded/STRS_HEK293_WT_sorb_{tp}_fwd.bw")
        rev_file <- glue("data/processed/rna/timecourse/output/mergeSignal/stranded/STRS_HEK293_WT_sorb_{tp}_rev.bw")
      }
      
      # Plot forward strand with fill
      plotSignal(
        fwd_file,
        params = p,
        x = 0.25,
        y = y_pos,
        height = 0.25,
        linecolor = color,
        fill = color,
        scale = TRUE,
        range = signalRange
      )
      
      # Plot reverse strand with fill
      plotSignal(
        rev_file,
        params = p,
        x = 0.25,
        y = y_pos,
        height = 0.25,
        linecolor = color,
        fill = color,
        scale = TRUE,
        range = signalRange
      )
      
      # Add timepoint label
      plotText(
        tp,
        fontcolor = color,
        rot = 90,
        y = y_pos + 0.1,
        x = 0.1
      )
    }
    
    # Plot genes with geneHighlights
    plotGenes(
      param = p,
      chrom = p$chrom,
      x = 0.25,
      y = 8.25,
      height = 1.0,
      fontsize = 8,
      geneHighlights = data.frame(
        gene = gene_symbol,
        color = "#DC3220"
      )
    )
    
    # Plot genome label
    plotGenomeLabel(
      params = p,
      x = 0.25,
      y = 9.35
    )
    
    # Highlight BOTH anchors of the main gained loop
    main_loop <- gained_gi[region$main_loop_idx]
    anchor1 <- anchors(main_loop, "first")
    anchor2 <- anchors(main_loop, "second")
    
    # Highlight anchor1
    annoHighlight(
      plot = control,
      fill = "lightgrey",
      chrom = as.character(seqnames(anchor1)),
      chromstart = start(anchor1),
      chromend = end(anchor1),
      default.units = "inches",
      y = 0.25,
      x = 0.25,
      height = 10.5,
      just = c("left", "top")
    )
    
    annoHighlight(
      plot = sorb,
      fill = "lightgrey",
      chrom = as.character(seqnames(anchor1)),
      chromstart = start(anchor1),
      chromend = end(anchor1),
      default.units = "inches",
      y = 0.25,
      x = 0.25,
      height = 10.5,
      just = c("left", "top")
    )
    
    # Highlight anchor2
    annoHighlight(
      plot = control,
      fill = "lightgrey",
      chrom = as.character(seqnames(anchor2)),
      chromstart = start(anchor2),
      chromend = end(anchor2),
      default.units = "inches",
      y = 0.25,
      x = 0.25,
      height = 10.5,
      just = c("left", "top")
    )
    
    annoHighlight(
      plot = sorb,
      fill = "lightgrey",
      chrom = as.character(seqnames(anchor2)),
      chromstart = start(anchor2),
      chromend = end(anchor2),
      default.units = "inches",
      y = 0.25,
      x = 0.25,
      height = 10.5,
      just = c("left", "top")
    )
  }
  
  dev.off()
  
  message("  Saved: ", output_file)
}

# Create survey plots -----------------------------------------------------

create_survey_plots(
  genes = genes_24h_peak,
  output_file = "plots/surveyPlots_24h_peaking.pdf"
)

# Summary -----------------------------------------------------------------

message("\n=== Summary ===")
message("Survey plots created for genes peaking at 24h at gained loop anchors:")
message("  Total genes: ", length(genes_24h_peak))
message("  Output: plots/surveyPlots_24h_peaking.pdf")

if(length(genes_24h_peak) > 0) {
  message("\nTop 10 genes (by 24h log2FC):")
  for(i in 1:min(10, length(genes_24h_peak))) {
    message(sprintf("  %2d. %s (log2FC@24h = %.2f)", 
                    i, 
                    genes_24h_peak$symbol[i], 
                    genes_24h_peak$wald_24h_log2FC[i]))
  }
}