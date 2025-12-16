################################################################################
# Survey Plots for Genes at Gained Loop Anchors
# 
# Creates plotgardener survey plots showing Hi-C contact maps with RNA-seq
# timecourse data for genes that overlap GAINED loop anchors and show
# significant temporal expression changes (LRT padj < 0.1)
# Ordered by LRT padj (most significant first)
# Centered on the gained loop to show both anchors
################################################################################

# Load packages -----------------------------------------------------------

library(plotgardener)
library(InteractionSet)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DESeq2)
library(RColorBrewer)
library(GenomicRanges)

# Load data ---------------------------------------------------------------

# Load DESeq2 object
dds <- readRDS("data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")

# Load differential loops
diff_loopCounts <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions()
diff_loopCounts <- keepStandardChromosomes(diff_loopCounts, pruning.mode = "coarse")

# Prepare gene data -------------------------------------------------------

# Get DGE results as GRanges with gene coordinates
dge_gr <- results(dds, format = "GRanges")
mcols(dge_gr) <- rowData(dds)

# Get promoter regions
dge_gr_prom <- promoters(dge_gr) |> 
  keepStandardChromosomes(pruning.mode = "coarse")
seqlevelsStyle(dge_gr_prom) <- "UCSC"

# Calculate LRT adjusted p-values
dge_gr_prom$padj_lrt <- p.adjust(dge_gr_prom$LRTPvalue, method = "BH")

# Define loop types -------------------------------------------------------

lostLoops <- diff_loopCounts[which(diff_loopCounts$padj < 0.1 & 
                                     diff_loopCounts$log2FoldChange < 0)] 

gainedLoops <- diff_loopCounts[which(diff_loopCounts$padj < 0.1 & 
                                       diff_loopCounts$log2FoldChange > 0)]

# Convert gained loops to GInteractions for overlap
as_gi <- function(x) x |> as.data.frame() |> as_ginteractions()
gained_gi <- as_gi(gainedLoops)

# Find genes overlapping GAINED loops only --------------------------------

gained_genes <- subsetByOverlaps(dge_gr_prom, gained_gi) |> 
  subset(padj_lrt < 0.1)

# Add loop type annotation
gained_genes$loop_type <- "gained"

# Sort genes by LRT padj (most significant first)
gained_genes <- gained_genes[order(gained_genes$padj_lrt)]

message("Found ", length(gained_genes), " genes at gained loop anchors with significant temporal expression changes")
message("Sorted by LRT padj (most significant first)")
message("Top 5 genes: ", paste(head(gained_genes$symbol, 5), collapse = ", "))

# For each gene, find which loop anchors it overlaps
# This will help us highlight the correct anchors later
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

# Store anchor information for each gene
anchor_info <- lapply(seq_along(gained_genes), function(i) {
  find_overlapping_anchors(gained_genes[i], gained_gi)
})

# Create regions for visualization centered on the loop ------------------

# Function to calculate loop center and create buffered region
create_loop_centered_region <- function(gene_idx, anchor_info, loops_gi, buffer = 200e3) {
  # Get the loop(s) this gene overlaps
  loop_indices <- anchor_info[[gene_idx]]$all_loops
  
  # If multiple loops, use the first one (could modify to choose differently)
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

# Create centered regions for each gene
centered_regions <- lapply(seq_along(gained_genes), function(i) {
  create_loop_centered_region(i, anchor_info, gained_gi, buffer = 200e3)
})

# Create Survey Plots -----------------------------------------------------

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

# Create PDF
pdf(
  file = "plots/surveyPlot_gainedLoopGenes_rnaTimecourse_10kb.pdf",
  width = 5.75,
  height = 9.1
)

# Loop through each gene region
for(i in seq_along(gained_genes)) {
  
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
  
  # Create page
  pageCreate(width = 5.75, height = 9.1, showGuides = FALSE)
  
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
  
  # Plot RNA-seq timecourse tracks
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
      fwd_file <- glue::glue("data/processed/rna/timecourse/output/mergeSignal/stranded/STRS_HEK293_WT_sorb_{tp}_fwd.bw")
      rev_file <- glue::glue("data/processed/rna/timecourse/output/mergeSignal/stranded/STRS_HEK293_WT_sorb_{tp}_rev.bw")
    }
    
    # Plot forward strand
    plotSignal(
      fwd_file,
      params = p,
      x = 0.25,
      y = y_pos,
      height = 0.25,
      linecolor = color,
      scale = TRUE,
      range = signalRange
    )
    
    # Plot reverse strand
    plotSignal(
      rev_file,
      params = p,
      x = 0.25,
      y = y_pos,
      height = 0.25,
      linecolor = color,
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
  
  # Plot genes
  plotGenes(
    param = p,
    chrom = p$chrom,
    x = 0.25,
    y = 8.25,
    height = 0.5
  )
  
  # Plot genome label
  plotGenomeLabel(
    params = p,
    x = 0.25,
    y = 8.85
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
    height = 9.1,
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
    height = 9.1,
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
    height = 9.1,
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
    height = 9.1,
    just = c("left", "top")
  )
  
  # Add gene symbol and padj as title
  gene_symbol <- gained_genes[i]$symbol
  gene_padj <- gained_genes[i]$padj_lrt
  
  plotText(
    label = glue::glue("{gene_symbol} (gained loop, LRT padj = {format(gene_padj, digits = 3, scientific = TRUE)})"),
    x = 2.875,
    y = 0.1,
    just = "center",
    fontface = "bold",
    fontsize = 9
  )
}

dev.off()

# Summary message ---------------------------------------------------------

message("\nSurvey plots created for ", length(gained_genes), " genes at gained loop anchors")
message("Ordered by LRT padj (most significant first)")
message("Plots centered on gained loop midpoint to show both anchors")
message("Output saved to: plots/surveyPlot_gainedLoopGenes_rnaTimecourse_10kb.pdf")